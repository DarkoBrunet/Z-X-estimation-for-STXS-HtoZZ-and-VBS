#!/usr/bin/env python3
"""Estimate the ZX composition systematic on MC SR yields.

Build two MC ZX estimates of the SR from the same CR events:

  N_3P1F^X = F^X * N_3P1F^{MC}   (ZZ subtracted)
  N_2P2F^X = (F^X)^2 * N_2P2F^{MC}

The mix is 2P2F-only. Two variants are drawn:

  noprompt: Z+1L and 2P2F mixes use HF+LF+conv (full MC stack)
  leftover: Z+1L drops WZ and keeps leftover prompt; 2P2F mix includes
            leftover prompt. ZZ is still subtracted from 3P1F yields.

The hybrid per m4l bin is the same as the SR estimator (2P2F if
3P1F <= 2*2P2F, else 3P1F-2P2F), frozen on the Z1L estimate.

  delta = (N_ZX^Z1L - N_ZX^{2P2F}) / N_ZX^Z1L

The inclusive calculation uses one effective fake-factor weight per flavour.
The per-pT calculation applies the weights to individual failing leptons.
Relative MC variations are transferred bin by bin to uncertainty bands
around the unchanged data-driven Z+X prediction.

Set VBS below independently of analyzeZX.VBS:
  True  -> ZpXVBS2024_results (nominal, ge1j, ge2j)
  False -> ZplusXpython_HtoZZ2024 (nominal only)
"""

from __future__ import print_function

import csv
import math
import os
import sys

import ROOT


_EOS = "/eos/user/b/brunetd"
# Independent of analyzeZX.VBS / estimateZX.VBS so Condor production can stay
# on VBS while this script is run on the H->ZZ archive.
VBS = True
VBS_DIR = os.path.join(_EOS, "ZpXVBS2024_results")
HTOZZ_DIR = os.path.join(_EOS, "ZplusXpython_HtoZZ2024")
ANALYSIS_LABEL = "VBS" if VBS else "HtoZZ"
if VBS:
    STAGE1_DIR = os.path.join(VBS_DIR, "ROOT_FILES")
    STAGE2_DIR = STAGE1_DIR
    TRANSFER_BASE = os.path.join(VBS_DIR, "FR_uncertainty_plots")
    FINAL_DIR = os.path.join(STAGE2_DIR, "final_numbers")
    PUBLIC_PLOT_DIR = os.path.join(VBS_DIR, "ZX_composition")
    PERPT_MC_FILE = os.path.join(FINAL_DIR, "composition_perpt_mc.root")
    DATA_FR_FILES = {
        "nominal": os.path.join(STAGE1_DIR, "Hist_Data_ptl3_Data.root"),
        "ge1j": os.path.join(STAGE1_DIR, "Hist_Data_ptl3_VBSFR_ge1j_Data.root"),
        "ge2j": os.path.join(STAGE1_DIR, "Hist_Data_ptl3_VBSFR_ge2j_Data.root"),
    }
    WZ_FR_FILES = {
        "nominal": os.path.join(STAGE1_DIR, "Hist_MC_ptl3_WZ.root"),
        "ge1j": os.path.join(STAGE1_DIR, "Hist_MC_ptl3_VBSFR_ge1j_WZ.root"),
        "ge2j": os.path.join(STAGE1_DIR, "Hist_MC_ptl3_VBSFR_ge2j_WZ.root"),
    }
    ESTIMATE_FILES = {
        "nominal": ("estimateZXData", "estimateZXZZ"),
        "ge1j": ("estimateZXData_ge1j", "estimateZXZZ_ge1j"),
        "ge2j": ("estimateZXData_ge2j", "estimateZXZZ_ge2j"),
    }
    DEFAULT_TAGS = ("nominal", "ge1j", "ge2j")
else:
    STAGE1_DIR = os.path.join(HTOZZ_DIR, "stage1")
    STAGE2_DIR = os.path.join(HTOZZ_DIR, "stage2")
    TRANSFER_BASE = os.path.join(HTOZZ_DIR, "FR_uncertainty_plots")
    FINAL_DIR = os.path.join(STAGE2_DIR, "final_numbers")
    PUBLIC_PLOT_DIR = os.path.join(HTOZZ_DIR, "ZX_composition")
    PERPT_MC_FILE = os.path.join(FINAL_DIR, "composition_perpt_mc.root")
    DATA_FR_FILES = {
        "nominal": os.path.join(STAGE1_DIR, "Hist_Data_ptl3_Data.root"),
    }
    WZ_FR_FILES = {
        "nominal": os.path.join(STAGE1_DIR, "Hist_MC_ptl3_WZ.root"),
    }
    ESTIMATE_FILES = {
        "nominal": ("estimateZXData", "estimateZXZZ"),
    }
    DEFAULT_TAGS = ("nominal",)
CHANNELS = (
    ("4e", "electron"),
    ("4mu", "muon"),
    ("2e2mu", "muon"),
    ("2mu2e", "electron"),
)
MC_2P2F_SAMPLES = ("DY50", "TT", "WZ", "ZZ")
MC_3P1F_SAMPLES = ("DY50", "TT", "WZ")  # ZZ subtracted from 3P1F, as in data
FLAVOUR_REGIONS = {
    "electron": ("EB", "EE"),
    "muon": ("MB", "ME"),
}
DENOM_NAMES = {
    "EB": "Data_FRel_EB_d",
    "EE": "Data_FRel_EE_d",
    "MB": "Data_FRmu_EB_d",
    "ME": "Data_FRmu_EE_d",
}
LUMI_LABEL = "109.82 fb^{-1} (13.6 TeV)"


def _resolve(basename):
    for name in (basename, basename + ".root"):
        path = os.path.join(STAGE2_DIR, name)
        if os.path.isfile(path):
            return path
    raise RuntimeError("Missing estimate file: %s" % basename)


def _clone(tfile, name, out_name):
    hist = tfile.Get(name)
    if not hist:
        raise RuntimeError("Missing %s in %s" % (name, tfile.GetName()))
    clone = hist.Clone(out_name)
    clone.SetDirectory(0)
    return clone


def _integral(hist):
    return float(hist.Integral()) if hist else 0.0


def _safe_scale(value, factor):
    return max(0.0, value * factor)


def _load_transfer(fr_tag, bin_id="inclusive"):
    path = os.path.join(TRANSFER_BASE, fr_tag, "FR_CR_transfer.csv")
    want = str(bin_id)
    with open(path, "r") as handle:
        rows = [
            row for row in csv.DictReader(handle)
            if str(row["bin"]) == want
        ]
    if not rows:
        raise RuntimeError("No transfer rows for bin %s in %s" % (want, path))
    table = {}
    for row in rows:
        table[(row["cr"], row["region"])] = {
            "F_z1l": float(row["F_z1l"]),
            "F_cr": float(row["F_cr"]),
        }
    return table


def _load_inclusive_transfer(fr_tag):
    return _load_transfer(fr_tag, "inclusive")


def _denominator_weights(fr_tag):
    data_path = DATA_FR_FILES[fr_tag]
    wz_path = WZ_FR_FILES[fr_tag]
    data_file = ROOT.TFile.Open(data_path, "READ")
    wz_file = ROOT.TFile.Open(wz_path, "READ")
    if not data_file or data_file.IsZombie():
        raise RuntimeError("Cannot open %s" % data_path)
    if not wz_file or wz_file.IsZombie():
        data_file.Close()
        raise RuntimeError("Cannot open %s" % wz_path)
    weights = {}
    try:
        for region, name in DENOM_NAMES.items():
            data_d = data_file.Get(name)
            wz_d = wz_file.Get(name)
            if not data_d:
                raise RuntimeError("Missing %s in %s" % (name, data_path))
            value = _integral(data_d)
            if wz_d:
                value -= _integral(wz_d)
            weights[region] = max(0.0, value)
    finally:
        data_file.Close()
        wz_file.Close()
    return weights


def _flavour_average(transfer, weights, cr_kind, flavour, field):
    regions = FLAVOUR_REGIONS[flavour]
    num = 0.0
    den = 0.0
    parts = []
    for region in regions:
        value = transfer[(cr_kind, region)][field]
        weight = weights[region]
        num += weight * value
        den += weight
        parts.append((region, value, weight))
    if den <= 0.0:
        raise RuntimeError("Empty FR denominator weights for %s" % flavour)
    return num / den, parts


MIX_VARIANTS = (
    {
        "key": "noprompt",
        "cr": "2P2F_noprompt",
        "label": "no prompt",
        "envelope": "2P2F mix, no prompt",
        "method": "inclusive FR, no prompt",
    },
    {
        "key": "leftover",
        "cr": "2P2F_leftover",
        "label": "leftover prompt (WZ out of Z+1L, ZZ out of 3P1F)",
        "envelope": "2P2F mix, leftover prompt",
        "method": "inclusive FR, leftover prompt",
    },
)
OBSOLETE_PLOT_STEMS = (
    "ZX_composition_band",
    "ZX_composition_band_4e",
    "ZX_composition_band_4mu",
    "ZX_composition_band_2e2mu",
    "ZX_composition_band_2mu2e",
    "ZX_composition_band_perpt",
    "ZX_composition_band_4e_perpt",
    "ZX_composition_band_4mu_perpt",
    "ZX_composition_band_2e2mu_perpt",
    "ZX_composition_band_2mu2e_perpt",
    "ZX_composition_delta_vs_pt",
    "ZX_composition_band_correlated",
    "ZX_composition_band_quadrature",
)


def _flavour_F(transfer, weights, cr_kind):
    flavour_F = {}
    flavour_report = []
    for flavour in ("electron", "muon"):
        f_z1l, parts_z = _flavour_average(
            transfer, weights, cr_kind, flavour, "F_z1l"
        )
        f_cr, _parts_c = _flavour_average(
            transfer, weights, cr_kind, flavour, "F_cr"
        )
        flavour_F[flavour] = {
            "F_z1l": f_z1l,
            "F_cr": f_cr,
        }
        flavour_report.append(
            "%s: F_Z1L=%.5f  F_2P2F=%.5f" % (flavour, f_z1l, f_cr)
        )
        flavour_report.append(
            "  Z1L regions: %s"
            % ", ".join("%s F=%.5f w=%.1f" % item for item in parts_z)
        )
    return flavour_F, flavour_report


def _remove_obsolete_plots(plot_dir):
    for stem in OBSOLETE_PLOT_STEMS:
        for extension in (".pdf", ".png", ".C", ".root"):
            path = os.path.join(plot_dir, stem + extension)
            if os.path.isfile(path):
                os.remove(path)


def _mc_cr_path(nickname):
    return os.path.join(STAGE1_DIR, "Hist_MC_ptl3_%s.root" % nickname)


def _sum_mc_cr(samples, hist_name, out_name):
    total = None
    opened = []
    try:
        for nickname in samples:
            path = _mc_cr_path(nickname)
            tfile = ROOT.TFile.Open(path, "READ")
            if not tfile or tfile.IsZombie():
                raise RuntimeError("Cannot open MC CR file: %s" % path)
            opened.append(tfile)
            hist = tfile.Get(hist_name)
            if not hist:
                raise RuntimeError("Missing %s in %s" % (hist_name, path))
            if total is None:
                total = hist.Clone(out_name)
                total.SetDirectory(0)
            else:
                total.Add(hist)
    finally:
        for tfile in opened:
            tfile.Close()
    if total is None:
        raise RuntimeError("No MC samples for %s" % hist_name)
    return total


def _scale_hist(hist, factor, out_name):
    out = hist.Clone(out_name)
    out.SetDirectory(0)
    out.Scale(factor)
    for ibin in range(1, out.GetNbinsX() + 1):
        if out.GetBinContent(ibin) < 0.0:
            out.SetBinContent(ibin, 0.0)
    return out


def _combine_channel(h_2p2f, h_3p1f, scale_2p2f, scale_3p1f, freeze_choice, out_name):
    out = h_2p2f.Clone(out_name)
    out.Reset("ICES")
    out.SetDirectory(0)
    n_2p2f_bins = 0
    n_diff_bins = 0
    for ibin in range(1, h_2p2f.GetNbinsX() + 1):
        y2 = h_2p2f.GetBinContent(ibin)
        y3 = h_3p1f.GetBinContent(ibin)
        e2 = h_2p2f.GetBinError(ibin)
        e3 = h_3p1f.GetBinError(ibin)
        use_2p2f = freeze_choice[ibin]
        if use_2p2f:
            value = _safe_scale(y2, scale_2p2f)
            error = e2 * scale_2p2f
            n_2p2f_bins += 1
        else:
            value = _safe_scale(y3, scale_3p1f) - _safe_scale(y2, scale_2p2f)
            error = math.sqrt(
                (e3 * scale_3p1f) ** 2 + (e2 * scale_2p2f) ** 2
            )
            n_diff_bins += 1
        out.SetBinContent(ibin, value)
        out.SetBinError(ibin, error)
    return out, n_2p2f_bins, n_diff_bins


def _nominal_choice(h_2p2f, h_3p1f):
    choice = {}
    n_2p2f = 0
    n_diff = 0
    for ibin in range(1, h_2p2f.GetNbinsX() + 1):
        use_2p2f = h_3p1f.GetBinContent(ibin) <= 2.0 * h_2p2f.GetBinContent(ibin)
        choice[ibin] = use_2p2f
        if use_2p2f:
            n_2p2f += 1
        else:
            n_diff += 1
    return choice, n_2p2f, n_diff


def _add_channels(hists, out_name):
    result = hists[0].Clone(out_name)
    result.SetDirectory(0)
    for hist in hists[1:]:
        result.Add(hist)
    return result


def _variation_envelope(nominal, varied, up_name, down_name):
    up = nominal.Clone(up_name)
    down = nominal.Clone(down_name)
    up.SetDirectory(0)
    down.SetDirectory(0)
    for ibin in range(1, nominal.GetNbinsX() + 1):
        y0 = nominal.GetBinContent(ibin)
        y1 = varied.GetBinContent(ibin)
        up.SetBinContent(ibin, max(y0, y1))
        down.SetBinContent(ibin, min(y0, y1))
    return up, down


def _transfer_envelope(data, mc_nom, mc_var, up_name, down_name):
    """Put the MC composition envelope on data, bin by bin, as a relative band.

    Central value stays the data ZX yield. If MC Z1L is empty the bin is
    left unscaled.
    """
    up = data.Clone(up_name)
    down = data.Clone(down_name)
    up.SetDirectory(0)
    down.SetDirectory(0)
    for ibin in range(1, data.GetNbinsX() + 1):
        y = data.GetBinContent(ibin)
        y0 = mc_nom.GetBinContent(ibin)
        y1 = mc_var.GetBinContent(ibin)
        if y0 > 0.0:
            up.SetBinContent(ibin, y * max(y0, y1) / y0)
            down.SetBinContent(ibin, y * min(y0, y1) / y0)
        else:
            up.SetBinContent(ibin, y)
            down.SetBinContent(ibin, y)
    return up, down


def _cms_prelim(canvas):
    lumi = ROOT.TLatex(0.95, 0.96, LUMI_LABEL)
    lumi.SetNDC()
    lumi.SetTextAlign(31)
    lumi.SetTextFont(42)
    lumi.SetTextSize(0.03)
    lumi.Draw()
    cms = ROOT.TLatex(0.15, 0.96, "CMS")
    cms.SetNDC()
    cms.SetTextFont(61)
    cms.SetTextSize(0.0375)
    cms.Draw()
    prelim = ROOT.TLatex(0.23, 0.96, "Preliminary")
    prelim.SetNDC()
    prelim.SetTextFont(52)
    prelim.SetTextSize(0.0285)
    prelim.Draw()
    canvas._cms_text = (lumi, cms, prelim)


def _fr_method_header(fr_tag, method, extra=""):
    """Legend header that states inclusive FR vs per-pT FR explicitly."""
    title = ANALYSIS_LABEL + " Z+X"
    if extra:
        title += " " + extra
    if fr_tag and fr_tag != "nominal":
        title += " [%s]" % fr_tag
    return "#splitline{%s}{%s}" % (title, method)


def _band_graph(nominal, up, down):
    graph = ROOT.TGraphAsymmErrors(nominal.GetNbinsX())
    ymax = 0.0
    for ibin in range(1, nominal.GetNbinsX() + 1):
        x = nominal.GetXaxis().GetBinCenter(ibin)
        ex = 0.5 * nominal.GetXaxis().GetBinWidth(ibin)
        y = nominal.GetBinContent(ibin)
        y_up = up.GetBinContent(ibin)
        y_down = down.GetBinContent(ibin)
        graph.SetPoint(ibin - 1, x, y)
        graph.SetPointError(
            ibin - 1,
            ex,
            ex,
            max(0.0, y - y_down),
            max(0.0, y_up - y),
        )
        ymax = max(ymax, y + nominal.GetBinError(ibin), y_up, y_down)
    return graph, ymax


def _draw_tot_band(
    nominal,
    up,
    down,
    out_base,
    header=None,
    envelope_label="2P2F mix envelope",
    method_label="inclusive FR",
    central_label="data ZX",
):
    canvas = ROOT.TCanvas("c_" + os.path.basename(out_base), out_base, 0, 0, 800, 600)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    canvas.SetBottomMargin(0.13)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    graph, ymax = _band_graph(nominal, up, down)

    frame = nominal.Clone(os.path.basename(out_base) + "_frame")
    frame.Reset("ICES")
    frame.SetDirectory(0)
    frame.GetXaxis().SetTitle("m_{4l} [GeV]")
    frame.GetYaxis().SetTitle("Z+X events / bin")
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.SetMinimum(0.0)
    frame.SetMaximum(max(1.45 * ymax, 0.05))
    frame.Draw("axis")

    graph.SetFillColorAlpha(ROOT.kAzure + 1, 0.45)
    graph.SetFillStyle(1001)
    graph.SetLineColor(ROOT.kAzure + 2)
    graph.SetLineWidth(2)
    graph.Draw("2 same")
    nominal.SetLineColor(ROOT.kBlack)
    nominal.SetMarkerColor(ROOT.kBlack)
    nominal.SetMarkerStyle(20)
    nominal.SetLineWidth(1)
    nominal.Draw("E1 X0 same")

    legend = ROOT.TLegend(0.38, 0.64, 0.92, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.028)
    legend.SetHeader(
        header or _fr_method_header("nominal", method_label)
    )
    legend.AddEntry(nominal, central_label, "LP")
    legend.AddEntry(graph, envelope_label, "F")
    legend.Draw()
    _cms_prelim(canvas)
    canvas._keep = (frame, graph, nominal, legend)
    canvas.SaveAs(out_base + ".pdf")
    canvas.SaveAs(out_base + ".png")


def _format_yield(value, error):
    return "%.4f +/- %.4f" % (value, error)


def _hist_error(hist):
    return math.sqrt(
        sum(
            hist.GetBinError(ibin) ** 2
            for ibin in range(1, hist.GetNbinsX() + 1)
        )
    )


def _write_perpt_bands(fr_tag, out_dir, public_dir, channel_data=None, tot_data=None):
    """Same m4l envelope plots as inclusive F, with event-by-event per-pT F."""
    if not os.path.isfile(PERPT_MC_FILE):
        print("Skipping per-pT composition: missing %s" % PERPT_MC_FILE)
        return
    tfile = ROOT.TFile.Open(PERPT_MC_FILE, "READ")
    if not tfile or tfile.IsZombie():
        print("Skipping per-pT composition: cannot open %s" % PERPT_MC_FILE)
        return
    probe = "mc_%s_z1l_noprompt_2P2F_4e" % fr_tag
    if not tfile.Get(probe):
        tfile.Close()
        print(
            "Skipping per-pT composition: %s not in %s. "
            "Refill composition_perpt_loop.py after the 2P2F mix change."
            % (probe, PERPT_MC_FILE)
        )
        return

    def load(kind, cr, channel):
        name = "mc_%s_%s_%s_%s" % (fr_tag, kind, cr, channel)
        hist = tfile.Get(name)
        if not hist:
            tfile.Close()
            raise RuntimeError("Missing %s in %s" % (name, PERPT_MC_FILE))
        clone = hist.Clone(name + "_use")
        clone.SetDirectory(0)
        return clone

    lines = [
        "%s Z+X composition per-pT FR [%s]" % (ANALYSIS_LABEL, fr_tag),
        "F treatment: per-pT FR (f(pT,eta) of each failing lepton).",
        "Mix: 2P2F-only. 3P1F is not used in the FR mix.",
        "3P1F: w*F_failing. 2P2F: w*F3*F4. Same hybrid as inclusive F.",
        "ZZ subtracted from 3P1F yields.",
        "",
    ]
    csv_rows = ["channel,variant,n_mc_z1l,n_mc_2p2f,delta"]
    plot_specs = []
    try:
        for variant in MIX_VARIANTS:
            z1l_kind = "z1l_" + variant["key"]
            cr_kind = "cr_" + variant["key"]
            channel_z1l = []
            channel_cr = []
            mc_n_z1l = []
            mc_n_cr = []
            lines.append("=== Mix: %s ===" % variant["label"])
            for channel, flavour in CHANNELS:
                h2_z1l = load(z1l_kind, "2P2F", channel)
                h3_z1l = load(z1l_kind, "3P1F", channel)
                h2_cr = load(cr_kind, "2P2F", channel)
                h3_cr = load(cr_kind, "3P1F", channel)
                choice, n_2p2f, n_diff = _nominal_choice(h2_z1l, h3_z1l)
                tot_z1l, _, _ = _combine_channel(
                    h2_z1l, h3_z1l, 1.0, 1.0, choice,
                    "perpt_%s_z1l_%s" % (variant["key"], channel),
                )
                tot_cr, _, _ = _combine_channel(
                    h2_cr, h3_cr, 1.0, 1.0, choice,
                    "perpt_%s_cr_%s" % (variant["key"], channel),
                )
                n_z1l = _integral(tot_z1l)
                n_cr = _integral(tot_cr)
                delta = (n_z1l - n_cr) / n_z1l if n_z1l > 0.0 else 0.0
                channel_z1l.append(tot_z1l)
                channel_cr.append(tot_cr)
                mc_n_z1l.append(n_z1l)
                mc_n_cr.append(n_cr)
                csv_rows.append(
                    "%s,%s,%.6f,%.6f,%.6f"
                    % (channel, variant["key"], n_z1l, n_cr, delta)
                )
                lines.append("=== %s (%s failing Z2) ===" % (channel, flavour))
                lines.append(
                    "MC hybrid bins (frozen on Z1L): %d use 2P2F, %d use 3P1F-2P2F"
                    % (n_2p2f, n_diff)
                )
                lines.append(
                    "N_Z1L=%.4f  N_2P2F=%.4f  delta=%.4f (%.1f%%)"
                    % (n_z1l, n_cr, delta, 100.0 * abs(delta))
                )
                lines.append("")
                ch_up, ch_down = _variation_envelope(
                    tot_z1l,
                    tot_cr,
                    "perpt_%s_up_%s" % (variant["key"], channel),
                    "perpt_%s_down_%s" % (variant["key"], channel),
                )
                plot_specs.append(
                    (
                        tot_z1l,
                        ch_up,
                        ch_down,
                        "ZX_composition_band_mc_%s_%s_perpt"
                        % (channel, variant["key"]),
                        _fr_method_header(
                            fr_tag,
                            "per-p_{T} FR, " + variant["key"],
                            extra=channel,
                        ),
                        variant["envelope"],
                        "MC ZX (F_{Z+1L})",
                    )
                )
                if channel_data is not None:
                    data_ch = channel_data[len(channel_z1l) - 1]
                    data_up, data_down = _transfer_envelope(
                        data_ch,
                        tot_z1l,
                        tot_cr,
                        "perpt_data_%s_up_%s" % (variant["key"], channel),
                        "perpt_data_%s_down_%s" % (variant["key"], channel),
                    )
                    plot_specs.append(
                        (
                            data_ch,
                            data_up,
                            data_down,
                            "ZX_composition_band_%s_%s_perpt"
                            % (channel, variant["key"]),
                            _fr_method_header(
                                fr_tag,
                                "per-p_{T} FR, " + variant["key"],
                                extra=channel,
                            ),
                            variant["envelope"],
                            "data ZX",
                        )
                    )
            tot_z1l = _add_channels(
                channel_z1l, "perpt_zx_z1l_" + variant["key"]
            )
            tot_cr = _add_channels(channel_cr, "perpt_zx_cr_" + variant["key"])
            tot_up, tot_down = _variation_envelope(
                tot_z1l,
                tot_cr,
                "perpt_zx_up_" + variant["key"],
                "perpt_zx_down_" + variant["key"],
            )
            n_z1l_inc = sum(mc_n_z1l)
            n_cr_inc = sum(mc_n_cr)
            delta_inc = (
                (n_z1l_inc - n_cr_inc) / n_z1l_inc if n_z1l_inc > 0.0 else 0.0
            )
            csv_rows.append(
                "inclusive,%s,%.6f,%.6f,%.6f"
                % (variant["key"], n_z1l_inc, n_cr_inc, delta_inc)
            )
            lines.append(
                "Inclusive %s: N_Z1L=%.4f  N_2P2F=%.4f  delta=%.4f (%.1f%%)"
                % (
                    variant["key"],
                    n_z1l_inc,
                    n_cr_inc,
                    delta_inc,
                    100.0 * abs(delta_inc),
                )
            )
            lines.append("")
            plot_specs.append(
                (
                    tot_z1l,
                    tot_up,
                    tot_down,
                    "ZX_composition_band_mc_%s_perpt" % variant["key"],
                    _fr_method_header(
                        fr_tag, "per-p_{T} FR, " + variant["key"]
                    ),
                    variant["envelope"],
                    "MC ZX (F_{Z+1L})",
                )
            )
            if tot_data is not None:
                data_up, data_down = _transfer_envelope(
                    tot_data,
                    tot_z1l,
                    tot_cr,
                    "perpt_data_%s_up" % variant["key"],
                    "perpt_data_%s_down" % variant["key"],
                )
                plot_specs.append(
                    (
                        tot_data,
                        data_up,
                        data_down,
                        "ZX_composition_band_%s_perpt" % variant["key"],
                        _fr_method_header(
                            fr_tag, "per-p_{T} FR, " + variant["key"]
                        ),
                        variant["envelope"],
                        "data ZX",
                    )
                )
    finally:
        tfile.Close()

    text = "\n".join(lines) + "\n"
    print(text)
    txt_path = os.path.join(out_dir, "composition_zx_yield_per_pt.txt")
    with open(txt_path, "w") as handle:
        handle.write(text)
    with open(
        os.path.join(public_dir, "composition_zx_yield_per_pt.txt"), "w"
    ) as handle:
        handle.write(text)
    csv_path = os.path.join(out_dir, "composition_zx_yield_per_pt.csv")
    with open(csv_path, "w") as handle:
        handle.write("\n".join(csv_rows) + "\n")
    for plot_dir in (out_dir, public_dir):
        _remove_obsolete_plots(plot_dir)
        for nominal, up, down, stem, header, envelope, central in plot_specs:
            _draw_tot_band(
                nominal,
                up,
                down,
                os.path.join(plot_dir, stem),
                header=header,
                envelope_label=envelope,
                method_label="per-p_{T} FR",
                central_label=central,
            )
    print("Wrote per-pT composition bands for %s" % fr_tag)


def run_tag(fr_tag):
    transfer = _load_inclusive_transfer(fr_tag)
    weights = _denominator_weights(fr_tag)
    data_name, zz_name = ESTIMATE_FILES[fr_tag]
    data_file = ROOT.TFile.Open(_resolve(data_name), "READ")
    zz_file = ROOT.TFile.Open(_resolve(zz_name), "READ")
    if not data_file or data_file.IsZombie():
        raise RuntimeError("Cannot open data estimate for %s" % fr_tag)
    if not zz_file or zz_file.IsZombie():
        data_file.Close()
        raise RuntimeError("Cannot open ZZ estimate for %s" % fr_tag)

    out_dir = os.path.join(FINAL_DIR, "composition_" + fr_tag)
    public_dir = os.path.join(PUBLIC_PLOT_DIR, fr_tag)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(public_dir, exist_ok=True)
    _remove_obsolete_plots(out_dir)
    _remove_obsolete_plots(public_dir)

    lines = [
        "%s Z+X composition systematic [%s]" % (ANALYSIS_LABEL, fr_tag),
        "F treatment: inclusive FR (one F per flavour = F_2P2F vs F_Z1L).",
        "Mix: 2P2F-only. 3P1F is not used in the FR mix.",
        "noprompt: HF+LF+conv from the full MC stack.",
        "leftover: WZ dropped from Z+1L; leftover prompt kept in the mix;",
        "          ZZ subtracted from 3P1F yields.",
        "MC SR ZX: CR event counts x F (3P1F) or F^2 (2P2F).",
        "3P1F MC samples: DY+TT+WZ (ZZ subtracted). 2P2F: DY+TT+WZ+ZZ.",
        "Hybrid: 2P2F if 3P1F <= 2*2P2F else 3P1F-2P2F, frozen on MC Z1L.",
        "delta = (N_MC_SR(F_Z1L) - N_MC_SR(F_2P2F)) / N_MC_SR(F_Z1L).",
        "Inclusive F: flavour-averaged with the WZ-subtracted Z+1L denominator.",
        "Inclusive and per-pT composition results are written separately.",
        "Relative MC variations define bands around the unchanged data ZX prediction.",
        "",
    ]

    stored = {}
    channel_data = []
    h2_mc_by_ch = {}
    h3_mc_by_ch = {}
    csv_rows = ["channel,variant,n_mc_z1l,n_mc_2p2f,delta"]
    try:
        for channel, flavour in CHANNELS:
            h2_mc_by_ch[channel] = _sum_mc_cr(
                MC_2P2F_SAMPLES,
                "h1D_m4l_2P2F_" + channel,
                "mc_unw_2P2F_" + channel,
            )
            h3_mc_by_ch[channel] = _sum_mc_cr(
                MC_3P1F_SAMPLES,
                "h1D_m4l_3P1F_" + channel,
                "mc_unw_3P1F_" + channel,
            )
            h2 = _clone(
                data_file,
                "h1D_m4l_SR_2P2F_" + channel,
                "nom_2P2F_" + channel,
            )
            h3 = _clone(
                data_file,
                "h1D_m4l_SR_3P1F_" + channel,
                "nom_3P1F_" + channel,
            )
            h3.Add(
                _clone(
                    zz_file,
                    "h1D_m4l_SR_3P1F_" + channel,
                    "zz_3P1F_" + channel,
                ),
                -1,
            )
            data_choice, n_data_2, n_data_diff = _nominal_choice(h2, h3)
            nom, _, _ = _combine_channel(
                h2, h3, 1.0, 1.0, data_choice, "tot_" + channel
            )
            channel_data.append(nom)
            stored.update(
                {
                    "h1D_m4l_SR_2P2F_" + channel: h2,
                    "h1D_m4l_SR_3P1F_ZZsub_" + channel: h3,
                    "h1D_m4l_SR_tot_" + channel: nom,
                }
            )
            lines.append(
                "data %s hybrid bins: %d use 2P2F, %d use 3P1F-2P2F; "
                "SR_tot=%.4f"
                % (channel, n_data_2, n_data_diff, _integral(nom))
            )
        tot_data = _add_channels(channel_data, "h1D_m4l_SR_tot")
        stored["h1D_m4l_SR_tot"] = tot_data
        lines.append(
            "nominal data SR_tot = %s (no composition scale)"
            % _format_yield(_integral(tot_data), _hist_error(tot_data))
        )
        lines.append("")

        for variant in MIX_VARIANTS:
            flavour_F, flavour_report = _flavour_F(
                transfer, weights, variant["cr"]
            )
            lines.append("=== Mix: %s ===" % variant["label"])
            lines.extend(flavour_report)
            lines.append("")
            channel_z1l = []
            channel_cr = []
            mc_n_z1l = []
            mc_n_cr = []
            plot_specs = []
            for channel, flavour in CHANNELS:
                f_z1l = flavour_F[flavour]["F_z1l"]
                f_cr = flavour_F[flavour]["F_cr"]
                h2_mc = h2_mc_by_ch[channel]
                h3_mc = h3_mc_by_ch[channel]
                h2_z1l = _scale_hist(
                    h2_mc,
                    f_z1l * f_z1l,
                    "mc_%s_z1l_2P2F_%s" % (variant["key"], channel),
                )
                h3_z1l = _scale_hist(
                    h3_mc,
                    f_z1l,
                    "mc_%s_z1l_3P1F_%s" % (variant["key"], channel),
                )
                h2_cr = _scale_hist(
                    h2_mc,
                    f_cr * f_cr,
                    "mc_%s_cr_2P2F_%s" % (variant["key"], channel),
                )
                h3_cr = _scale_hist(
                    h3_mc,
                    f_cr,
                    "mc_%s_cr_3P1F_%s" % (variant["key"], channel),
                )
                choice, n_2p2f, n_diff = _nominal_choice(h2_z1l, h3_z1l)
                tot_z1l, _, _ = _combine_channel(
                    h2_z1l,
                    h3_z1l,
                    1.0,
                    1.0,
                    choice,
                    "mc_zx_%s_z1l_%s" % (variant["key"], channel),
                )
                tot_cr, _, _ = _combine_channel(
                    h2_cr,
                    h3_cr,
                    1.0,
                    1.0,
                    choice,
                    "mc_zx_%s_cr_%s" % (variant["key"], channel),
                )
                n_z1l = _integral(tot_z1l)
                n_cr = _integral(tot_cr)
                if n_z1l <= 0.0:
                    raise RuntimeError(
                        "Empty MC Z1L ZX yield for %s (%s)"
                        % (channel, variant["key"])
                    )
                delta = (n_z1l - n_cr) / n_z1l
                channel_z1l.append(tot_z1l)
                channel_cr.append(tot_cr)
                mc_n_z1l.append(n_z1l)
                mc_n_cr.append(n_cr)
                csv_rows.append(
                    "%s,%s,%.6f,%.6f,%.6f"
                    % (channel, variant["key"], n_z1l, n_cr, delta)
                )
                stored["mc_zx_%s_z1l_%s" % (variant["key"], channel)] = tot_z1l
                stored["mc_zx_%s_cr_%s" % (variant["key"], channel)] = tot_cr
                ch_up, ch_down = _variation_envelope(
                    tot_z1l,
                    tot_cr,
                    "mc_zx_%s_up_%s" % (variant["key"], channel),
                    "mc_zx_%s_down_%s" % (variant["key"], channel),
                )
                stored["mc_zx_%s_up_%s" % (variant["key"], channel)] = ch_up
                stored["mc_zx_%s_down_%s" % (variant["key"], channel)] = ch_down
                data_ch = channel_data[len(channel_z1l) - 1]
                data_up, data_down = _transfer_envelope(
                    data_ch,
                    tot_z1l,
                    tot_cr,
                    "data_zx_%s_up_%s" % (variant["key"], channel),
                    "data_zx_%s_down_%s" % (variant["key"], channel),
                )
                stored["data_zx_%s_up_%s" % (variant["key"], channel)] = data_up
                stored["data_zx_%s_down_%s" % (variant["key"], channel)] = data_down
                plot_specs.append(
                    (
                        tot_z1l,
                        ch_up,
                        ch_down,
                        "ZX_composition_band_mc_%s_%s"
                        % (channel, variant["key"]),
                        _fr_method_header(
                            fr_tag, variant["method"], extra=channel
                        ),
                        variant["envelope"],
                        "MC ZX (F_{Z+1L})",
                    )
                )
                plot_specs.append(
                    (
                        data_ch,
                        data_up,
                        data_down,
                        "ZX_composition_band_%s_%s"
                        % (channel, variant["key"]),
                        _fr_method_header(
                            fr_tag, variant["method"], extra=channel
                        ),
                        variant["envelope"],
                        "data ZX",
                    )
                )
                lines.append("=== %s (%s failing Z2) ===" % (channel, flavour))
                lines.append(
                    "MC hybrid bins (frozen on Z1L): %d use 2P2F, %d use 3P1F-2P2F"
                    % (n_2p2f, n_diff)
                )
                lines.append(
                    "F_Z1L=%.5f  F_2P2F=%.5f" % (f_z1l, f_cr)
                )
                lines.append(
                    "N_Z1L=%.4f  N_2P2F=%.4f  delta=%.4f (%.1f%%)"
                    % (n_z1l, n_cr, delta, 100.0 * abs(delta))
                )
                lines.append("")

            tot_z1l = _add_channels(
                channel_z1l, "mc_zx_%s_z1l" % variant["key"]
            )
            tot_cr = _add_channels(channel_cr, "mc_zx_%s_cr" % variant["key"])
            tot_up, tot_down = _variation_envelope(
                tot_z1l,
                tot_cr,
                "mc_zx_%s_up" % variant["key"],
                "mc_zx_%s_down" % variant["key"],
            )
            stored["mc_zx_%s_z1l" % variant["key"]] = tot_z1l
            stored["mc_zx_%s_cr" % variant["key"]] = tot_cr
            stored["mc_zx_%s_up" % variant["key"]] = tot_up
            stored["mc_zx_%s_down" % variant["key"]] = tot_down
            n_z1l_inc = sum(mc_n_z1l)
            n_cr_inc = sum(mc_n_cr)
            delta_inc = (
                (n_z1l_inc - n_cr_inc) / n_z1l_inc if n_z1l_inc > 0.0 else 0.0
            )
            csv_rows.append(
                "inclusive,%s,%.6f,%.6f,%.6f"
                % (variant["key"], n_z1l_inc, n_cr_inc, delta_inc)
            )
            lines.append(
                "Inclusive %s: N_Z1L=%.4f  N_2P2F=%.4f  delta=%.4f (%.1f%%)"
                % (
                    variant["key"],
                    n_z1l_inc,
                    n_cr_inc,
                    delta_inc,
                    100.0 * abs(delta_inc),
                )
            )
            lines.append("")
            plot_specs.append(
                (
                    tot_z1l,
                    tot_up,
                    tot_down,
                    "ZX_composition_band_mc_" + variant["key"],
                    _fr_method_header(fr_tag, variant["method"]),
                    variant["envelope"],
                    "MC ZX (F_{Z+1L})",
                )
            )
            data_up, data_down = _transfer_envelope(
                tot_data,
                tot_z1l,
                tot_cr,
                "data_zx_%s_up" % variant["key"],
                "data_zx_%s_down" % variant["key"],
            )
            stored["data_zx_%s_up" % variant["key"]] = data_up
            stored["data_zx_%s_down" % variant["key"]] = data_down
            plot_specs.append(
                (
                    tot_data,
                    data_up,
                    data_down,
                    "ZX_composition_band_" + variant["key"],
                    _fr_method_header(fr_tag, variant["method"]),
                    variant["envelope"],
                    "data ZX",
                )
            )
            for plot_dir in (out_dir, public_dir):
                for nominal, up, down, stem, header, envelope, central in plot_specs:
                    _draw_tot_band(
                        nominal,
                        up,
                        down,
                        os.path.join(plot_dir, stem),
                        header=header,
                        envelope_label=envelope,
                        method_label=variant["method"],
                        central_label=central,
                    )
    finally:
        data_file.Close()
        zz_file.Close()

    text = "\n".join(lines) + "\n"
    text_path = os.path.join(out_dir, "composition_zx_yield.txt")
    with open(text_path, "w") as handle:
        handle.write(text)
    with open(os.path.join(public_dir, "composition_zx_yield.txt"), "w") as handle:
        handle.write(text)
    print(text)
    print("Wrote %s" % text_path)

    csv_path = os.path.join(out_dir, "composition_zx_yield.csv")
    with open(csv_path, "w") as handle:
        handle.write("\n".join(csv_rows) + "\n")

    root_path = os.path.join(out_dir, "composition_zx_yield.root")
    root_file = ROOT.TFile(root_path, "RECREATE")
    for name, hist in stored.items():
        written = hist.Clone(name)
        written.SetDirectory(root_file)
        written.Write()
    root_file.Close()
    print("Wrote %s" % root_path)
    print("Wrote %s" % csv_path)
    _write_perpt_bands(
        fr_tag, out_dir, public_dir, channel_data=channel_data, tot_data=tot_data
    )


def main():
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    os.makedirs(FINAL_DIR, exist_ok=True)
    tags = sys.argv[1:] or list(DEFAULT_TAGS)
    for tag in tags:
        print("\n======== %s ========" % tag)
        run_tag(tag)
    return 0


if __name__ == "__main__":
    sys.exit(main())
