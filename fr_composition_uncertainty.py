"""
Fake-rate composition uncertainties.

  - Z+1L truth-origin composition (denominator) with Clopper-Pearson errors
  - 2P2F/3P1F CR composition fractions with Clopper-Pearson errors
  - composition-weighted average fake rates per EB/EE/MB/ME

Inputs: Hist_ID_ptl3_{ZZ,WZ,TT,DY50}.root from stage 1 (NanoAOD analyzeZX).
  Set VBS below to choose the input directory and output folder:
    VBS=True  -> VBS Hist_ID (nominal, ge1j and ge2j)
    VBS=False -> H->ZZ Hist_ID on EOS (latest HtoZZ stage1)
  Both modes build a 2P2F-only CR-mix fake rate (3P1F is not used in
  the mix). Two variants are written: noprompt (HF+LF+conv) and leftover
  prompt (WZ dropped from Z+1L; leftover prompt kept in the mix).
  ZZ is subtracted on the 3P1F yield, not in this mix. This flag is
  independent of analyzeZX.VBS.

Z+1L MC measures origin fake rates f_i (HF, LF, conversion, and prompt).
noprompt uses HF/LF/conversion only. leftover keeps leftover prompt after
dropping the WZ sample from Z+1L, matching the leftover bias in data FR.

Residual-prompt SF (Z+1L denominator):
  SF = prompt(MC without WZ) / prompt(all MC samples)
Diagnostic only; it is not applied as a correction or a separate systematic.

"""

from __future__ import print_function

import csv
import math
import os
import sys

import ROOT
from ROOT import TColor

# Analysis mode for this script only (independent of analyzeZX.VBS).
# False = H->ZZ (inclusive CR, m4l from ~70 GeV); True = VBS-gated CR.
# Production jobs keep analyzeZX.VBS / estimateZX.VBS as they are.
VBS = True

_EOS = "/eos/user/b/brunetd"
HTOZZ_HIST_ID_DIR = os.path.join(_EOS, "ZplusXpython_HtoZZ2024", "stage1")
HTOZZ_DATA_FR_FILE = os.path.join(
    HTOZZ_HIST_ID_DIR, "Hist_Data_ptl3_Data.root"
)
HTOZZ_WZ_REMOVED_FR_FILE = os.path.join(
    HTOZZ_HIST_ID_DIR, "Hist_Data_ptl3_WZremoved.root"
)
VBS_HIST_ID_DIR = os.path.join(_EOS, "ZpXVBS2024_results", "ROOT_FILES")
HTOZZ_OUTPUT_DIR = os.path.join(
    _EOS, "ZplusXpython_HtoZZ2024", "FR_uncertainty_plots"
)
VBS_OUTPUT_DIR = os.path.join(
    _EOS, "ZpXVBS2024_results", "FR_uncertainty_plots"
)

VAR_NAME = "ptl3"
MC_SAMPLES = ("ZZ", "WZ", "TT", "DY50")


def _resolve_input_dir():
    if not VBS:
        return HTOZZ_HIST_ID_DIR
    return VBS_HIST_ID_DIR


INPUT_DIR = _resolve_input_dir()
BASE_OUT_DIR = VBS_OUTPUT_DIR if VBS else HTOZZ_OUTPUT_DIR
OUT_DIR = BASE_OUT_DIR
FR_TAGS = (None, "ge1j", "ge2j") if VBS else (None,)
ANALYSIS_MODE = "VBS" if VBS else "HtoZZ"

CR_SECTIONS = [
    "2P2F",
    "2P2F4e",
    "2P2F4mu",
    "2P2F2e2mu",
    "2P2F2mu2e",
    "3P1F",
    "3P1F4e",
    "3P1F4mu",
    "3P1F2e2mu",
    "3P1F2mu2e",
]

FR_REGIONS = ("EB", "EE", "MB", "ME")
FR_REGION_LABELS = {
    "EB": "electron barrel",
    "EE": "electron endcap",
    "MB": "muon barrel",
    "ME": "muon endcap",
}
DATA_FR_NAMES = {
    "EB": "Data_FRel_EB",
    "EE": "Data_FRel_EE",
    "MB": "Data_FRmu_EB",
    "ME": "Data_FRmu_EE",
}
CR_FINAL_STATE_PAIRS = {
    "total": ("2P2F", "3P1F"),
    "4e": ("2P2F4e", "3P1F4e"),
    "4mu": ("2P2F4mu", "3P1F4mu"),
    "2e2mu": ("2P2F2e2mu", "3P1F2e2mu"),
    "2mu2e": ("2P2F2mu2e", "3P1F2mu2e"),
}
CR_COMBINED_SECTIONS = {
    # Z2 is the fake-rate pair: 4e and 2mu2e have electron Z2 candidates;
    # 4mu and 2e2mu have muon Z2 candidates.
    "2P2Fe": ("2P2F4e", "2P2F2mu2e"),
    "3P1Fe": ("3P1F4e", "3P1F2mu2e"),
    "2P2Fmu": ("2P2F4mu", "2P2F2e2mu"),
    "3P1Fmu": ("3P1F4mu", "3P1F2e2mu"),
    # One CR mix per flavour: 2P2F and 3P1F added before the fake rate.
    "electron": ("2P2F4e", "2P2F2mu2e", "3P1F4e", "3P1F2mu2e"),
    "muon": ("2P2F4mu", "2P2F2e2mu", "3P1F4mu", "3P1F2e2mu"),
}
CR_FLAVOUR_PAIRS = {
    "electron": ("2P2Fe", "3P1Fe"),
    "muon": ("2P2Fmu", "3P1Fmu"),
}
CR_SUMMARY_PAIRS = dict(CR_FINAL_STATE_PAIRS)
CR_SUMMARY_PAIRS.update(CR_FLAVOUR_PAIRS)

# Fill colours for ID stacks
COLOR_PROMPT = ROOT.kBlue - 7
COLOR_FAKE = ROOT.kRed - 7
COLOR_CONV = ROOT.kGreen - 7
COLOR_BDFAKE = ROOT.kRed - 10

LUMI_LABEL = "109.82 fb^{-1} (13.6 TeV)"


def _safe_integral(hist):
    if not hist:
        return 0.0
    return float(hist.Integral())


def _integral_variance(hist):
    """Variance of an integral from the histogram's Sumw2 errors."""
    if not hist:
        return 0.0
    return sum(
        float(hist.GetBinError(ibin)) ** 2
        for ibin in range(1, hist.GetNbinsX() + 1)
    )


def _weighted_efficiency(numerator, denominator):
    """Weighted pass/total ratio and its Sumw2-based standard error."""
    den = _safe_integral(denominator)
    num = _safe_integral(numerator)
    if den <= 0:
        return 0.0, 0.0
    value = num / den
    num_variance = _integral_variance(numerator)
    den_variance = _integral_variance(denominator)
    # The numerator is a subset of the denominator, so Cov(num, den)=Var(num).
    variance = (
        (1.0 - 2.0 * value) * num_variance
        + value * value * den_variance
    ) / (den * den)
    return value, math.sqrt(max(0.0, variance))


def _weighted_composition(histograms):
    """Fractions and covariance matrix for independent weighted yields."""
    keys = tuple(histograms)
    yields = {key: _safe_integral(histograms[key]) for key in keys}
    variances = {
        key: _integral_variance(histograms[key]) for key in keys
    }
    total = sum(yields.values())
    if total <= 0:
        return None, None
    fractions = {key: yields[key] / total for key in keys}
    covariance = {}
    for key_i in keys:
        covariance[key_i] = {}
        for key_j in keys:
            value = 0.0
            for source in keys:
                deriv_i = (
                    (1.0 if source == key_i else 0.0) - fractions[key_i]
                ) / total
                deriv_j = (
                    (1.0 if source == key_j else 0.0) - fractions[key_j]
                ) / total
                value += deriv_i * deriv_j * variances[source]
            covariance[key_i][key_j] = value
    return fractions, covariance


def _clone_sum(files, hist_name, out_name):
    first = None
    for tfile in files:
        hist = tfile.Get(hist_name)
        if not hist:
            continue
        if first is None:
            first = hist.Clone(out_name)
            first.SetDirectory(0)
        else:
            first.Add(hist)
    return first


def _open_mc_files(
    var_name=VAR_NAME, samples=MC_SAMPLES, input_dir=None, fr_tag=None
):
    if input_dir is None:
        input_dir = INPUT_DIR
    files = []
    for sample in samples:
        if fr_tag is None:
            filename = "Hist_ID_%s_%s.root" % (var_name, sample)
        else:
            filename = "Hist_ID_%s_VBSFR_%s_%s.root" % (
                var_name, fr_tag, sample
            )
        path = os.path.join(input_dir, filename)
        if not os.path.isfile(path):
            print("WARNING: missing %s — skipping" % path)
            continue
        tfile = ROOT.TFile.Open(path, "READ")
        if not tfile or tfile.IsZombie():
            print("WARNING: cannot open %s — skipping" % path)
            continue
        files.append(tfile)
        print("Opened %s" % path)
    if not files:
        raise RuntimeError(
            "No Hist_ID files found in %s for samples %s (FR tag: %s)"
            % (input_dir, samples, fr_tag or "nominal")
        )
    return files


def _open_required_file(path):
    if not os.path.isfile(path):
        raise RuntimeError("Missing required ROOT file: %s" % path)
    tfile = ROOT.TFile.Open(path, "READ")
    if not tfile or tfile.IsZombie():
        raise RuntimeError("Cannot open ROOT file: %s" % path)
    return tfile


def _clone_required(tfile, name, out_name):
    hist = tfile.Get(name)
    if not hist:
        raise RuntimeError("Missing %s in %s" % (name, tfile.GetName()))
    clone = hist.Clone(out_name)
    clone.SetDirectory(0)
    return clone


def _add_histograms(histograms, out_name):
    valid = [hist for hist in histograms if hist]
    if not valid:
        return None
    result = valid[0].Clone(out_name)
    result.SetDirectory(0)
    for hist in valid[1:]:
        result.Add(hist)
    return result


def _ratio_hist(numerator, denominator, out_name):
    ratio = numerator.Clone(out_name)
    ratio.SetDirectory(0)
    ratio.Divide(numerator, denominator, 1.0, 1.0, "B")
    return ratio


def _files_without_wz(files):
    return _files_excluding(files, ("WZ",))


def _files_excluding(files, nicknames):
    """Drop Hist_ID files whose basename ends with _<nickname>.root."""
    nicknames = tuple(nicknames)
    kept = []
    for tfile in files:
        base = os.path.basename(tfile.GetName())
        if any(base.endswith("_%s.root" % nick) for nick in nicknames):
            continue
        kept.append(tfile)
    return kept


def _origin_denominator_yield(files, region, origin_prefix, out_name):
    """Integral (+ variance) of one Z+1L origin denominator histogram."""
    hist = _clone_sum(
        files, "%s%s_d" % (origin_prefix, region), out_name
    )
    return _safe_integral(hist), _integral_variance(hist)


def compute_residual_prompt_sf(files, tex_handle):
    """
    Residual-prompt SF in Z+1L (denominator):

      SF = prompt(MC without WZ) / prompt(all MC)

    Numerator: prompt yield from ZZ+TT+DY (WZ excluded).
    Denominator: prompt yield from all MC samples (ZZ+WZ+TT+DY).
    Printed, written to FR_residual_prompt_sf.txt, and appended as LaTeX.
    Diagnostic only. The CR-transfer drops prompt from the mix instead
    of subtracting WZ, so this SF is not applied separately.
    """
    files_no_wz = _files_without_wz(files)
    if not files_no_wz:
        print(
            "WARNING: no non-WZ Hist_ID files — skipping residual-prompt SF"
        )
        return {}
    if len(files_no_wz) == len(files):
        print(
            "WARNING: WZ Hist_ID file missing — SF still uses available prompt(all)"
        )

    results = {}
    print("\n=== Residual-prompt scale factor (Z+1L) ===")
    print(
        "Definition: SF = prompt(MC, WZ excluded) / prompt(all MC samples)"
    )
    print("Diagnostic only: not applied as a correction or separate systematic.")

    for region in FR_REGIONS:
        num, var_num = _origin_denominator_yield(
            files_no_wz, region, "Hist_prompt", "prompt_d_nowz_" + region
        )
        den, var_den = _origin_denominator_yield(
            files, region, "Hist_prompt", "prompt_d_all_" + region
        )
        if den <= 0:
            print("  %s: prompt(all MC) denominator empty — SF undefined" % region)
            continue
        sf = num / den
        rel2 = 0.0
        if num > 0:
            rel2 += var_num / (num * num)
        rel2 += var_den / (den * den)
        sf_err = sf * math.sqrt(max(0.0, rel2))
        results[region] = {
            "prompt_no_wz": num,
            "prompt_all": den,
            "var_num": var_num,
            "var_den": var_den,
            "sf": sf,
            "sf_err": sf_err,
        }
        print(
            "  %s (%s): prompt(no WZ)=%.5f, prompt(all)=%.5f, "
            "SF=%.4f ± %.4f (%.2f%%)"
            % (
                region,
                FR_REGION_LABELS[region],
                num,
                den,
                sf,
                sf_err,
                100.0 * sf,
            )
        )

    for flavour, regions in (
        ("electron", ("EB", "EE")),
        ("muon", ("MB", "ME")),
    ):
        num = sum(
            results[region]["prompt_no_wz"]
            for region in regions
            if region in results
        )
        den = sum(
            results[region]["prompt_all"]
            for region in regions
            if region in results
        )
        var_num = sum(
            results[region]["var_num"]
            for region in regions
            if region in results
        )
        var_den = sum(
            results[region]["var_den"]
            for region in regions
            if region in results
        )
        if den <= 0:
            continue
        sf = num / den
        rel2 = 0.0
        if num > 0:
            rel2 += var_num / (num * num)
        rel2 += var_den / (den * den)
        sf_err = sf * math.sqrt(max(0.0, rel2))
        key = "combined_" + flavour
        results[key] = {
            "prompt_no_wz": num,
            "prompt_all": den,
            "var_num": var_num,
            "var_den": var_den,
            "sf": sf,
            "sf_err": sf_err,
        }
        print(
            "  combined %s: prompt(no WZ)=%.5f, prompt(all)=%.5f, "
            "SF=%.4f ± %.4f (%.2f%%)"
            % (flavour, num, den, sf, sf_err, 100.0 * sf)
        )

    txt_path = os.path.join(OUT_DIR, "FR_residual_prompt_sf.txt")
    with open(txt_path, "w") as text_handle:
        text_handle.write(
            "Residual-prompt scale factor (Z+1L denominator)\n"
            "SF = prompt(MC without WZ) / prompt(all MC samples)\n"
            "Diagnostic only; not applied as a correction or separate "
            "systematic.\n\n"
        )
        for key in list(FR_REGIONS) + [
            "combined_electron",
            "combined_muon",
        ]:
            if key not in results:
                continue
            entry = results[key]
            label = FR_REGION_LABELS.get(key, key)
            text_handle.write(
                "%s: prompt_no_wz=%.6f  prompt_all=%.6f  "
                "SF=%.6f  SF_err=%.6f  SF_percent=%.3f%%\n"
                % (
                    label,
                    entry["prompt_no_wz"],
                    entry["prompt_all"],
                    entry["sf"],
                    entry["sf_err"],
                    100.0 * entry["sf"],
                )
            )
    print("Wrote %s" % txt_path)

    if tex_handle is not None:
        tex_handle.write("\\begin{table}\n")
        tex_handle.write("\\centering\n")
        tex_handle.write(
            "\\caption{Residual-prompt scale factor in Z+1L: "
            "$\\mathrm{SF}=N_{\\mathrm{prompt}}(\\mathrm{MC\\ without\\ WZ})"
            "/N_{\\mathrm{prompt}}(\\mathrm{all\\ MC})$. "
            "This is a diagnostic only and is not applied as a correction "
            "or separate systematic.}\n"
        )
        tex_handle.write("\\label{tab:residual_prompt_sf}\n")
        tex_handle.write("\\begin{tabular}{|c|c|c|c|}\n")
        tex_handle.write("\\hline\n")
        tex_handle.write(
            "Region & Prompt (no WZ) & Prompt (all MC) & "
            "SF $\\pm$ (stat) \\\\ \\hline\n"
        )
        row_keys = list(FR_REGIONS) + [
            "combined_electron",
            "combined_muon",
        ]
        for key in row_keys:
            if key not in results:
                continue
            entry = results[key]
            label = FR_REGION_LABELS.get(
                key,
                key.replace("combined_", "combined ").replace("_", " "),
            )
            tex_handle.write(
                "%s & %.5f & %.5f & $%.4f \\pm %.4f$ \\\\ \\hline\n"
                % (
                    label,
                    entry["prompt_no_wz"],
                    entry["prompt_all"],
                    entry["sf"],
                    entry["sf_err"],
                )
            )
        tex_handle.write("\\end{tabular}\n")
        tex_handle.write("\\end{table}\n\n")

    return results


def _per_bin_composition_weighted_fr(n_by_origin, d_by_origin, keys, out_name):
    """Per-bin FR = sum_o w_o FR_o with w_o = D_o / sum D (exact identity with N/D)."""
    template = None
    for key in keys:
        if d_by_origin.get(key) is not None:
            template = d_by_origin[key]
            break
    if template is None:
        return None
    out = template.Clone(out_name)
    out.Reset("ICES")
    out.SetDirectory(0)
    for ibin in range(1, out.GetNbinsX() + 1):
        d_sum = 0.0
        n_sum = 0.0
        for key in keys:
            d_hist = d_by_origin.get(key)
            n_hist = n_by_origin.get(key)
            if not d_hist or not n_hist:
                continue
            d_bin = d_hist.GetBinContent(ibin)
            n_bin = n_hist.GetBinContent(ibin)
            if d_bin > 0:
                d_sum += d_bin
                n_sum += n_bin
        if d_sum <= 0:
            out.SetBinContent(ibin, 0.0)
            out.SetBinError(ibin, 0.0)
            continue
        # Equivalent to sum_o (D_o/D) * (N_o/D_o); keep the explicit form.
        value = 0.0
        for key in keys:
            d_hist = d_by_origin.get(key)
            n_hist = n_by_origin.get(key)
            if not d_hist or not n_hist:
                continue
            d_bin = d_hist.GetBinContent(ibin)
            if d_bin <= 0:
                continue
            weight = d_bin / d_sum
            value += weight * (n_hist.GetBinContent(ibin) / d_bin)
        # Binomial-like error from the equivalent total counts in this bin.
        err = math.sqrt(max(0.0, value * (1.0 - value) / d_sum))
        out.SetBinContent(ibin, value)
        out.SetBinError(ibin, err)
    return out


def _integrated_ratio(numerator, denominator):
    den = _safe_integral(denominator)
    if den <= 0:
        return 0.0
    return _safe_integral(numerator) / den


def _mc_origin_histograms(files, region, prefix):
    """Return summed numerator/denominator histograms for each MC origin."""
    output = {}
    for key, hist_prefix in (
        ("Fake", "Hist_fakes"),
        ("BDfake", "Hist_BDfakes"),
        ("Conv", "Hist_conv"),
        ("Prompt", "Hist_prompt"),
    ):
        output[key] = {
            "n": _clone_sum(
                files,
                "%s%s_n" % (hist_prefix, region),
                "%s_%s_%s_n" % (prefix, key, region),
            ),
            "d": _clone_sum(
                files,
                "%s%s_d" % (hist_prefix, region),
                "%s_%s_%s_d" % (prefix, key, region),
            ),
        }
    return output


def _combined_origin_fr(origin_hists, origins, out_name):
    numerator = _add_histograms(
        [origin_hists[key]["n"] for key in origins], out_name + "_n"
    )
    denominator = _add_histograms(
        [origin_hists[key]["d"] for key in origins], out_name + "_d"
    )
    if not numerator or not denominator:
        raise RuntimeError("Cannot build %s: missing origin histograms" % out_name)
    return (
        _ratio_hist(numerator, denominator, out_name),
        _integrated_ratio(numerator, denominator),
    )


NONPROMPT_ORIGIN_KEYS = ("Fake", "BDfake", "Conv")
ALL_ORIGIN_KEYS = ("Fake", "BDfake", "Conv", "Prompt")


_ODDS_FR_CAP = 0.999999


def _origin_odds(fr):
    """f/(1-f) of one origin FR (capped just below 1)."""
    fr = min(max(fr, 0.0), _ODDS_FR_CAP)
    return fr / (1.0 - fr)


def _origin_odds_error(fr, fr_error):
    """Error of f/(1-f): df / (1-f)^2."""
    fr = min(max(fr, 0.0), _ODDS_FR_CAP)
    return fr_error / ((1.0 - fr) ** 2)


def _odds_to_effective_fr(odds):
    """Effective FR f_eff such that f_eff/(1-f_eff) = odds."""
    return odds / (1.0 + odds)


# CR composition fractions c_i count the FAILING leptons of each origin.
# Each failing lepton is weighted by F_i = f_i/(1-f_i), so the CR-mix weight
# is the composition average of the per-origin F_i, not f/(1-f) of the
# average f. (Z+1L F = N/(D-N) already is the failing-lepton-weighted mean
# of F_i.) The three helpers below return the effective FR
# f_eff = F/(1+F) with F = sum_i c_i F_i, so that f_eff/(1-f_eff) is the
# correct CR-mix F wherever the code converts f to F afterwards.
def _effective_cr_fr(origin_hists, composition):
    """CR-mix effective FR using only HF/LF/conversion (prompt dropped)."""
    raw = {key: composition[key] / 100.0 for key in NONPROMPT_ORIGIN_KEYS}
    total = sum(raw.values())
    if total <= 0:
        return 0.0
    odds = 0.0
    for key in NONPROMPT_ORIGIN_KEYS:
        origin_fr = _integrated_ratio(
            origin_hists[key]["n"], origin_hists[key]["d"]
        )
        odds += (raw[key] / total) * _origin_odds(origin_fr)
    return _odds_to_effective_fr(odds)


def _composition_weighted_fr(
    perc, covariance, origin_frs, origin_fr_variances, keys
):
    """CR-composition effective FR for the given origin keys.

    Fractions are taken from ``perc`` (percent, failing leptons) and
    renormalized to sum to 1 over ``keys``. The CR-transfer uses
    Fake+BD+Conv only. The mix is averaged in F=f/(1-f) space and returned
    as f_eff = F/(1+F) (see the note above ``_effective_cr_fr``).
    """
    raw = {key: perc[key] / 100.0 for key in keys}
    total = sum(raw.values())
    if total <= 0:
        return 0.0, 0.0
    fractions = {key: raw[key] / total for key in keys}
    origin_odds = {key: _origin_odds(origin_frs[key]) for key in keys}
    odds = sum(fractions[key] * origin_odds[key] for key in keys)
    odds_variance = sum(
        fractions[key] ** 2
        * _origin_odds_error(
            origin_frs[key], math.sqrt(max(0.0, origin_fr_variances[key]))
        ) ** 2
        for key in keys
    )
    # Constant-total approximation: Cov(F_i/S, F_j/S) ≈ Cov(F_i, F_j) / S^2.
    composition_variance = sum(
        origin_odds[key_i]
        * origin_odds[key_j]
        * covariance[key_i][key_j]
        / (total * total)
        for key_i in keys
        for key_j in keys
    )
    odds_error = math.sqrt(max(0.0, odds_variance + composition_variance))
    value = _odds_to_effective_fr(odds)
    return value, odds_error / ((1.0 + odds) ** 2)


def _cms_prelim(canvas):
    tex1 = ROOT.TLatex(0.95, 0.96, LUMI_LABEL)
    tex1.SetNDC()
    tex1.SetTextAlign(31)
    tex1.SetTextFont(42)
    tex1.SetTextSize(0.03)
    tex1.Draw("goff")
    tex2 = ROOT.TLatex(0.15, 0.96, "CMS")
    tex2.SetNDC()
    tex2.SetTextFont(61)
    tex2.SetTextSize(0.0375)
    tex2.Draw("goff")
    tex3 = ROOT.TLatex(0.23, 0.96, "Preliminary")
    tex3.SetNDC()
    tex3.SetTextFont(52)
    tex3.SetTextSize(0.0285)
    tex3.Draw("goff")
    # Keep references alive until SaveAs
    canvas._cms_tex = (tex1, tex2, tex3)


def _style_canvas(name):
    c = ROOT.TCanvas(name, "myPlots", 0, 0, 800, 600)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetBorderSize(2)
    c.SetTickx(1)
    c.SetTicky(1)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.13)
    return c


def write_latex_ul_err(
    handle, caption, titles, values, uerrors, lerrors,
    value_label="Percentage"
):
    handle.write("\\begin{table}\n")
    handle.write("\\centering\n")
    handle.write("\\caption{%s}\n" % caption)
    handle.write("\\label{Labela}\n")
    handle.write("\\begin{tabular}{|c|c|}\n")
    handle.write("\\hline\n")
    handle.write(
        "\\textbf{Type} & \\textbf{%s} \\\\ \\hline\n" % value_label
    )
    for title, value, uerr, lerr in zip(titles, values, uerrors, lerrors):
        handle.write(
            "\\textbf{%s} & \\textbf{$%.5f^{+%.5f}_{-%.5f}$} \\\\ \\hline\n"
            % (title, value, uerr, lerr)
        )
    handle.write("\\end{tabular}\n")
    handle.write("\\end{table}\n\n")
    print("\n=== %s ===" % caption)
    for title, value, uerr, lerr in zip(titles, values, uerrors, lerrors):
        print("  %s: %.5f +%.5f -%.5f" % (title, value, uerr, lerr))


def _draw_id_stack(h_prompt, h_fake, h_conv, h_bdfake, out_base, header, xtitle, ytitle):
    if not all((h_prompt, h_fake, h_conv, h_bdfake)):
        print("Skipping stack %s: missing histogram" % out_base)
        return

    c = _style_canvas("c_" + out_base)
    leg = ROOT.TLegend(0.50, 0.65, 0.90, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetHeader(header)

    for hist, color, label in (
        (h_prompt, COLOR_PROMPT, "Prompt"),
        (h_fake, COLOR_FAKE, "Fakes"),
        (h_conv, COLOR_CONV, "Conversion"),
        (h_bdfake, COLOR_BDFAKE, "BD Fakes"),
    ):
        hist.SetFillColor(color)
        hist.SetLineColor(ROOT.kBlack)
        hist.SetFillStyle(1001)
        hist.SetLineWidth(1)
        leg.AddEntry(hist, label, "F")

    stack = ROOT.THStack("stack_" + out_base, "")
    stack.Add(h_conv)
    stack.Add(h_prompt)
    stack.Add(h_fake)
    stack.Add(h_bdfake)
    stack.Draw("hist")
    stack.GetXaxis().SetTitle(xtitle)
    stack.GetYaxis().SetTitle(ytitle)
    stack.GetXaxis().SetTitleSize(0.05)
    stack.GetYaxis().SetTitleSize(0.05)
    stack.GetXaxis().SetLabelSize(0.05)
    stack.GetYaxis().SetLabelSize(0.05)
    stack.GetXaxis().SetTitleOffset(1.1)
    stack.GetYaxis().SetTitleOffset(1.35)
    ymax = max(stack.GetMaximum(), 1.0)
    stack.SetMaximum(1.35 * ymax)
    leg.Draw("goff")
    _cms_prelim(c)
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".pdf"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".C"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".root"))


def _draw_fr_origin_overlay(fr_prompt, fr_fake, fr_conv, fr_bdfake, region, var_name):
    if not all((fr_prompt, fr_fake, fr_conv, fr_bdfake)):
        print("Skipping FR overlay %s: missing FR hist" % region)
        return

    c = _style_canvas("c_FR_origin_" + region)
    leg = ROOT.TLegend(0.50, 0.65, 0.90, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetHeader(
        "pp #rightarrow Z + l_{loose}, " + FR_REGION_LABELS[region]
    )

    frame = None
    for hist, color, style, label in (
        (fr_prompt, COLOR_PROMPT, 20, "FR Prompt"),
        (fr_fake, COLOR_FAKE, 21, "FR Fakes"),
        (fr_conv, COLOR_CONV, 22, "FR Conversion"),
        (fr_bdfake, COLOR_BDFAKE, 23, "FR BD fakes"),
    ):
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetMarkerStyle(style)
        hist.SetLineWidth(2)
        hist.SetFillStyle(0)
        hist.GetXaxis().SetTitle("p_{T}(l_{3}) [GeV]")
        hist.GetYaxis().SetTitle("Fake Rate")
        hist.GetXaxis().SetTitleSize(0.05)
        hist.GetYaxis().SetTitleSize(0.05)
        hist.GetXaxis().SetLabelSize(0.05)
        hist.GetYaxis().SetLabelSize(0.05)
        hist.SetMinimum(0.0)
        hist.SetMaximum(1.05)
        if frame is None:
            frame = hist
            hist.Draw("E1")
        else:
            hist.Draw("E1 same")
        leg.AddEntry(hist, label, "LP")

    leg.Draw("goff")
    _cms_prelim(c)
    out_base = "FR_%s_%s" % (var_name, region)
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".pdf"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".C"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".root"))


def _style_fr_curve(hist, color, marker, label, legend):
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(marker)
    hist.SetLineWidth(2)
    hist.SetFillStyle(0)
    hist.GetXaxis().SetTitle("p_{T}(l_{3}) [GeV]")
    hist.GetYaxis().SetTitle("Fake Rate")
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetLabelSize(0.05)
    hist.GetYaxis().SetLabelSize(0.05)
    legend.AddEntry(hist, label, "LP")


def _draw_z1l_prompt_subtraction_study(region, curves):
    """Z+1L validation: total FR vs weighted FR, and WZ-subtracted vs no-prompt."""
    required = ("mc_total", "weighted", "mc_no_wz", "weighted_no_prompt")
    if not all(curves.get(key) for key in required):
        print(
            "Skipping Z+1L prompt-subtraction study for %s: missing curve"
            % region
        )
        return

    c = _style_canvas("c_FR_Z1L_prompt_study_" + region)
    leg = ROOT.TLegend(0.36, 0.58, 0.92, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.026)
    leg.SetHeader("Z+1L, " + FR_REGION_LABELS[region])

    styled = (
        (curves["mc_total"], ROOT.kBlack, 20, "MC total FR"),
        (curves["weighted"], ROOT.kBlue + 1, 24, "Weighted FR"),
        (curves["mc_no_wz"], ROOT.kRed + 1, 21, "MC total FR, WZ removed"),
        (
            curves["weighted_no_prompt"],
            ROOT.kGreen + 2,
            25,
            "Weighted FR, no prompt",
        ),
    )
    ymax = 0.0
    for hist, _color, _marker, _label in styled:
        for ibin in range(1, hist.GetNbinsX() + 1):
            ymax = max(
                ymax,
                hist.GetBinContent(ibin) + hist.GetBinError(ibin),
            )
    for index, (hist, color, marker, label) in enumerate(styled):
        _style_fr_curve(hist, color, marker, label, leg)
        hist.SetMinimum(0.0)
        hist.SetMaximum(max(1.35 * ymax, 0.05))
        # Markers only (no connecting line between bins).
        hist.SetLineWidth(1)
        hist.Draw("E1 X0" if index == 0 else "E1 X0 same")

    leg.Draw("goff")
    _cms_prelim(c)
    out_base = "FR_Z1L_prompt_study_" + region
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".pdf"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".C"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".root"))
    c._z1l_prompt_study_keep = styled + (leg,)


def _cr_flavour_section(cr_kind, region):
    """Map detector region to flavour-combined CR section (2P2Fe/…, 3P1Fmu/…)."""
    flavour = "e" if region in ("EB", "EE") else "mu"
    return "%s%s" % (cr_kind, flavour)


def _cr_combined_flavour_section(region):
    """2P2F+3P1F origin mix for the fake-rate flavour of this region."""
    return "electron" if region in ("EB", "EE") else "muon"


def _build_cr_mix(
    cr_files,
    source_frs,
    source_integrated,
    source_variances,
    cr_kind,
    region,
    tag_label,
    origin_keys=None,
    cr_treatment=None,
):
    """Composition-weighted CR fake rate.

    Adopted mix is 2P2F-only. origin_keys is Fake+BD+Conv (noprompt) or
    all origins including leftover prompt.
    """
    origin_keys = origin_keys or NONPROMPT_ORIGIN_KEYS
    if cr_kind == "combined":
        section = _cr_combined_flavour_section(region)
        default_treatment = "2P2F+3P1F, all processes, prompt dropped"
    elif cr_kind.startswith("2P2F"):
        section = _cr_flavour_section("2P2F", region)
        default_treatment = "2P2F only"
    elif cr_kind.startswith("3P1F"):
        section = _cr_flavour_section("3P1F", region)
        default_treatment = "3P1F only"
    else:
        section = _cr_flavour_section(cr_kind, region)
        default_treatment = "all processes"
    if cr_treatment is None:
        cr_treatment = default_treatment
    composition, covariance = _cr_composition_with_covariance(
        cr_files,
        section,
        "transfer_%s_%s_%s" % (tag_label, cr_kind, region),
    )
    if not composition or not covariance:
        raise RuntimeError("Missing %s composition for %s" % (cr_kind, region))
    cr_fr = _cr_weighted_fr_curve(
        source_frs,
        composition,
        origin_keys,
        "FR_%s_transfer_%s_%s" % (cr_kind, tag_label, region),
        covariance=covariance,
    )
    cr_value, cr_error = _composition_weighted_fr(
        composition,
        covariance,
        source_integrated,
        source_variances,
        origin_keys,
    )
    F_cr = _fr_to_odds(cr_value)
    F_cr_error = _odds_error(cr_value, cr_error)
    if F_cr is None:
        raise RuntimeError("Invalid integrated %s FR for %s" % (cr_kind, region))
    return {
        "composition": composition,
        "covariance": covariance,
        "fr": cr_fr,
        "f": cr_value,
        "f_err": cr_error,
        "F": F_cr,
        "F_err": F_cr_error,
        "treatment": cr_treatment,
    }


def _origin_fr_histograms(files, region, tag):
    """Per-origin Z+1L FR(pt) histograms from the given MC file set."""
    output = {}
    for key, hist_prefix in (
        ("Fake", "Hist_fakes"),
        ("BDfake", "Hist_BDfakes"),
        ("Conv", "Hist_conv"),
        ("Prompt", "Hist_prompt"),
    ):
        numerator = _clone_sum(
            files,
            "%s%s_n" % (hist_prefix, region),
            "n_%s_%s_%s" % (tag, key, region),
        )
        denominator = _clone_sum(
            files,
            "%s%s_d" % (hist_prefix, region),
            "d_%s_%s_%s" % (tag, key, region),
        )
        if not numerator or not denominator:
            return None
        output[key] = _ratio_hist(
            numerator,
            denominator,
            "FR_%s_%s_%s" % (tag, key, region),
        )
    return output


def _cr_composition_with_covariance(files, section, tag):
    """CR origin percentages and covariance (supports combined sections)."""
    if section in CR_COMBINED_SECTIONS:
        histograms = {
            "Fake": _clone_sum_cr_sections(
                files, "Fake", CR_COMBINED_SECTIONS[section],
                "crfake_%s_%s" % (tag, section),
            ),
            "BDfake": _clone_sum_cr_sections(
                files, "BDfake", CR_COMBINED_SECTIONS[section],
                "crbd_%s_%s" % (tag, section),
            ),
            "Conv": _clone_sum_cr_sections(
                files, "Conv", CR_COMBINED_SECTIONS[section],
                "crconv_%s_%s" % (tag, section),
            ),
            "Prompt": _clone_sum_cr_sections(
                files, "Prom", CR_COMBINED_SECTIONS[section],
                "crprom_%s_%s" % (tag, section),
            ),
        }
    else:
        histograms = {
            "Fake": _clone_sum(
                files, "h1D_m4l_Fake_" + section, "crfake_%s_%s" % (tag, section)
            ),
            "BDfake": _clone_sum(
                files, "h1D_m4l_BDfake_" + section, "crbd_%s_%s" % (tag, section)
            ),
            "Conv": _clone_sum(
                files, "h1D_m4l_Conv_" + section, "crconv_%s_%s" % (tag, section)
            ),
            "Prompt": _clone_sum(
                files, "h1D_m4l_Prom_" + section, "crprom_%s_%s" % (tag, section)
            ),
        }
    if not all(histograms.values()):
        return None, None
    fractions, covariance = _weighted_composition(histograms)
    if fractions is None:
        return None, None
    percentages = {key: 100.0 * fractions[key] for key in fractions}
    return percentages, covariance


def _cr_composition_percentages(files, section, tag):
    """CR origin percentages for a section (supports combined 2P2Fe/…)."""
    percentages, _covariance = _cr_composition_with_covariance(
        files, section, tag
    )
    return percentages


def _cr_weighted_fr_curve(
    fr_by_origin, perc, keys, out_name, covariance=None
):
    """Apply constant CR composition weights to per-bin origin FRs.

    The CR fractions count failing leptons, so each pT bin is averaged in
    F=f/(1-f) space and stored as the effective FR f_eff = F/(1+F)
    (f_eff/(1-f_eff) is then the correct CR-mix F).
    The CR fractions are integrated over m4l and therefore constant in pT.
    When supplied, ``covariance`` propagates their common statistical
    uncertainty into each pT bin.
    """
    raw = {key: perc[key] / 100.0 for key in keys}
    total = sum(raw.values())
    if total <= 0:
        return None
    weights = {key: raw[key] / total for key in keys}
    template = None
    for key in keys:
        if fr_by_origin.get(key) is not None:
            template = fr_by_origin[key]
            break
    if template is None:
        return None
    out = template.Clone(out_name)
    out.Reset("ICES")
    out.SetDirectory(0)
    for ibin in range(1, out.GetNbinsX() + 1):
        odds = 0.0
        err2 = 0.0
        for key in keys:
            hist = fr_by_origin.get(key)
            if not hist:
                continue
            odds += weights[key] * _origin_odds(hist.GetBinContent(ibin))
            err2 += (
                weights[key]
                * _origin_odds_error(
                    hist.GetBinContent(ibin), hist.GetBinError(ibin)
                )
            ) ** 2
        if covariance is not None:
            err2 += sum(
                _origin_odds(fr_by_origin[key_i].GetBinContent(ibin))
                * _origin_odds(fr_by_origin[key_j].GetBinContent(ibin))
                * covariance[key_i][key_j]
                / (total * total)
                for key_i in keys
                for key_j in keys
            )
        out.SetBinContent(ibin, _odds_to_effective_fr(odds))
        out.SetBinError(
            ibin, math.sqrt(max(0.0, err2)) / ((1.0 + odds) ** 2)
        )
    return out


def _wz_removed_data_path(fr_tag):
    if not VBS:
        return HTOZZ_WZ_REMOVED_FR_FILE
    if fr_tag is None:
        filename = "Hist_Data_ptl3_WZremoved.root"
    else:
        filename = "Hist_Data_ptl3_VBSFR_%s_WZremoved.root" % fr_tag
    return os.path.join(VBS_HIST_ID_DIR, filename)


def _fr_to_odds(fr):
    if fr < 0.0 or fr >= 1.0:
        return None
    return fr / (1.0 - fr)


def _odds_error(fr, fr_error):
    if fr < 0.0 or fr >= 1.0:
        return 0.0
    return fr_error / ((1.0 - fr) ** 2)


def _odds_ratio_hist(cr_fr, z1l_fr, out_name):
    """Return kappa=F_CR/F_Z1L and |kappa-1| per pT bin."""
    kappa = cr_fr.Clone(out_name)
    kappa.Reset("ICES")
    kappa.SetDirectory(0)
    delta = cr_fr.Clone(out_name.replace("kappa", "delta"))
    delta.Reset("ICES")
    delta.SetDirectory(0)

    for ibin in range(1, cr_fr.GetNbinsX() + 1):
        f_cr = cr_fr.GetBinContent(ibin)
        f_z1l = z1l_fr.GetBinContent(ibin)
        F_cr = _fr_to_odds(f_cr)
        F_z1l = _fr_to_odds(f_z1l)
        if F_cr is None or F_z1l is None or F_z1l <= 0.0:
            continue
        value = F_cr / F_z1l
        F_cr_error = _odds_error(f_cr, cr_fr.GetBinError(ibin))
        F_z1l_error = _odds_error(f_z1l, z1l_fr.GetBinError(ibin))
        # Approximate MC-stat error: common Z+1L origin inputs are treated
        # as independent here. The central kappa is the systematic size;
        # this error is stored as a diagnostic.
        variance = (
            (F_cr_error / F_z1l) ** 2
            + (F_cr * F_z1l_error / (F_z1l * F_z1l)) ** 2
        )
        error = math.sqrt(max(0.0, variance))
        kappa.SetBinContent(ibin, value)
        kappa.SetBinError(ibin, error)
        delta.SetBinContent(ibin, abs(value - 1.0))
        delta.SetBinError(ibin, error)
    return kappa, delta


def _data_band_graph(data_fr, delta, out_name):
    """Map a relative odds-space uncertainty into an asymmetric FR band."""
    graph = ROOT.TGraphAsymmErrors(data_fr.GetNbinsX())
    graph.SetName(out_name)
    axis = data_fr.GetXaxis()
    for ibin in range(1, data_fr.GetNbinsX() + 1):
        point = ibin - 1
        fr = data_fr.GetBinContent(ibin)
        rel = max(0.0, delta.GetBinContent(ibin))
        odds = _fr_to_odds(fr)
        x = axis.GetBinCenter(ibin)
        ex = 0.5 * axis.GetBinWidth(ibin)
        if odds is None:
            graph.SetPoint(point, x, fr)
            graph.SetPointError(point, ex, ex, 0.0, 0.0)
            continue
        odds_down = max(0.0, odds * (1.0 - rel))
        odds_up = odds * (1.0 + rel)
        fr_down = odds_down / (1.0 + odds_down)
        fr_up = odds_up / (1.0 + odds_up)
        graph.SetPoint(point, x, fr)
        graph.SetPointError(
            point,
            ex,
            ex,
            max(0.0, fr - fr_down),
            max(0.0, fr_up - fr),
        )
    return graph


def _draw_cr_prompt_subtraction_study(
    cr_kind, region, process_label, curves, reference_label=None
):
    """MC-only CR study: reference weighted FR vs no-prompt weighted FR."""
    required = ("mc_reference", "weighted_no_prompt")
    if not all(curves.get(key) for key in required):
        print(
            "Skipping %s prompt-subtraction study for %s: missing curve"
            % (cr_kind, region)
        )
        return

    if reference_label is None:
        if process_label:
            reference_label = "CR-weighted FR, %s removed" % process_label
        else:
            reference_label = "CR-weighted FR (WZ included)"

    c = _style_canvas("c_FR_%s_prompt_study_%s" % (cr_kind, region))
    leg = ROOT.TLegend(0.36, 0.64, 0.92, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.026)
    leg.SetHeader("%s, %s" % (cr_kind, FR_REGION_LABELS[region]))

    styled = (
        (
            curves["mc_reference"],
            ROOT.kRed + 1,
            21,
            reference_label,
        ),
        (
            curves["weighted_no_prompt"],
            ROOT.kGreen + 2,
            25,
            "CR-weighted FR, no prompt",
        ),
    )
    ymax = 0.0
    for hist, _color, _marker, _label in styled:
        for ibin in range(1, hist.GetNbinsX() + 1):
            ymax = max(
                ymax,
                hist.GetBinContent(ibin) + hist.GetBinError(ibin),
            )
    for index, (hist, color, marker, label) in enumerate(styled):
        _style_fr_curve(hist, color, marker, label, leg)
        hist.SetMinimum(0.0)
        hist.SetMaximum(max(1.35 * ymax, 0.05))
        hist.SetLineWidth(1)
        hist.Draw("E1 X0" if index == 0 else "E1 X0 same")

    leg.Draw("goff")
    _cms_prelim(c)
    out_base = "FR_%s_prompt_study_%s" % (cr_kind, region)
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".pdf"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".C"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".root"))
    c._cr_prompt_study_keep = styled + (leg,)


def _draw_vbs_transfer_curves(region, styled):
    """Draw Z+1L vs 2P2F-mix fake rates for noprompt and leftover-prompt."""
    c = _style_canvas("c_FR_CR_transfer_" + region)
    leg = ROOT.TLegend(0.28, 0.58, 0.92, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.024)
    leg.SetHeader(ANALYSIS_MODE + ", " + FR_REGION_LABELS[region])
    ymax = max(
        hist.GetBinContent(ibin) + hist.GetBinError(ibin)
        for hist, _color, _marker, _label in styled
        for ibin in range(1, hist.GetNbinsX() + 1)
    )
    for index, (hist, color, marker, label) in enumerate(styled):
        _style_fr_curve(hist, color, marker, label, leg)
        hist.SetMinimum(0.0)
        hist.SetMaximum(max(1.40 * ymax, 0.02))
        hist.SetLineWidth(1)
        hist.Draw("E1 X0" if index == 0 else "E1 X0 same")
    leg.Draw("goff")
    _cms_prelim(c)
    c._vbs_transfer_keep = tuple(styled) + (leg,)
    out_base = "FR_CR_transfer_" + region
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".pdf"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".C"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".root"))


def _draw_vbs_transfer_ratio(region, kappa_2p2f, kappa_3p1f):
    """Draw odds-space CR/Z+1L transfer ratios."""
    c = _style_canvas("c_FR_CR_transfer_ratio_" + region)
    frame = kappa_2p2f.Clone("frame_FR_CR_transfer_ratio_" + region)
    frame.Reset("ICES")
    frame.SetDirectory(0)
    frame.GetXaxis().SetTitle("p_{T}(l) [GeV]")
    frame.GetYaxis().SetTitle("#kappa_{R} = F_{R}^{MC}/F_{Z+1L}^{MC}")
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetYaxis().SetTitleOffset(1.35)

    nonzero = [
        hist.GetBinContent(ibin)
        for hist in (kappa_2p2f, kappa_3p1f)
        for ibin in range(1, hist.GetNbinsX() + 1)
        if hist.GetBinContent(ibin) > 0.0
    ]
    low = min(nonzero + [1.0])
    high = max(nonzero + [1.0])
    span = max(0.15, high - low)
    frame.SetMinimum(max(0.0, low - 0.30 * span))
    frame.SetMaximum(high + 0.45 * span)
    frame.Draw("axis")

    leg = ROOT.TLegend(0.42, 0.72, 0.92, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetHeader(ANALYSIS_MODE + ", " + FR_REGION_LABELS[region])
    for hist, color, marker, label in (
        (kappa_2p2f, ROOT.kBlue + 1, 21, "2P2F"),
        (kappa_3p1f, ROOT.kRed + 1, 22, "3P1F"),
    ):
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetMarkerStyle(marker)
        hist.SetLineWidth(1)
        hist.Draw("E1 X0 same")
        leg.AddEntry(hist, label, "LP")

    unity = ROOT.TLine(
        frame.GetXaxis().GetXmin(),
        1.0,
        frame.GetXaxis().GetXmax(),
        1.0,
    )
    unity.SetLineColor(ROOT.kGray + 2)
    unity.SetLineStyle(2)
    unity.SetLineWidth(2)
    unity.Draw("same")
    leg.AddEntry(unity, "#kappa=1", "L")
    leg.Draw("goff")
    _cms_prelim(c)
    c._vbs_ratio_keep = (
        frame,
        kappa_2p2f,
        kappa_3p1f,
        unity,
        leg,
    )
    out_base = "FR_CR_transfer_ratio_" + region
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".pdf"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".C"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".root"))


def _draw_vbs_data_band(region, cr_kind, data_fr, band):
    """Draw the unchanged data FR with the CR-transfer systematic envelope."""
    c = _style_canvas("c_FR_data_band_%s_%s" % (cr_kind, region))
    frame = data_fr.Clone("frame_FR_data_band_%s_%s" % (cr_kind, region))
    frame.Reset("ICES")
    frame.SetDirectory(0)
    frame.GetXaxis().SetTitle("p_{T}(l) [GeV]")
    frame.GetYaxis().SetTitle("Fake rate")
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetYaxis().SetTitleOffset(1.35)

    ymax = 0.0
    for ibin in range(1, data_fr.GetNbinsX() + 1):
        point = ibin - 1
        ymax = max(
            ymax,
            data_fr.GetBinContent(ibin) + data_fr.GetBinError(ibin),
            band.GetPointY(point) + band.GetErrorYhigh(point),
        )
    frame.SetMinimum(0.0)
    frame.SetMaximum(max(1.45 * ymax, 0.02))
    frame.Draw("axis")

    color = ROOT.kBlue + 1 if cr_kind == "2P2F" else ROOT.kRed + 1
    band.SetFillColorAlpha(color, 0.30)
    band.SetLineColor(color)
    band.SetMarkerSize(0)
    band.Draw("2 same")
    data_fr.SetLineColor(ROOT.kBlack)
    data_fr.SetMarkerColor(ROOT.kBlack)
    data_fr.SetMarkerStyle(20)
    data_fr.SetLineWidth(1)
    data_fr.Draw("E1 X0 same")

    leg = ROOT.TLegend(0.40, 0.72, 0.92, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetHeader(
        "%s %s, %s" % (ANALYSIS_MODE, cr_kind, FR_REGION_LABELS[region])
    )
    leg.AddEntry(data_fr, "Data FR (WZ removed)", "LP")
    leg.AddEntry(band, "MC CR-transfer band", "F")
    leg.Draw("goff")
    _cms_prelim(c)
    c._vbs_band_keep = (frame, data_fr, band, leg)
    out_base = "FR_data_band_%s_%s" % (cr_kind, region)
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".pdf"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".C"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".root"))


def _draw_z1l_closure(region, curves):
    """Compare data and MC Z+1L FRs before/after WZ removal."""
    c = _style_canvas("c_FR_Z1L_closure_" + region)
    leg = ROOT.TLegend(0.42, 0.58, 0.90, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetHeader("Z+1L, " + FR_REGION_LABELS[region])

    styled = (
        (curves["data_raw"], ROOT.kBlack, 20, "Data, before WZ removal"),
        (curves["data_wz_removed"], ROOT.kRed + 1, 21, "Data, WZ removed"),
        (curves["mc_all"], ROOT.kBlue + 1, 24, "MC all origins, incl. WZ"),
        (curves["mc_all_no_wz"], ROOT.kAzure + 7, 25, "MC all origins, no WZ"),
        (
            curves["mc_nonprompt"],
            ROOT.kGreen + 2,
            26,
            "MC nonprompt, incl. WZ",
        ),
        (
            curves["mc_nonprompt_no_wz"],
            ROOT.kMagenta + 1,
            32,
            "MC nonprompt, no WZ",
        ),
    )
    ymax = max(
        max(hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1))
        for hist, _color, _marker, _label in styled
    )
    for index, (hist, color, marker, label) in enumerate(styled):
        _style_fr_curve(hist, color, marker, label, leg)
        hist.SetMinimum(0.0)
        hist.SetMaximum(max(1.35 * ymax, 0.02))
        hist.Draw("E1" if index == 0 else "E1 same")

    leg.Draw("goff")
    _cms_prelim(c)
    out_base = "FR_Z1L_closure_" + region
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".pdf"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".C"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".root"))


def _draw_composition_fr_summary(region, titles, values, uerrs, lerrs):
    """Compare CR-weighted nonprompt FRs with the Z+1L nonprompt FR."""
    entries = {
        title: (value, uerr, lerr)
        for title, value, uerr, lerr in zip(titles, values, uerrs, lerrs)
        if not title.endswith(" + upper err")
    }
    if region in ("EB", "EE"):
        categories = [
            ("4e", "2P2F4e", "3P1F4e"),
            ("2#mu2e", "2P2F2mu2e", "3P1F2mu2e"),
        ]
    else:
        categories = [
            ("4#mu", "2P2F4mu", "3P1F4mu"),
            ("2e2#mu", "2P2F2e2mu", "3P1F2e2mu"),
        ]
    required = ["Z+1L nonprompt FR"] + [
        key
        for _label, key2, key3 in categories
        for key in (key2 + " (no prompt)", key3 + " (no prompt)")
    ]
    missing = [key for key in required if key not in entries]
    if missing:
        raise RuntimeError(
            "Missing FR summary entries for %s: %s"
            % (region, ", ".join(missing))
        )

    n = len(categories)
    c = _style_canvas("c_FR_composition_" + region)
    c.SetBottomMargin(0.17)
    graph_2p2f_np = ROOT.TGraphAsymmErrors(n)
    graph_3p1f_np = ROOT.TGraphAsymmErrors(n)
    ymax = entries["Z+1L nonprompt FR"][0]

    for i, (_label, key_2p2f, key_3p1f) in enumerate(categories):
        for graph, key, x_offset in (
            (graph_2p2f_np, key_2p2f + " (no prompt)", -0.10),
            (graph_3p1f_np, key_3p1f + " (no prompt)", 0.10),
        ):
            value, uerr, lerr = entries[key]
            graph.SetPoint(i, i + 0.5 + x_offset, value)
            graph.SetPointError(
                i, 0.0, 0.0, max(0.0, lerr), max(0.0, uerr)
            )
            ymax = max(ymax, value + max(0.0, uerr))

    frame = ROOT.TH1D(
        "frame_FR_composition_" + region, "", n, 0.0, float(n)
    )
    frame.SetDirectory(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(max(1.35 * ymax, 0.05))
    frame.SetStats(0)
    frame.GetYaxis().SetTitle("Composition-weighted fake rate")
    frame.GetXaxis().SetTitle("Control-region final state")
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetLabelSize(0.042)
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.GetXaxis().SetTitleOffset(1.35)
    for i, (label, _key2, _key3) in enumerate(categories):
        frame.GetXaxis().SetBinLabel(i + 1, label)

    for graph, color, marker in (
        (graph_2p2f_np, ROOT.kBlue + 1, 20),
        (graph_3p1f_np, ROOT.kRed + 1, 21),
    ):
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(1.2)
        graph.SetMarkerColor(color)
        graph.SetLineColor(color)
        graph.SetLineWidth(2)

    nonprompt_fr = entries["Z+1L nonprompt FR"][0]
    nonprompt_line = ROOT.TLine(0.0, nonprompt_fr, float(n), nonprompt_fr)
    nonprompt_line.SetLineColor(ROOT.kGray + 2)
    nonprompt_line.SetLineStyle(2)
    nonprompt_line.SetLineWidth(2)

    frame.Draw("AXIS")
    graph_2p2f_np.Draw("P Z SAME")
    graph_3p1f_np.Draw("P Z SAME")
    nonprompt_line.Draw("SAME")

    leg = ROOT.TLegend(0.42, 0.68, 0.92, 0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.028)
    leg.SetHeader(FR_REGION_LABELS[region].capitalize())
    leg.AddEntry(graph_2p2f_np, "2P2F weighted FR (no prompt)", "LP")
    leg.AddEntry(graph_3p1f_np, "3P1F weighted FR (no prompt)", "LP")
    leg.AddEntry(nonprompt_line, "Z+1L nonprompt FR", "L")
    leg.Draw("goff")

    c._fr_composition_keep = (
        frame,
        graph_2p2f_np,
        graph_3p1f_np,
        nonprompt_line,
        leg,
    )
    _cms_prelim(c)
    out_base = "FR_composition_" + region
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".pdf"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".C"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".root"))


def _draw_flavour_combined_fr_summary(region_entries):
    """Compare electron- and muon-combined nonprompt CR-weighted FRs."""
    regions = ("EB", "EE", "MB", "ME")
    labels = ("EB (e)", "EE (e)", "MB (#mu)", "ME (#mu)")
    n = len(regions)
    c = _style_canvas("c_FR_composition_flavour_combined")
    c.SetBottomMargin(0.17)
    graph_z1l_np = ROOT.TGraphAsymmErrors(n)
    graph_2p2f_np = ROOT.TGraphAsymmErrors(n)
    graph_3p1f_np = ROOT.TGraphAsymmErrors(n)
    ymax = 0.0

    for i, region in enumerate(regions):
        flavour = "e" if region in ("EB", "EE") else "mu"
        entries = region_entries[region]
        for graph, key, x_offset in (
            (graph_z1l_np, "Z+1L nonprompt FR", -0.15),
            (graph_2p2f_np, "2P2F" + flavour + " (no prompt)", 0.0),
            (graph_3p1f_np, "3P1F" + flavour + " (no prompt)", 0.15),
        ):
            if key not in entries:
                raise RuntimeError(
                    "Missing combined FR summary entry %s for %s"
                    % (key, region)
                )
            value, uerr, lerr = entries[key]
            graph.SetPoint(i, i + 0.5 + x_offset, value)
            graph.SetPointError(
                i, 0.0, 0.0, max(0.0, lerr), max(0.0, uerr)
            )
            ymax = max(ymax, value + max(0.0, uerr))

    frame = ROOT.TH1D(
        "frame_FR_composition_flavour_combined", "", n, 0.0, float(n)
    )
    frame.SetDirectory(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(max(1.45 * ymax, 0.05))
    frame.SetStats(0)
    frame.GetYaxis().SetTitle("Composition-weighted fake rate")
    frame.GetXaxis().SetTitle("Lepton flavour and detector region")
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.GetXaxis().SetTitleOffset(1.35)
    for i, label in enumerate(labels):
        frame.GetXaxis().SetBinLabel(i + 1, label)

    for graph, color, marker in (
        (graph_z1l_np, ROOT.kBlack, 20),
        (graph_2p2f_np, ROOT.kBlue + 1, 21),
        (graph_3p1f_np, ROOT.kRed + 1, 22),
    ):
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(1.15)
        graph.SetMarkerColor(color)
        graph.SetLineColor(color)
        graph.SetLineWidth(2)

    frame.Draw("AXIS")
    graph_z1l_np.Draw("P Z SAME")
    graph_2p2f_np.Draw("P Z SAME")
    graph_3p1f_np.Draw("P Z SAME")

    leg = ROOT.TLegend(0.36, 0.62, 0.92, 0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.026)
    leg.SetHeader("Combined failing-Z2 flavour categories")
    leg.AddEntry(graph_z1l_np, "Z+1L nonprompt FR", "LP")
    leg.AddEntry(graph_2p2f_np, "2P2F combined FR (no prompt)", "LP")
    leg.AddEntry(graph_3p1f_np, "3P1F combined FR (no prompt)", "LP")
    leg.AddEntry(0, "e combined: CR origins from 4e + 2#mu2e", "")
    leg.AddEntry(0, "#mu combined: CR origins from 4#mu + 2e2#mu", "")
    leg.Draw("goff")

    c._fr_flavour_combined_keep = (
        frame,
        graph_z1l_np,
        graph_2p2f_np,
        graph_3p1f_np,
        leg,
    )
    _cms_prelim(c)
    out_base = "FR_composition_flavour_combined"
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".pdf"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".C"))
    c.SaveAs(os.path.join(OUT_DIR, out_base + ".root"))


def _write_composition_summary_latex(handle, summaries):
    """Write absolute average-FR and CR-weighted FR values."""
    handle.write("\\begin{table}\n")
    handle.write("\\centering\n")
    handle.write(
        "\\caption{H$\\to$ZZ fake-rate composition with WZ excluded. "
        "The Z+1L reference and CR-weighted FRs include Fake/BD/Conv/Prompt "
        "origins, computed from the current stage-1 histograms.}\n"
    )
    handle.write("\\begin{tabular}{|c|c|c|c|}\n")
    handle.write("\\hline\n")
    handle.write(
        "Region/FS & Z+1L average FR & $f_{2\\mathrm{P}2\\mathrm{F}}$ & "
        "$f_{3\\mathrm{P}1\\mathrm{F}}$ \\\\ \\hline\n"
    )
    for region in FR_REGIONS:
        final_states = (
            ("4e", "2mu2e", "electron")
            if region in ("EB", "EE")
            else ("4mu", "2e2mu", "muon")
        )
        for final_state in final_states:
            entry = summaries[region][final_state]
            handle.write(
                "%s/%s & %.5f & %.5f & %.5f \\\\ \\hline\n"
                % (
                    region,
                    final_state,
                    summaries[region]["z1l_mc_all_no_wz"],
                    entry["f_2p2f"],
                    entry["f_3p1f"],
                )
            )
    handle.write("\\end{tabular}\n")
    handle.write("\\end{table}\n\n")


def compute_z1l_composition(files, tex_handle):
    """ denominator composition in Z+1L by region."""
    for region in FR_REGIONS:
        d_conv = _clone_sum(files, "Hist_conv%s_d" % region, "sum_conv_d_" + region)
        d_fake = _clone_sum(files, "Hist_fakes%s_d" % region, "sum_fake_d_" + region)
        d_bd = _clone_sum(files, "Hist_BDfakes%s_d" % region, "sum_bd_d_" + region)
        d_prom = _clone_sum(files, "Hist_prompt%s_d" % region, "sum_prom_d_" + region)
        if not all((d_conv, d_fake, d_bd, d_prom)):
            print("Skipping Z+1L composition for %s: missing hist" % region)
            continue

        histograms = {
            "Fake": d_fake,
            "BDfake": d_bd,
            "Conv": d_conv,
            "Prompt": d_prom,
        }
        fractions, covariance = _weighted_composition(histograms)
        if fractions is None:
            print("Skipping Z+1L composition for %s: empty" % region)
            continue

        titles, values, uerrs, lerrs = [], [], [], []
        for label, key in (
            ("Fakes", "Fake"),
            ("BD fakes", "BDfake"),
            ("Conversions", "Conv"),
            ("Prompt", "Prompt"),
        ):
            perc = 100.0 * fractions[key]
            error = 100.0 * math.sqrt(max(0.0, covariance[key][key]))
            titles.append(label)
            values.append(perc)
            uerrs.append(error)
            lerrs.append(error)

        write_latex_ul_err(
            tex_handle,
            "Z+1L denominator composition (%s)" % FR_REGION_LABELS[region],
            titles,
            values,
            uerrs,
            lerrs,
        )


def _clone_sum_cr_sections(files, origin, sections, out_name):
    histograms = []
    for section in sections:
        histogram = _clone_sum(
            files,
            "h1D_m4l_%s_%s" % (origin, section),
            "%s_%s" % (out_name, section),
        )
        if histogram:
            histograms.append(histogram)
    return _add_histograms(histograms, out_name)


def _store_cr_composition(
    section, h_fake, h_bd, h_conv, h_prompt, categor_perc,
    composition_cov, tex_handle
):
    if not all((h_fake, h_bd, h_conv, h_prompt)):
        print("Skipping CR composition for %s: missing hist" % section)
        return

    histograms = {
        "Fake": h_fake,
        "BDfake": h_bd,
        "Conv": h_conv,
        "Prompt": h_prompt,
    }
    fractions, covariance = _weighted_composition(histograms)
    if fractions is None:
        print("Skipping CR composition for %s: empty yield" % section)
        return

    titles, values, uerrs, lerrs = [], [], [], []
    categor_perc[section] = {}
    composition_cov[section] = covariance

    for label, key in (
        ("Fakes", "Fake"),
        ("BD fakes", "BDfake"),
        ("Conversions", "Conv"),
        ("Prompt", "Prompt"),
    ):
        perc = 100.0 * fractions[key]
        error = 100.0 * math.sqrt(max(0.0, covariance[key][key]))
        categor_perc[section][key] = perc
        titles.append(label)
        values.append(perc)
        uerrs.append(error)
        lerrs.append(error)

    write_latex_ul_err(
        tex_handle,
        section + " categorisation",
        titles,
        values,
        uerrs,
        lerrs,
    )


def compute_cr_composition(files, tex_handle):
    """CR composition of Fake/BD/Conv/Prompt origins."""
    categor_perc = {}
    composition_cov = {}

    for section in CR_SECTIONS:
        h_fake = _clone_sum(
            files, "h1D_m4l_Fake_" + section, "sum_fake_" + section
        )
        h_bd = _clone_sum(
            files, "h1D_m4l_BDfake_" + section, "sum_bd_" + section
        )
        h_conv = _clone_sum(
            files, "h1D_m4l_Conv_" + section, "sum_conv_" + section
        )
        h_prompt = _clone_sum(
            files, "h1D_m4l_Prom_" + section, "sum_prompt_" + section
        )
        _store_cr_composition(
            section,
            h_fake,
            h_bd,
            h_conv,
            h_prompt,
            categor_perc,
            composition_cov,
            tex_handle,
        )

    for section, source_sections in CR_COMBINED_SECTIONS.items():
        _store_cr_composition(
            section,
            _clone_sum_cr_sections(
                files, "Fake", source_sections, "sum_fake_" + section
            ),
            _clone_sum_cr_sections(
                files, "BDfake", source_sections, "sum_bd_" + section
            ),
            _clone_sum_cr_sections(
                files, "Conv", source_sections, "sum_conv_" + section
            ),
            _clone_sum_cr_sections(
                files, "Prom", source_sections, "sum_prompt_" + section
            ),
            categor_perc,
            composition_cov,
            tex_handle,
        )

    return categor_perc, composition_cov


def compute_weighted_fr(
    files, categor_perc, composition_cov, tex_handle, cr_files=None
):
    """Origin-averaged FR and composition-weighted FR.

    ``files`` provide Z+1L origin FR histograms (jet-tagged when applicable).
    ``cr_files`` provide CR origin yields for the process-removal study; defaults
    to ``files``. For VBS ge1j/ge2j, pass nominal CR files here.
    """
    if cr_files is None:
        cr_files = files
    region_entries = {}
    for region in FR_REGIONS:
        n_fake = _clone_sum(files, "Hist_fakes%s_n" % region, "n_fake_" + region)
        d_fake = _clone_sum(files, "Hist_fakes%s_d" % region, "d_fake_" + region)
        n_bd = _clone_sum(files, "Hist_BDfakes%s_n" % region, "n_bd_" + region)
        d_bd = _clone_sum(files, "Hist_BDfakes%s_d" % region, "d_bd_" + region)
        n_conv = _clone_sum(files, "Hist_conv%s_n" % region, "n_conv_" + region)
        d_conv = _clone_sum(files, "Hist_conv%s_d" % region, "d_conv_" + region)
        n_prom = _clone_sum(files, "Hist_prompt%s_n" % region, "n_prom_" + region)
        d_prom = _clone_sum(files, "Hist_prompt%s_d" % region, "d_prom_" + region)

        if not all((n_fake, d_fake, n_bd, d_bd, n_conv, d_conv, n_prom, d_prom)):
            print("Skipping weighted FR for %s: missing hist" % region)
            continue

        # ID stacks (numerator / denominator) 
        _draw_id_stack(
            n_prom.Clone("plot_n_prom_" + region),
            n_fake.Clone("plot_n_fake_" + region),
            n_conv.Clone("plot_n_conv_" + region),
            n_bd.Clone("plot_n_bd_" + region),
            "ID_%s%s_n" % (VAR_NAME, region),
            "pp #rightarrow Z + l_{loose} (numerator)",
            "p_{T}(l_{3}) [GeV]",
            "Events",
        )
        _draw_id_stack(
            d_prom.Clone("plot_d_prom_" + region),
            d_fake.Clone("plot_d_fake_" + region),
            d_conv.Clone("plot_d_conv_" + region),
            d_bd.Clone("plot_d_bd_" + region),
            "ID_%s%s_d" % (VAR_NAME, region),
            "pp #rightarrow Z + l_{loose} (denominator)",
            "p_{T}(l_{3}) [GeV]",
            "Events",
        )

        # Origin-dependent FR histograms for overlay plots
        fr_fake = _ratio_hist(n_fake, d_fake, "FR_fake_" + region)
        fr_bd = _ratio_hist(n_bd, d_bd, "FR_bd_" + region)
        fr_conv = _ratio_hist(n_conv, d_conv, "FR_conv_" + region)
        fr_prom = _ratio_hist(n_prom, d_prom, "FR_prom_" + region)
        _draw_fr_origin_overlay(fr_prom, fr_fake, fr_conv, fr_bd, region, VAR_NAME)

        n_by_origin = {
            "Fake": n_fake,
            "BDfake": n_bd,
            "Conv": n_conv,
            "Prompt": n_prom,
        }
        d_by_origin = {
            "Fake": d_fake,
            "BDfake": d_bd,
            "Conv": d_conv,
            "Prompt": d_prom,
        }
        n_all = _add_histograms(
            (n_fake, n_bd, n_conv, n_prom), "n_all_" + region
        )
        d_all = _add_histograms(
            (d_fake, d_bd, d_conv, d_prom), "d_all_" + region
        )
        fr_mc_total = _ratio_hist(n_all, d_all, "FR_mc_total_" + region)
        fr_weighted = _per_bin_composition_weighted_fr(
            n_by_origin,
            d_by_origin,
            ALL_ORIGIN_KEYS,
            "FR_weighted_" + region,
        )
        fr_weighted_no_prompt = _per_bin_composition_weighted_fr(
            n_by_origin,
            d_by_origin,
            NONPROMPT_ORIGIN_KEYS,
            "FR_weighted_noprompt_" + region,
        )

        files_no_wz = _files_without_wz(files)
        n_fake_nw = _clone_sum(
            files_no_wz, "Hist_fakes%s_n" % region, "n_fake_nowz_" + region
        )
        d_fake_nw = _clone_sum(
            files_no_wz, "Hist_fakes%s_d" % region, "d_fake_nowz_" + region
        )
        n_bd_nw = _clone_sum(
            files_no_wz, "Hist_BDfakes%s_n" % region, "n_bd_nowz_" + region
        )
        d_bd_nw = _clone_sum(
            files_no_wz, "Hist_BDfakes%s_d" % region, "d_bd_nowz_" + region
        )
        n_conv_nw = _clone_sum(
            files_no_wz, "Hist_conv%s_n" % region, "n_conv_nowz_" + region
        )
        d_conv_nw = _clone_sum(
            files_no_wz, "Hist_conv%s_d" % region, "d_conv_nowz_" + region
        )
        n_prom_nw = _clone_sum(
            files_no_wz, "Hist_prompt%s_n" % region, "n_prom_nowz_" + region
        )
        d_prom_nw = _clone_sum(
            files_no_wz, "Hist_prompt%s_d" % region, "d_prom_nowz_" + region
        )
        fr_mc_no_wz = None
        if all(
            (
                n_fake_nw,
                d_fake_nw,
                n_bd_nw,
                d_bd_nw,
                n_conv_nw,
                d_conv_nw,
                n_prom_nw,
                d_prom_nw,
            )
        ):
            n_all_nw = _add_histograms(
                (n_fake_nw, n_bd_nw, n_conv_nw, n_prom_nw),
                "n_all_nowz_" + region,
            )
            d_all_nw = _add_histograms(
                (d_fake_nw, d_bd_nw, d_conv_nw, d_prom_nw),
                "d_all_nowz_" + region,
            )
            fr_mc_no_wz = _ratio_hist(
                n_all_nw, d_all_nw, "FR_mc_nowz_" + region
            )
        _draw_z1l_prompt_subtraction_study(
            region,
            {
                "mc_total": fr_mc_total,
                "weighted": fr_weighted,
                "mc_no_wz": fr_mc_no_wz,
                "weighted_no_prompt": fr_weighted_no_prompt,
            },
        )

        # CR studies (MC only). Both CRs use the full MC mix; prompt is
        # dropped from the weighted rate, not by removing a process.
        fr_source_no_wz = _origin_fr_histograms(
            files_no_wz, region, "source_no_wz"
        )
        if not fr_source_no_wz:
            print(
                "Skipping CR prompt studies for %s: missing no-WZ source FRs"
                % region
            )
            continue
        for cr_kind in ("2P2F", "3P1F"):
            section = _cr_flavour_section(cr_kind, region)
            if section not in categor_perc:
                print(
                    "Skipping %s prompt study for %s: missing %s composition"
                    % (cr_kind, region, section)
                )
                continue
            perc_all = categor_perc[section]
            weighted_no_prompt = _cr_weighted_fr_curve(
                fr_source_no_wz,
                perc_all,
                NONPROMPT_ORIGIN_KEYS,
                "FR_%s_noprompt_%s" % (cr_kind, region),
            )
            mc_reference = _cr_weighted_fr_curve(
                fr_source_no_wz,
                perc_all,
                ALL_ORIGIN_KEYS,
                "FR_%s_inclusive_%s" % (cr_kind, region),
            )
            _draw_cr_prompt_subtraction_study(
                cr_kind,
                region,
                None,
                {
                    "mc_reference": mc_reference,
                    "weighted_no_prompt": weighted_no_prompt,
                },
                reference_label=(
                    "CR-weighted FR, all processes (prompt included)"
                ),
            )

        n_nonprompt = _add_histograms(
            (n_fake, n_bd, n_conv), "n_nonprompt_" + region
        )
        d_nonprompt = _add_histograms(
            (d_fake, d_bd, d_conv), "d_nonprompt_" + region
        )
        if _safe_integral(d_all) <= 0:
            print("Skipping weighted FR for %s: empty denominator" % region)
            continue

        avr_fr_mc, err_mc = _weighted_efficiency(n_all, d_all)
        nonprompt_fr, nonprompt_err = _weighted_efficiency(
            n_nonprompt, d_nonprompt
        )
        avr_fake, err_fake = _weighted_efficiency(n_fake, d_fake)
        avr_bd, err_bd = _weighted_efficiency(n_bd, d_bd)
        avr_conv, err_conv = _weighted_efficiency(n_conv, d_conv)
        avr_prom, err_prom = _weighted_efficiency(n_prom, d_prom)

        titles = [
            "Average FR",
            "Z+1L nonprompt FR",
            "Avr Fake FR",
            "Avr BD fake FR",
            "Avr Conv FR",
            "Avr Prompt FR",
        ]
        values = [
            avr_fr_mc,
            nonprompt_fr,
            avr_fake,
            avr_bd,
            avr_conv,
            avr_prom,
        ]
        uerrs = [
            err_mc,
            nonprompt_err,
            err_fake,
            err_bd,
            err_conv,
            err_prom,
        ]
        lerrs = list(uerrs)

        origin_frs = {
            "Fake": avr_fake,
            "BDfake": avr_bd,
            "Conv": avr_conv,
            "Prompt": avr_prom,
        }
        origin_fr_variances = {
            "Fake": err_fake * err_fake,
            "BDfake": err_bd * err_bd,
            "Conv": err_conv * err_conv,
            "Prompt": err_prom * err_prom,
        }

        for section in CR_SECTIONS + list(CR_COMBINED_SECTIONS):
            if section not in categor_perc:
                continue
            perc = categor_perc[section]
            covariance = composition_cov[section]

            # Inclusive composition-weighted FR (Fake + BD + Conv + Prompt)
            avr_cat, avr_cat_error = _composition_weighted_fr(
                perc,
                covariance,
                origin_frs,
                origin_fr_variances,
                ALL_ORIGIN_KEYS,
            )
            titles.append(section)
            values.append(avr_cat)
            uerrs.append(avr_cat_error)
            lerrs.append(avr_cat_error)

            # Same estimate with Prompt removed and Fake/BD/Conv renormalized
            avr_np, avr_np_error = _composition_weighted_fr(
                perc,
                covariance,
                origin_frs,
                origin_fr_variances,
                NONPROMPT_ORIGIN_KEYS,
            )
            titles.append(section + " (no prompt)")
            values.append(avr_np)
            uerrs.append(avr_np_error)
            lerrs.append(avr_np_error)

        valid_sections = (
            ("2P2F4e", "2P2F2mu2e", "3P1F4e", "3P1F2mu2e",
             "2P2Fe", "3P1Fe")
            if region in ("EB", "EE")
            else
            ("2P2F4mu", "2P2F2e2mu", "3P1F4mu", "3P1F2e2mu",
             "2P2Fmu", "3P1Fmu")
        )
        valid_titles = set(valid_sections)
        for section in valid_sections:
            valid_titles.add(section + " (no prompt)")
        table_indices = [
            index
            for index, title in enumerate(titles)
            if title.startswith("Avr ")
            or title == "Average FR"
            or title == "Z+1L nonprompt FR"
            or title in valid_titles
        ]
        write_latex_ul_err(
            tex_handle,
            "FR " + region,
            [titles[index] for index in table_indices],
            [values[index] for index in table_indices],
            [uerrs[index] for index in table_indices],
            [lerrs[index] for index in table_indices],
            value_label="Fake rate",
        )
        _draw_composition_fr_summary(
            region,
            titles,
            values,
            uerrs,
            lerrs,
        )
        region_entries[region] = {
            title: (value, uerr, lerr)
            for title, value, uerr, lerr in zip(
                titles, values, uerrs, lerrs
            )
            if not title.endswith(" + upper err")
        }

    if len(region_entries) == len(FR_REGIONS):
        _draw_flavour_combined_fr_summary(region_entries)
    else:
        print(
            "Skipping flavour-combined FR summary: missing one or more regions"
        )


def compute_htozz_composition_summary(files, categor_perc, tex_handle):
    """Build H->ZZ closure plots and diagnostic composition ratios."""
    files_no_wz = [
        tfile
        for tfile in files
        if not os.path.basename(tfile.GetName()).endswith("_WZ.root")
    ]
    if len(files_no_wz) != len(files) - 1:
        raise RuntimeError(
            "Expected exactly one WZ Hist_ID file, found %d removable files"
            % (len(files) - len(files_no_wz))
        )

    data_file = _open_required_file(HTOZZ_DATA_FR_FILE)
    wz_removed_file = _open_required_file(HTOZZ_WZ_REMOVED_FR_FILE)
    summaries = {}

    root_path = os.path.join(OUT_DIR, "FR_composition_summary.root")
    root_file = ROOT.TFile(root_path, "RECREATE")

    try:
        for region in FR_REGIONS:
            data_name = DATA_FR_NAMES[region]
            # Stage 1 already stores both ratios; do not repeat WZ subtraction.
            data_raw_fr = _clone_required(
                data_file, data_name, "FR_data_raw_" + region
            )
            data_wz_removed_fr = _clone_required(
                wz_removed_file,
                data_name,
                "FR_data_wz_removed_" + region,
            )

            origins_all = _mc_origin_histograms(
                files, region, "mc_all_processes"
            )
            origins_no_wz = _mc_origin_histograms(
                files_no_wz, region, "mc_no_wz"
            )
            mc_all_fr, mc_all_integrated = _combined_origin_fr(
                origins_all,
                ("Fake", "BDfake", "Conv", "Prompt"),
                "FR_mc_all_" + region,
            )
            mc_all_no_wz_fr, mc_all_no_wz_integrated = _combined_origin_fr(
                origins_no_wz,
                ("Fake", "BDfake", "Conv", "Prompt"),
                "FR_mc_all_no_wz_" + region,
            )
            mc_np_fr, mc_np_integrated = _combined_origin_fr(
                origins_all,
                ("Fake", "BDfake", "Conv"),
                "FR_mc_nonprompt_" + region,
            )
            mc_np_no_wz_fr, mc_np_no_wz_integrated = _combined_origin_fr(
                origins_no_wz,
                ("Fake", "BDfake", "Conv"),
                "FR_mc_nonprompt_no_wz_" + region,
            )

            region_summary = {
                "z1l_mc_all": mc_all_integrated,
                "z1l_mc_all_no_wz": mc_all_no_wz_integrated,
                "z1l_mc_nonprompt": mc_np_integrated,
                "z1l_mc_nonprompt_no_wz": mc_np_no_wz_integrated,
            }

            baseline = mc_all_no_wz_integrated
            if baseline <= 0:
                raise RuntimeError(
                    "Empty H->ZZ all-origin MC Z+1L baseline for %s" % region
                )
            for final_state, (section_2p2f, section_3p1f) in (
                CR_SUMMARY_PAIRS.items()
            ):
                if (
                    section_2p2f not in categor_perc
                    or section_3p1f not in categor_perc
                ):
                    continue
                f_2p2f = _effective_cr_fr(
                    origins_no_wz, categor_perc[section_2p2f]
                )
                f_3p1f = _effective_cr_fr(
                    origins_no_wz, categor_perc[section_3p1f]
                )
                r_2p2f = f_2p2f / baseline
                r_3p1f = f_3p1f / baseline
                delta = max(abs(r_2p2f - 1.0), abs(r_3p1f - 1.0))
                region_summary[final_state] = {
                    "f_2p2f": f_2p2f,
                    "f_3p1f": f_3p1f,
                    "r_2p2f": r_2p2f,
                    "r_3p1f": r_3p1f,
                    "delta_comp": delta,
                }
            summaries[region] = region_summary

            closure_curves = {
                "data_raw": data_raw_fr,
                "data_wz_removed": data_wz_removed_fr,
                "mc_all": mc_all_fr,
                "mc_all_no_wz": mc_all_no_wz_fr,
                "mc_nonprompt": mc_np_fr,
                "mc_nonprompt_no_wz": mc_np_no_wz_fr,
            }
            _draw_z1l_closure(region, closure_curves)

            root_file.cd()
            for hist in closure_curves.values():
                hist.Write()
            summary_hist = ROOT.TH1D(
                "FR_composition_values_" + region,
                "FR composition summary " + region,
                7,
                0.5,
                7.5,
            )
            flavour_key = "electron" if region in ("EB", "EE") else "muon"
            scalar_values = (
                ("MC all", region_summary["z1l_mc_all"]),
                ("MC all no WZ", region_summary["z1l_mc_all_no_wz"]),
                ("MC nonprompt", region_summary["z1l_mc_nonprompt"]),
                (
                    "MC nonprompt no WZ",
                    region_summary["z1l_mc_nonprompt_no_wz"],
                ),
                ("2P2F", region_summary[flavour_key]["f_2p2f"]),
                ("3P1F", region_summary[flavour_key]["f_3p1f"]),
                (
                    "Diagnostic delta",
                    region_summary[flavour_key]["delta_comp"],
                ),
            )
            for ibin, (label, value) in enumerate(scalar_values, 1):
                summary_hist.GetXaxis().SetBinLabel(ibin, label)
                summary_hist.SetBinContent(ibin, value)
                # These bins mix different scalar diagnostics and do not have
                # a common statistical-error definition. Explicitly suppress
                # TH1's default sqrt(content) error.
                summary_hist.SetBinError(ibin, 0.0)
            summary_hist.Write()
    finally:
        root_file.Close()
        data_file.Close()
        wz_removed_file.Close()

    _write_composition_summary_latex(tex_handle, summaries)

    text_path = os.path.join(OUT_DIR, "FR_composition_summary.txt")
    with open(text_path, "w") as text_handle:
        text_handle.write(
            "CR truth composition counts the actual failing Z2 lepton(s).\n"
            "R and delta are diagnostic comparisons to the all-origin "
            "(Fake/BD/Conv/Prompt), no-WZ Z+1L MC baseline; they are not "
            "an adopted systematic. CR f is the effective FR F/(1+F) with "
            "F = sum_i c_i f_i/(1-f_i) over failing-lepton fractions c_i.\n"
        )
        for region in FR_REGIONS:
            text_handle.write("\n=== %s (%s) ===\n" % (
                region, FR_REGION_LABELS[region]
            ))
            text_handle.write(
                "MC Z+1L: all=%.6f all(no WZ)=%.6f "
                "nonprompt=%.6f nonprompt(no WZ)=%.6f\n"
                % (
                    summaries[region]["z1l_mc_all"],
                    summaries[region]["z1l_mc_all_no_wz"],
                    summaries[region]["z1l_mc_nonprompt"],
                    summaries[region]["z1l_mc_nonprompt_no_wz"],
                )
            )
            final_states = (
                ("4e", "2mu2e", "electron")
                if region in ("EB", "EE")
                else ("4mu", "2e2mu", "muon")
            )
            for final_state in final_states:
                entry = summaries[region][final_state]
                text_handle.write(
                    "%-9s f2P2F=%.6f R2=%.4f "
                    "f3P1F=%.6f R3=%.4f diagnostic_delta=%.4f\n"
                    % (
                        final_state,
                        entry["f_2p2f"],
                        entry["r_2p2f"],
                        entry["f_3p1f"],
                        entry["r_3p1f"],
                        entry["delta_comp"],
                    )
                )

    print("\n=== H->ZZ diagnostic composition-ratio summary ===")
    for region in FR_REGIONS:
        flavour_key = "electron" if region in ("EB", "EE") else "muon"
        entry = summaries[region][flavour_key]
        print(
            "  %s (%s CR): Z1L(all,noWZ)=%.5f, 2P2F=%.5f (R=%.3f), "
            "3P1F=%.5f (R=%.3f), diagnostic delta=%.1f%%"
            % (
                region,
                flavour_key,
                summaries[region]["z1l_mc_all_no_wz"],
                entry["f_2p2f"],
                entry["r_2p2f"],
                entry["f_3p1f"],
                entry["r_3p1f"],
                100.0 * entry["delta_comp"],
            )
        )
    print("Wrote %s" % root_path)
    print("Wrote %s" % text_path)
    return summaries


def _z1l_effective(files, origin_keys, region, tag):
    """Integrated and per-pT Z+1L FR from the given samples and origin keys."""
    origins = _mc_origin_histograms(files, region, tag)
    frs = {}
    integrated = {}
    variances = {}
    for origin in ALL_ORIGIN_KEYS:
        frs[origin] = _ratio_hist(
            origins[origin]["n"],
            origins[origin]["d"],
            "FR_%s_%s_%s" % (tag, origin, region),
        )
        value, error = _weighted_efficiency(
            origins[origin]["n"], origins[origin]["d"]
        )
        integrated[origin] = value
        variances[origin] = error * error
    z1l_fr, _unused = _combined_origin_fr(
        origins, origin_keys, "FR_Z1L_%s_%s" % (tag, region)
    )
    z1l_n = _add_histograms(
        [origins[origin]["n"] for origin in origin_keys],
        "n_Z1L_%s_%s" % (tag, region),
    )
    z1l_d = _add_histograms(
        [origins[origin]["d"] for origin in origin_keys],
        "d_Z1L_%s_%s" % (tag, region),
    )
    z1l_value, z1l_error = _weighted_efficiency(z1l_n, z1l_d)
    F_z1l = _fr_to_odds(z1l_value)
    F_z1l_error = _odds_error(z1l_value, z1l_error)
    if F_z1l is None or F_z1l <= 0.0:
        raise RuntimeError("Invalid Z+1L FR for %s (%s)" % (region, tag))
    return {
        "origins": origins,
        "frs": frs,
        "integrated": integrated,
        "variances": variances,
        "z1l_fr": z1l_fr,
        "f": z1l_value,
        "f_err": z1l_error,
        "F": F_z1l,
        "F_err": F_z1l_error,
    }


TRANSFER_VARIANTS = (
    {
        "cr": "2P2F_noprompt",
        "z1l_tag": "noprompt",
        "drop_wz_z1l": False,
        "origin_keys": NONPROMPT_ORIGIN_KEYS,
        "source_treatment": "Z1L all MC, HF+LF+conv",
        "cr_treatment": "2P2F only, prompt dropped",
        "z1l_label": "Z+1L no prompt",
        "cr_label": "2P2F mix, no prompt",
    },
    {
        "cr": "2P2F_leftover",
        "z1l_tag": "leftover",
        "drop_wz_z1l": True,
        "origin_keys": ALL_ORIGIN_KEYS,
        "source_treatment": "Z1L no WZ, all origins (leftover prompt)",
        "cr_treatment": "2P2F only, leftover prompt included",
        "z1l_label": "Z+1L leftover prompt (WZ out)",
        "cr_label": "2P2F mix, leftover prompt",
    },
)


def compute_vbs_cr_transfer(origin_files, cr_files, fr_tag):
    """Build 2P2F-only CR-mix fake rates, noprompt vs leftover-prompt.

    noprompt: Z+1L and 2P2F mixes use Fake+BD+Conv only (full MC stack).
    leftover: Z+1L drops the WZ sample and keeps leftover prompt; 2P2F mix
    includes prompt. ZZ subtraction stays on the 3P1F yield, not this mix.
    """
    tag_label = fr_tag or "nominal"
    source_files = origin_files

    data_path = _wz_removed_data_path(fr_tag)
    data_file = _open_required_file(data_path)
    root_path = os.path.join(OUT_DIR, "FR_CR_transfer.root")
    csv_path = os.path.join(OUT_DIR, "FR_CR_transfer.csv")
    text_path = os.path.join(OUT_DIR, "FR_CR_transfer_summary.txt")
    root_file = ROOT.TFile(root_path, "RECREATE")
    if not root_file or root_file.IsZombie():
        data_file.Close()
        raise RuntimeError("Cannot create ROOT output: %s" % root_path)

    fieldnames = [
        "analysis",
        "fr_tag",
        "cr",
        "region",
        "bin",
        "pt_low",
        "pt_high",
        "source_treatment",
        "cr_treatment",
        "c_fake",
        "c_bdfake",
        "c_conv",
        "c_prompt",
        "f_data",
        "f_band_down",
        "f_band_up",
        "f_z1l",
        "f_z1l_stat",
        "f_cr",
        "f_cr_stat",
        "F_z1l",
        "F_z1l_stat",
        "F_cr",
        "F_cr_stat",
        "kappa",
        "kappa_stat_approx",
        "delta",
    ]
    rows = []
    summaries = {}
    report = [
        "%s CR/Z+1L fake-rate transfer, 2P2F-only [%s]"
        % (ANALYSIS_MODE, tag_label),
        "noprompt: Z+1L and 2P2F mixes use HF+LF+conv (full MC stack).",
        "leftover: Z+1L drops WZ and keeps leftover prompt; 2P2F mix includes prompt.",
        "3P1F is not used in the FR mix. ZZ is subtracted on the 3P1F yield.",
        "CR fractions are inclusive in m4l and constant across pT bins.",
        "CR-mix F = sum_i c_i f_i/(1-f_i) over failing-lepton fractions c_i;",
        "f_cr is the effective FR F/(1+F) so that f_cr/(1-f_cr) = F_cr.",
        "",
    ]
    leftover_files = _files_without_wz(source_files)
    if not leftover_files:
        raise RuntimeError("No non-WZ Hist_ID files for leftover-prompt Z+1L")

    def _append_rows(
        cr_kind, region, z1l, mix, data_fr, source_treatment, cr_treatment
    ):
        composition = mix["composition"]
        z1l_fr = z1l["z1l_fr"]
        cr_fr = mix["fr"]
        axis = z1l_fr.GetXaxis()
        for ibin in range(1, z1l_fr.GetNbinsX() + 1):
            f_z1l_bin = z1l_fr.GetBinContent(ibin)
            f_cr_bin = cr_fr.GetBinContent(ibin)
            F_z1l_bin = _fr_to_odds(f_z1l_bin)
            F_cr_bin = _fr_to_odds(f_cr_bin)
            rows.append(
                {
                    "analysis": ANALYSIS_MODE,
                    "fr_tag": tag_label,
                    "cr": cr_kind,
                    "region": region,
                    "bin": ibin,
                    "pt_low": axis.GetBinLowEdge(ibin),
                    "pt_high": axis.GetBinUpEdge(ibin),
                    "source_treatment": source_treatment,
                    "cr_treatment": cr_treatment,
                    "c_fake": composition["Fake"] / 100.0,
                    "c_bdfake": composition["BDfake"] / 100.0,
                    "c_conv": composition["Conv"] / 100.0,
                    "c_prompt": composition["Prompt"] / 100.0,
                    "f_data": data_fr.GetBinContent(ibin),
                    "f_band_down": "",
                    "f_band_up": "",
                    "f_z1l": f_z1l_bin,
                    "f_z1l_stat": z1l_fr.GetBinError(ibin),
                    "f_cr": f_cr_bin,
                    "f_cr_stat": cr_fr.GetBinError(ibin),
                    "F_z1l": F_z1l_bin if F_z1l_bin is not None else "",
                    "F_z1l_stat": _odds_error(
                        f_z1l_bin, z1l_fr.GetBinError(ibin)
                    ),
                    "F_cr": F_cr_bin if F_cr_bin is not None else "",
                    "F_cr_stat": _odds_error(
                        f_cr_bin, cr_fr.GetBinError(ibin)
                    ),
                    "kappa": "",
                    "kappa_stat_approx": "",
                    "delta": "",
                }
            )
        rows.append(
            {
                "analysis": ANALYSIS_MODE,
                "fr_tag": tag_label,
                "cr": cr_kind,
                "region": region,
                "bin": "inclusive",
                "pt_low": "",
                "pt_high": "",
                "source_treatment": source_treatment,
                "cr_treatment": cr_treatment,
                "c_fake": composition["Fake"] / 100.0,
                "c_bdfake": composition["BDfake"] / 100.0,
                "c_conv": composition["Conv"] / 100.0,
                "c_prompt": composition["Prompt"] / 100.0,
                "f_data": "",
                "f_band_down": "",
                "f_band_up": "",
                "f_z1l": z1l["f"],
                "f_z1l_stat": z1l["f_err"],
                "f_cr": mix["f"],
                "f_cr_stat": mix["f_err"],
                "F_z1l": z1l["F"],
                "F_z1l_stat": z1l["F_err"],
                "F_cr": mix["F"],
                "F_cr_stat": mix["F_err"],
                "kappa": "",
                "kappa_stat_approx": "",
                "delta": "",
            }
        )

    try:
        for region in FR_REGIONS:
            data_fr = _clone_required(
                data_file,
                DATA_FR_NAMES[region],
                "FR_data_WZremoved_%s_%s" % (tag_label, region),
            )
            region_result = {"data_fr": data_fr}
            report.append(
                "=== %s (%s) ===" % (region, FR_REGION_LABELS[region])
            )
            draw_curves = []
            colors = {
                "2P2F_noprompt": (ROOT.kBlack, ROOT.kBlue + 1, 20, 21),
                "2P2F_leftover": (ROOT.kGreen + 2, ROOT.kRed + 1, 22, 23),
            }
            for variant in TRANSFER_VARIANTS:
                z1l_files = leftover_files if variant["drop_wz_z1l"] else source_files
                z1l = _z1l_effective(
                    z1l_files,
                    variant["origin_keys"],
                    region,
                    "transfer_%s_%s" % (tag_label, variant["z1l_tag"]),
                )
                mix = _build_cr_mix(
                    cr_files,
                    z1l["frs"],
                    z1l["integrated"],
                    z1l["variances"],
                    variant["cr"],
                    region,
                    tag_label,
                    origin_keys=variant["origin_keys"],
                    cr_treatment=variant["cr_treatment"],
                )
                composition = mix["composition"]
                region_result[variant["cr"]] = {
                    "z1l": z1l,
                    "mix": mix,
                    "F_z1l": z1l["F"],
                    "F_cr": mix["F"],
                }
                report.append(
                    "%s Z+1L: f=%.6f +/- %.6f, F=%.6f +/- %.6f"
                    % (
                        variant["z1l_tag"],
                        z1l["f"],
                        z1l["f_err"],
                        z1l["F"],
                        z1l["F_err"],
                    )
                )
                report.append(
                    "%s 2P2F mix Fake=%.4f BDfake=%.4f Conv=%.4f Prompt=%.4f"
                    % (
                        variant["z1l_tag"],
                        composition["Fake"] / 100.0,
                        composition["BDfake"] / 100.0,
                        composition["Conv"] / 100.0,
                        composition["Prompt"] / 100.0,
                    )
                )
                report.append(
                    "%s 2P2F: f=%.6f +/- %.6f, F=%.6f +/- %.6f"
                    % (
                        variant["z1l_tag"],
                        mix["f"],
                        mix["f_err"],
                        mix["F"],
                        mix["F_err"],
                    )
                )
                _append_rows(
                    variant["cr"],
                    region,
                    z1l,
                    mix,
                    data_fr,
                    variant["source_treatment"],
                    variant["cr_treatment"],
                )
                z_color, c_color, z_m, c_m = colors[variant["cr"]]
                draw_curves.append(
                    (z1l["z1l_fr"], z_color, z_m, variant["z1l_label"])
                )
                draw_curves.append(
                    (mix["fr"], c_color, c_m, variant["cr_label"])
                )

            report.append("")
            summaries[region] = region_result
            _draw_vbs_transfer_curves(region, draw_curves)
            root_file.cd()
            data_fr.Write()
            for variant in TRANSFER_VARIANTS:
                entry = region_result[variant["cr"]]
                entry["z1l"]["z1l_fr"].Write()
                entry["mix"]["fr"].Write()
    finally:
        root_file.Close()
        data_file.Close()

    with open(csv_path, "w") as csv_handle:
        writer = csv.DictWriter(csv_handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    with open(text_path, "w") as text_handle:
        text_handle.write("\n".join(report) + "\n")

    print("\n=== %s 2P2F-mix fake-rate summary [%s] ===" % (ANALYSIS_MODE, tag_label))
    for region in FR_REGIONS:
        nop = summaries[region]["2P2F_noprompt"]
        left = summaries[region]["2P2F_leftover"]
        print(
            "  %s: noprompt F_Z1L=%.5f F_2P2F=%.5f | leftover F_Z1L=%.5f F_2P2F=%.5f"
            % (
                region,
                nop["F_z1l"],
                nop["F_cr"],
                left["F_z1l"],
                left["F_cr"],
            )
        )
    print("Wrote %s" % root_path)
    print("Wrote %s" % csv_path)
    print("Wrote %s" % text_path)
    return summaries


def _remove_obsolete_plot_files():
    """Remove plot formats superseded by the current composition summary."""
    prefixes = (
        "FR_composition_envelope_",
        "FR_average_",
        "FR_CR_transfer_ratio_",
    )
    for prefix in prefixes:
        for region in FR_REGIONS:
            for extension in (".pdf", ".png", ".C", ".root"):
                path = os.path.join(OUT_DIR, prefix + region + extension)
                if os.path.exists(path):
                    os.remove(path)
                    print("Removed obsolete output %s" % path)
    for cr_kind in ("2P2F", "3P1F"):
        for region in FR_REGIONS:
            for extension in (".pdf", ".png", ".C", ".root"):
                path = os.path.join(
                    OUT_DIR, "FR_data_band_%s_%s%s" % (cr_kind, region, extension)
                )
                if os.path.exists(path):
                    os.remove(path)
                    print("Removed obsolete output %s" % path)


def main():
    global OUT_DIR
    ROOT.gROOT.SetBatch(True)

    mode = "VBS" if VBS else "HtoZZ"
    print("Mode: %s (VBS=%s)" % (mode, VBS))
    print("Input Hist_ID directory: %s" % INPUT_DIR)
    print("Base output directory: %s" % BASE_OUT_DIR)
    print("Reading MC Hist_ID files and computing FR composition uncertainties.")

    for fr_tag in FR_TAGS:
        tag_label = fr_tag or "nominal"
        OUT_DIR = os.path.join(BASE_OUT_DIR, tag_label)
        os.makedirs(OUT_DIR, exist_ok=True)
        _remove_obsolete_plot_files()
        print("\n=== FR composition mode: %s ===" % tag_label)
        print("Output directory: %s" % OUT_DIR)

        origin_files = _open_mc_files(fr_tag=fr_tag)
        # CR composition is defined by the analysis selection and is common to
        # nominal/ge1j/ge2j. Jet-tagged files contain only Z+1L origin hists.
        cr_files = origin_files if fr_tag is None else _open_mc_files()
        # Full MC stack, including WZ. Prompt is dropped later from the
        # CR-weighted transfer, not by removing a process from the sample list.
        analysis_files = origin_files
        composition_files = cr_files

        tex_path = os.path.join(
            OUT_DIR, "FR_composition_uncertainty_tables.tex"
        )
        with open(tex_path, "w") as tex_handle:
            tex_handle.write(
                "%% Auto-generated by fr_composition_uncertainty.py "
                "(%s, %s)\n\n" % (mode, tag_label)
            )
            # Full MC set: numerator is prompt(no WZ); denominator is prompt(all MC).
            compute_residual_prompt_sf(origin_files, tex_handle)
            compute_z1l_composition(analysis_files, tex_handle)
            categor_perc, composition_cov = compute_cr_composition(
                composition_files, tex_handle
            )
            compute_weighted_fr(
                analysis_files,
                categor_perc,
                composition_cov,
                tex_handle,
                cr_files=composition_files,
            )
            compute_vbs_cr_transfer(origin_files, cr_files, fr_tag)
            if not VBS and fr_tag is None:
                compute_htozz_composition_summary(
                    origin_files, categor_perc, tex_handle
                )

        for tfile in origin_files:
            tfile.Close()
        if cr_files is not origin_files:
            for tfile in cr_files:
                tfile.Close()

        print("Wrote LaTeX tables to %s" % tex_path)
        print("Wrote plots to %s/" % OUT_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
