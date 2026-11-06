import ROOT
from ROOT import TLorentzVector
import glob
import math
import os
import re
import numpy as np
from MC_composition import PartOrigin, PartOriginFromFlav


# VBS mode
# ----------------------------------------------------------------------------
# When VBS is True, the CR composition fills (2P2F/3P1F/4P m4l histograms) are
# additionally gated by the ZZjj-inclusive VBS selection (on-shell Zs, m4l>180,
# and two tagging jets with mjj>100). The 2P2F/3P1F classification is left
# untouched. IMPORTANT: the fake-rate measurement itself (the Z+1l block) is
# NOT gated, so the measured fake rates are identical to the nominal analysis.
# When VBS is True the m4l histograms are also re-booked to the VBS range.
# Keep this flag in sync with the VBS flag in estimateZX.py.
VBS = False

# VBS selection cut values (ZZjj inclusive, following the VBS analysis note).
VBS_MZ_LOW = 60.0
VBS_MZ_HIGH = 120.0
VBS_M4L_MIN = 180.0
VBS_JET_PT_MIN = 30.0
VBS_JET_ETA_MAX = 4.7
VBS_JET_LEP_DR_MIN = 0.4
VBS_MJJ_MIN = 100.0
VBS_NJET_MIN = 2
# tight jet ID + tight lepton veto (jetId bitmask). Set to None to disable.
VBS_JET_ID = 6

VBS_FR_JET_THRESHOLDS = [
    ("ge1j", 1),
    ("ge2j", 2),
]
_FR_RAW_KEYS = [
    "h1D_FRel_EB",
    "h1D_FRel_EE",
    "h1D_FRmu_EB",
    "h1D_FRmu_EE",
    "h1D_FRel_EB_d",
    "h1D_FRel_EE_d",
    "h1D_FRmu_EB_d",
    "h1D_FRmu_EE_d",
]


def _delta_r(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    # |delta phi| folded into [0, pi] via acos(cos(.)) to avoid 2pi wrapping
    dphi = math.acos(max(-1.0, min(1.0, math.cos(phi1 - phi2))))
    return math.sqrt(deta * deta + dphi * dphi)


def _chunk_dirname(chunk_path):
    path = chunk_path
    if path.startswith("root://"):
        path = re.sub(r"^root://[^/]+/*", "/", path)
    return os.path.basename(os.path.dirname(path))


def _cleaned_vbs_jets(event, lep_etas, lep_phis):
    """Jets passing VBS kinematics/ID and cleaned against the provided leptons."""
    cleaned = []  # (pt, TLorentzVector) of jets surviving cleaning
    for ij in range(event.nJet):
        jpt = event.Jet_pt[ij]
        if jpt <= VBS_JET_PT_MIN:
            continue
        jeta = event.Jet_eta[ij]
        if abs(jeta) >= VBS_JET_ETA_MAX:
            continue
        if VBS_JET_ID is not None and event.Jet_jetId[ij] != VBS_JET_ID:
            continue
        jphi = event.Jet_phi[ij]
        # jet-lepton cleaning against all four ZLL-candidate leptons
        too_close = False
        for leta, lphi in zip(lep_etas, lep_phis):
            if _delta_r(leta, lphi, jeta, jphi) < VBS_JET_LEP_DR_MIN:
                too_close = True
                break
        if too_close:
            continue
        jet = ROOT.TLorentzVector()
        jet.SetPtEtaPhiM(jpt, jeta, jphi, event.Jet_mass[ij])
        cleaned.append((jpt, jet))

    cleaned.sort(key=lambda x: x[0], reverse=True)
    return cleaned


def _pass_vbs_jets(event, lep_etas, lep_phis):
    """ZZjj-inclusive jet requirement for the VBS selection.

    Jets are cleaned by hand against the four ZLL-candidate leptons: the stored
    nCleanedJets*/JetLeadingIdx branches are cleaned only against ZZFullSel
    leptons and are therefore not valid for the CR fake leptons (which fail
    ZZFullSel by construction). Returns True if there are at least two cleaned
    jets passing (pT, |eta|, jet-ID) whose two highest-pT members have
    mjj > VBS_MJJ_MIN.
    """
    cleaned = _cleaned_vbs_jets(event, lep_etas, lep_phis)
    if len(cleaned) < VBS_NJET_MIN:
        return False

    mjj = (cleaned[0][1] + cleaned[1][1]).M()
    return mjj > VBS_MJJ_MIN


def _pass_vbs_selection(event, zllIdx, lep_etas, lep_phis, m4l):
    """Full ZZjj-inclusive VBS gate: on-shell Zs, m4l > 180, and tagging jets.

    Applied on top of the existing 2P2F/3P1F CR composition logic.
    """
    if not (VBS_MZ_LOW < event.ZLLCand_Z1mass[zllIdx] < VBS_MZ_HIGH):
        return False
    if not (VBS_MZ_LOW < event.ZLLCand_Z2mass[zllIdx] < VBS_MZ_HIGH):
        return False
    if not (m4l > VBS_M4L_MIN):
        return False
    return _pass_vbs_jets(event, lep_etas, lep_phis)


def _pass_vbs_selection_zzcand(event, zzIdx, lep_etas, lep_phis, m4l):
    """VBS gate for the 4P signal region (bestCandIdx / ZZCand)."""
    if not (VBS_MZ_LOW < event.ZZCand_Z1mass[zzIdx] < VBS_MZ_HIGH):
        return False
    if not (VBS_MZ_LOW < event.ZZCand_Z2mass[zzIdx] < VBS_MZ_HIGH):
        return False
    if not (m4l > VBS_M4L_MIN):
        return False
    return _pass_vbs_jets(event, lep_etas, lep_phis)


def _n_failed_leptons(event, lep_Hindex):
    """Count leptons failing the tight 4l ID used in CR composition."""
    n_failed = 0
    for k in range(4):
        lep_tight = event.Lepton_ZZFullSel[lep_Hindex[k]]
        lep_iso = getLeptonIso(event, lep_Hindex[k])
        idL = event.Lepton_pdgId[lep_Hindex[k]]
        if not (
            lep_tight
            and ((abs(idL) == 11) or (abs(idL) == 13 and lep_iso < 0.35))
        ):
            n_failed += 1
    return n_failed


def _fill_4p_m4l_histograms(h1D_m4l_4P, h1D_m4l_4P_4e, h1D_m4l_4P_4mu, h1D_m4l_4P_2e2mu, h1D_m4l_4P_2mu2e, m4l, weight, idL):
    """Fill inclusive and flavor 4P m4l histograms."""
    h1D_m4l_4P.Fill(m4l, weight)
    pdg_sum = abs(idL[0]) + abs(idL[1]) + abs(idL[2]) + abs(idL[3])
    if pdg_sum == 44:
        h1D_m4l_4P_4e.Fill(m4l, weight)
    elif pdg_sum == 52:
        h1D_m4l_4P_4mu.Fill(m4l, weight)
    elif pdg_sum == 48:
        if abs(idL[2]) + abs(idL[3]) == 26:
            h1D_m4l_4P_2e2mu.Fill(m4l, weight)
        if abs(idL[2]) + abs(idL[3]) == 22:
            h1D_m4l_4P_2mu2e.Fill(m4l, weight)


def m4l_plot_config():
    """Binning/range for m4l plots; keep in sync with analyzeZX / estimateZX."""
    if VBS:
        return {"nbins": 51, "xmin": 180.0, "xmax": 1200.0, "bin_width": 20.0}
    return {"nbins": 40, "xmin": 70.0, "xmax": 870.0, "bin_width": 20.0}


def _cr_plot_config(varName):
    return m4l_plot_config()


def _book_fr_histograms(tag, ptl_bins_el, ptl_bins_mu):
    hists = {
        "h1D_FRel_EB": ROOT.TH1D("h1D_FRel_EB_" + tag, "h1D_FRel_EB_" + tag, 6, ptl_bins_el),
        "h1D_FRel_EE": ROOT.TH1D("h1D_FRel_EE_" + tag, "h1D_FRel_EE_" + tag, 6, ptl_bins_el),
        "h1D_FRmu_EB": ROOT.TH1D("h1D_FRmu_EB_" + tag, "h1D_FRmu_EB_" + tag, 7, ptl_bins_mu),
        "h1D_FRmu_EE": ROOT.TH1D("h1D_FRmu_EE_" + tag, "h1D_FRmu_EE_" + tag, 7, ptl_bins_mu),
        "h1D_FRel_EB_d": ROOT.TH1D("h1D_FRel_EB_d_" + tag, "h1D_FRel_EB_d_" + tag, 6, ptl_bins_el),
        "h1D_FRel_EE_d": ROOT.TH1D("h1D_FRel_EE_d_" + tag, "h1D_FRel_EE_d_" + tag, 6, ptl_bins_el),
        "h1D_FRmu_EB_d": ROOT.TH1D("h1D_FRmu_EB_d_" + tag, "h1D_FRmu_EB_d_" + tag, 7, ptl_bins_mu),
        "h1D_FRmu_EE_d": ROOT.TH1D("h1D_FRmu_EE_d_" + tag, "h1D_FRmu_EE_d_" + tag, 7, ptl_bins_mu),
    }
    for hist in hists.values():
        hist.Sumw2()
    return hists


_ORIGIN_HIST_PREFIXES = ("Hist_conv", "Hist_prompt", "Hist_fakes", "Hist_BDfakes")


def _book_origin_histograms(tag, ptl_bins_el, ptl_bins_mu):
    """Book truth-origin numerator/denominator histograms for one FR jet tag."""
    hists = {}
    for prefix in _ORIGIN_HIST_PREFIXES:
        for index, region in enumerate(_ANALYZEZX_FR_ORDER):
            bins = ptl_bins_el if index < barrel_endcap_region.MB_n else ptl_bins_mu
            nbins = len(bins) - 1
            name = "%s%s_%s" % (prefix, region, tag)
            hists["%s%s" % (prefix, region)] = ROOT.TH1D(
                name, name, nbins, bins
            )
            hists["%s%s" % (prefix, region)].Sumw2()
    return hists


def _flatten_jet_fr_hists(jet_fr_hists):
    flat = {}
    for tag, _min_jets in VBS_FR_JET_THRESHOLDS:
        for key in _FR_RAW_KEYS:
            flat["%s_%s" % (key, tag)] = _detach_histogram(jet_fr_hists[tag][key])
    return flat


def _flatten_jet_origin_hists(jet_origin_hists):
    flat = {}
    for tag, _min_jets in VBS_FR_JET_THRESHOLDS:
        for key, hist in jet_origin_hists[tag].items():
            flat["%s_%s" % (key, tag)] = _detach_histogram(hist)
    return flat


def _fill_fr_hist_set(hists, idL3, etaL3, FillVar, weight, lep_tight, lep_iso):
    if ((abs(idL3) == 11) and (math.fabs(etaL3) < 1.497)):
        hists["h1D_FRel_EB_d"].Fill(FillVar, weight)
        if lep_tight:
            hists["h1D_FRel_EB"].Fill(FillVar, weight)

    if ((abs(idL3) == 11) and (math.fabs(etaL3) > 1.497)):
        hists["h1D_FRel_EE_d"].Fill(FillVar, weight)
        if lep_tight:
            hists["h1D_FRel_EE"].Fill(FillVar, weight)

    if ((abs(idL3) == 13) and (math.fabs(etaL3) < 1.2)):
        hists["h1D_FRmu_EB_d"].Fill(FillVar, weight)
        if lep_tight and (lep_iso < 0.35):
            hists["h1D_FRmu_EB"].Fill(FillVar, weight)

    if ((abs(idL3) == 13) and (math.fabs(etaL3) > 1.2)):
        hists["h1D_FRmu_EE_d"].Fill(FillVar, weight)
        if lep_tight and (lep_iso < 0.35):
            hists["h1D_FRmu_EE"].Fill(FillVar, weight)


def _fill_origin_hist_set(
    hists, isData, event, lep_Hindex, region_d, region_n,
    FillVar, weight, passes_numerator
):
    """Fill one jet-tagged truth-origin FR denominator and numerator."""
    region_d_name = _ANALYZEZX_FR_ORDER[region_d]
    region_n_name = _ANALYZEZX_FR_ORDER[region_n]
    PartOriginFromFlav(
        isData, event, lep_Hindex, FillVar, weight,
        hists["Hist_prompt" + region_d_name],
        hists["Hist_fakes" + region_d_name],
        hists["Hist_BDfakes" + region_d_name],
        hists["Hist_conv" + region_d_name],
        False,
    )
    if passes_numerator:
        PartOriginFromFlav(
            isData, event, lep_Hindex, FillVar, weight,
            hists["Hist_prompt" + region_n_name],
            hists["Hist_fakes" + region_n_name],
            hists["Hist_BDfakes" + region_n_name],
            hists["Hist_conv" + region_n_name],
            False,
        )


def _fill_cr_origins(
    cr_hists, cr_region, isData, event, lep_Hindex, fill_value, weight,
    lepton_positions
):
    """Fill one CR origin entry for each failing Z2 lepton.

    Z2 lepton 3 is not guaranteed to be the failing lepton: NanoAOD orders
    the two Z2 leptons by charge, not by tight-ID status. Therefore, 3P1F
    must classify whichever Z2 lepton actually fails, while 2P2F must count
    both failing leptons. Always using lepton 3 would sometimes classify the
    passing lepton in 3P1F and would discard half of the 2P2F information.

    Each failing lepton is filled with weight/N (N = number of failing Z2
    leptons in this event) so a 2P2F event with two failing leptons still
    contributes exactly one event's worth of weight in total -- matching the
    single per-event fill of h1D_m4l_2P2F/3P1F used elsewhere. Without this,
    2P2F origin histograms (used by CRMatchPlot and the composition
    uncertainty) were exactly double the true yield, since both failing
    leptons were filled at full weight. Splitting the weight leaves the
    Fake:BDfake:Conv *ratio* unchanged (it is a per-event normalization),
    so it does not affect the composition fractions -- only the absolute
    yield used for plotting/closure checks.
    """
    per_lepton_weight = weight / len(lepton_positions)
    for lepton_position in lepton_positions:
        PartOriginFromFlav(
            isData,
            event,
            lep_Hindex,
            fill_value,
            per_lepton_weight,
            cr_hists[type_of_fake.Prom][cr_region],
            cr_hists[type_of_fake.Fake][cr_region],
            cr_hists[type_of_fake.BDfake][cr_region],
            cr_hists[type_of_fake.Conv][cr_region],
            False,
            lepton_position=lepton_position,
        )


_ANALYZEZX_CR_SECTION = [
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
    "4P",
    "4P4e",
    "4P4mu",
    "4P2e2mu",
    "4P2mu2e",
]
_ANALYZEZX_FR_ORDER = ["EB_n", "EB_d", "EE_n", "EE_d", "MB_n", "MB_d", "ME_n", "ME_d"]

_ANALYZEZX_SINGLE_KEYS = [
    "h1D_FRel_EB",
    "h1D_FRel_EE",
    "h1D_FRmu_EB",
    "h1D_FRmu_EE",
    "h1D_FRel_EB_d",
    "h1D_FRel_EE_d",
    "h1D_FRmu_EB_d",
    "h1D_FRmu_EE_d",
    "h1D_m4l_2P2F",
    "h1D_m4l_3P1F",
    "h1D_m4l_4P",
    "h1D_m4l_2P2F_2e2mu",
    "h1D_m4l_3P1F_2e2mu",
    "h1D_m4l_2P2F_2mu2e",
    "h1D_m4l_3P1F_2mu2e",
    "h1D_m4l_2P2F_4e",
    "h1D_m4l_3P1F_4e",
    "h1D_m4l_2P2F_4mu",
    "h1D_m4l_3P1F_4mu",
    "h1D_m4l_4P_2e2mu",
    "h1D_m4l_4P_4e",
    "h1D_m4l_4P_4mu",
    "h1D_m4l_4P_2mu2e",
]
_ANALYZEZX_JET_FR_KEYS = [
    "%s_%s" % (key, tag)
    for tag, _min_jets in VBS_FR_JET_THRESHOLDS
    for key in _FR_RAW_KEYS
]
_ANALYZEZX_SINGLE_KEYS += _ANALYZEZX_JET_FR_KEYS
_ANALYZEZX_JET_ORIGIN_KEYS = [
    "%s%s_%s" % (prefix, region, tag)
    for tag, _min_jets in VBS_FR_JET_THRESHOLDS
    for prefix in _ORIGIN_HIST_PREFIXES
    for region in _ANALYZEZX_FR_ORDER
]
_ANALYZEZX_SINGLE_KEYS += _ANALYZEZX_JET_ORIGIN_KEYS

_ANALYZEZX_LIST_KEYS = [
    "CR_Prom",
    "CR_Fake",
    "CR_Conv",
    "CR_BDfake",
    "Hist_conv",
    "Hist_prompt",
    "Hist_fakes",
    "Hist_BDfakes",
]



def setCavasAndStyles(canvasName, c, stat):
    #setup canvas
    c = ROOT.TCanvas(canvasName,"myPlots",0,0,800,600)
    c.cd(1)
    c.SetLogy(0)
    ROOT.gStyle.SetOptStat(stat)
    ROOT.gStyle.SetPalette(1)

class barrel_endcap_region():
        EB_n = 0
        EB_d = 1
        EE_n = 2
        EE_d = 3
        MB_n = 4
        MB_d = 5
        ME_n = 6
        ME_d = 7

class type_of_fake():
    Prom = 0
    Fake = 1
    Conv = 2
    BDfake = 3

class CRregion():
    _2P2F = 0
    _2P2F_4e = 1
    _2P2F_4mu = 2
    _2P2F_2e2mu = 3
    _2P2F_2mu2e = 4
    _3P1F = 5
    _3P1F_4e = 6
    _3P1F_4mu = 7
    _3P1F_2e2mu = 8
    _3P1F_2mu2e = 9
    _4P = 10
    _4P_4e = 11
    _4P_4mu = 12
    _4P_2e2mu = 13
    _4P_2mu2e = 14

def setHistProperties(hist, lineWidth, lineStyle, lineColor, fillStyle, fillColor, xAxisTitle, yAxisTitle):
    if not(hist):
        return -1
    hist.SetLineWidth(lineWidth)
    hist.SetLineStyle(lineStyle)
    hist.SetLineColor(lineColor)
    # fill
    hist.SetFillStyle(fillStyle)
    hist.SetFillColor(fillColor)
    # divisions, offsets, sizes
    hist.GetXaxis().SetNdivisions(510)
    hist.GetYaxis().SetNdivisions(510)
    hist.GetXaxis().SetLabelSize(0.05)
    hist.GetYaxis().SetLabelSize(0.05)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleOffset(1.2)
    # titles
    if (xAxisTitle!="skip"):
        hist.GetXaxis().SetTitle(xAxisTitle)
    if (yAxisTitle!="skip"):
        hist.GetYaxis().SetTitle(yAxisTitle)
    # return
    return 0

#Reduntant now, after we check that the new version is good
#REMOVE THIS
def setNEvents(processFileName):
    lNEvents = 1
    if (processFileName=="Data"):
        lNEvents = 1
    if (processFileName=="DY10"):
        lNEvents = 37951928.0
    if (processFileName=="DY50"): 
        lNEvents = 99795992.0
    if (processFileName=="TT"): 
        lNEvents = 63667448.0
    if (processFileName=="WZ"):
        lNEvents = 6739437.0
    if (processFileName=="ZZ"):
        lNEvents = 97457264.0
    return lNEvents


def getLeptonIso(event, idx):
    """Return FSR-corrected relIso for merged Lepton index."""
    if idx < event.nElectron:
        return event.Electron_pfRelIso03FsrCorr[idx]
    return event.Muon_pfRelIso03FsrCorr[idx - event.nElectron]


def read_genEventSumw_from_events_tree(events_tree):
    """Sum genEventSumw from the Runs tree (per-file MC normalization)."""
    f = events_tree.GetCurrentFile()
    if f is None:
        raise RuntimeError("Cannot read Runs.genEventSumw: Events tree has no associated TFile")
    return read_genEventSumw_from_file(f)


def read_genEventSumw_from_file(tfile):
    """Sum genEventSumw from the Runs tree of one ROOT file."""
    runs = tfile.Get("Runs")
    if runs is None:
        raise RuntimeError("Runs tree not found in %s" % tfile.GetName())
    genEventSumw = 0.0
    for iRun in range(runs.GetEntries()):
        runs.GetEntry(iRun)
        genEventSumw += runs.genEventSumw
    return genEventSumw


def sum_genEventSumw_from_chunk_paths(chunk_paths):
    """Campaign-level genEventSumw for chunked MC (sum over all chunk Runs trees).

    Chunked ntuples store per-chunk Runs.genEventSumw. Histograms merged across
    chunks must all use this full-sample denominator; using each chunk's local
    genEventSumw inflates MC yields by ~the number of chunks.
    """
    total = 0.0
    n_chunks = len(chunk_paths)
    for i, fpath in enumerate(chunk_paths):
        in_file = ROOT.TFile.Open(fpath, "READ")
        if in_file is None or in_file.IsZombie():
            raise RuntimeError("Cannot open chunk file: %s" % fpath)
        total += read_genEventSumw_from_file(in_file)
        in_file.Close()
        if n_chunks > 50 and (i + 1) % 500 == 0:
            print(
                "--- campaign genEventSumw: %d/%d chunks, running sum = %.6g"
                % (i + 1, n_chunks, total)
            )
    return total


def _detach_histogram(hist):
    if hist is None:
        return None
    cloned = hist.Clone()
    cloned.SetDirectory(0)
    return cloned


def _analyzeZX_list_length(list_key):
    if list_key.startswith("CR_"):
        return len(_ANALYZEZX_CR_SECTION)
    return barrel_endcap_region.ME_d + 1


def _read_histogram_from_file(tfile, name):
    key = tfile.GetKey(name)
    if key is None:
        return None
    obj = key.ReadObj()
    if obj is None:
        return None
    if not obj.InheritsFrom(ROOT.TH1.Class()):
        raise RuntimeError(
            "Object %s in %s is %s, expected TH1"
            % (name, tfile.GetName(), obj.ClassName())
        )
    cloned = obj.Clone()
    cloned.SetDirectory(0)
    return cloned


def _pack_analyzeZX_histograms(
    h1D_FRel_EB,
    h1D_FRel_EE,
    h1D_FRmu_EB,
    h1D_FRmu_EE,
    h1D_FRel_EB_d,
    h1D_FRel_EE_d,
    h1D_FRmu_EB_d,
    h1D_FRmu_EE_d,
    h1D_m4l_2P2F,
    h1D_m4l_3P1F,
    h1D_m4l_4P,
    h1D_m4l_2P2F_2e2mu,
    h1D_m4l_3P1F_2e2mu,
    h1D_m4l_2P2F_2mu2e,
    h1D_m4l_3P1F_2mu2e,
    h1D_m4l_2P2F_4e,
    h1D_m4l_3P1F_4e,
    h1D_m4l_2P2F_4mu,
    h1D_m4l_3P1F_4mu,
    h1D_m4l_4P_2e2mu,
    h1D_m4l_4P_4e,
    h1D_m4l_4P_4mu,
    h1D_m4l_4P_2mu2e,
    CR_Prom,
    CR_Fake,
    CR_Conv,
    CR_BDfake,
    Hist_conv,
    Hist_prompt,
    Hist_fakes,
    Hist_BDfakes,
    jet_fr_hists=None,
    jet_origin_hists=None,
):
    hists = {
        "h1D_FRel_EB": _detach_histogram(h1D_FRel_EB),
        "h1D_FRel_EE": _detach_histogram(h1D_FRel_EE),
        "h1D_FRmu_EB": _detach_histogram(h1D_FRmu_EB),
        "h1D_FRmu_EE": _detach_histogram(h1D_FRmu_EE),
        "h1D_FRel_EB_d": _detach_histogram(h1D_FRel_EB_d),
        "h1D_FRel_EE_d": _detach_histogram(h1D_FRel_EE_d),
        "h1D_FRmu_EB_d": _detach_histogram(h1D_FRmu_EB_d),
        "h1D_FRmu_EE_d": _detach_histogram(h1D_FRmu_EE_d),
        "h1D_m4l_2P2F": _detach_histogram(h1D_m4l_2P2F),
        "h1D_m4l_3P1F": _detach_histogram(h1D_m4l_3P1F),
        "h1D_m4l_4P": _detach_histogram(h1D_m4l_4P),
        "h1D_m4l_2P2F_2e2mu": _detach_histogram(h1D_m4l_2P2F_2e2mu),
        "h1D_m4l_3P1F_2e2mu": _detach_histogram(h1D_m4l_3P1F_2e2mu),
        "h1D_m4l_2P2F_2mu2e": _detach_histogram(h1D_m4l_2P2F_2mu2e),
        "h1D_m4l_3P1F_2mu2e": _detach_histogram(h1D_m4l_3P1F_2mu2e),
        "h1D_m4l_2P2F_4e": _detach_histogram(h1D_m4l_2P2F_4e),
        "h1D_m4l_3P1F_4e": _detach_histogram(h1D_m4l_3P1F_4e),
        "h1D_m4l_2P2F_4mu": _detach_histogram(h1D_m4l_2P2F_4mu),
        "h1D_m4l_3P1F_4mu": _detach_histogram(h1D_m4l_3P1F_4mu),
        "h1D_m4l_4P_2e2mu": _detach_histogram(h1D_m4l_4P_2e2mu),
        "h1D_m4l_4P_4e": _detach_histogram(h1D_m4l_4P_4e),
        "h1D_m4l_4P_4mu": _detach_histogram(h1D_m4l_4P_4mu),
        "h1D_m4l_4P_2mu2e": _detach_histogram(h1D_m4l_4P_2mu2e),
        "CR_Prom": [_detach_histogram(h) for h in CR_Prom],
        "CR_Fake": [_detach_histogram(h) for h in CR_Fake],
        "CR_Conv": [_detach_histogram(h) for h in CR_Conv],
        "CR_BDfake": [_detach_histogram(h) for h in CR_BDfake],
        "Hist_conv": [_detach_histogram(h) for h in Hist_conv],
        "Hist_prompt": [_detach_histogram(h) for h in Hist_prompt],
        "Hist_fakes": [_detach_histogram(h) for h in Hist_fakes],
        "Hist_BDfakes": [_detach_histogram(h) for h in Hist_BDfakes],
    }
    if jet_fr_hists is not None:
        hists.update(_flatten_jet_fr_hists(jet_fr_hists))
    if jet_origin_hists is not None:
        hists.update(_flatten_jet_origin_hists(jet_origin_hists))
    return hists


def _merge_analyzeZX_histograms(accumulated, chunk_hists):
    for key in _ANALYZEZX_SINGLE_KEYS:
        accumulated[key].Add(chunk_hists[key])

    for list_key in _ANALYZEZX_LIST_KEYS:
        for i in range(len(accumulated[list_key])):
            accumulated[list_key][i].Add(chunk_hists[list_key][i])


def _is_data_nickname(Nickname):
    return "Data" in Nickname


def _hist_data_mc_label(isData):
    return "Data" if isData else "MC"


def _partial_id_path_for_main(main_path):
    base = os.path.basename(main_path)
    match = re.match(r"^Hist_(?:Data|MC)_([^_]+)_([^_]+)_(.+)\.root$", base)
    if match is None:
        raise RuntimeError("Cannot derive Hist_ID path from partial file: %s" % main_path)
    varName, Nickname, suffix = match.groups()
    return os.path.join(os.path.dirname(main_path), "Hist_ID_%s_%s_%s.root" % (varName, Nickname, suffix))


def _write_analyzeZX_partial(hists, Nickname, varName, isData, out_dir, suffix):
    """Write raw (pre-ratio) histograms for one Condor batch; suffix identifies the chunk range."""
    os.makedirs(out_dir, exist_ok=True)
    isData_string = _hist_data_mc_label(isData)
    main_path = os.path.join(
        out_dir,
        "Hist_%s_%s_%s_%s.root" % (isData_string, varName, Nickname, suffix),
    )
    id_path = os.path.join(
        out_dir,
        "Hist_ID_%s_%s_%s.root" % (varName, Nickname, suffix),
    )

    main_file = ROOT.TFile(main_path, "RECREATE")
    for key in _ANALYZEZX_SINGLE_KEYS:
        hist = hists[key]
        hist.SetName(key)
        hist.Write()
    main_file.Close()

    id_file = ROOT.TFile(id_path, "RECREATE")
    for list_key in _ANALYZEZX_LIST_KEYS:
        for i, hist in enumerate(hists[list_key]):
            hist.SetName("%s_%d" % (list_key, i))
            hist.Write()
    id_file.Close()

    return main_path, id_path


def _read_analyzeZX_partial(main_path):
    """Load raw histogram dict written by _write_analyzeZX_partial."""
    id_path = _partial_id_path_for_main(main_path)
    if not os.path.isfile(id_path):
        raise RuntimeError("Missing paired Hist_ID partial: %s" % id_path)

    hists = {}
    main_file = ROOT.TFile.Open(main_path, "READ")
    if main_file is None or main_file.IsZombie():
        raise RuntimeError("Cannot open partial Hist file: %s" % main_path)
    for key in _ANALYZEZX_SINGLE_KEYS:
        hist = _read_histogram_from_file(main_file, key)
        if hist is None:
            main_file.Close()
            raise RuntimeError("Missing histogram %s in %s" % (key, main_path))
        hists[key] = hist
    main_file.Close()

    id_file = ROOT.TFile.Open(id_path, "READ")
    if id_file is None or id_file.IsZombie():
        raise RuntimeError("Cannot open partial Hist_ID file: %s" % id_path)
    for list_key in _ANALYZEZX_LIST_KEYS:
        hists[list_key] = []
        for i in range(_analyzeZX_list_length(list_key)):
            hist = _read_histogram_from_file(id_file, "%s_%d" % (list_key, i))
            if hist is None:
                id_file.Close()
                raise RuntimeError(
                    "Missing histogram %s_%d in %s" % (list_key, i, id_path)
                )
            hists[list_key].append(hist)
    id_file.Close()
    return hists


def merge_analyzeZX_partials(
    partial_dir, Nickname, varName, isData, final_dir=".", partial_paths=None
):
    """Sum partial Hist_* files for one sample and write the standard stage-1 outputs.

    Pass partial_paths from the current campaign manifest. A bare glob of
    partial_dir would also pick up leftover files from older chunks-per-job
    groupings and double-count events.
    """
    if partial_paths is None:
        isData_string = _hist_data_mc_label(isData)
        pattern = os.path.join(
            partial_dir,
            "Hist_%s_%s_%s_*.root" % (isData_string, varName, Nickname),
        )
        partial_paths = sorted(glob.glob(pattern))
        if not partial_paths:
            raise RuntimeError("No partial files found: %s" % pattern)
    else:
        partial_paths = list(partial_paths)
        if not partial_paths:
            raise RuntimeError(
                "No manifest partials for %s / %s in %s"
                % (varName, Nickname, partial_dir)
            )
        missing = [path for path in partial_paths if not os.path.isfile(path)]
        if missing:
            raise RuntimeError(
                "Missing %d stage-1 partials for %s / %s, e.g. %s"
                % (len(missing), varName, Nickname, missing[0])
            )

    accumulated = None
    for main_path in partial_paths:
        chunk_hists = _read_analyzeZX_partial(main_path)
        if accumulated is None:
            accumulated = chunk_hists
        else:
            _merge_analyzeZX_histograms(accumulated, chunk_hists)

    _write_analyzeZX_histograms(accumulated, Nickname, varName, isData, out_dir=final_dir)
    return partial_paths


def _write_fr_histograms_to_file(out_path, fr_hists):
    out_file = ROOT.TFile(out_path, "RECREATE")

    ratios = {
        "Data_FRel_EB": fr_hists["h1D_FRel_EB"].Clone(),
        "Data_FRel_EE": fr_hists["h1D_FRel_EE"].Clone(),
        "Data_FRmu_EB": fr_hists["h1D_FRmu_EB"].Clone(),
        "Data_FRmu_EE": fr_hists["h1D_FRmu_EE"].Clone(),
    }
    ratios["Data_FRel_EB"].Divide(fr_hists["h1D_FRel_EB_d"])
    ratios["Data_FRel_EE"].Divide(fr_hists["h1D_FRel_EE_d"])
    ratios["Data_FRmu_EB"].Divide(fr_hists["h1D_FRmu_EB_d"])
    ratios["Data_FRmu_EE"].Divide(fr_hists["h1D_FRmu_EE_d"])

    for key, out_name in (
        ("h1D_FRel_EB", "Data_FRel_EB_n"),
        ("h1D_FRel_EE", "Data_FRel_EE_n"),
        ("h1D_FRmu_EB", "Data_FRmu_EB_n"),
        ("h1D_FRmu_EE", "Data_FRmu_EE_n"),
        ("h1D_FRel_EB_d", "Data_FRel_EB_d"),
        ("h1D_FRel_EE_d", "Data_FRel_EE_d"),
        ("h1D_FRmu_EB_d", "Data_FRmu_EB_d"),
        ("h1D_FRmu_EE_d", "Data_FRmu_EE_d"),
    ):
        hist = fr_hists[key].Clone(out_name)
        hist.SetName(out_name)
        hist.Write()

    for out_name, hist in ratios.items():
        hist.SetName(out_name)
        hist.Write()

    out_file.Close()


def _write_analyzeZX_histograms(hists, Nickname, varName, isData, out_dir="."):
    h1D_FRel_EB = hists["h1D_FRel_EB"]
    h1D_FRel_EE = hists["h1D_FRel_EE"]
    h1D_FRmu_EB = hists["h1D_FRmu_EB"]
    h1D_FRmu_EE = hists["h1D_FRmu_EE"]
    h1D_FRel_EB_d = hists["h1D_FRel_EB_d"]
    h1D_FRel_EE_d = hists["h1D_FRel_EE_d"]
    h1D_FRmu_EB_d = hists["h1D_FRmu_EB_d"]
    h1D_FRmu_EE_d = hists["h1D_FRmu_EE_d"]
    h1D_m4l_2P2F = hists["h1D_m4l_2P2F"]
    h1D_m4l_3P1F = hists["h1D_m4l_3P1F"]
    h1D_m4l_4P = hists["h1D_m4l_4P"]
    h1D_m4l_2P2F_2e2mu = hists["h1D_m4l_2P2F_2e2mu"]
    h1D_m4l_3P1F_2e2mu = hists["h1D_m4l_3P1F_2e2mu"]
    h1D_m4l_2P2F_2mu2e = hists["h1D_m4l_2P2F_2mu2e"]
    h1D_m4l_3P1F_2mu2e = hists["h1D_m4l_3P1F_2mu2e"]
    h1D_m4l_2P2F_4e = hists["h1D_m4l_2P2F_4e"]
    h1D_m4l_3P1F_4e = hists["h1D_m4l_3P1F_4e"]
    h1D_m4l_2P2F_4mu = hists["h1D_m4l_2P2F_4mu"]
    h1D_m4l_3P1F_4mu = hists["h1D_m4l_3P1F_4mu"]
    CR_Prom = hists["CR_Prom"]
    CR_Fake = hists["CR_Fake"]
    CR_Conv = hists["CR_Conv"]
    CR_BDfake = hists["CR_BDfake"]
    Hist_conv = hists["Hist_conv"]
    Hist_prompt = hists["Hist_prompt"]
    Hist_fakes = hists["Hist_fakes"]
    Hist_BDfakes = hists["Hist_BDfakes"]

    isData_string = _hist_data_mc_label(isData)
    os.makedirs(out_dir, exist_ok=True)

    SaveRootFile = ROOT.TFile(
        os.path.join(out_dir, "Hist_%s_%s_%s.root" % (isData_string, varName, Nickname)),
        "RECREATE",
    )

    h1D_FRel_EE_n = h1D_FRel_EE.Clone()
    h1D_FRmu_EE_n = h1D_FRmu_EE.Clone()
    h1D_FRel_EB_n = h1D_FRel_EB.Clone()
    h1D_FRmu_EB_n = h1D_FRmu_EB.Clone()

    h1D_FRel_EB_n.Divide(h1D_FRel_EB_d)
    h1D_FRel_EE_n.Divide(h1D_FRel_EE_d)
    h1D_FRmu_EB_n.Divide(h1D_FRmu_EB_d)
    h1D_FRmu_EE_n.Divide(h1D_FRmu_EE_d)

    h1D_FRel_EB_n.GetYaxis().SetRangeUser(0.01,0.35)
    h1D_FRel_EE_n.GetYaxis().SetRangeUser(0.01,0.35)
    h1D_FRmu_EB_n.GetYaxis().SetRangeUser(0.04,0.35)
    h1D_FRmu_EE_n.GetYaxis().SetRangeUser(0.04,0.35)

    h1D_FRel_EB.SetName("Data_FRel_EB_n")
    h1D_FRel_EB.Write()
    h1D_FRel_EE.SetName("Data_FRel_EE_n")
    h1D_FRel_EE.Write()
    h1D_FRmu_EB.SetName("Data_FRmu_EB_n")
    h1D_FRmu_EB.Write()
    h1D_FRmu_EE.SetName("Data_FRmu_EE_n")
    h1D_FRmu_EE.Write()

    h1D_FRel_EB_d.SetName("Data_FRel_EB_d")
    h1D_FRel_EB_d.Write()
    h1D_FRel_EE_d.SetName("Data_FRel_EE_d")
    h1D_FRel_EE_d.Write()
    h1D_FRmu_EB_d.SetName("Data_FRmu_EB_d")
    h1D_FRmu_EB_d.Write()
    h1D_FRmu_EE_d.SetName("Data_FRmu_EE_d")
    h1D_FRmu_EE_d.Write()

    h1D_FRel_EB_n.SetName("Data_FRel_EB")
    h1D_FRel_EB_n.Write()
    h1D_FRel_EE_n.SetName("Data_FRel_EE")
    h1D_FRel_EE_n.Write()
    h1D_FRmu_EB_n.SetName("Data_FRmu_EB")
    h1D_FRmu_EB_n.Write()
    h1D_FRmu_EE_n.SetName("Data_FRmu_EE")
    h1D_FRmu_EE_n.Write()

    h1D_m4l_2P2F.SetName("h1D_m4l_2P2F")
    h1D_m4l_2P2F.Write()
    h1D_m4l_2P2F_2e2mu.SetName("h1D_m4l_2P2F_2e2mu")
    h1D_m4l_2P2F_2e2mu.Write()
    h1D_m4l_2P2F_2mu2e.SetName("h1D_m4l_2P2F_2mu2e")
    h1D_m4l_2P2F_2mu2e.Write()
    h1D_m4l_2P2F_4e.SetName("h1D_m4l_2P2F_4e")
    h1D_m4l_2P2F_4e.Write()
    h1D_m4l_2P2F_4mu.SetName("h1D_m4l_2P2F_4mu")
    h1D_m4l_2P2F_4mu.Write()

    h1D_m4l_3P1F.SetName("h1D_m4l_3P1F")
    h1D_m4l_3P1F.Write()
    h1D_m4l_3P1F_2e2mu.SetName("h1D_m4l_3P1F_2e2mu")
    h1D_m4l_3P1F_2e2mu.Write()
    h1D_m4l_3P1F_2mu2e.SetName("h1D_m4l_3P1F_2mu2e")
    h1D_m4l_3P1F_2mu2e.Write()
    h1D_m4l_3P1F_4e.SetName("h1D_m4l_3P1F_4e")
    h1D_m4l_3P1F_4e.Write()
    h1D_m4l_3P1F_4mu.SetName("h1D_m4l_3P1F_4mu")
    h1D_m4l_3P1F_4mu.Write()

    h1D_m4l_4P.SetName("h1D_m4l_4P")
    h1D_m4l_4P.Write()
    h1D_m4l_4P_2e2mu = hists["h1D_m4l_4P_2e2mu"]
    h1D_m4l_4P_2e2mu.SetName("h1D_m4l_4P_2e2mu")
    h1D_m4l_4P_2e2mu.Write()
    h1D_m4l_4P_4e = hists["h1D_m4l_4P_4e"]
    h1D_m4l_4P_4e.SetName("h1D_m4l_4P_4e")
    h1D_m4l_4P_4e.Write()
    h1D_m4l_4P_4mu = hists["h1D_m4l_4P_4mu"]
    h1D_m4l_4P_4mu.SetName("h1D_m4l_4P_4mu")
    h1D_m4l_4P_4mu.Write()
    h1D_m4l_4P_2mu2e = hists["h1D_m4l_4P_2mu2e"]
    h1D_m4l_4P_2mu2e.SetName("h1D_m4l_4P_2mu2e")
    h1D_m4l_4P_2mu2e.Write()

    SaveRootFile.Close()

    if VBS and varName == "ptl3":
        for tag, _min_jets in VBS_FR_JET_THRESHOLDS:
            fr_hists = {}
            for key in _FR_RAW_KEYS:
                fr_hists[key] = hists["%s_%s" % (key, tag)]
            _write_fr_histograms_to_file(
                os.path.join(
                    out_dir,
                    "Hist_%s_%s_VBSFR_%s_%s.root"
                    % (isData_string, varName, tag, Nickname),
                ),
                fr_hists,
            )

    SaveRootFile_Types_of_fake = ROOT.TFile(
        os.path.join(out_dir, "Hist_ID_%s_%s.root" % (varName, Nickname)),
        "RECREATE",
    )

    for count in range(len(_ANALYZEZX_CR_SECTION)):
        section = _ANALYZEZX_CR_SECTION[count]
        CR_Prom[count].SetName("h1D_m4l_Prom_" + section)
        CR_Prom[count].Write()
        CR_Fake[count].SetName("h1D_m4l_Fake_" + section)
        CR_Fake[count].Write()
        CR_Conv[count].SetName("h1D_m4l_Conv_" + section)
        CR_Conv[count].Write()
        CR_BDfake[count].SetName("h1D_m4l_BDfake_" + section)
        CR_BDfake[count].Write()

    for count in range(barrel_endcap_region.ME_d + 1):
        region = _ANALYZEZX_FR_ORDER[count]
        Hist_conv[count].SetName("Hist_conv" + region)
        Hist_conv[count].Write()
        Hist_prompt[count].SetName("Hist_prompt" + region)
        Hist_prompt[count].Write()
        Hist_fakes[count].SetName("Hist_fakes" + region)
        Hist_fakes[count].Write()
        Hist_BDfakes[count].SetName("Hist_BDfakes" + region)
        Hist_BDfakes[count].Write()

    SaveRootFile_Types_of_fake.Close()

    if VBS and varName == "ptl3":
        for tag, _min_jets in VBS_FR_JET_THRESHOLDS:
            tagged_id_file = ROOT.TFile(
                os.path.join(
                    out_dir,
                    "Hist_ID_%s_VBSFR_%s_%s.root" % (varName, tag, Nickname),
                ),
                "RECREATE",
            )
            for prefix in _ORIGIN_HIST_PREFIXES:
                for region in _ANALYZEZX_FR_ORDER:
                    key = "%s%s_%s" % (prefix, region, tag)
                    hist = hists[key].Clone("%s%s" % (prefix, region))
                    hist.SetName("%s%s" % (prefix, region))
                    hist.Write()
            tagged_id_file.Close()


def analyzeZXChunkBatch(
    file_paths, Nickname, varName, out_dir, suffix, mc_genEventSumw=None
):
    """Process a Condor batch of ntuples and write one partial Hist_* pair."""
    if not file_paths:
        raise RuntimeError("No chunk files provided for sample %s" % Nickname)

    accumulated = None
    for i, fpath in enumerate(file_paths):
        chunk_genEventSumw = mc_genEventSumw
        if isinstance(mc_genEventSumw, dict):
            norm_key = _chunk_dirname(fpath).split("_Chunk", 1)[0]
            chunk_genEventSumw = mc_genEventSumw.get(norm_key)
            if chunk_genEventSumw is None:
                raise RuntimeError(
                    "Missing genEventSumw for component %s in sample %s"
                    % (norm_key, Nickname)
                )
        inFile = ROOT.TFile.Open(fpath, "READ")
        if inFile is None or inFile.IsZombie():
            raise RuntimeError("Cannot open chunk file: %s" % fpath)
        tree = inFile.Get("Events")
        if tree is None:
            inFile.Close()
            raise RuntimeError("Events tree not found in: %s" % fpath)
        print(
            "--- Chunk %d/%d: \"%s\" (%d events)"
            % (i + 1, len(file_paths), fpath, tree.GetEntries())
        )
        chunk_hists = analyzeZX(
            tree,
            Nickname,
            varName,
            write_output=False,
            mc_genEventSumw=chunk_genEventSumw,
        )
        inFile.Close()
        if accumulated is None:
            accumulated = chunk_hists
        else:
            _merge_analyzeZX_histograms(accumulated, chunk_hists)

    isData = "Data" in Nickname
    return _write_analyzeZX_partial(accumulated, Nickname, varName, isData, out_dir, suffix)


def analyzeZXFromChunks(file_paths, Nickname, varName="ptl3", mc_genEventSumw=None):
    if not file_paths:
        raise RuntimeError("No chunk files provided for sample %s" % Nickname)

    isData = "Data" in Nickname
    if not isData and mc_genEventSumw is None:
        print(
            "--- Computing campaign genEventSumw for %s (%d chunks)"
            % (Nickname, len(file_paths))
        )
        mc_genEventSumw = sum_genEventSumw_from_chunk_paths(file_paths)
        print("--- Campaign genEventSumw for %s = %.6g" % (Nickname, mc_genEventSumw))

    accumulated = None
    for i, fpath in enumerate(file_paths):
        inFile = ROOT.TFile.Open(fpath, "READ")
        if inFile is None or inFile.IsZombie():
            raise RuntimeError("Cannot open chunk file: %s" % fpath)
        tree = inFile.Get("Events")
        if tree is None:
            inFile.Close()
            raise RuntimeError("Events tree not found in: %s" % fpath)
        print(
            "--- Chunk %d/%d: \"%s\" (%d events)"
            % (i + 1, len(file_paths), fpath, tree.GetEntries())
        )
        chunk_hists = analyzeZX(
            tree,
            Nickname,
            varName,
            write_output=False,
            mc_genEventSumw=mc_genEventSumw,
        )
        inFile.Close()
        if accumulated is None:
            accumulated = chunk_hists
        else:
            _merge_analyzeZX_histograms(accumulated, chunk_hists)

    _write_analyzeZX_histograms(accumulated, Nickname, varName, isData)

def getMET_pt(event):
    """MET pt: PFMET_pt (Nano v13+/2024) or MET_pt (older ntuples)."""
    if hasattr(event, "PFMET_pt"):
        return event.PFMET_pt
    return event.MET_pt


def analyzeZX(
    fTemplateTree, Nickname, varName="ptl3", write_output=True, mc_genEventSumw=None
):
    print ("--- Initiating the analyzeZX procedure for file nicknamed as: "+ Nickname +".")

    LUMI_INT = 109820
    lineWidth = 2
    leg_xl = 0.50
    leg_xr = 0.90
    leg_yb = 0.72 
    leg_yt = 0.90
    

    varAxLabel = ""
    
    var_plotHigh = 120
    var_plotLow = 60
    var_nBins = 20
    varAxLabel = "m_{Z1}"
    
    if (varName=="mEt"):
        var_plotHigh = 50
        var_plotLow = 0
        var_nBins = 10
        varAxLabel = "E_{T,miss}"
    

    PtlBins = [5.,10.,20.,30.,40.,50.,80.]
    PtlBins = np.array(PtlBins)
    PtlBinsMu = [5.,7.,10.,20.,30.,40.,50.,80.]
    PtlBinsMu = np.array(PtlBinsMu)
    
    #initiate numerator and denominator histograms for FR computation

    binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.
    sUnit = "GeV"    
    
    
    if (varName=="ptl3"):
        h1D_dummy = ROOT.TH1D("dummy", "dummy",6, PtlBins)
        setHistProperties(h1D_dummy,1,1,1,0,0,varAxLabel,"Misidentification rate")
        #define FR numerator histograms

        h1D_FRel_EB   = ROOT.TH1D("h1D_FRel_EB","h1D_FRel_EB",6, PtlBins) 
        h1D_FRel_EB.Sumw2()
        h1D_FRel_EE   = ROOT.TH1D("h1D_FRel_EE","h1D_FRel_EE",6, PtlBins) 
        h1D_FRel_EE.Sumw2()
        
        h1D_FRmu_EB   = ROOT.TH1D("h1D_FRmu_EB","h1D_FRmu_EB",7, PtlBinsMu) 
        h1D_FRmu_EB.Sumw2()
        h1D_FRmu_EE   = ROOT.TH1D("h1D_FRmu_EE","h1D_FRmu_EE",7, PtlBinsMu) 
        h1D_FRmu_EE.Sumw2()
        
        # define FR denominator histograms
        h1D_FRel_EB_d   = ROOT.TH1D("h1D_FRel_EB_d","h1D_FRel_EB_d",6, PtlBins) 
        h1D_FRel_EB_d.Sumw2()
        h1D_FRel_EE_d   = ROOT.TH1D("h1D_FRel_EE_d","h1D_FRel_EE_d",6, PtlBins) 
        h1D_FRel_EE_d.Sumw2()
        
        h1D_FRmu_EB_d   = ROOT.TH1D("h1D_FRmu_EB_d","h1D_FRmu_EB_d",7, PtlBinsMu) 
        h1D_FRmu_EB_d.Sumw2()
        h1D_FRmu_EE_d   = ROOT.TH1D("h1D_FRmu_EE_d","h1D_FRmu_EE_d",7, PtlBinsMu) 
        h1D_FRmu_EE_d.Sumw2()    
    else:

        h1D_dummy = ROOT.TH1D("dummy", "dummy", var_nBins, var_plotLow, var_plotHigh)
        setHistProperties(h1D_dummy,1,1,1,0,0,varAxLabel,"Misidentification rate")

        # define FR numerator histograms
        h1D_FRel_EB   = ROOT.TH1D("h1D_FRel_EB","h1D_FRel_EB",var_nBins, var_plotLow, var_plotHigh)
        h1D_FRel_EB.Sumw2()
        h1D_FRel_EE   = ROOT.TH1D("h1D_FRel_EE","h1D_FRel_EE",var_nBins, var_plotLow, var_plotHigh)
        h1D_FRel_EE.Sumw2()
        h1D_FRmu_EB   = ROOT.TH1D("h1D_FRmu_EB","h1D_FRmu_EB",var_nBins, var_plotLow, var_plotHigh) 
        h1D_FRmu_EB.Sumw2()
        h1D_FRmu_EE   = ROOT.TH1D("h1D_FRmu_EE","h1D_FRmu_EE",var_nBins, var_plotLow, var_plotHigh) 
        h1D_FRmu_EE.Sumw2()
    
        # define FR denominator histograms
        h1D_FRel_EB_d   = ROOT.TH1D("h1D_FRel_EB_d","h1D_FRel_EB_d",var_nBins, var_plotLow, var_plotHigh) 
        h1D_FRel_EB_d.Sumw2()
        h1D_FRel_EE_d   = ROOT.TH1D("h1D_FRel_EE_d","h1D_FRel_EE_d",var_nBins, var_plotLow, var_plotHigh) 
        h1D_FRel_EE_d.Sumw2()
        h1D_FRmu_EB_d   = ROOT.TH1D("h1D_FRmu_EB_d","h1D_FRmu_EB_d",var_nBins, var_plotLow, var_plotHigh) 
        h1D_FRmu_EB_d.Sumw2()
        h1D_FRmu_EE_d   = ROOT.TH1D("h1D_FRmu_EE_d","h1D_FRmu_EE_d",var_nBins, var_plotLow, var_plotHigh) 
        h1D_FRmu_EE_d.Sumw2()


    jet_fr_hists = {}
    jet_origin_hists = {}
    for tag, _min_jets in VBS_FR_JET_THRESHOLDS:
        jet_fr_hists[tag] = _book_fr_histograms(tag, PtlBins, PtlBinsMu)
        jet_origin_hists[tag] = _book_origin_histograms(
            tag, PtlBins, PtlBinsMu
        )

    cr_cfg = _cr_plot_config(varName)
    CR_var_plotHigh = cr_cfg["xmax"]
    CR_var_plotLow = cr_cfg["xmin"]
    CR_var_nBins = cr_cfg["nbins"]

   # define CR histograms
    h1D_m4l_2P2F = ROOT.TH1D("h1D_m4l_2P2F","h1D_m4l_2P2F",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_2P2F.Sumw2()
    h1D_m4l_3P1F = ROOT.TH1D("h1D_m4l_3P1F","h1D_m4l_3P1F",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_3P1F.Sumw2()
    h1D_m4l_4P   = ROOT.TH1D("h1D_m4l_4P","h1D_m4l_4P",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_4P.Sumw2()
    
    # define CR histograms
    h1D_m4l_2P2F_2e2mu = ROOT.TH1D("h1D_m4l_2P2F_2e2mu","h1D_m4l_2P2F_2e2mu",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_2P2F_2e2mu.Sumw2()
    h1D_m4l_3P1F_2e2mu = ROOT.TH1D("h1D_m4l_3P1F_2e2mu","h1D_m4l_3P1F_2e2mu",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_3P1F_2e2mu.Sumw2()
    h1D_m4l_4P_2e2mu   = ROOT.TH1D("h1D_m4l_4P_2e2mu","h1D_m4l_4P_2e2mu",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_4P_2e2mu.Sumw2()
    
    # define CR histograms
    h1D_m4l_2P2F_2mu2e = ROOT.TH1D("h1D_m4l_2P2F_2mu2e","h1D_m4l_2P2F_2mu2e",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_2P2F_2mu2e.Sumw2()
    h1D_m4l_3P1F_2mu2e = ROOT.TH1D("h1D_m4l_3P1F_2mu2e","h1D_m4l_3P1F_2mu2e",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_3P1F_2mu2e.Sumw2()
    h1D_m4l_4P_2mu2e   = ROOT.TH1D("h1D_m4l_4P_2mu2e","h1D_m4l_4P_2mu2e",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_4P_2mu2e.Sumw2()
    
    # define CR histograms
    h1D_m4l_2P2F_4e = ROOT.TH1D("h1D_m4l_2P2F_4e","h1D_m4l_2P2F_4e",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_2P2F_4e.Sumw2()
    h1D_m4l_3P1F_4e = ROOT.TH1D("h1D_m4l_3P1F_4e","h1D_m4l_3P1F_4e",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_3P1F_4e.Sumw2()
    h1D_m4l_4P_4e   = ROOT.TH1D("h1D_m4l_4P_4e","h1D_m4l_4P_4e",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_4P_4e.Sumw2()
    
    # define CR histograms
    h1D_m4l_2P2F_4mu = ROOT.TH1D("h1D_m4l_2P2F_4mu","h1D_m4l_2P2F_4mu",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_2P2F_4mu.Sumw2()
    h1D_m4l_3P1F_4mu = ROOT.TH1D("h1D_m4l_3P1F_4mu","h1D_m4l_3P1F_4mu",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_3P1F_4mu.Sumw2()
    h1D_m4l_4P_4mu   = ROOT.TH1D("h1D_m4l_4P_4mu","h1D_m4l_4P_4mu",CR_var_nBins, CR_var_plotLow, CR_var_plotHigh) 
    h1D_m4l_4P_4mu.Sumw2()






    isData = ("Data" in Nickname)
    iEvt = -1
    nentries = fTemplateTree.GetEntries()
    lNEvents = setNEvents(Nickname)
     
    # Begin: Define histograms for uncertainty studies (separating fakes per type of fake) 

    order = ["EB_n","EB_d","EE_n","EE_d","MB_n","MB_d","ME_n","ME_d"]
    CRSection = _ANALYZEZX_CR_SECTION

    ID ={"Prom","Fake","Conv","BDfake"};

    ## Define CR histograms for Prompt/Fake/Converson and BDfakes    
    CR = [] 
    CR_Prom = []
    CR_Fake = []
    CR_Conv = []
    CR_BDfake = []

    for count in range(len(CRSection)):     
    
        CR_Prom.append(ROOT.TH1D("h1D_m4l_Prom_"+CRSection[count],"h1D_m4l_Prom_"+CRSection[count],CR_var_nBins, CR_var_plotLow, CR_var_plotHigh))
        CR_Prom[count].Sumw2()

        CR_Fake.append(ROOT.TH1D("h1D_m4l_Fake_"+CRSection[count],"h1D_m4l_Fake_"+CRSection[count],CR_var_nBins, CR_var_plotLow, CR_var_plotHigh))
        CR_Fake[count].Sumw2()

        CR_Conv.append(ROOT.TH1D("h1D_m4l_Conv_"+CRSection[count],"h1D_m4l_Conv_"+CRSection[count],CR_var_nBins, CR_var_plotLow, CR_var_plotHigh))
        CR_Conv[count].Sumw2()

        CR_BDfake.append(ROOT.TH1D("h1D_m4l_BDfake_"+CRSection[count],"h1D_m4l_BDfake_"+CRSection[count],CR_var_nBins, CR_var_plotLow, CR_var_plotHigh))
        CR_BDfake[count].Sumw2()

    CR.append(CR_Prom)
    CR.append(CR_Fake)
    CR.append(CR_Conv)
    CR.append(CR_BDfake)

    ## Define numerator/denominator histograms for Prompt/Fake/Converson and BDfakes in Z+1l region   

    Hist_conv = []
    Hist_prompt = []
    Hist_fakes = []
    Hist_BDfakes = []

    for count in range(barrel_endcap_region.ME_d+1):
    
        if (varName=="ptl3") and (count < barrel_endcap_region.MB_n):
            Hist_conv.append(ROOT.TH1D("Hist_conv"+order[count] ,"Hist_conv"+order[count],6,PtlBins))
            Hist_conv[count].Sumw2()

            Hist_prompt.append(ROOT.TH1D("Hist_prompt"+order[count],"Hist_prompt"+order[count],6,PtlBins))
            Hist_prompt[count].Sumw2()

            Hist_fakes.append(ROOT.TH1D("Hist_fakes"+order[count],"Hist_fakes"+order[count],6,PtlBins))
            Hist_fakes[count].Sumw2()

            Hist_BDfakes.append(ROOT.TH1D("Hist_BDfakes"+order[count],"Hist_BDfakes"+order[count],6,PtlBins))
            Hist_BDfakes[count].Sumw2()

        elif (varName=="ptl3") and (count >= barrel_endcap_region.MB_n):
            Hist_conv.append(ROOT.TH1D("Hist_conv"+order[count] ,"Hist_conv"+order[count],7,PtlBinsMu))
            Hist_conv[count].Sumw2()

            Hist_prompt.append(ROOT.TH1D("Hist_prompt"+order[count],"Hist_prompt"+order[count],7,PtlBinsMu))
            Hist_prompt[count].Sumw2()

            Hist_fakes.append(ROOT.TH1D("Hist_fakes"+order[count],"Hist_fakes"+order[count],7,PtlBinsMu))
            Hist_fakes[count].Sumw2()

            Hist_BDfakes.append(ROOT.TH1D("Hist_BDfakes"+order[count],"Hist_BDfakes"+order[count],7,PtlBinsMu))
            Hist_BDfakes[count].Sumw2()
        
        else:
            Hist_conv.append(ROOT.TH1D("Hist_conv"+order[count] ,"Hist_conv"+order[count],var_nBins, var_plotLow, var_plotHigh))
            Hist_conv[count].Sumw2()

            Hist_prompt.append(ROOT.TH1D("Hist_prompt"+order[count],"Hist_prompt"+order[count],var_nBins, var_plotLow, var_plotHigh))
            Hist_prompt[count].Sumw2()

            Hist_fakes.append(ROOT.TH1D("Hist_fakes"+order[count],"Hist_fakes"+order[count],var_nBins, var_plotLow, var_plotHigh))
            Hist_fakes[count].Sumw2()

            Hist_BDfakes.append(ROOT.TH1D("Hist_BDfakes"+order[count],"Hist_BDfakes"+order[count],var_nBins, var_plotLow, var_plotHigh))
            Hist_BDfakes[count].Sumw2()
        

    # End: Define histograms for uncertainty studies (separating fakes per type of fake) 
 
    genEventSumw = 1.0
    if not isData:
        if mc_genEventSumw is not None:
            genEventSumw = mc_genEventSumw
            print(
                "--- MC normalization: campaign genEventSumw = " + str(genEventSumw)
            )
        else:
            genEventSumw = read_genEventSumw_from_events_tree(fTemplateTree)
            print("--- MC normalization: genEventSumw = " + str(genEventSumw))
        if genEventSumw <= 0:
            raise RuntimeError(
                "Invalid genEventSumw=%s for sample %s" % (genEventSumw, Nickname)
            )

    has_zl_candidate = (
        fTemplateTree.GetBranch("bestZIdx") is not None
        and fTemplateTree.GetBranch("ZLCand_lepIdx") is not None
    )
    has_zll_candidate = (
        fTemplateTree.GetBranch("ZLLbest2P2FIdx") is not None
        and fTemplateTree.GetBranch("ZLLbest3P1FIdx") is not None
    )
    is_sr_only_sample = Nickname in ("HZZ", "ggZZ")
    if not is_sr_only_sample and (
        not has_zl_candidate or not has_zll_candidate
    ):
        raise RuntimeError(
            "Sample %s is missing Z+X control-region branches "
            "(ZLCand/ZLLCand). Check the NanoAOD production configuration."
            % Nickname
        )
    if is_sr_only_sample:
        # HZZ signal and ggZZ continuum samples are only used for nominal 4P
        # H->ZZ SR plots, not as Z+X CR/FR inputs.
        has_zl_candidate = False
        has_zll_candidate = False

    for event in fTemplateTree:
        iEvt+=1
        if(iEvt%50000==0):
            print ("---- Processing event: " + str(iEvt) + "/" + str(nentries))     
        # weight
            weight = 1.0
        if not isData:
            weight = LUMI_INT * event.overallEventWeight / genEventSumw

        if (
            varName in ("ptl3", "mZ1", "mEt")
            and has_zl_candidate
            and (event.bestZIdx >= 0)
            and (event.ZLCand_lepIdx >= 0)
        ):
            lep_Hindex = [
                event.ZCand_l1Idx[event.bestZIdx],
                event.ZCand_l2Idx[event.bestZIdx],
                event.ZLCand_lepIdx,
            ]

            lep_tight = event.Lepton_ZZFullSel[lep_Hindex[2]]
            lep_iso = getLeptonIso(event, lep_Hindex[2])
            idL3 = event.Lepton_pdgId[lep_Hindex[2]]
            lep = ROOT.TLorentzVector()
            lep.SetPtEtaPhiM(event.Lepton_pt[lep_Hindex[2]],event.Lepton_eta[lep_Hindex[2]],event.Lepton_phi[lep_Hindex[2]],event.Lepton_mass[lep_Hindex[2]])
            pTL3  = lep.Pt()
            etaL3 = lep.Eta()
            phiL3 = lep.Phi()
            
            lep_id = event.Lepton_pdgId

            lep_1 = ROOT.TLorentzVector()
            lep_2 = ROOT.TLorentzVector()
            lep_1.SetPtEtaPhiM(event.Lepton_pt[lep_Hindex[0]],event.Lepton_eta[lep_Hindex[0]],event.Lepton_phi[lep_Hindex[0]],event.Lepton_mass[lep_Hindex[0]])
            lep_2.SetPtEtaPhiM(event.Lepton_pt[lep_Hindex[1]],event.Lepton_eta[lep_Hindex[1]],event.Lepton_phi[lep_Hindex[1]],event.Lepton_mass[lep_Hindex[1]])
            
            massZ1 = (lep_1+lep_2).M()


            TestVar=False
            FillVar=0.
            
            if (varName=="mZ1"):
                TestVar= 1
                FillVar= massZ1
            
            if (varName=="mEt"):
                TestVar= 1
                FillVar = getMET_pt(event)

            
            if (varName=="ptl3"):
                TestVar = (math.fabs(massZ1 - 91.188) < 7) and (getMET_pt(event) < 25) 
                FillVar= pTL3

            if ((abs(idL3) == 11) and (math.fabs(etaL3) < 1.497) and TestVar):
                h1D_FRel_EB_d.Fill(FillVar, weight)
                PartOriginFromFlav(isData, event, lep_Hindex, FillVar, weight, Hist_prompt[barrel_endcap_region.EB_d], Hist_fakes[barrel_endcap_region.EB_d], Hist_BDfakes[barrel_endcap_region.EB_d], Hist_conv[barrel_endcap_region.EB_d], False)
                if (lep_tight and TestVar):
                    PartOriginFromFlav(isData, event, lep_Hindex, FillVar, weight, Hist_prompt[barrel_endcap_region.EB_n], Hist_fakes[barrel_endcap_region.EB_n],Hist_BDfakes[barrel_endcap_region.EB_n],Hist_conv[barrel_endcap_region.EB_n], False)
                    h1D_FRel_EB.Fill(FillVar, weight)

            if ((abs(idL3) == 11) and (math.fabs(etaL3) > 1.497) and TestVar):
                h1D_FRel_EE_d.Fill(FillVar, weight)
                PartOriginFromFlav(isData, event, lep_Hindex, FillVar, weight, Hist_prompt[barrel_endcap_region.EE_d], Hist_fakes[barrel_endcap_region.EE_d],Hist_BDfakes[barrel_endcap_region.EE_d], Hist_conv[barrel_endcap_region.EE_d], False)
                
                if (lep_tight and TestVar):
                    PartOriginFromFlav(isData, event, lep_Hindex, FillVar, weight, Hist_prompt[barrel_endcap_region.EE_n], Hist_fakes[barrel_endcap_region.EE_n], Hist_BDfakes[barrel_endcap_region.EE_n], Hist_conv[barrel_endcap_region.EE_n], False)
                    h1D_FRel_EE.Fill(FillVar, weight)
           
            if ((abs(idL3) == 13) and (math.fabs(etaL3) < 1.2) and TestVar):
                h1D_FRmu_EB_d.Fill(FillVar, weight)
                PartOriginFromFlav(isData, event, lep_Hindex, FillVar, weight, Hist_prompt[barrel_endcap_region.MB_d], Hist_fakes[barrel_endcap_region.MB_d], Hist_BDfakes[barrel_endcap_region.MB_d], Hist_conv[barrel_endcap_region.MB_d], False)
                
                if (lep_tight and (lep_iso < 0.35) and TestVar):
                    PartOriginFromFlav(isData, event, lep_Hindex, FillVar, weight, Hist_prompt[barrel_endcap_region.MB_n], Hist_fakes[barrel_endcap_region.MB_n],Hist_BDfakes[barrel_endcap_region.MB_n], Hist_conv[barrel_endcap_region.MB_n], False)
                    h1D_FRmu_EB.Fill(FillVar, weight)

            if ((abs(idL3) == 13) and (math.fabs(etaL3) > 1.2) and TestVar):
                h1D_FRmu_EE_d.Fill(FillVar, weight)
                PartOriginFromFlav(isData, event, lep_Hindex, FillVar, weight, Hist_prompt[barrel_endcap_region.ME_d], Hist_fakes[barrel_endcap_region.ME_d],Hist_BDfakes[barrel_endcap_region.ME_d], Hist_conv[barrel_endcap_region.ME_d], False)
                
                if (lep_tight and (lep_iso<0.35) and TestVar):
                    PartOriginFromFlav(isData, event, lep_Hindex, FillVar, weight, Hist_prompt[barrel_endcap_region.ME_n], Hist_fakes[barrel_endcap_region.ME_n],Hist_BDfakes[barrel_endcap_region.ME_n], Hist_conv[barrel_endcap_region.ME_n], False)
                    h1D_FRmu_EE.Fill(FillVar, weight)
  
            if VBS and varName == "ptl3" and TestVar:
                cleaned_jets = _cleaned_vbs_jets(
                    event,
                    [lep_1.Eta(), lep_2.Eta(), etaL3],
                    [lep_1.Phi(), lep_2.Phi(), phiL3],
                )
                for tag, min_jets in VBS_FR_JET_THRESHOLDS:
                    if len(cleaned_jets) >= min_jets:
                        _fill_fr_hist_set(
                            jet_fr_hists[tag],
                            idL3,
                            etaL3,
                            FillVar,
                            weight,
                            lep_tight,
                            lep_iso,
                        )
                        if abs(idL3) == 11 and math.fabs(etaL3) < 1.497:
                            _fill_origin_hist_set(
                                jet_origin_hists[tag], isData, event, lep_Hindex,
                                barrel_endcap_region.EB_d,
                                barrel_endcap_region.EB_n,
                                FillVar, weight, lep_tight,
                            )
                        elif abs(idL3) == 11 and math.fabs(etaL3) > 1.497:
                            _fill_origin_hist_set(
                                jet_origin_hists[tag], isData, event, lep_Hindex,
                                barrel_endcap_region.EE_d,
                                barrel_endcap_region.EE_n,
                                FillVar, weight, lep_tight,
                            )
                        elif abs(idL3) == 13 and math.fabs(etaL3) < 1.2:
                            _fill_origin_hist_set(
                                jet_origin_hists[tag], isData, event, lep_Hindex,
                                barrel_endcap_region.MB_d,
                                barrel_endcap_region.MB_n,
                                FillVar, weight,
                                lep_tight and lep_iso < 0.35,
                            )
                        elif abs(idL3) == 13 and math.fabs(etaL3) > 1.2:
                            _fill_origin_hist_set(
                                jet_origin_hists[tag], isData, event, lep_Hindex,
                                barrel_endcap_region.ME_d,
                                barrel_endcap_region.ME_n,
                                FillVar, weight,
                                lep_tight and lep_iso < 0.35,
                            )
  

        elif (
            has_zll_candidate
            and ((event.ZLLbest2P2FIdx >= 0) or (event.ZLLbest3P1FIdx >= 0))
            and varName == "ptl3"
        ):
            for zllIdx in (event.ZLLbest3P1FIdx, event.ZLLbest2P2FIdx):
                if zllIdx < 0:
                    continue

                lep_Hindex = [
                    event.ZLLCand_Z1l1Idx[zllIdx],
                    event.ZLLCand_Z1l2Idx[zllIdx],
                    event.ZLLCand_Z2l1Idx[zllIdx],
                    event.ZLLCand_Z2l2Idx[zllIdx],
                ]

                lep_tight = []
                lep_iso = []
                idL = []
                for k in range(4):
                    lep_tight.append(event.Lepton_ZZFullSel[lep_Hindex[k]])
                    lep_iso.append(getLeptonIso(event, lep_Hindex[k]))
                    idL.append(event.Lepton_pdgId[lep_Hindex[k]])

                lep_3 = ROOT.TLorentzVector()
                lep_3.SetPtEtaPhiM(event.Lepton_pt[lep_Hindex[2]], event.Lepton_eta[lep_Hindex[2]], event.Lepton_phi[lep_Hindex[2]], event.Lepton_mass[lep_Hindex[2]])
                lep_4 = ROOT.TLorentzVector()
                lep_4.SetPtEtaPhiM(event.Lepton_pt[lep_Hindex[3]], event.Lepton_eta[lep_Hindex[3]], event.Lepton_phi[lep_Hindex[3]], event.Lepton_mass[lep_Hindex[3]])
                lep_1 = ROOT.TLorentzVector()
                lep_2 = ROOT.TLorentzVector()
                lep_1.SetPtEtaPhiM(event.Lepton_pt[lep_Hindex[0]], event.Lepton_eta[lep_Hindex[0]], event.Lepton_phi[lep_Hindex[0]], event.Lepton_mass[lep_Hindex[0]])
                lep_2.SetPtEtaPhiM(event.Lepton_pt[lep_Hindex[1]], event.Lepton_eta[lep_Hindex[1]], event.Lepton_phi[lep_Hindex[1]], event.Lepton_mass[lep_Hindex[1]])

                nFailedLeptons3 = not(lep_tight[2] and ((abs(idL[2])==11) or (abs(idL[2])==13 and lep_iso[2]<0.35)))
                nFailedLeptons4 = not(lep_tight[3] and ((abs(idL[3])==11) or (abs(idL[3])==13 and lep_iso[3]<0.35)))
                # Do not assume that lepton 3 is the fake candidate: Z2 is
                # charge-ordered. Select the actual failing lepton(s) instead.
                failed_z2_positions = [
                    position
                    for position, failed in (
                        (2, nFailedLeptons3),
                        (3, nFailedLeptons4),
                    )
                    if failed
                ]
                m4l = event.ZLLCand_mass[zllIdx]


                if VBS:
                    etaL = [lep_1.Eta(), lep_2.Eta(), lep_3.Eta(), lep_4.Eta()]
                    phiL = [lep_1.Phi(), lep_2.Phi(), lep_3.Phi(), lep_4.Phi()]
                    if not _pass_vbs_selection(event, zllIdx, etaL, phiL, m4l):
                        continue
                        
                if len(failed_z2_positions) == 1:
                    h1D_m4l_3P1F.Fill(m4l, weight)
                    _fill_cr_origins(CR, CRregion._3P1F, isData, event, lep_Hindex, m4l, weight, failed_z2_positions)
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44):
                        _fill_cr_origins(CR, CRregion._3P1F_4e, isData, event, lep_Hindex, m4l, weight, failed_z2_positions)
                        h1D_m4l_3P1F_4e.Fill(m4l, weight)
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52):
                        _fill_cr_origins(CR, CRregion._3P1F_4mu, isData, event, lep_Hindex, m4l, weight, failed_z2_positions)
                        h1D_m4l_3P1F_4mu.Fill(m4l, weight)
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48):
                        if (abs(idL[2])+abs(idL[3])==22):
                            h1D_m4l_3P1F_2mu2e.Fill(m4l, weight)
                            _fill_cr_origins(CR, CRregion._3P1F_2mu2e, isData, event, lep_Hindex, m4l, weight, failed_z2_positions)
                        if (abs(idL[2])+abs(idL[3])==26):
                            h1D_m4l_3P1F_2e2mu.Fill(m4l, weight)
                            _fill_cr_origins(CR, CRregion._3P1F_2e2mu, isData, event, lep_Hindex, m4l, weight, failed_z2_positions)

                if len(failed_z2_positions) == 2:
                    h1D_m4l_2P2F.Fill(m4l, weight)
                    _fill_cr_origins(CR, CRregion._2P2F, isData, event, lep_Hindex, m4l, weight, failed_z2_positions)
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44):
                        _fill_cr_origins(CR, CRregion._2P2F_4e, isData, event, lep_Hindex, m4l, weight, failed_z2_positions)
                        h1D_m4l_2P2F_4e.Fill(m4l, weight)
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52):
                        _fill_cr_origins(CR, CRregion._2P2F_4mu, isData, event, lep_Hindex, m4l, weight, failed_z2_positions)
                        h1D_m4l_2P2F_4mu.Fill(m4l, weight)
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48):
                        if (abs(idL[2])+abs(idL[3])==22):
                            h1D_m4l_2P2F_2mu2e.Fill(m4l, weight)
                            _fill_cr_origins(CR, CRregion._2P2F_2mu2e, isData, event, lep_Hindex, m4l, weight, failed_z2_positions)
                        if (abs(idL[2])+abs(idL[3])==26):
                            h1D_m4l_2P2F_2e2mu.Fill(m4l, weight)
                            _fill_cr_origins(CR, CRregion._2P2F_2e2mu, isData, event, lep_Hindex, m4l, weight, failed_z2_positions)
                
        elif (event.bestCandIdx >= 0) and (varName == "ptl3"):
            # 4P signal region: full H4l candidates (not ZLL CR indices).
            zzIdx = event.bestCandIdx
            lep_Hindex = [
                event.ZZCand_Z1l1Idx[zzIdx],
                event.ZZCand_Z1l2Idx[zzIdx],
                event.ZZCand_Z2l1Idx[zzIdx],
                event.ZZCand_Z2l2Idx[zzIdx],
            ]
            if _n_failed_leptons(event, lep_Hindex) != 0:
                continue

            idL = [event.Lepton_pdgId[lep_Hindex[k]] for k in range(4)]
            lep_etas = [event.Lepton_eta[lep_Hindex[k]] for k in range(4)]
            lep_phis = [event.Lepton_phi[lep_Hindex[k]] for k in range(4)]
            m4l = event.ZZCand_mass[zzIdx]
            fill_var = m4l

            if VBS and not _pass_vbs_selection_zzcand(
                event, zzIdx, lep_etas, lep_phis, m4l
            ):
                continue

            _fill_4p_m4l_histograms(
                h1D_m4l_4P,
                h1D_m4l_4P_4e,
                h1D_m4l_4P_4mu,
                h1D_m4l_4P_2e2mu,
                h1D_m4l_4P_2mu2e,
                fill_var,
                weight,
                idL,
            )
            PartOriginFromFlav(isData, event, lep_Hindex, fill_var, weight, CR[type_of_fake.Prom][CRregion._4P], CR[type_of_fake.Fake][CRregion._4P], CR[type_of_fake.BDfake][CRregion._4P], CR[type_of_fake.Conv][CRregion._4P], False)
            if ((abs(idL[0]) + abs(idL[1]) + abs(idL[2]) + abs(idL[3])) == 44):
                PartOriginFromFlav(isData, event, lep_Hindex, fill_var, weight, CR[type_of_fake.Prom][CRregion._4P_4e], CR[type_of_fake.Fake][CRregion._4P_4e], CR[type_of_fake.BDfake][CRregion._4P_4e], CR[type_of_fake.Conv][CRregion._4P_4e], False)
            if ((abs(idL[0]) + abs(idL[1]) + abs(idL[2]) + abs(idL[3])) == 52):
                PartOriginFromFlav(isData, event, lep_Hindex, fill_var, weight, CR[type_of_fake.Prom][CRregion._4P_4mu], CR[type_of_fake.Fake][CRregion._4P_4mu], CR[type_of_fake.BDfake][CRregion._4P_4mu], CR[type_of_fake.Conv][CRregion._4P_4mu], False)
            if ((abs(idL[0]) + abs(idL[1]) + abs(idL[2]) + abs(idL[3])) == 48):
                if abs(idL[2]) + abs(idL[3]) == 26:
                    PartOriginFromFlav(isData, event, lep_Hindex, fill_var, weight, CR[type_of_fake.Prom][CRregion._4P_2e2mu], CR[type_of_fake.Fake][CRregion._4P_2e2mu], CR[type_of_fake.BDfake][CRregion._4P_2e2mu], CR[type_of_fake.Conv][CRregion._4P_2e2mu], False)
                if abs(idL[2]) + abs(idL[3]) == 22:
                    PartOriginFromFlav(isData, event, lep_Hindex, fill_var, weight, CR[type_of_fake.Prom][CRregion._4P_2mu2e], CR[type_of_fake.Fake][CRregion._4P_2mu2e], CR[type_of_fake.BDfake][CRregion._4P_2mu2e], CR[type_of_fake.Conv][CRregion._4P_2mu2e], False)

      


    hists = _pack_analyzeZX_histograms(
        h1D_FRel_EB,
        h1D_FRel_EE,
        h1D_FRmu_EB,
        h1D_FRmu_EE,
        h1D_FRel_EB_d,
        h1D_FRel_EE_d,
        h1D_FRmu_EB_d,
        h1D_FRmu_EE_d,
        h1D_m4l_2P2F,
        h1D_m4l_3P1F,
        h1D_m4l_4P,
        h1D_m4l_2P2F_2e2mu,
        h1D_m4l_3P1F_2e2mu,
        h1D_m4l_2P2F_2mu2e,
        h1D_m4l_3P1F_2mu2e,
        h1D_m4l_2P2F_4e,
        h1D_m4l_3P1F_4e,
        h1D_m4l_2P2F_4mu,
        h1D_m4l_3P1F_4mu,
        h1D_m4l_4P_2e2mu,
        h1D_m4l_4P_4e,
        h1D_m4l_4P_4mu,
        h1D_m4l_4P_2mu2e,
        CR_Prom,
        CR_Fake,
        CR_Conv,
        CR_BDfake,
        Hist_conv,
        Hist_prompt,
        Hist_fakes,
        Hist_BDfakes,
        jet_fr_hists,
        jet_origin_hists,
    )

    if write_output:
        _write_analyzeZX_histograms(hists, Nickname, varName, isData)
    else:
        return hists
    
