from analyzeZX import setNEvents, setHistProperties, read_genEventSumw_from_events_tree, sum_genEventSumw_from_chunk_paths
import glob
import math
import os
import ROOT


# ---------------------------------------------------------------------------
# VBS mode
# ---------------------------------------------------------------------------
# When VBS is True, every CR/SR histogram fill is additionally gated by the
# ZZjj-inclusive VBS selection:
#   * on-shell Z bosons : 60 < mZ1, mZ2 < 120 GeV
#   * four-lepton mass   : m4l > 180 GeV
#   * >= 2 tagging jets  : pT > 30 GeV, |eta| < 4.7, cleaned by dR(l,j) > 0.4
#                          against the four ZLL-candidate leptons, mjj > 100 GeV
# The fake-rate measurement and the 2P2F/3P1F CR classification are NOT changed;
# the VBS cuts are simply required on top of the existing logic.
# Flip this single flag by hand: set it to True to run the VBS selection, False
# for the nominal H->ZZ behaviour. Condor jobs import this same module from AFS,
# so the value set here is exactly what the workers use (no env var needed).
# Keep this flag in sync with analyzeZX.py.
# False selects the inclusive H->ZZ workflow and disables VBS jet/m4l cuts.
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
# Require tight jet ID (jetId==6, i.e. tight + tightLepVeto) as used by the
# ZZAnalysis jetFiller. Set to None to disable the jet-ID requirement.
VBS_JET_ID = 6


def _detach_histogram(hist):
    if hist is None:
        return None
    cloned = hist.Clone()
    cloned.SetDirectory(0)
    return cloned


_ESTIMATE_SR_KEYS = [
    "h1D_m4l_SR_2P2F",
    "h1D_m4l_SR_3P1F",
    "h1D_m4l_Add_2P2F",
    "h1D_m4l_SR_2P2F_4e",
    "h1D_m4l_SR_3P1F_4e",
    "h1D_m4l_Add_2P2F_4e",
    "h1D_m4l_SR_2P2F_4mu",
    "h1D_m4l_SR_3P1F_4mu",
    "h1D_m4l_Add_2P2F_4mu",
    "h1D_m4l_SR_2P2F_2e2mu",
    "h1D_m4l_SR_3P1F_2e2mu",
    "h1D_m4l_Add_2P2F_2e2mu",
    "h1D_m4l_SR_2P2F_2mu2e",
    "h1D_m4l_SR_3P1F_2mu2e",
    "h1D_m4l_Add_2P2F_2mu2e",
]


def _pack_estimateZX_histograms(hists_dict):
    return {key: _detach_histogram(hists_dict[key]) for key in _ESTIMATE_SR_KEYS}


def _merge_estimateZX_histograms(accumulated, chunk_hists):
    for key in _ESTIMATE_SR_KEYS:
        accumulated[key].Add(chunk_hists[key])


def _estimate_output_name(Nickname, varName="ptl3", fr_tag=None):
    if varName == "ptl3":
        base = "estimateZX" + Nickname
    else:
        base = "estimateZX%s_%s" % (Nickname, varName)
    if fr_tag and fr_tag != "nominal":
        return "%s_%s" % (base, fr_tag)
    return base


def _estimate_partial_name(Nickname, varName, suffix, fr_tag=None):
    return "%s_%s.root" % (_estimate_output_name(Nickname, varName, fr_tag), suffix)


def _estimate_plot_config(varName):
    if VBS:
        return {"nbins": 51, "xmin": 180.0, "xmax": 1200.0, "bin_width": 20.0}
    return {"nbins": 40, "xmin": 70.0, "xmax": 870.0, "bin_width": 20.0}


def _load_fake_rate_histograms(fake_rate_file):
    fr_file = ROOT.TFile(fake_rate_file, "READ")
    if not fr_file or fr_file.IsZombie():
        raise RuntimeError("Cannot open fake rate file: %s" % fake_rate_file)
    fr_hists = {}
    for key, hname in (
        ("h1D_FRel_EB", "Data_FRel_EB"),
        ("h1D_FRel_EE", "Data_FRel_EE"),
        ("h1D_FRmu_EB", "Data_FRmu_EB"),
        ("h1D_FRmu_EE", "Data_FRmu_EE"),
    ):
        hist = fr_file.Get(hname)
        # Use "not hist" rather than "is None": a missing key yields a null
        # ROOT pointer (CPyCppyy_NoneType), which is not Python None.
        if not hist:
            fr_file.Close()
            raise RuntimeError("Missing %s in %s" % (hname, fake_rate_file))
        # Detach the clone from the file directory BEFORE closing the file,
        # otherwise ROOT deletes it on Close() and we get a null pointer.
        cloned = hist.Clone()
        cloned.SetDirectory(0)
        fr_hists[key] = cloned
    fr_file.Close()
    return fr_hists


def _write_estimateZX_partial(hists, Nickname, varName, out_dir, suffix, fr_tag=None):
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, _estimate_partial_name(Nickname, varName, suffix, fr_tag))
    out_file = ROOT.TFile(out_path, "RECREATE")
    for key in _ESTIMATE_SR_KEYS:
        hist = hists[key]
        hist.SetName(key)
        hist.Write()
    out_file.Close()
    return out_path


def _read_estimateZX_partial(partial_path):
    hists = {}
    in_file = ROOT.TFile.Open(partial_path, "READ")
    if in_file is None or in_file.IsZombie():
        raise RuntimeError("Cannot open partial estimateZX file: %s" % partial_path)
    for key in _ESTIMATE_SR_KEYS:
        key_obj = in_file.GetKey(key)
        if key_obj is None:
            in_file.Close()
            raise RuntimeError("Missing histogram %s in %s" % (key, partial_path))
        obj = key_obj.ReadObj()
        if obj is None or not obj.InheritsFrom(ROOT.TH1.Class()):
            in_file.Close()
            raise RuntimeError("Object %s in %s is not a TH1" % (key, partial_path))
        cloned = obj.Clone()
        cloned.SetDirectory(0)
        hists[key] = cloned
    in_file.Close()
    return hists


def merge_estimateZX_partials(
    partial_dir,
    Nickname,
    fake_rate_file,
    final_dir=".",
    varName="ptl3",
    fr_tag=None,
    partial_paths=None,
):
    if partial_paths is None:
        pattern = os.path.join(
            partial_dir, "%s_*.root" % _estimate_output_name(Nickname, varName, fr_tag)
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
                "Missing %d stage-2 partials for %s / %s, e.g. %s"
                % (len(missing), varName, Nickname, missing[0])
            )

    accumulated = None
    for partial_path in partial_paths:
        chunk_hists = _read_estimateZX_partial(partial_path)
        if accumulated is None:
            accumulated = chunk_hists
        else:
            _merge_estimateZX_histograms(accumulated, chunk_hists)

    fr_hists = _load_fake_rate_histograms(fake_rate_file)
    _write_estimateZX_histograms(
        accumulated, Nickname, fr_hists, out_dir=final_dir, varName=varName, fr_tag=fr_tag
    )
    return partial_paths


def estimateZXChunkBatch(
    fake_rate_file,
    file_paths,
    Nickname,
    out_dir,
    suffix,
    mc_genEventSumw=None,
    varName="ptl3",
    fr_tag=None,
):
    if not file_paths:
        raise RuntimeError("No chunk files provided for sample %s" % Nickname)

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
        chunk_hists = estimateZX(
            fake_rate_file,
            tree,
            Nickname,
            write_output=False,
            mc_genEventSumw=mc_genEventSumw,
            varName=varName,
        )
        inFile.Close()
        if accumulated is None:
            accumulated = chunk_hists
        else:
            _merge_estimateZX_histograms(accumulated, chunk_hists)

    return _write_estimateZX_partial(
        accumulated, Nickname, varName, out_dir, suffix, fr_tag=fr_tag
    )


def _write_estimateZX_histograms(hists, Nickname, fr_hists, out_dir=".", varName="ptl3", fr_tag=None):
    os.makedirs(out_dir, exist_ok=True)
    SaveRootFile = ROOT.TFile(
        os.path.join(out_dir, _estimate_output_name(Nickname, varName, fr_tag)), "RECREATE"
    )

    fr_hists["h1D_FRel_EB"].SetName("h1D_FRel_EB")
    fr_hists["h1D_FRel_EB"].Write()
    fr_hists["h1D_FRel_EE"].SetName("h1D_FRel_EE")
    fr_hists["h1D_FRel_EE"].Write()
    fr_hists["h1D_FRmu_EB"].SetName("h1D_FRmu_EB")
    fr_hists["h1D_FRmu_EB"].Write()
    fr_hists["h1D_FRmu_EE"].SetName("h1D_FRmu_EE")
    fr_hists["h1D_FRmu_EE"].Write()

    for key in _ESTIMATE_SR_KEYS:
        hists[key].SetName(key)
        hists[key].Write()

    SaveRootFile.Close()


def estimateZXFromChunks(
    FakeRateFile, file_paths, Nickname, mc_genEventSumw=None, varName="ptl3", fr_tag=None
):
    if not file_paths:
        raise RuntimeError("No chunk files provided for sample %s" % Nickname)

    fr_hists = _load_fake_rate_histograms(FakeRateFile)

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
        chunk_hists = estimateZX(
            FakeRateFile,
            tree,
            Nickname,
            write_output=False,
            mc_genEventSumw=mc_genEventSumw,
            varName=varName,
        )
        inFile.Close()
        if accumulated is None:
            accumulated = chunk_hists
        else:
            _merge_estimateZX_histograms(accumulated, chunk_hists)

    _write_estimateZX_histograms(
        accumulated, Nickname, fr_hists, out_dir=".", varName=varName, fr_tag=fr_tag
    )


def getFR(lep_id, lep_pt, lep_eta, h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE):

    if ((math.fabs(lep_id) == 11) and (math.fabs(lep_eta) < 1.497)):
        return h1D_FRel_EB.GetBinContent(h1D_FRel_EB.FindBin(lep_pt))

    if ((math.fabs(lep_id) == 11) and (math.fabs(lep_eta) > 1.497)):
        return h1D_FRel_EE.GetBinContent(h1D_FRel_EE.FindBin(lep_pt))

    if ((math.fabs(lep_id) == 13) and (math.fabs(lep_eta) < 1.2)):
        return h1D_FRmu_EB.GetBinContent(h1D_FRmu_EB.FindBin(lep_pt))

    if ((math.fabs(lep_id) == 13) and (math.fabs(lep_eta) > 1.2)):
        return h1D_FRmu_EE.GetBinContent(h1D_FRmu_EE.FindBin(lep_pt))

    return 0


def _delta_r(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    # |delta phi| folded into [0, pi] via acos(cos(.)) to avoid 2pi wrapping
    dphi = math.acos(max(-1.0, min(1.0, math.cos(phi1 - phi2))))
    return math.sqrt(deta * deta + dphi * dphi)


def _pass_vbs_jets(event, lep_etas, lep_phis):
    """ZZjj-inclusive jet requirement for the VBS selection.

    Jets are cleaned by hand against the four ZLL-candidate leptons: the stored
    nCleanedJets*/JetLeadingIdx branches are cleaned only against ZZFullSel
    leptons and are therefore not valid for the CR fake leptons (which fail
    ZZFullSel by construction). Returns True if there are at least two cleaned
    jets passing (pT, |eta|, jet-ID) whose two highest-pT members have
    mjj > VBS_MJJ_MIN.
    """
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

    if len(cleaned) < VBS_NJET_MIN:
        return False

    cleaned.sort(key=lambda x: x[0], reverse=True)
    mjj = (cleaned[0][1] + cleaned[1][1]).M()
    return mjj > VBS_MJJ_MIN


def _pass_vbs_selection(event, zllIdx, lep_etas, lep_phis, m4l):
    """Full ZZjj-inclusive VBS gate: on-shell Zs, m4l > 180, and tagging jets.

    Applied on top of the existing 2P2F/3P1F CR logic.
    """
    if not (VBS_MZ_LOW < event.ZLLCand_Z1mass[zllIdx] < VBS_MZ_HIGH):
        return False
    if not (VBS_MZ_LOW < event.ZLLCand_Z2mass[zllIdx] < VBS_MZ_HIGH):
        return False
    if not (m4l > VBS_M4L_MIN):
        return False
    return _pass_vbs_jets(event, lep_etas, lep_phis)





def estimateZX(
    FakeRateFile,
    tree,
    Nickname,
    write_output=True,
    mc_genEventSumw=None,
    varName="ptl3",
    fr_tag=None,
):
    LUMI_INT = 109820

    # define dummy histogram for CRs
    plot_cfg = _estimate_plot_config(varName)
    var_plotHigh = plot_cfg["xmax"]
    var_plotLow = plot_cfg["xmin"]
    var_nBins = plot_cfg["nbins"]

    h1D_dummyCR = ROOT.TH1D("dummyCR", "dummyCR", var_nBins, var_plotLow, var_plotHigh)

    varAxLabel = "m_{4l} [GeV]"

    binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.
    setHistProperties(h1D_dummyCR,1,1,1,0,0,varAxLabel,"Events/"+str(binWidth)+" [GeV]")
    
    # define CR histograms
    h1D_m4l_SR_2P2F = ROOT.TH1D("h1D_m4l_SR_2P2F","h1D_m4l_SR_2P2F",var_nBins, var_plotLow, var_plotHigh)
    h1D_m4l_SR_2P2F.Sumw2()
    h1D_m4l_SR_3P1F = ROOT.TH1D("h1D_m4l_SR_3P1F","h1D_m4l_SR_3P1F",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_SR_3P1F.Sumw2()
    h1D_m4l_Add_2P2F = ROOT.TH1D("h1D_m4l_Add_2P2F","h1D_m4l_Add_2P2F",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_Add_2P2F.Sumw2()
    
    # define CR histograms_4e
    h1D_m4l_SR_2P2F_4e = ROOT.TH1D("h1D_m4l_SR_2P2F_4e","h1D_m4l_SR_2P2F_4e",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_SR_2P2F_4e.Sumw2()
    h1D_m4l_SR_3P1F_4e= ROOT.TH1D("h1D_m4l_SR_3P1F_4e","h1D_m4l_SR_3P1F_4e",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_SR_3P1F_4e.Sumw2()
    h1D_m4l_Add_2P2F_4e = ROOT.TH1D("h1D_m4l_Add_2P2F_4e","h1D_m4l_Add_2P2F_4e",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_Add_2P2F_4e.Sumw2()
    
    # define CR histograms_4mu
    h1D_m4l_SR_2P2F_4mu = ROOT.TH1D("h1D_m4l_SR_2P2F_4mu","h1D_m4l_SR_2P2F_4mu",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_SR_2P2F_4mu.Sumw2()
    h1D_m4l_SR_3P1F_4mu = ROOT.TH1D("h1D_m4l_SR_3P1F_4mu","h1D_m4l_SR_3P1F_4mu",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_SR_3P1F_4mu.Sumw2()
    h1D_m4l_Add_2P2F_4mu = ROOT.TH1D("h1D_m4l_Add_2P2F_4mu","h1D_m4l_Add_2P2F_4mu",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_Add_2P2F_4mu.Sumw2()
    
    
    # define CR histograms_2e2mu
    h1D_m4l_SR_2P2F_2e2mu = ROOT.TH1D("h1D_m4l_SR_2P2F_2e2mu","h1D_m4l_SR_2P2F_2e2mu",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_SR_2P2F_2e2mu.Sumw2()
    h1D_m4l_SR_3P1F_2e2mu = ROOT.TH1D("h1D_m4l_SR_3P1F_2e2mu","h1D_m4l_SR_3P1F_2e2mu",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_SR_3P1F_2e2mu.Sumw2()
    h1D_m4l_Add_2P2F_2e2mu = ROOT.TH1D("h1D_m4l_Add_2P2F_2e2mu","h1D_m4l_Add_2P2F_2e2mu",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_Add_2P2F_2e2mu.Sumw2()
    
    # define CR histograms_2mu2e
    h1D_m4l_SR_2P2F_2mu2e = ROOT.TH1D("h1D_m4l_SR_2P2F_2mu2e","h1D_m4l_SR_2P2F_2mu2e",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_SR_2P2F_2mu2e.Sumw2()
    h1D_m4l_SR_3P1F_2mu2e = ROOT.TH1D("h1D_m4l_SR_3P1F_2mu2e","h1D_m4l_SR_3P1F_2mu2e",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_SR_3P1F_2mu2e.Sumw2()
    h1D_m4l_Add_2P2F_2mu2e = ROOT.TH1D("h1D_m4l_Add_2P2F_2mu2e","h1D_m4l_Add_2P2F_2mu2e",var_nBins, var_plotLow, var_plotHigh) 
    h1D_m4l_Add_2P2F_2mu2e.Sumw2()


    # Read in the FRs 

    FR_file = ROOT.TFile(FakeRateFile, "READ")
 
    h1D_FRel_EB = FR_file.Get("Data_FRel_EB")
    h1D_FRel_EE = FR_file.Get("Data_FRel_EE")
    h1D_FRmu_EB = FR_file.Get("Data_FRmu_EB")
    h1D_FRmu_EE = FR_file.Get("Data_FRmu_EE")

    isData = ("Data" in Nickname)
    iEvt = -1
    nentries = tree.GetEntries()
    lNEvents = setNEvents(Nickname)

    genEventSumw = 1.0
    if not isData:
        if mc_genEventSumw is not None:
            genEventSumw = mc_genEventSumw
            print(
                "--- MC normalization: campaign genEventSumw = " + str(genEventSumw)
            )
        else:
            genEventSumw = read_genEventSumw_from_events_tree(tree)
            print("--- MC normalization: genEventSumw = " + str(genEventSumw))
        if genEventSumw <= 0:
            raise RuntimeError(
                "Invalid genEventSumw=%s for sample %s" % (genEventSumw, Nickname)
            )

    for event in tree:
        iEvt+=1
        if(iEvt%50000==0):
            print ("---- Processing event: " + str(iEvt) + "/" + str(nentries))     
        # weight
        weight = 1.0
        if not isData:
            weight = LUMI_INT * event.overallEventWeight / genEventSumw

        if ((event.ZLLbest3P1FIdx >= 0) or (event.ZLLbest2P2FIdx >= 0)):
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
                pTL = []
                etaL = []
                phiL = []
                lep_p4 = []
                for k in range(4):
                    lep_tight.append(event.Lepton_ZZFullSel[lep_Hindex[k]])
                    if lep_Hindex[k] < event.nElectron:
                        lep_iso.append(event.Electron_pfRelIso03FsrCorr[lep_Hindex[k]])
                    else:
                        lep_iso.append(event.Muon_pfRelIso03FsrCorr[lep_Hindex[k] - event.nElectron])
                    idL.append(event.Lepton_pdgId[lep_Hindex[k]])
                    lep = ROOT.TLorentzVector()
                    lep.SetPtEtaPhiM(
                        event.Lepton_pt[lep_Hindex[k]],
                        event.Lepton_eta[lep_Hindex[k]],
                        event.Lepton_phi[lep_Hindex[k]],
                        event.Lepton_mass[lep_Hindex[k]],
                    )
                    pTL.append(lep.Pt())
                    etaL.append(lep.Eta())
                    phiL.append(lep.Phi())
                    lep_p4.append(lep)

                fr3 = getFR(idL[2], pTL[2], etaL[2], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE)
                fr4 = getFR(idL[3], pTL[3], etaL[3], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE)
                m4l = event.ZLLCand_mass[zllIdx]
                fill_var = m4l


                if VBS and not _pass_vbs_selection(event, zllIdx, etaL, phiL, m4l):
                    continue

                nFailedLeptons3 = not(lep_tight[2] and ((abs(idL[2])==11) or (abs(idL[2])==13 and lep_iso[2]<0.35)))
                nFailedLeptons4 = not(lep_tight[3] and ((abs(idL[3])==11) or (abs(idL[3])==13 and lep_iso[3]<0.35)))
                nFailedLeptonsZ2 = nFailedLeptons3 + nFailedLeptons4
                
                if (nFailedLeptonsZ2 == 1):
                    fr_bool_3 = (not(lep_tight[2] and ((abs(idL[2])==11) or (abs(idL[2])==13 and lep_iso[2]<0.35))))
                    fr_3 = fr_bool_3 * (fr3/(1-fr3))

                    fr_bool_4 = (not(lep_tight[3] and ((abs(idL[3])==11) or (abs(idL[3])==13 and lep_iso[3]<0.35))))
                    fr_4 = fr_bool_4 * (fr4/(1-fr4))

                    fr = fr_3 + fr_4
                    h1D_m4l_SR_3P1F.Fill(fill_var, weight * fr)
                    
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44):
                        h1D_m4l_SR_3P1F_4e.Fill(fill_var, weight* fr)
                        
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52):
                        h1D_m4l_SR_3P1F_4mu.Fill(fill_var, weight* fr)
                    
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48):
                        if ((abs(idL[2])+abs(idL[3]))== 26):
                            h1D_m4l_SR_3P1F_2e2mu.Fill(fill_var, weight* fr)

                        if ((abs(idL[2])+abs(idL[3]))== 22):
                            h1D_m4l_SR_3P1F_2mu2e.Fill(fill_var, weight* fr)
                
                elif (nFailedLeptonsZ2 == 2):
                    fr= (fr3/(1-fr3)) * (fr4/(1-fr4))
                    fr2 =(fr3/(1-fr3)) + (fr4/(1-fr4))
                    
                    # *Add_2P2F* represents 2P2F contribution in 3P1F
                    h1D_m4l_Add_2P2F.Fill(fill_var, weight * fr2)
                    h1D_m4l_SR_2P2F.Fill(fill_var, weight * fr)

                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44):
                        h1D_m4l_SR_2P2F_4e.Fill(fill_var, weight* fr)
                        h1D_m4l_Add_2P2F_4e.Fill(fill_var, weight * fr2)
                        
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52):
                        h1D_m4l_SR_2P2F_4mu.Fill(fill_var, weight* fr)
                        h1D_m4l_Add_2P2F_4mu.Fill(fill_var, weight * fr2)
                        
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48):
                        if ((abs(idL[2])+abs(idL[3]))== 26):
                            h1D_m4l_SR_2P2F_2e2mu.Fill(fill_var, weight* fr)
                            h1D_m4l_Add_2P2F_2e2mu.Fill(fill_var, weight * fr2)
                            
                        if ((abs(idL[2])+abs(idL[3]))== 22):
                            h1D_m4l_SR_2P2F_2mu2e.Fill(fill_var, weight* fr)
                            h1D_m4l_Add_2P2F_2mu2e.Fill(fill_var, weight * fr2)
                    
          
    sr_hists = {
        "h1D_m4l_SR_2P2F": h1D_m4l_SR_2P2F,
        "h1D_m4l_SR_3P1F": h1D_m4l_SR_3P1F,
        "h1D_m4l_Add_2P2F": h1D_m4l_Add_2P2F,
        "h1D_m4l_SR_2P2F_4e": h1D_m4l_SR_2P2F_4e,
        "h1D_m4l_SR_3P1F_4e": h1D_m4l_SR_3P1F_4e,
        "h1D_m4l_Add_2P2F_4e": h1D_m4l_Add_2P2F_4e,
        "h1D_m4l_SR_2P2F_4mu": h1D_m4l_SR_2P2F_4mu,
        "h1D_m4l_SR_3P1F_4mu": h1D_m4l_SR_3P1F_4mu,
        "h1D_m4l_Add_2P2F_4mu": h1D_m4l_Add_2P2F_4mu,
        "h1D_m4l_SR_2P2F_2e2mu": h1D_m4l_SR_2P2F_2e2mu,
        "h1D_m4l_SR_3P1F_2e2mu": h1D_m4l_SR_3P1F_2e2mu,
        "h1D_m4l_Add_2P2F_2e2mu": h1D_m4l_Add_2P2F_2e2mu,
        "h1D_m4l_SR_2P2F_2mu2e": h1D_m4l_SR_2P2F_2mu2e,
        "h1D_m4l_SR_3P1F_2mu2e": h1D_m4l_SR_3P1F_2mu2e,
        "h1D_m4l_Add_2P2F_2mu2e": h1D_m4l_Add_2P2F_2mu2e,
    }

    if write_output:
        _write_estimateZX_histograms(
            _pack_estimateZX_histograms(sr_hists),
            Nickname,
            {
                "h1D_FRel_EB": h1D_FRel_EB,
                "h1D_FRel_EE": h1D_FRel_EE,
                "h1D_FRmu_EB": h1D_FRmu_EB,
                "h1D_FRmu_EE": h1D_FRmu_EE,
            },
            varName=varName,
            fr_tag=fr_tag,
        )
        FR_file.Close()
    else:
        FR_file.Close()
        return _pack_estimateZX_histograms(sr_hists)

    
    