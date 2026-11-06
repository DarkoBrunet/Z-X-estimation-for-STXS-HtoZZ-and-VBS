#!/usr/bin/env python3
"""Event-by-event MC ZX using per-pT 2P2F-mix composition FRs.

Origin FRs vs pT are already in FR_CR_transfer.csv. This loop applies the
f(pT, eta) of each failing lepton, then F=f/(1-f), the same way
estimateZX applies data fake rates:

  3P1F: w * F_failing
  2P2F: w * F3 * F4

Two mixes: noprompt (HF+LF+conv) and leftover prompt (WZ out of Z+1L,
prompt kept in the 2P2F mix). 3P1F is not used in the mix. ZZ is omitted
from 3P1F yields.

Set VBS below independently of analyzeZX.VBS:
  True  -> ZpXVBS2024_results, VBS CR/SR gate, ge1j/ge2j
  False -> ZplusXpython_HtoZZ2024, inclusive CR, nominal only

Run:
  python3 composition_perpt_loop.py --test-batches 1
  python3 composition_perpt_loop.py --write-condor-submit --workers 16
  condor_submit condor_htozz/composition_perpt/composition_perpt.sub
  python3 composition_perpt_loop.py --merge --workers 16
"""

from __future__ import print_function

import argparse
import csv
import json
import math
import os
import re
import sys

import ROOT


_EOS = "/eos/user/b/brunetd"
VBS = True
VBS_DIR = os.path.join(_EOS, "ZpXVBS2024_results")
HTOZZ_DIR = os.path.join(_EOS, "ZplusXpython_HtoZZ2024")
if VBS:
    STAGE1_DIR = os.path.join(VBS_DIR, "ROOT_FILES")
    TRANSFER_BASE = os.path.join(VBS_DIR, "FR_uncertainty_plots")
    MANIFEST_PATH = os.path.join(
        STAGE1_DIR, "partials_stage1", "stage1_batches_ptl3.json"
    )
    OUT_DIR = os.path.join(STAGE1_DIR, "final_numbers")
    DATA_FR_FILES = {
        "nominal": os.path.join(STAGE1_DIR, "Hist_Data_ptl3_Data.root"),
        "ge1j": os.path.join(STAGE1_DIR, "Hist_Data_ptl3_VBSFR_ge1j_Data.root"),
        "ge2j": os.path.join(STAGE1_DIR, "Hist_Data_ptl3_VBSFR_ge2j_Data.root"),
    }
    FR_TAGS = ("nominal", "ge1j", "ge2j")
    NBINS = 51
    M4L_MIN = 180.0
    M4L_MAX = 1200.0
    CONDOR_SUBDIR = "composition_perpt"
else:
    STAGE1_DIR = os.path.join(HTOZZ_DIR, "stage1")
    TRANSFER_BASE = os.path.join(HTOZZ_DIR, "FR_uncertainty_plots")
    MANIFEST_PATH = os.path.join(
        STAGE1_DIR, "partials", "stage1_batches_ptl3.json"
    )
    OUT_DIR = os.path.join(HTOZZ_DIR, "stage2", "final_numbers")
    DATA_FR_FILES = {
        "nominal": os.path.join(STAGE1_DIR, "Hist_Data_ptl3_Data.root"),
    }
    FR_TAGS = ("nominal",)
    NBINS = 40
    M4L_MIN = 70.0
    M4L_MAX = 870.0
    CONDOR_SUBDIR = "composition_perpt_htozz"
PARTIAL_DIR = os.path.join(OUT_DIR, "composition_perpt_partials")
MERGED_PATH = os.path.join(OUT_DIR, "composition_perpt_mc.root")
MC_NICKNAMES = ("DY50", "TT", "WZ", "ZZ")
CHANNELS = ("4e", "4mu", "2e2mu", "2mu2e")
KINDS = ("z1l_noprompt", "cr_noprompt", "z1l_leftover", "cr_leftover")
CRS = ("2P2F", "3P1F")
REGION_HIST = {
    "EB": "Data_FRel_EB",
    "EE": "Data_FRel_EE",
    "MB": "Data_FRmu_EB",
    "ME": "Data_FRmu_EE",
}
GETFR_KEYS = ("el_EB", "el_EE", "mu_EB", "mu_EE")
REGION_TO_GETFR = {
    "EB": "el_EB",
    "EE": "el_EE",
    "MB": "mu_EB",
    "ME": "mu_EE",
}
LUMI_INT = 109820.0
VBS_MZ_LOW = 60.0
VBS_MZ_HIGH = 120.0
VBS_M4L_MIN = 180.0
VBS_JET_PT_MIN = 30.0
VBS_JET_ETA_MAX = 4.7
VBS_JET_LEP_DR_MIN = 0.4
VBS_MJJ_MIN = 100.0
VBS_NJET_MIN = 2
VBS_JET_ID = 6


def _odds(fake):
    if fake <= 0.0 or fake >= 1.0:
        return 0.0
    return fake / (1.0 - fake)


def getFR(lep_id, lep_pt, lep_eta, hists):
    abs_id = math.fabs(lep_id)
    abs_eta = math.fabs(lep_eta)
    if abs_id == 11 and abs_eta < 1.497:
        hist = hists["el_EB"]
    elif abs_id == 11:
        hist = hists["el_EE"]
    elif abs_id == 13 and abs_eta < 1.2:
        hist = hists["mu_EB"]
    elif abs_id == 13:
        hist = hists["mu_EE"]
    else:
        return 0.0
    return float(hist.GetBinContent(hist.FindBin(lep_pt)))


def _delta_r(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = math.acos(max(-1.0, min(1.0, math.cos(phi1 - phi2))))
    return math.sqrt(deta * deta + dphi * dphi)


def _pass_vbs_jets(event, lep_etas, lep_phis):
    cleaned = []
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
    cleaned.sort(key=lambda item: item[0], reverse=True)
    return (cleaned[0][1] + cleaned[1][1]).M() > VBS_MJJ_MIN


def _pass_vbs_selection(event, zllIdx, lep_etas, lep_phis, m4l):
    if not (VBS_MZ_LOW < event.ZLLCand_Z1mass[zllIdx] < VBS_MZ_HIGH):
        return False
    if not (VBS_MZ_LOW < event.ZLLCand_Z2mass[zllIdx] < VBS_MZ_HIGH):
        return False
    if not (m4l > VBS_M4L_MIN):
        return False
    return _pass_vbs_jets(event, lep_etas, lep_phis)


def _failed_tight(tight, pdg, iso):
    return not (tight and ((abs(pdg) == 11) or (abs(pdg) == 13 and iso < 0.35)))


def _channel(idL):
    total = abs(idL[0]) + abs(idL[1]) + abs(idL[2]) + abs(idL[3])
    if total == 44:
        return "4e"
    if total == 52:
        return "4mu"
    if total == 48:
        z2 = abs(idL[2]) + abs(idL[3])
        if z2 == 26:
            return "2e2mu"
        if z2 == 22:
            return "2mu2e"
    return None


def _hist_name(fr_tag, kind, cr, channel):
    return "mc_%s_%s_%s_%s" % (fr_tag, kind, cr, channel)


def _book_hists():
    hists = {}
    for fr_tag in FR_TAGS:
        for kind in KINDS:
            for cr in CRS:
                for channel in CHANNELS:
                    name = _hist_name(fr_tag, kind, cr, channel)
                    hist = ROOT.TH1D(name, name, NBINS, M4L_MIN, M4L_MAX)
                    hist.Sumw2()
                    hist.SetDirectory(0)
                    hists[name] = hist
    return hists


def _add_hists(total, extra):
    if total is None:
        total = {}
        for name, hist in extra.items():
            clone = hist.Clone(name)
            clone.SetDirectory(0)
            total[name] = clone
        return total
    for name, hist in extra.items():
        total[name].Add(hist)
    return total


def _build_fr_hists():
    """TH1 f(pT) per region from the transfer CSV, data-FR binning."""
    out = {}
    for fr_tag in FR_TAGS:
        data_file = ROOT.TFile.Open(DATA_FR_FILES[fr_tag], "READ")
        if not data_file or data_file.IsZombie():
            raise RuntimeError("Cannot open %s" % DATA_FR_FILES[fr_tag])
        csv_path = os.path.join(TRANSFER_BASE, fr_tag, "FR_CR_transfer.csv")
        rows = []
        with open(csv_path, "r") as handle:
            rows = [
                row for row in csv.DictReader(handle)
                if row["bin"] != "inclusive"
            ]
        tag_hists = {}
        for kind, cr_key, field in (
            ("z1l_noprompt", "2P2F_noprompt", "f_z1l"),
            ("cr_noprompt", "2P2F_noprompt", "f_cr"),
            ("z1l_leftover", "2P2F_leftover", "f_z1l"),
            ("cr_leftover", "2P2F_leftover", "f_cr"),
        ):
            packed = {}
            for region, hname in REGION_HIST.items():
                proto = data_file.Get(hname)
                if not proto:
                    data_file.Close()
                    raise RuntimeError("Missing %s in %s" % (hname, DATA_FR_FILES[fr_tag]))
                hist = proto.Clone("fr_%s_%s_%s" % (fr_tag, kind, region))
                hist.Reset("ICES")
                hist.SetDirectory(0)
                packed[REGION_TO_GETFR[region]] = hist
            for row in rows:
                if row["cr"] != cr_key:
                    continue
                key = REGION_TO_GETFR.get(row["region"])
                if key not in packed:
                    continue
                hist = packed[key]
                pt_lo = float(row["pt_low"])
                ibin = hist.GetXaxis().FindFixBin(pt_lo + 1e-4)
                if ibin < 1 or ibin > hist.GetNbinsX():
                    continue
                hist.SetBinContent(ibin, float(row[field]))
            tag_hists[kind] = packed
        data_file.Close()
        out[fr_tag] = tag_hists
    return out


def _chunk_dirname(chunk_path):
    path = chunk_path
    if path.startswith("root://"):
        path = re.sub(r"^root://[^/]+/*", "/", path)
    return os.path.basename(os.path.dirname(path))


def _sumw_for_chunk(sumw_map, nickname, chunk_path):
    value = sumw_map[nickname]
    if isinstance(value, dict):
        component = _chunk_dirname(chunk_path).split("_Chunk", 1)[0]
        if component not in value:
            raise RuntimeError(
                "Missing genEventSumw for %s / %s" % (nickname, component)
            )
        return float(value[component])
    return float(value)


def _fill_file(fpath, nickname, genEventSumw, fr_hists, hists):
    tfile = ROOT.TFile.Open(fpath, "READ")
    if not tfile or tfile.IsZombie():
        print("WARNING: cannot open %s" % fpath)
        return 0
    tree = tfile.Get("Events")
    if not tree:
        tfile.Close()
        print("WARNING: no Events in %s" % fpath)
        return 0
    tree.SetBranchStatus("*", 0)
    for branch in (
        "ZLLbest3P1FIdx",
        "ZLLbest2P2FIdx",
        "ZLLCand_Z1l1Idx",
        "ZLLCand_Z1l2Idx",
        "ZLLCand_Z2l1Idx",
        "ZLLCand_Z2l2Idx",
        "ZLLCand_mass",
        "ZLLCand_Z1mass",
        "ZLLCand_Z2mass",
        "Lepton_ZZFullSel",
        "Lepton_pdgId",
        "Lepton_pt",
        "Lepton_eta",
        "Lepton_phi",
        "nLepton",
        "nElectron",
        "Electron_pfRelIso03FsrCorr",
        "Muon_pfRelIso03FsrCorr",
        "overallEventWeight",
        "nJet",
        "Jet_pt",
        "Jet_eta",
        "Jet_phi",
        "Jet_mass",
        "Jet_jetId",
    ):
        tree.SetBranchStatus(branch, 1)
    n_filled = 0
    include_3p1f = nickname != "ZZ"
    for event in tree:
        if event.ZLLbest3P1FIdx < 0 and event.ZLLbest2P2FIdx < 0:
            continue
        weight = LUMI_INT * event.overallEventWeight / genEventSumw
        for zllIdx in (event.ZLLbest3P1FIdx, event.ZLLbest2P2FIdx):
            if zllIdx < 0:
                continue
            lep_Hindex = [
                event.ZLLCand_Z1l1Idx[zllIdx],
                event.ZLLCand_Z1l2Idx[zllIdx],
                event.ZLLCand_Z2l1Idx[zllIdx],
                event.ZLLCand_Z2l2Idx[zllIdx],
            ]
            tight = []
            iso = []
            idL = []
            pTL = []
            etaL = []
            phiL = []
            for k in range(4):
                idx = lep_Hindex[k]
                tight.append(event.Lepton_ZZFullSel[idx])
                if idx < event.nElectron:
                    iso.append(event.Electron_pfRelIso03FsrCorr[idx])
                else:
                    iso.append(event.Muon_pfRelIso03FsrCorr[idx - event.nElectron])
                idL.append(event.Lepton_pdgId[idx])
                pTL.append(event.Lepton_pt[idx])
                etaL.append(event.Lepton_eta[idx])
                phiL.append(event.Lepton_phi[idx])
            m4l = event.ZLLCand_mass[zllIdx]
            if VBS and not _pass_vbs_selection(event, zllIdx, etaL, phiL, m4l):
                continue
            fail3 = _failed_tight(tight[2], idL[2], iso[2])
            fail4 = _failed_tight(tight[3], idL[3], iso[3])
            nfail = int(fail3) + int(fail4)
            channel = _channel(idL)
            if channel is None or nfail not in (1, 2):
                continue
            n_filled += 1
            cr = "3P1F" if nfail == 1 else "2P2F"
            if cr == "3P1F" and not include_3p1f:
                continue
            for fr_tag in FR_TAGS:
                f3 = {}
                f4 = {}
                for kind in KINDS:
                    f3[kind] = getFR(idL[2], pTL[2], etaL[2], fr_hists[fr_tag][kind])
                    f4[kind] = getFR(idL[3], pTL[3], etaL[3], fr_hists[fr_tag][kind])
                F3 = {kind: _odds(val) for kind, val in f3.items()}
                F4 = {kind: _odds(val) for kind, val in f4.items()}
                for kind in KINDS:
                    if cr == "3P1F":
                        wgt = weight * (fail3 * F3[kind] + fail4 * F4[kind])
                    else:
                        wgt = weight * F3[kind] * F4[kind]
                    hists[_hist_name(fr_tag, kind, cr, channel)].Fill(m4l, wgt)
    tfile.Close()
    return n_filled


def _load_jobs(max_batches=0, nicknames=None):
    with open(MANIFEST_PATH, "r") as handle:
        manifest = json.load(handle)
    sumw_map = manifest["campaign_genEventSumw"]
    jobs = []
    for batch in manifest["batches"]:
        nickname = batch["nickname"]
        if nickname not in MC_NICKNAMES:
            continue
        if nicknames and nickname not in nicknames:
            continue
        jobs.append(
            {
                "batch_index": batch["batch_index"],
                "nickname": nickname,
                "chunks": batch["chunks"],
                "sumw_map": sumw_map,
            }
        )
        if max_batches and len(jobs) >= max_batches:
            break
    return jobs


def _worker(payload):
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    fr_hists = _build_fr_hists()
    hists = _book_hists()
    n_filled = 0
    n_files = 0
    for ijob, job in enumerate(payload["jobs"]):
        nickname = job["nickname"]
        for fpath in job["chunks"]:
            genEventSumw = _sumw_for_chunk(job["sumw_map"], nickname, fpath)
            n_filled += _fill_file(fpath, nickname, genEventSumw, fr_hists, hists)
            n_files += 1
        if (ijob + 1) % 10 == 0 or ijob + 1 == len(payload["jobs"]):
            print(
                "  worker %d: %d/%d batches, %d CR fills"
                % (payload["worker_id"], ijob + 1, len(payload["jobs"]), n_filled),
                flush=True,
            )
    out_path = payload["out_path"]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tfile = ROOT.TFile(out_path, "RECREATE")
    for hist in hists.values():
        hist.SetDirectory(tfile)
        hist.Write()
    tfile.Close()
    return payload["worker_id"], n_files, n_filled, out_path


def merge_partials(partial_paths, out_path):
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    total = None
    for path in partial_paths:
        tfile = ROOT.TFile.Open(path, "READ")
        if not tfile or tfile.IsZombie():
            raise RuntimeError("Cannot open partial %s" % path)
        extra = {}
        for key in tfile.GetListOfKeys():
            hist = tfile.Get(key.GetName())
            clone = hist.Clone(key.GetName())
            clone.SetDirectory(0)
            extra[key.GetName()] = clone
        tfile.Close()
        total = _add_hists(total, extra)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    out = ROOT.TFile(out_path, "RECREATE")
    for hist in total.values():
        hist.SetDirectory(out)
        hist.Write()
    out.Close()
    print("Wrote %s from %d partials" % (out_path, len(partial_paths)))


def run_fill(workers, max_batches=0, nicknames=None, worker_id=None, merge_only=False):
    jobs = _load_jobs(max_batches=max_batches, nicknames=nicknames)
    if not jobs:
        raise RuntimeError("No MC batches to process")
    os.makedirs(PARTIAL_DIR, exist_ok=True)
    n_workers = max(1, min(int(workers), len(jobs)))
    chunks = [[] for _ in range(n_workers)]
    for i, job in enumerate(jobs):
        chunks[i % n_workers].append(job)
    payloads = []
    for i, group in enumerate(chunks):
        if not group:
            continue
        payloads.append(
            {
                "worker_id": i,
                "jobs": group,
                "out_path": os.path.join(PARTIAL_DIR, "worker_%d.root" % i),
            }
        )
    if merge_only:
        partials = [item["out_path"] for item in payloads]
        missing = [path for path in partials if not os.path.isfile(path)]
        if missing:
            raise RuntimeError("Missing %d partials, e.g. %s" % (len(missing), missing[0]))
        merge_partials(partials, MERGED_PATH)
        return MERGED_PATH
    if worker_id is not None:
        wanted = [item for item in payloads if item["worker_id"] == int(worker_id)]
        if not wanted:
            raise RuntimeError("No jobs for worker_id %s" % worker_id)
        payloads = wanted
        n_workers = 1
    print(
        "Per-pT MC ZX fill [%s]: %d batches, %d workers%s"
        % (
            "VBS" if VBS else "HtoZZ",
            len(jobs),
            len(payloads) if worker_id is None else n_workers,
            "" if worker_id is None else " (worker_id=%s)" % worker_id,
        )
    )
    if len(payloads) == 1:
        results = [_worker(payloads[0])]
    else:
        raise RuntimeError(
            "Refusing to multiprocessing-Pool on lxplus. "
            "Run one worker: --fill --workers N --worker-id K "
            "or submit condor/composition_perpt.sub"
        )
    partials = []
    for wid, n_files, n_filled, out_path in results:
        print(
            "  worker %d: %d files, %d CR fills -> %s"
            % (wid, n_files, n_filled, out_path)
        )
        partials.append(out_path)
    if worker_id is None:
        merge_partials(partials, MERGED_PATH)
    return MERGED_PATH


def write_condor_submit(n_workers=16):
    condor_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "condor_htozz",
        CONDOR_SUBDIR,
    )
    os.makedirs(condor_dir, exist_ok=True)
    wrapper = os.path.join(condor_dir, "run_composition_perpt_batch.sh")
    sub_path = os.path.join(condor_dir, "composition_perpt.sub")
    batch_list = os.path.join(condor_dir, "worker_ids.txt")
    with open(batch_list, "w") as handle:
        for i in range(int(n_workers)):
            handle.write("%d\n" % i)
    with open(wrapper, "w") as handle:
        handle.write(
            """#!/bin/bash
set -e
WORKER_ID="$1"
cd /afs/cern.ch/user/b/brunetd/CMSSW_14_1_6/src
cmsenv
cd ZplusXpython
exec python3 composition_perpt_loop.py --fill --workers %d --worker-id "$WORKER_ID"
"""
            % int(n_workers)
        )
    os.chmod(wrapper, 0o755)
    with open(sub_path, "w") as handle:
        handle.write(
            """universe        = vanilla
executable      = %(wrapper)s
arguments       = $(worker_id)
output          = /dev/null
error           = /dev/null
log             = %(condor_dir)s/composition_perpt.$(ClusterId).log
should_transfer_files = NO
getenv          = True
request_cpus    = 1
request_memory  = 6000
request_disk    = 2000000
+JobFlavour     = "tomorrow"
queue worker_id from %(batch_list)s
"""
            % {
                "wrapper": wrapper,
                "condor_dir": condor_dir,
                "batch_list": batch_list,
            }
        )
    print("Condor submit file: %s" % sub_path)
    print("Submit with: condor_submit %s" % sub_path)
    return sub_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fill", action="store_true")
    parser.add_argument("--test-batches", type=int, default=0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--worker-id", type=int, default=None)
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--write-condor-submit", action="store_true")
    parser.add_argument("--nickname", action="append", default=None)
    args = parser.parse_args()
    if args.write_condor_submit:
        write_condor_submit(n_workers=args.workers if args.workers > 1 else 16)
        return 0
    max_batches = args.test_batches
    if not args.fill and not max_batches and not args.merge:
        parser.error("use --fill, --merge, --test-batches N, or --write-condor-submit")
    run_fill(
        workers=1 if max_batches else args.workers,
        max_batches=max_batches,
        nicknames=args.nickname,
        worker_id=args.worker_id,
        merge_only=args.merge,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
