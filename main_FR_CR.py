#import stage
import argparse
import json
import numpy as np
import ROOT
from ROOT import TCanvas, TH1F, TF1, TLegend, gPad, THStack, TColor
#import tdrstyle
import glob
import fnmatch
import os.path
import os
import re
import subprocess
from array import array
import yaml
import sys
import time
import math

from analyzeZX import (
    VBS as _VBS_MODE,
    analyzeZX,
    analyzeZXChunkBatch,
    analyzeZXFromChunks,
    merge_analyzeZX_partials,
    sum_genEventSumw_from_chunk_paths,
    _hist_data_mc_label,
    _is_data_nickname,
)
# Fill variables produced by one stage-1 campaign. Each name is a separate
# Condor submit (Hist_*_<var>_*.root). There is no sideband variable.
STAGE_VARIABLES = ("ptl3", "mZ1", "mEt")


def _eos_url(path):
    if path.startswith("root://"):
        return path
    if not path.startswith("/"):
        path = "/" + path
    return "root://eosuser.cern.ch/" + path


def _chunk_sort_key(path):
    match = re.search(r"_Chunk(\d+)", path)
    if match is None:
        return path
    return int(match.group(1))


def _discover_chunk_files_via_eos(chunk_dir, chunk_glob):
    """List ZZ4lAnalysis.root files via the eosuser MGM (bypasses FUSE)."""
    env = os.environ.copy()
    env["EOS_MGM_URL"] = "root://eosuser.cern.ch"
    try:
        proc = subprocess.run(
            ["eos", "find", "-f", chunk_dir, "--name", "ZZ4lAnalysis.root"],
            check=False,
            capture_output=True,
            text=True,
            env=env,
            timeout=600,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        print("eos find fallback failed for %s: %s" % (chunk_dir, exc))
        return []
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        print("eos find fallback failed for %s: %s" % (chunk_dir, err))
        return []

    matched = []
    for line in proc.stdout.splitlines():
        path = line.strip()
        if not path:
            continue
        rel = os.path.relpath(path, chunk_dir)
        if fnmatch.fnmatch(rel, chunk_glob):
            matched.append(path)
    return sorted(matched, key=_chunk_sort_key)


def discover_chunk_files(chunk_dir, chunk_glob="*_Chunk*/ZZ4lAnalysis.root"):
    pattern = os.path.join(chunk_dir, chunk_glob)
    # EOS FUSE mounts (/eos/user) occasionally return an empty listing on a
    # transient stall even though the files exist. glob.glob swallows the
    # underlying OSError and hands back [], so retry a few times before giving
    # up to avoid aborting a whole submit on a momentary hiccup.
    files = []
    attempts = 5
    for attempt in range(attempts):
        files = sorted(glob.glob(pattern), key=_chunk_sort_key)
        if files:
            break
        if attempt < attempts - 1:
            print(
                "No chunk files yet for pattern (attempt %d/%d): %s -- retrying..."
                % (attempt + 1, attempts, pattern)
            )
            time.sleep(3 * (attempt + 1))
    if not files and chunk_dir.startswith("/eos/"):
        print(
            "FUSE listing empty for %s — falling back to eos find (eosuser MGM)"
            % pattern
        )
        files = _discover_chunk_files_via_eos(chunk_dir, chunk_glob)
    if not files:
        raise RuntimeError("No chunk files found for pattern: %s" % pattern)
    return [_eos_url(f) for f in files]


def _chunk_dirname(chunk_path):
    path = chunk_path
    if path.startswith("root://"):
        path = re.sub(r"^root://[^/]+/*", "/", path)
    return os.path.basename(os.path.dirname(path))


def partial_hist_suffix(chunk_paths):
    """Suffix appended to Hist_* outputs, built from chunk folder names."""
    first = _chunk_dirname(chunk_paths[0])
    last = _chunk_dirname(chunk_paths[-1])
    if first == last:
        return first
    return "%s_to_%s" % (first, last)


# False = one hadded ZZ4lAnalysis.root per sample).
# True  = 2024 EOS chunk ntuples: fill per chunk, sum histograms, write one Hist_* per sample.
use_chunks = True

# EOS archives: H->ZZ uses stage1/; VBS overwrites ZpXVBS2024_results/ROOT_FILES/.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

_EOS_USER = "/eos/user/b/brunetd"
HTOZZ_EOS_DIR = os.path.join(_EOS_USER, "ZplusXpython_HtoZZ2024")
VBS_EOS_DIR = os.path.join(_EOS_USER, "ZpXVBS2024_results")
ANALYSIS_EOS_DIR = VBS_EOS_DIR if _VBS_MODE else HTOZZ_EOS_DIR
# Smaller batches = more Condor jobs and shorter wall time. Override with
# --chunks-per-job (1 is the maximum split).
CHUNKS_PER_JOB = 2


def set_chunks_per_job(n_chunks):
    global CHUNKS_PER_JOB
    n_chunks = int(n_chunks)
    if n_chunks < 1:
        raise ValueError("chunks-per-job must be >= 1")
    CHUNKS_PER_JOB = n_chunks


def chunks_per_job():
    return CHUNKS_PER_JOB
if _VBS_MODE:
    FINAL_HIST_DIR = os.path.join(ANALYSIS_EOS_DIR, "ROOT_FILES")
    PARTIAL_HIST_DIR = os.path.join(FINAL_HIST_DIR, "partials_stage1")
else:
    FINAL_HIST_DIR = os.path.join(ANALYSIS_EOS_DIR, "stage1")
    PARTIAL_HIST_DIR = os.path.join(FINAL_HIST_DIR, "partials")
# Standard schedds cannot read/write submit files, batch lists, or
# out/err/log files directly on /eos, so those stay on AFS. Only the
# histogram outputs written by the python job itself (via xrootd/TFile
# root://eosuser.cern.ch// URLs) live on EOS.
CONDOR_DIR = os.path.join(SCRIPT_DIR, "condor_htozz", "stage1")


# CJLST NanoAOD-based file list (2022, pre+post EE combined MC)
fileList = [
    "root://eosuser.cern.ch//eos/user/b/brunetd/ZXvariable_synch/STXS_samples/Data2022/ZZ4lAnalysis_AllData2022.root",
    "root://eosuser.cern.ch//eos/user/b/brunetd/MC2022_COMBINED_SAMPLES/DYJetsToLL/ZZ4lAnalysis.root",
    "root://eosuser.cern.ch//eos/user/b/brunetd/MC2022_COMBINED_SAMPLES/TTto2L2Nu/ZZ4lAnalysis.root",
    "root://eosuser.cern.ch//eos/user/b/brunetd/MC2022_COMBINED_SAMPLES/WZto3LNu/ZZ4lAnalysis.root",
    "root://eosuser.cern.ch//eos/user/b/brunetd/ZXvariable_synch/ZZto4l_MC2022_CJLST/ZZ4lAnalysis.root",
]
# Nicknames label output ROOT files (Hist_MC_ptl3_<nickname>.root). DY50 = DYJetsToLL M(ll)>50 GeV (not 50 fb-1).
RootNickNames = ["Data", "DY50", "TT", "WZ", "ZZ"]

varName = "all"  # iterates STAGE_VARIABLES: ptl3, mZ1, mEt


def stage_variables(requested):
    """Return the variables to run. 'all' is ptl3, mZ1, mEt."""
    if requested in (None, "", "all"):
        return list(STAGE_VARIABLES)
    if requested not in STAGE_VARIABLES:
        raise ValueError(
            "Unknown variable %s (allowed: %s or all)"
            % (requested, ", ".join(STAGE_VARIABLES))
        )
    return [requested]

# 2024 chunk-based ntuples (use when use_chunks = True; mirrors fileList / RootNickNames above)
chunkList = [

    {
        "nickname": "Data",
        "chunk_dir": "/eos/user/b/brunetd/DATA2024/DATA2024_Instructions_CONDOR",
        "chunk_glob": "*_Chunk*/ZZ4lAnalysis.root",
    },
    {
        "nickname": "DY50",
        "chunk_dir": "/eos/user/b/brunetd/DYJets2024",
        "chunk_glob": "DYJetsTo*_Chunk*/ZZ4lAnalysis.root",
        # 2024 DY is three exclusive samples (2E, 2Mu, 2Tau). Normalize each
        # component by its own Runs.genEventSumw, then merge them as DY50.
        "normalize_by_chunk_prefix": True,
    },
    {
        "nickname": "TT",
        "chunk_dir": "/eos/user/b/brunetd/TT2L2Nu2024/TT2L2Nu2024_Instructions_CONDOR",
        "chunk_glob": "TTto2L2Nu_Chunk*/ZZ4lAnalysis.root",
    },
    {
        "nickname": "WZ",
        "chunk_dir": "/eos/user/b/brunetd/WZ2024/WZ2024_Instructions_CONDOR",
        "chunk_glob": "WZto3LNu_Chunk*/ZZ4lAnalysis.root",
    },
    {
        "nickname": "ZZ",
        "chunk_dir": "/eos/user/b/brunetd/ZZTo4l2024MC/ZZ2024_CONDOR_Instructions",
        "chunk_glob": "ZZTo4l_Chunk*/ZZ4lAnalysis.root",
    },
    # HZZ and ggZZ are required by the nominal H->ZZ SR plotting macros.
    {
        "nickname": "HZZ",
        "chunk_dir": "/eos/user/b/brunetd/HZZ1252024/HZZ125_2024_CONDOR_Instructions",
        "chunk_glob": "*H125_Chunk*/ZZ4lAnalysis.root",
    },
    {
        "nickname": "ggZZ",
        "chunk_dir": "/eos/user/b/brunetd/ggZZ2024/ggZZ_2024_CONDOR_Instructions",
        "chunk_glob": "ggTo*_Contin_MCFM701_Chunk*/ZZ4lAnalysis.root",
    },
]


def _samples_for_var(var_name):
    """HZZ/ggZZ are needed only for the H->ZZ 4P histograms, filled with ptl3."""
    if var_name == "ptl3":
        return chunkList
    return [
        sample
        for sample in chunkList
        if sample["nickname"] not in ("HZZ", "ggZZ")
    ]


def _manifest_path(var_name):
    return os.path.join(PARTIAL_HIST_DIR, "stage1_batches_%s.json" % var_name)


def compute_campaign_genEventSumw_by_nickname(samples=None):
    """Sum Runs.genEventSumw over all chunks for each MC sample."""
    if samples is None:
        samples = chunkList
    totals = {}
    for sample in samples:
        nickname = sample["nickname"]
        if _is_data_nickname(nickname):
            continue
        chunk_files = discover_chunk_files(
            sample["chunk_dir"],
            sample.get("chunk_glob", "*_Chunk*/ZZ4lAnalysis.root"),
        )
        print(
            "Computing campaign genEventSumw for %s (%d chunks)..."
            % (nickname, len(chunk_files))
        )
        if sample.get("normalize_by_chunk_prefix", False):
            component_totals = {}
            for chunk_file in chunk_files:
                component = _chunk_dirname(chunk_file).split("_Chunk", 1)[0]
                component_totals.setdefault(component, []).append(chunk_file)
            totals[nickname] = {}
            for component, component_files in sorted(component_totals.items()):
                totals[nickname][component] = sum_genEventSumw_from_chunk_paths(
                    component_files
                )
                print(
                    "  %s/%s campaign genEventSumw = %.6g"
                    % (nickname, component, totals[nickname][component])
                )
        else:
            totals[nickname] = sum_genEventSumw_from_chunk_paths(chunk_files)
            print("  %s campaign genEventSumw = %.6g" % (nickname, totals[nickname]))
    return totals


def build_stage1_batches(var_name):
    batches = []
    for sample in _samples_for_var(var_name):
        nickname = sample["nickname"]
        chunk_files = discover_chunk_files(
            sample["chunk_dir"],
            sample.get("chunk_glob", "*_Chunk*/ZZ4lAnalysis.root"),
        )
        for start in range(0, len(chunk_files), chunks_per_job()):
            batch_chunks = chunk_files[start : start + chunks_per_job()]
            suffix = partial_hist_suffix(batch_chunks)
            is_data = _is_data_nickname(nickname)
            data_mc = _hist_data_mc_label(is_data)
            batches.append(
                {
                    "batch_index": len(batches),
                    "nickname": nickname,
                    "varName": var_name,
                    "suffix": suffix,
                    "chunks": batch_chunks,
                    "hist_main": "Hist_%s_%s_%s_%s.root"
                    % (data_mc, var_name, nickname, suffix),
                    "hist_id": "Hist_ID_%s_%s_%s.root"
                    % (var_name, nickname, suffix),
                }
            )
    return batches


def write_stage1_manifest(var_name):
    os.makedirs(PARTIAL_HIST_DIR, exist_ok=True)
    batches = build_stage1_batches(var_name)
    campaign_genEventSumw = compute_campaign_genEventSumw_by_nickname(
        _samples_for_var(var_name)
    )
    manifest = {
        "varName": var_name,
        "chunks_per_job": chunks_per_job(),
        "partial_hist_dir": PARTIAL_HIST_DIR,
        "final_hist_dir": FINAL_HIST_DIR,
        "campaign_genEventSumw": campaign_genEventSumw,
        "batches": batches,
    }
    path = _manifest_path(var_name)
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2)
    print("Wrote %d batches to %s" % (len(batches), path))
    return manifest


def load_stage1_manifest(var_name):
    path = _manifest_path(var_name)
    with open(path) as handle:
        return json.load(handle)


def run_condor_batch(batch_index, var_name):
    manifest = load_stage1_manifest(var_name)
    if batch_index < 0 or batch_index >= len(manifest["batches"]):
        raise RuntimeError(
            "batch_index %d out of range (0..%d)"
            % (batch_index, len(manifest["batches"]) - 1)
        )
    batch = manifest["batches"][batch_index]
    mc_genEventSumw = None
    if not _is_data_nickname(batch["nickname"]):
        mc_genEventSumw = manifest.get("campaign_genEventSumw", {}).get(
            batch["nickname"]
        )
        if mc_genEventSumw is None:
            raise RuntimeError(
                "Manifest missing campaign_genEventSumw for %s — re-run "
                "--write-condor-submit --var %s"
                % (batch["nickname"], var_name)
            )
    print(
        "Condor batch %d/%d: sample %s, suffix %s, %d chunks"
        % (
            batch_index,
            len(manifest["batches"]) - 1,
            batch["nickname"],
            batch["suffix"],
            len(batch["chunks"]),
        )
    )
    main_path, id_path = analyzeZXChunkBatch(
        batch["chunks"],
        batch["nickname"],
        var_name,
        PARTIAL_HIST_DIR,
        batch["suffix"],
        mc_genEventSumw=mc_genEventSumw,
    )
    print("Wrote partial outputs:\n  %s\n  %s" % (main_path, id_path))


def merge_stage1_partials(var_name):
    os.makedirs(FINAL_HIST_DIR, exist_ok=True)
    manifest = load_stage1_manifest(var_name)
    for sample in _samples_for_var(var_name):
        nickname = sample["nickname"]
        is_data = _is_data_nickname(nickname)
        expected = [
            os.path.join(PARTIAL_HIST_DIR, batch["hist_main"])
            for batch in manifest["batches"]
            if batch["nickname"] == nickname
        ]
        if not expected:
            raise RuntimeError(
                "Manifest %s has no batches for sample %s"
                % (_manifest_path(var_name), nickname)
            )
        glob_pattern = os.path.join(
            PARTIAL_HIST_DIR,
            "Hist_%s_%s_%s_*.root" % (_hist_data_mc_label(is_data), var_name, nickname),
        )
        leftover_n = len(glob.glob(glob_pattern)) - len(expected)
        partial_paths = merge_analyzeZX_partials(
            PARTIAL_HIST_DIR,
            nickname,
            var_name,
            is_data,
            final_dir=FINAL_HIST_DIR,
            partial_paths=expected,
        )
        extra = ""
        if leftover_n > 0:
            extra = " (ignored %d leftover files not in the current manifest)" % leftover_n
        print(
            "Merged %d partial files for %s -> Hist_*_%s_%s.root in %s%s"
            % (
                len(partial_paths),
                nickname,
                var_name,
                nickname,
                FINAL_HIST_DIR,
                extra,
            )
        )


def write_condor_submit_files(var_name):
    manifest = write_stage1_manifest(var_name)
    condor_dir = CONDOR_DIR
    os.makedirs(condor_dir, exist_ok=True)
    os.makedirs(PARTIAL_HIST_DIR, exist_ok=True)

    batch_list_path = os.path.join(condor_dir, "stage1_batch_list_%s.txt" % var_name)
    with open(batch_list_path, "w") as handle:
        for batch in manifest["batches"]:
            handle.write("%d\n" % batch["batch_index"])

    sub_path = os.path.join(condor_dir, "stage1_batches_%s.sub" % var_name)
    with open(sub_path, "w") as handle:
        handle.write(
            """universe        = vanilla
executable      = %(executable)s
arguments       = %(var)s $(batch_index)
# Do not write one .out/.err per job on AFS (errno 27 / schedd load).
output          = /dev/null
error           = /dev/null
log             = %(condor_dir)s/stage1_batches_%(var)s.$(ClusterId).log
should_transfer_files = NO
getenv          = True
request_cpus    = 2
request_memory  = 6000
request_disk    = 4000000
+JobFlavour     = "tomorrow"
max_idle        = 400
queue batch_index from %(batch_list)s
"""
            % {
                "var": var_name,
                "executable": os.path.join(SCRIPT_DIR, "run_stage1_batch.sh"),
                "condor_dir": condor_dir,
                "batch_list": batch_list_path,
            }
        )

    print("Condor submit file: %s" % sub_path)
    print("Batch list (%d jobs): %s" % (len(manifest["batches"]), batch_list_path))
    print("Submit with: condor_submit %s" % sub_path)
    return sub_path, batch_list_path


def run_sequential_chunks(var_name):
    print("Stage 1 (chunks): FR + Z+1L/CR histograms, varName = %s\n" % var_name)
    for sample in _samples_for_var(var_name):
        nickname = sample["nickname"]
        chunk_files = discover_chunk_files(
            sample["chunk_dir"],
            sample.get("chunk_glob", "*_Chunk*/ZZ4lAnalysis.root"),
        )
        print(
            "\n=== Sample %s: processing %d chunks from %s ===\n"
            % (nickname, len(chunk_files), sample["chunk_dir"])
        )
        analyzeZXFromChunks(chunk_files, nickname, var_name)


def run_sequential_legacy(var_name):
    print("Stage 1: FR + Z+1L/CR histograms, varName = %s\n" % var_name)
    for i in range(len(fileList)):
        inFile = ROOT.TFile.Open(fileList[i], "READ")
        tree = inFile.Get("Events")
        print(
            '- File: "%s" has been opened.\n-- It has %d events.\n'
            % (fileList[i], tree.GetEntries())
        )
        analyzeZX(tree, RootNickNames[i], var_name)
        inFile.Close()


def _parse_args():
    parser = argparse.ArgumentParser(description="Stage 1: FR + CR histograms (main_FR_CR)")
    parser.add_argument(
        "--var",
        default=varName,
        help="One of %s, or 'all' (default: all three)"
        % ", ".join(STAGE_VARIABLES),
    )
    parser.add_argument(
        "--chunks-per-job",
        type=int,
        default=CHUNKS_PER_JOB,
        metavar="N",
        help="Chunk files per Condor job (default: %(default)s)",
    )
    parser.add_argument(
        "--condor-batch",
        type=int,
        metavar="N",
        help="Process Condor batch N from stage1_batches_<var>.json",
    )
    parser.add_argument(
        "--merge-partials",
        action="store_true",
        help="Merge partial Hist_* files into final Hist_* outputs",
    )
    parser.add_argument(
        "--write-condor-submit",
        action="store_true",
        help="Write manifest + condor/stage1_batches_<var>.sub",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    set_chunks_per_job(args.chunks_per_job)
    active_vars = stage_variables(args.var)
    print("Stage 1 variables: %s" % ", ".join(active_vars))

    if args.write_condor_submit:
        for active_var in active_vars:
            print("\n=== Stage 1 Condor submit: %s ===" % active_var)
            write_condor_submit_files(active_var)
    elif args.condor_batch is not None:
        if len(active_vars) != 1:
            raise RuntimeError("--condor-batch needs a single --var, not all")
        run_condor_batch(args.condor_batch, active_vars[0])
    elif args.merge_partials:
        for active_var in active_vars:
            print("\n=== Merge stage-1 partials: %s ===" % active_var)
            merge_stage1_partials(active_var)
    elif use_chunks:
        for active_var in active_vars:
            run_sequential_chunks(active_var)
    else:
        for active_var in active_vars:
            run_sequential_legacy(active_var)
