#import stage
import argparse
import json
import os
import sys

import numpy as np
import ROOT
from ROOT import TCanvas, TH1F, TF1, TLegend, gPad, THStack, TColor

from estimateZX import (
    VBS,
    _estimate_output_name,
    estimateZX,
    estimateZXChunkBatch,
    estimateZXFromChunks,
    merge_estimateZX_partials,
)
from main_FR_CR import (
    ANALYSIS_EOS_DIR,
    CHUNKS_PER_JOB,
    FINAL_HIST_DIR as STAGE1_HIST_DIR,
    chunks_per_job,
    compute_campaign_genEventSumw_by_nickname,
    discover_chunk_files,
    partial_hist_suffix,
    set_chunks_per_job,
    use_chunks,
)

FR_TAG_CHOICES = ("nominal", "ge1j", "ge2j") if VBS else ("nominal",)
FAKE_RATE_FILES = {
    "nominal": "Hist_Data_ptl3_WZremoved.root",
    "ge1j": "Hist_Data_ptl3_VBSFR_ge1j_WZremoved.root",
    "ge2j": "Hist_Data_ptl3_VBSFR_ge2j_WZremoved.root",
}
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# VBS: stage-2 estimates live next to stage-1 hists in ROOT_FILES/ (overwrite).
# H->ZZ: stage2/ subdirectory.
FINAL_ESTIMATE_DIR = (
    STAGE1_HIST_DIR if VBS else os.path.join(ANALYSIS_EOS_DIR, "stage2")
)
PARTIAL_ESTIMATE_DIR = os.path.join(
    STAGE1_HIST_DIR if VBS else FINAL_ESTIMATE_DIR,
    "partials_stage2" if VBS else "partials",
)
# Submit file, batch list, and out/err/log files must live on AFS since
# standard schedds cannot read/write them directly on /eos (see main_FR_CR.py).
CONDOR_DIR = os.path.join(SCRIPT_DIR, "condor_htozz", "stage2")
varName = "ptl3"

fileList = [
    "root://eosuser.cern.ch//eos/user/b/brunetd/ZXvariable_synch/STXS_samples/Data2022/ZZ4lAnalysis_AllData2022.root",
    "root://eosuser.cern.ch//eos/user/b/brunetd/ZXvariable_synch/ZZto4l_MC2022_CJLST/ZZ4lAnalysis.root",
]
RootNickNames = ["Data", "ZZ"]

chunkList = [
    {
        "nickname": "Data",
        "chunk_dir": "/eos/user/b/brunetd/DATA2024/DATA2024_Instructions_CONDOR",
        "chunk_glob": "*_Chunk*/ZZ4lAnalysis.root",
    },
    {
        "nickname": "ZZ",
        "chunk_dir": "/eos/user/b/brunetd/ZZTo4l2024MC/ZZ2024_CONDOR_Instructions",
        "chunk_glob": "ZZTo4l_Chunk*/ZZ4lAnalysis.root",
    },
]


def _manifest_path(var_name):
    return os.path.join(PARTIAL_ESTIMATE_DIR, "stage2_batches_%s.json" % var_name)


def _fake_rate_path(fr_tag):
    return os.path.join(STAGE1_HIST_DIR, FAKE_RATE_FILES[fr_tag])


def build_stage2_batches(var_name):
    batches = []
    for sample in chunkList:
        nickname = sample["nickname"]
        chunk_files = discover_chunk_files(
            sample["chunk_dir"],
            sample.get("chunk_glob", "*_Chunk*/ZZ4lAnalysis.root"),
        )
        for start in range(0, len(chunk_files), chunks_per_job()):
            batch_chunks = chunk_files[start : start + chunks_per_job()]
            suffix = partial_hist_suffix(batch_chunks)
            batches.append(
                {
                    "batch_index": len(batches),
                    "nickname": nickname,
                    "varName": var_name,
                    "suffix": suffix,
                    "chunks": batch_chunks,
                    "partial_files": {
                        fr_tag: "%s_%s.root"
                        % (_estimate_output_name(nickname, var_name, fr_tag), suffix)
                        for fr_tag in FR_TAG_CHOICES
                    },
                }
            )
    return batches


def write_stage2_manifest(var_name):
    os.makedirs(PARTIAL_ESTIMATE_DIR, exist_ok=True)
    batches = build_stage2_batches(var_name)
    campaign_genEventSumw = compute_campaign_genEventSumw_by_nickname(chunkList)
    manifest = {
        "varName": var_name,
        "fr_tags": list(FR_TAG_CHOICES),
        "fake_rate_files": FAKE_RATE_FILES,
        "chunks_per_job": chunks_per_job(),
        "partial_dir": PARTIAL_ESTIMATE_DIR,
        "final_dir": FINAL_ESTIMATE_DIR,
        "campaign_genEventSumw": campaign_genEventSumw,
        "batches": batches,
    }
    path = _manifest_path(var_name)
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2)
    print("Wrote %d batches to %s" % (len(batches), path))
    return manifest


def load_stage2_manifest(var_name):
    with open(_manifest_path(var_name)) as handle:
        return json.load(handle)


def run_condor_batch(batch_index, var_name):
    manifest = load_stage2_manifest(var_name)
    if batch_index < 0 or batch_index >= len(manifest["batches"]):
        raise RuntimeError(
            "batch_index %d out of range (0..%d)"
            % (batch_index, len(manifest["batches"]) - 1)
        )
    batch = manifest["batches"][batch_index]
    mc_genEventSumw = None
    if batch["nickname"] != "Data":
        mc_genEventSumw = manifest.get("campaign_genEventSumw", {}).get(
            batch["nickname"]
        )
        if mc_genEventSumw is None:
            raise RuntimeError(
                "Manifest missing campaign_genEventSumw for %s — re-run "
                "--write-condor-submit for stage 2"
                % batch["nickname"]
            )

    for fr_tag in FR_TAG_CHOICES:
        fake_rate_path = _fake_rate_path(fr_tag)
        if not os.path.isfile(fake_rate_path):
            raise RuntimeError(
                "Missing %s — run stage-1 ptl3 merge + WZremoval_from_FR_comp.py first"
                % fake_rate_path
            )
        print(
            "Condor batch %d/%d: sample %s, fr_tag %s, suffix %s, %d chunks"
            % (
                batch_index,
                len(manifest["batches"]) - 1,
                batch["nickname"],
                fr_tag,
                batch["suffix"],
                len(batch["chunks"]),
            )
        )
        out_path = estimateZXChunkBatch(
            fake_rate_path,
            batch["chunks"],
            batch["nickname"],
            PARTIAL_ESTIMATE_DIR,
            batch["suffix"],
            mc_genEventSumw=mc_genEventSumw,
            varName=var_name,
            fr_tag=fr_tag,
        )
        print("Wrote partial output: %s" % out_path)


def merge_stage2_partials(var_name):
    os.makedirs(FINAL_ESTIMATE_DIR, exist_ok=True)
    manifest = load_stage2_manifest(var_name)
    for fr_tag in FR_TAG_CHOICES:
        fake_rate_path = _fake_rate_path(fr_tag)
        if not os.path.isfile(fake_rate_path):
            raise RuntimeError("Missing %s" % fake_rate_path)
        for sample in chunkList:
            nickname = sample["nickname"]
            expected = [
                os.path.join(PARTIAL_ESTIMATE_DIR, batch["partial_files"][fr_tag])
                for batch in manifest["batches"]
                if batch["nickname"] == nickname
            ]
            if not expected:
                raise RuntimeError(
                    "Manifest %s has no batches for sample %s"
                    % (_manifest_path(var_name), nickname)
                )
            partial_paths = merge_estimateZX_partials(
                PARTIAL_ESTIMATE_DIR,
                nickname,
                fake_rate_path,
                final_dir=FINAL_ESTIMATE_DIR,
                varName=var_name,
                fr_tag=fr_tag,
                partial_paths=expected,
            )
            print(
                "Merged %d partial files for %s -> %s in %s"
                % (
                    len(partial_paths),
                    nickname,
                    _estimate_output_name(nickname, var_name, fr_tag),
                    FINAL_ESTIMATE_DIR,
                )
            )


def write_condor_submit_files(var_name):
    manifest = write_stage2_manifest(var_name)
    condor_dir = CONDOR_DIR
    os.makedirs(condor_dir, exist_ok=True)

    batch_list_path = os.path.join(condor_dir, "stage2_batch_list_%s.txt" % var_name)
    with open(batch_list_path, "w") as handle:
        for batch in manifest["batches"]:
            handle.write("%d\n" % batch["batch_index"])

    sub_path = os.path.join(condor_dir, "stage2_batches_%s.sub" % var_name)
    with open(sub_path, "w") as handle:
        handle.write(
            """universe        = vanilla
executable      = %(executable)s
arguments       = %(var)s $(batch_index)
# Do not write one .out/.err per job on AFS (errno 27 / schedd load).
output          = /dev/null
error           = /dev/null
log             = %(condor_dir)s/stage2_batches_%(var)s.$(ClusterId).log
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
                "executable": os.path.join(SCRIPT_DIR, "run_stage2_batch.sh"),
                "condor_dir": condor_dir,
                "batch_list": batch_list_path,
            }
        )

    print("Condor submit file: %s" % sub_path)
    print("Batch list (%d jobs): %s" % (len(manifest["batches"]), batch_list_path))
    print("Submit with: condor_submit %s" % sub_path)
    print(
        "Each job writes partials for all FR variants: %s"
        % ", ".join(FR_TAG_CHOICES)
    )


def run_sequential_chunks(var_name):
    print("Stage 2: ZX SR estimate for Data and ZZ, varName = %s.\n" % var_name)
    for fr_tag in FR_TAG_CHOICES:
        if not os.path.isfile(_fake_rate_path(fr_tag)):
            raise RuntimeError("Missing %s" % _fake_rate_path(fr_tag))
    for sample in chunkList:
        nickname = sample["nickname"]
        chunk_files = discover_chunk_files(
            sample["chunk_dir"],
            sample.get("chunk_glob", "*_Chunk*/ZZ4lAnalysis.root"),
        )
        print(
            "\n=== Sample %s: processing %d chunks from %s ===\n"
            % (nickname, len(chunk_files), sample["chunk_dir"])
        )
        for fr_tag in FR_TAG_CHOICES:
            print("--- FR variant: %s ---" % fr_tag)
            estimateZXFromChunks(
                _fake_rate_path(fr_tag),
                chunk_files,
                nickname,
                varName=var_name,
                fr_tag=fr_tag,
            )


def run_sequential_legacy(var_name):
    print("Stage 2: ZX SR estimate for Data and ZZ, varName = %s.\n" % var_name)
    for fr_tag in FR_TAG_CHOICES:
        if not os.path.isfile(_fake_rate_path(fr_tag)):
            raise RuntimeError("Missing %s" % _fake_rate_path(fr_tag))
    for i in range(len(fileList)):
        inFile = ROOT.TFile.Open(fileList[i], "READ")
        tree = inFile.Get("Events")
        print(
            '- File: "%s" has been opened.\n-- It has %d events.\n'
            % (fileList[i], tree.GetEntries())
        )
        for fr_tag in FR_TAG_CHOICES:
            print("--- FR variant: %s ---" % fr_tag)
            estimateZX(
                _fake_rate_path(fr_tag),
                tree,
                RootNickNames[i],
                varName=var_name,
                fr_tag=fr_tag,
            )
        inFile.Close()


def _parse_args():
    parser = argparse.ArgumentParser(description="Stage 2: estimateZX SR histograms")
    parser.add_argument(
        "--var",
        default=varName,
        help="Stage-2 estimate uses the ptl3 fake rates (default: %(default)s)",
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
        help="Process Condor batch N from stage2_batches_<var>.json",
    )
    parser.add_argument(
        "--merge-partials",
        action="store_true",
        help=(
            "Merge partial estimateZX* files into final outputs for all FR "
            "variants (estimateZXData, estimateZXData_ge1j, estimateZXData_ge2j, ...)"
        ),
    )
    parser.add_argument(
        "--write-condor-submit",
        action="store_true",
        help="Write manifest + condor/stage2_batches_<var>.sub",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    set_chunks_per_job(args.chunks_per_job)
    active_vars = [args.var]

    if args.write_condor_submit:
        for active_var in active_vars:
            print("\n=== Stage 2 Condor submit: %s ===" % active_var)
            write_condor_submit_files(active_var)
    elif args.condor_batch is not None:
        if len(active_vars) != 1:
            raise RuntimeError("--condor-batch needs a single --var, not all")
        run_condor_batch(args.condor_batch, active_vars[0])
    elif args.merge_partials:
        for active_var in active_vars:
            print("\n=== Merge stage-2 partials: %s ===" % active_var)
            merge_stage2_partials(active_var)
    elif use_chunks:
        for active_var in active_vars:
            run_sequential_chunks(active_var)
    else:
        for active_var in active_vars:
            run_sequential_legacy(active_var)
