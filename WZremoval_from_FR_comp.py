import os

import ROOT

from analyzeZX import VBS
from main_FR_CR import FINAL_HIST_DIR as STAGE1_DIR


def _clone_required(tfile, name):
    hist = tfile.Get(name)
    if not hist:
        raise RuntimeError("Missing %s in %s" % (name, tfile.GetName()))
    cloned = hist.Clone()
    cloned.SetDirectory(0)
    return cloned


def _write_wz_removed(data_path, wz_path, out_path, required=False):
    if not os.path.isfile(data_path) or not os.path.isfile(wz_path):
        if required:
            raise RuntimeError("Missing input(s): %s, %s" % (data_path, wz_path))
        print("Skipping WZ removal, missing %s or %s" % (data_path, wz_path))
        return

    data_file = ROOT.TFile(data_path, "READ")
    wz_file = ROOT.TFile(wz_path, "READ")

    n_EB = _clone_required(data_file, "Data_FRel_EB_n")
    n_EE = _clone_required(data_file, "Data_FRel_EE_n")
    n_MB = _clone_required(data_file, "Data_FRmu_EB_n")
    n_ME = _clone_required(data_file, "Data_FRmu_EE_n")

    d_EB = _clone_required(data_file, "Data_FRel_EB_d")
    d_EE = _clone_required(data_file, "Data_FRel_EE_d")
    d_MB = _clone_required(data_file, "Data_FRmu_EB_d")
    d_ME = _clone_required(data_file, "Data_FRmu_EE_d")

    n_EB.Add(_clone_required(wz_file, "Data_FRel_EB_n"), -1)
    d_EB.Add(_clone_required(wz_file, "Data_FRel_EB_d"), -1)
    n_EB.Divide(d_EB)

    n_EE.Add(_clone_required(wz_file, "Data_FRel_EE_n"), -1)
    d_EE.Add(_clone_required(wz_file, "Data_FRel_EE_d"), -1)
    n_EE.Divide(d_EE)

    n_MB.Add(_clone_required(wz_file, "Data_FRmu_EB_n"), -1)
    d_MB.Add(_clone_required(wz_file, "Data_FRmu_EB_d"), -1)
    n_MB.Divide(d_MB)

    n_ME.Add(_clone_required(wz_file, "Data_FRmu_EE_n"), -1)
    d_ME.Add(_clone_required(wz_file, "Data_FRmu_EE_d"), -1)
    n_ME.Divide(d_ME)

    out_file = ROOT.TFile(out_path, "RECREATE")
    n_EB.SetName("Data_FRel_EB")
    n_EB.Write()
    n_EE.SetName("Data_FRel_EE")
    n_EE.Write()
    n_MB.SetName("Data_FRmu_EB")
    n_MB.Write()
    n_ME.SetName("Data_FRmu_EE")
    n_ME.Write()

    data_file.Close()
    wz_file.Close()
    out_file.Close()
    print("Wrote %s" % out_path)


print("WZ removal stage-1 directory: %s (VBS=%s)" % (STAGE1_DIR, VBS))

_write_wz_removed(
    os.path.join(STAGE1_DIR, "Hist_Data_ptl3_Data.root"),
    os.path.join(STAGE1_DIR, "Hist_MC_ptl3_WZ.root"),
    os.path.join(STAGE1_DIR, "Hist_Data_ptl3_WZremoved.root"),
    required=True,
)

for tag in (("ge1j", "ge2j") if VBS else ()):
    _write_wz_removed(
        os.path.join(
            STAGE1_DIR, "Hist_Data_ptl3_VBSFR_%s_Data.root" % tag
        ),
        os.path.join(
            STAGE1_DIR, "Hist_MC_ptl3_VBSFR_%s_WZ.root" % tag
        ),
        os.path.join(
            STAGE1_DIR,
            "Hist_Data_ptl3_VBSFR_%s_WZremoved.root" % tag,
        ),
    )
