import ROOT
import math
import os
import sys

from analyzeZX import VBS

_EOS_USER = "/eos/user/b/brunetd"
HTOZZ_EOS_DIR = os.path.join(_EOS_USER, "ZplusXpython_HtoZZ2024")
VBS_EOS_DIR = os.path.join(_EOS_USER, "ZpXVBS2024_results")
ANALYSIS_EOS_DIR = VBS_EOS_DIR if VBS else HTOZZ_EOS_DIR
# VBS estimates sit in ROOT_FILES/; H->ZZ uses stage2/.
STAGE2_DIR = (
    os.path.join(ANALYSIS_EOS_DIR, "ROOT_FILES")
    if VBS
    else os.path.join(ANALYSIS_EOS_DIR, "stage2")
)
# Final Z+X bin-by-bin combination (SR_tot = SR_3P1F - SR_2P2F). Not a
# separate numbered stage; its inputs/outputs live under stage 2.
FINAL_NUMBERS_DIR = os.path.join(STAGE2_DIR, "final_numbers")
os.makedirs(FINAL_NUMBERS_DIR, exist_ok=True)


class _Tee(object):
    """Mirror terminal output to a text file."""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, message):
        for stream in self.streams:
            stream.write(message)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


_original_stdout = sys.stdout
_output_txt_path = os.path.join(FINAL_NUMBERS_DIR, "estimate_final_numbers.txt")
_output_txt = open(_output_txt_path, "w")
sys.stdout = _Tee(_original_stdout, _output_txt)

def _resolve_estimate_file(basename):
    for filename in (basename, basename + ".root"):
        path = os.path.join(STAGE2_DIR, filename)
        if os.path.isfile(path):
            return path
    raise RuntimeError("Missing estimate file: %s" % basename)


_estimate_path = _resolve_estimate_file("estimateZXData")
_zz_path = _resolve_estimate_file("estimateZXZZ")

_fileData = ROOT.TFile(_estimate_path, "READ")



h1D_m4l_SR_2P2F = _fileData.Get("h1D_m4l_SR_2P2F")
h1D_m4l_SR_2P2F_4e = _fileData.Get("h1D_m4l_SR_2P2F_4e")
h1D_m4l_SR_2P2F_4mu = _fileData.Get("h1D_m4l_SR_2P2F_4mu")
h1D_m4l_SR_2P2F_2e2mu = _fileData.Get("h1D_m4l_SR_2P2F_2e2mu")
h1D_m4l_SR_2P2F_2mu2e = _fileData.Get("h1D_m4l_SR_2P2F_2mu2e")

h1D_m4l_SR_3P1F = _fileData.Get("h1D_m4l_SR_3P1F")
h1D_m4l_SR_3P1F_4e = _fileData.Get("h1D_m4l_SR_3P1F_4e")
h1D_m4l_SR_3P1F_4mu = _fileData.Get("h1D_m4l_SR_3P1F_4mu")
h1D_m4l_SR_3P1F_2e2mu = _fileData.Get("h1D_m4l_SR_3P1F_2e2mu")
h1D_m4l_SR_3P1F_2mu2e = _fileData.Get("h1D_m4l_SR_3P1F_2mu2e")


print("------------3P1F Data-----------------")
print("h1D_m4l_SR_3P1F = " + str(h1D_m4l_SR_3P1F.Integral()))
print("h1D_m4l_SR_3P1F_4e = " + str(h1D_m4l_SR_3P1F_4e.Integral()))
print("h1D_m4l_SR_3P1F_4mu = " + str(h1D_m4l_SR_3P1F_4mu.Integral()))
print("h1D_m4l_SR_3P1F_2e2mu = " + str(h1D_m4l_SR_3P1F_2e2mu.Integral()))
print("h1D_m4l_SR_3P1F_2mu2e = " + str(h1D_m4l_SR_3P1F_2mu2e.Integral()))
print("--------------------------------")


#Just to initialize them with correct binning, will do a bin by bin estimate in the following lines
h1D_m4l_SR_tot = h1D_m4l_SR_2P2F.Clone()
h1D_m4l_SR_tot_4e = h1D_m4l_SR_2P2F_4e.Clone()
h1D_m4l_SR_tot_4mu = h1D_m4l_SR_2P2F_4mu.Clone()
h1D_m4l_SR_tot_2e2mu = h1D_m4l_SR_2P2F_2e2mu.Clone()
h1D_m4l_SR_tot_2mu2e = h1D_m4l_SR_2P2F_2mu2e.Clone()

_fileZZ = ROOT.TFile(_zz_path, "READ")


ZZ_4e = _fileZZ.Get("h1D_m4l_SR_3P1F_4e")
ZZ_4mu = _fileZZ.Get("h1D_m4l_SR_3P1F_4mu")
ZZ_2e2mu = _fileZZ.Get("h1D_m4l_SR_3P1F_2e2mu")
ZZ_2mu2e = _fileZZ.Get("h1D_m4l_SR_3P1F_2mu2e")
ZZ_tot = _fileZZ.Get("h1D_m4l_SR_3P1F")

print("------------3P1F ZZ-----------------")
print("ZZ_SR_3P1F = " + str(ZZ_tot.Integral()))
print("ZZ_SR_3P1F_4e = " + str(ZZ_4e.Integral()))
print("ZZ_SR_3P1F_4mu = " + str(ZZ_4mu.Integral()))
print("ZZ_SR_3P1F_2e2mu = " + str(ZZ_2e2mu.Integral()))
print("ZZ_SR_3P1F_2mu2e = " + str(ZZ_2mu2e.Integral()))
print("--------------------------------")


h1D_m4l_SR_3P1F.Add(ZZ_tot,-1)
h1D_m4l_SR_3P1F_4e.Add(ZZ_4e,-1)
h1D_m4l_SR_3P1F_4mu.Add(ZZ_4mu,-1)
h1D_m4l_SR_3P1F_2e2mu.Add(ZZ_2e2mu,-1)
h1D_m4l_SR_3P1F_2mu2e.Add(ZZ_2mu2e,-1)

var_nBins = h1D_m4l_SR_2P2F.GetNbinsX()

for c in range (1, var_nBins + 1):

    if (h1D_m4l_SR_3P1F_4e.GetBinContent(c) <= 2*h1D_m4l_SR_2P2F_4e.GetBinContent(c)):
        h1D_m4l_SR_tot_4e.SetBinContent(c,h1D_m4l_SR_2P2F_4e.GetBinContent(c))
        h1D_m4l_SR_tot_4e.SetBinError(c,h1D_m4l_SR_2P2F_4e.GetBinError(c))
    
    else:
        h1D_m4l_SR_tot_4e.SetBinContent(c,h1D_m4l_SR_3P1F_4e.GetBinContent(c)-h1D_m4l_SR_2P2F_4e.GetBinContent(c))
        h1D_m4l_SR_tot_4e.SetBinError(c,math.sqrt(h1D_m4l_SR_2P2F_4e.GetBinError(c)*h1D_m4l_SR_2P2F_4e.GetBinError(c)+h1D_m4l_SR_3P1F_4e.GetBinError(c)*h1D_m4l_SR_3P1F_4e.GetBinError(c)))


    if (h1D_m4l_SR_3P1F_4mu.GetBinContent(c) <= 2*h1D_m4l_SR_2P2F_4mu.GetBinContent(c)):
    
        h1D_m4l_SR_tot_4mu.SetBinContent(c,h1D_m4l_SR_2P2F_4mu.GetBinContent(c))
        h1D_m4l_SR_tot_4mu.SetBinError(c,h1D_m4l_SR_2P2F_4mu.GetBinError(c))

    else:
        h1D_m4l_SR_tot_4mu.SetBinContent(c,h1D_m4l_SR_3P1F_4mu.GetBinContent(c)-h1D_m4l_SR_2P2F_4mu.GetBinContent(c))
        h1D_m4l_SR_tot_4mu.SetBinError(c,math.sqrt(h1D_m4l_SR_2P2F_4mu.GetBinError(c)*h1D_m4l_SR_2P2F_4mu.GetBinError(c)+h1D_m4l_SR_3P1F_4mu.GetBinError(c)*h1D_m4l_SR_3P1F_4mu.GetBinError(c)))



    if (h1D_m4l_SR_3P1F_2e2mu.GetBinContent(c) <= 2*h1D_m4l_SR_2P2F_2e2mu.GetBinContent(c)):
        h1D_m4l_SR_tot_2e2mu.SetBinContent(c,h1D_m4l_SR_2P2F_2e2mu.GetBinContent(c))
        h1D_m4l_SR_tot_2e2mu.SetBinError(c,h1D_m4l_SR_2P2F_2e2mu.GetBinError(c))

    else:
        h1D_m4l_SR_tot_2e2mu.SetBinContent(c,h1D_m4l_SR_3P1F_2e2mu.GetBinContent(c)-h1D_m4l_SR_2P2F_2e2mu.GetBinContent(c))
        h1D_m4l_SR_tot_2e2mu.SetBinError(c,math.sqrt(h1D_m4l_SR_2P2F_2e2mu.GetBinError(c)*h1D_m4l_SR_2P2F_2e2mu.GetBinError(c)+h1D_m4l_SR_3P1F_2e2mu.GetBinError(c)*h1D_m4l_SR_3P1F_2e2mu.GetBinError(c)))


    if (h1D_m4l_SR_3P1F_2mu2e.GetBinContent(c) <= 2*h1D_m4l_SR_2P2F_2mu2e.GetBinContent(c)):
        h1D_m4l_SR_tot_2mu2e.SetBinContent(c,h1D_m4l_SR_2P2F_2mu2e.GetBinContent(c))
        h1D_m4l_SR_tot_2mu2e.SetBinError(c,h1D_m4l_SR_2P2F_2mu2e.GetBinError(c))
    
    else:
        h1D_m4l_SR_tot_2mu2e.SetBinContent(c,h1D_m4l_SR_3P1F_2mu2e.GetBinContent(c)-h1D_m4l_SR_2P2F_2mu2e.GetBinContent(c))
        h1D_m4l_SR_tot_2mu2e.SetBinError(c,math.sqrt(h1D_m4l_SR_2P2F_2mu2e.GetBinError(c)*h1D_m4l_SR_2P2F_2mu2e.GetBinError(c)+h1D_m4l_SR_3P1F_2mu2e.GetBinError(c)*h1D_m4l_SR_3P1F_2mu2e.GetBinError(c)))

h1D_m4l_SR_tot.Reset()
h1D_m4l_SR_tot.Add(h1D_m4l_SR_tot_4e)
h1D_m4l_SR_tot.Add(h1D_m4l_SR_tot_4mu)
h1D_m4l_SR_tot.Add(h1D_m4l_SR_tot_2e2mu)
h1D_m4l_SR_tot.Add(h1D_m4l_SR_tot_2mu2e)

print("--------Final Z+X estimates----------")
print("h1D_m4l_SR_2P2F = " + str(h1D_m4l_SR_2P2F.Integral()))
print("h1D_m4l_SR_2P2F_4e = " + str(h1D_m4l_SR_2P2F_4e.Integral()))
print("h1D_m4l_SR_2P2F_4mu = " + str(h1D_m4l_SR_2P2F_4mu.Integral()))
print("h1D_m4l_SR_2P2F_2e2mu = " + str(h1D_m4l_SR_2P2F_2e2mu.Integral()))
print("h1D_m4l_SR_2P2F_2mu2e = " + str(h1D_m4l_SR_2P2F_2mu2e.Integral()))
print("--------------------------------")

print("--------------------------------")
print("h1D_m4l_SR_3P1F = " + str(h1D_m4l_SR_3P1F.Integral()))
print("h1D_m4l_SR_3P1F_4e = " + str(h1D_m4l_SR_3P1F_4e.Integral()))
print("h1D_m4l_SR_3P1F_4mu = " + str(h1D_m4l_SR_3P1F_4mu.Integral()))
print("h1D_m4l_SR_3P1F_2e2mu = " + str(h1D_m4l_SR_3P1F_2e2mu.Integral()))
print("h1D_m4l_SR_3P1F_2mu2e = " + str(h1D_m4l_SR_3P1F_2mu2e.Integral()))
print("--------------------------------")

print("--------------------------------")
print("h1D_m4l_SR_tot = " + str(h1D_m4l_SR_tot.Integral()))
print("h1D_m4l_SR_tot_4e = " + str(h1D_m4l_SR_tot_4e.Integral()))
print("h1D_m4l_SR_tot_4mu = " + str(h1D_m4l_SR_tot_4mu.Integral()))
print("h1D_m4l_SR_tot_2e2mu = " + str(h1D_m4l_SR_tot_2e2mu.Integral()))
print("h1D_m4l_SR_tot_2mu2e = " + str(h1D_m4l_SR_tot_2mu2e.Integral()))
print("--------------------------------")


_output_root_path = os.path.join(FINAL_NUMBERS_DIR, "estimate_final_numbers.root")
_output_root = ROOT.TFile(_output_root_path, "RECREATE")
for name, hist in [
    ("h1D_m4l_SR_2P2F", h1D_m4l_SR_2P2F),
    ("h1D_m4l_SR_2P2F_4e", h1D_m4l_SR_2P2F_4e),
    ("h1D_m4l_SR_2P2F_4mu", h1D_m4l_SR_2P2F_4mu),
    ("h1D_m4l_SR_2P2F_2e2mu", h1D_m4l_SR_2P2F_2e2mu),
    ("h1D_m4l_SR_2P2F_2mu2e", h1D_m4l_SR_2P2F_2mu2e),
    ("h1D_m4l_SR_3P1F_ZZsub", h1D_m4l_SR_3P1F),
    ("h1D_m4l_SR_3P1F_ZZsub_4e", h1D_m4l_SR_3P1F_4e),
    ("h1D_m4l_SR_3P1F_ZZsub_4mu", h1D_m4l_SR_3P1F_4mu),
    ("h1D_m4l_SR_3P1F_ZZsub_2e2mu", h1D_m4l_SR_3P1F_2e2mu),
    ("h1D_m4l_SR_3P1F_ZZsub_2mu2e", h1D_m4l_SR_3P1F_2mu2e),
    ("h1D_m4l_SR_tot", h1D_m4l_SR_tot),
    ("h1D_m4l_SR_tot_4e", h1D_m4l_SR_tot_4e),
    ("h1D_m4l_SR_tot_4mu", h1D_m4l_SR_tot_4mu),
    ("h1D_m4l_SR_tot_2e2mu", h1D_m4l_SR_tot_2e2mu),
    ("h1D_m4l_SR_tot_2mu2e", h1D_m4l_SR_tot_2mu2e),
]:
    out_hist = hist.Clone(name)
    out_hist.SetDirectory(_output_root)
    out_hist.Write()
_output_root.Close()

print("Wrote %s" % _output_txt_path)
print("Wrote %s" % _output_root_path)

sys.stdout = _original_stdout
_output_txt.close()