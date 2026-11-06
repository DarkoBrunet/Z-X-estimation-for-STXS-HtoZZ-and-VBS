# Z + X estimation framework


This repository was developed using the methods developed by Vukašin Milošević and adjusted for running on NanoAOD files, adding some additional studies along the way (mainly for a detailed origin composition study).
This repository estimates the reducible Z+X background for  four-lepton
analyses, H->ZZ and VBS ZZjj. It uses CJLST NanoAOD 
ntuples, data-driven fake rates measurement, and uses MC samples
for prompt subtraction, and composition studies in different Control regions.

**Please keep in mind that the repository is under active development. Many of the analysis scripts can and will be cleaned up and optimised later, significantly altering the contents inside them. Diagnostic scripts, Condor submission wrappers, and other supporting files may be revised, consolidated, or removed as workflows are streamlined and the repository is reorganized.**

The main estimate has two event-processing stages. Stage 1 fills the fake-rate, control-region, and signal-region input histograms. After merging Stage 1, WZ is subtracted from the data fake-rate numerator and denominator. Stage 2 applies these corrected fake rates to data and ZZ control-region events. After merging Stage 2, the weighted ZZ contribution is subtracted from 3P1F, and the final Z+X yields and plots are produced.

## Environment and configuration

The production scripts expect the CMSSW 14_1_6 environment, PyROOT, access to the
configured EOS input and output areas, and HTCondor for full campaigns.

Input locations, EOS output locations, and analysis mode are configured in the
Python files. In the current checkout, `analyzeZX.py` and `estimateZX.py` have
`VBS = False`, selecting the H->ZZ paths. Composition-study scripts have
independent mode flags and currently default to `VBS = True`. Check the flags,
input files, output archive, and selection before running any workflow. 


### 1. Stage 1: fake-rate and control-region histograms

Create the batch manifest and Condor submit file, submit it, then merge only
after all jobs have completed successfully:

```bash
python3 main_FR_CR.py --write-condor-submit --var ptl3
condor_submit condor_htozz/stage1/stage1_batches_ptl3.sub
# Check that every expected batch succeeded, then:
python3 main_FR_CR.py --merge-partials --var ptl3
```

Repeat with `--var mZ1` or `--var mEt` if needed.

### 2. Remove WZ from the data fake-rate measurement

```bash
python3 WZremoval_from_FR_comp.py
```

This subtracts WZ from the data numerator and denominator and writes the
fake-rate ROOT file consumed by Stage 2. It requires the merged Stage 1 `ptl3`
data and WZ histograms.

### 3. Stage 2: apply fake rates to data and ZZ

```bash
python3 main_estimateZX_ntuples.py --write-condor-submit --var ptl3
condor_submit condor_htozz/stage2/stage2_batches_ptl3.sub
# Check that every expected batch succeeded, then:
python3 main_estimateZX_ntuples.py --merge-partials --var ptl3
```

Stage 2 requires the merged Stage 1 inputs and WZ-subtracted fake-rate file.
Its Condor jobs process data and ZZ chunks and write the weighted 2P2F/3P1F
estimates.

### 4. Final yields and plots

```bash
python3 estimate_final_numbers_macro.py
python3 plotting_macros.py
```

The first script subtracts the weighted ZZ contribution from 3P1F and forms the
final per-bin Z+X yields. The plotting script reads the configured EOS inputs
and writes the standard CR, SR, Z+1L, and fake-rate plots.

These scripts write to the configured EOS archive and may replace existing
outputs.

## Composition studies

These scripts are used for determining wheter the origin composition of fakes in the Control regions are different in the Z1L region, where we measure the FR.

1. `fr_composition_uncertainty.py` reads Stage 1 MC histograms and produces
   origin/composition tables and fake-rate transfer plots.
2. `composition_perpt_loop.py` applies per-pT composition rates to MC events.

3. `propagate_composition_zx.py` combines the transfer and estimate inputs and
   writes composition yield tables and plots.

Align their independent `VBS` flags and EOS inputs before running.

