# Pial collateral network stability model

MATLAB code for the master's thesis *How do vascular loops in cerebral circulation
survive: endothelial electrical coupling as a loop stabilising mechanism*
(Computational Science, University of Amsterdam).

The model builds a two dimensional pial arterial network on the cortical surface,
lets vessel radii adapt to wall shear stress until they reach a steady state, and
compares four coupling variants (no coupling, arterial electrical coupling, arterial
plus fixed capillary current, and arterial plus a flow regulated capillary current).
The main structural readout is the number of surviving collateral loops.

## Layout

```
src/               all MATLAB source
  ElectricModel/   endothelial electrical network code
  penetrators.csv  anatomical input (see Data below)
  FinalResults/    generated .mat and .png outputs
ReportResults/     final figures used in the thesis (PDF)
REPO_GUIDE.md      detailed guide: model equations and which script makes each figure
```

## Running it

In MATLAB:

```matlab
cd src
smoke_test          % quick end-to-end check (< 1 min)
```

Each script adds its own folder to the path, so once you are in `src/` you can run
any experiment by name, for example `RunGeometryEnsemble`, `TestAutoregulationCurve`,
or `MultistabilityDemo`. Outputs are written to `src/FinalResults/`. `REPO_GUIDE.md`
lists which script produces each thesis figure and table.

Requires MATLAB (tested on R2023+). The paired significance test in the geometry
ensemble uses `signrank` from the Statistics Toolbox; nothing else needs a toolbox.

## Data

`src/penetrators.csv` is the catalogue of penetrating arterioles used as the
downstream boundary of the model: their radii, lengths, and areal density. These were
extracted from the H01 human cerebral cortex electron microscopy reconstruction
(Shapson-Coe et al., 2024) by the penetrator extraction step described in the thesis
Methods, which identifies penetrating vessels in the H01 blood vessel segmentation by
pial proximity, cortical depth, and orientation. Twenty-three penetrators were found,
with radii of 1.7 to 52.6 um and an areal density of 3.35 per mm2; the model samples
this catalogue with replacement to reach the physiological penetrator count.

The H01 dataset is released under a Creative Commons Attribution 4.0 licence
(https://h01-release.storage.googleapis.com). The values derived from it in
`penetrators.csv` carry the same attribution requirement.

## Licence

Code in this repository is released under the MIT licence (see `LICENSE`). The
`penetrators.csv` data is derived from H01 and is CC BY 4.0.
