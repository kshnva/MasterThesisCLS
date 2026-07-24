# Pial collateral network stability model

MATLAB code for the thesis *How do vascular loops in cerebral circulation survive:
endothelial electrical coupling as a loop stabilising mechanism*.

The model builds a two dimensional pial arterial network, lets vessel radii adapt to
wall shear stress until they reach a steady state, and compares four coupling variants.
The structural readout is the number of surviving collateral loops (cycle rank of the
pial subgraph). Anatomical inputs (penetrator radii, lengths, areal density) come from
the H01 human cortex reconstruction and are stored in `penetrators.csv`.

## Requirements

MATLAB (tested with R2023+). No toolboxes beyond base MATLAB are required; the paired
significance test in the geometry ensemble uses `signrank` (Statistics Toolbox). All
source is under `src/`; each script adds its own folder to the path, so from MATLAB
`cd src` and run any script by name (e.g. `RunGeometryEnsemble`).

## The four coupling models

| Name | Arterial EC coupling | Capillary current | Flow-regulated capillary |
|------|:---:|:---:|:---:|
| `NoCoupling` | no | no | no |
| `ArtCoupling` | yes | no | no |
| `FullCoupling` | yes | fixed | no |
| `FullCouplingFlowReg` | yes | yes | yes (flow-deficit driven) |

## Canonical network

`ncombsx = 2`, `branchlevel = 4` (so `ncombsy = 4`, eight hexagons, about 83 mm2),
penetrator density 3.35 /mm2 giving roughly 279 penetrators and 22 collateral loops.
Operating pressure 86 mmHg, venous pressure 15 mmHg.

## Which script produces which result

Each script is self-contained: it builds the network, calibrates the flow-regulation
target `Qreg` from a NoCoupling steady state, runs the models, prints a table and saves
a `.mat` plus a figure.

| Script | Thesis section | Output files |
|--------|----------------|--------------|
| `RunGeometryEnsemble.m` | 4.1 loop survival, Table 4.1 | `FinalResults/GeometryEnsemble_Results.mat`, `.png` |
| `ExpFlow.m` | 4.2 perfusion and occlusion, Tables 4.2 to 4.3 | `FinalResults/Exp_Flow_Results.mat`, `.png` |
| `TestAutoregulationCurve.m` | 4.3 autoregulation, Table 4.3 | `AutoregulationCurve_Results.mat`, `.png` |
| `PlotAI.m` | 4.3 (plot only) | `AI_plot.png` (reads the .mat above) |
| `TestFlowDeficitNearP0.m` | 4.3.1 fine pressure sweep, Table 4.4 | `FlowDeficitNearP0_Results.mat`, `.png` |
| `SensitivityLambda.m` | 4.4 coupling range, Table 4.5 | `FinalResults/SensitivityLambda_Results.mat`, `.png` |
| `SensitivityLengthLambda.m` | 4.4 length vs lambda | `SensitivityLengthLambda.png` |
| `SensitivityKchannel.m` | 4.5 K+ channel gain, Table 4.6 | `FinalResults/SensitivityKchannel_Results.mat`, `.png` |
| `MultistabilityDemo.m` | 4.6 multistability, Table 4.7 | `FinalResults/MultistabilityDemo_Results.mat`, `.png` |
| `TestPhysiologicalPenetrators.m` | 3.2 penetrator calibration | `PhysiologicalPenetrators_Results.mat` |

Only hard dependency: run `TestAutoregulationCurve.m` before `PlotAI.m`, which reads its
saved `.mat`.

## Note on the flow-sensing gain

The structural and occlusion experiments (`sens_ss.m`, used by the ensemble, ExpFlow,
sensitivity and multistability scripts) run with `kSensQ = 0.6`. The autoregulation
scripts (`TestAutoregulationCurve.m`, `TestFlowDeficitNearP0.m`) set `KSENSQ = 3.0`, the
value that maximises loop survival. This difference is discussed in the thesis
limitations.

## Core files (shared by every script)

Geometry: `defaultPars.m`, `leptomeningeal_geometry.m`, `honeycomb_geometry.m`,
`tree_geometry.m`, `jointrees.m`, `MakeNodeTable.m`, `LengthFromPosition.m`.

Penetrators: `load_penetrator_data.m`, `sample_penetrators.m`,
`integrate_sampled_penetrators.m`, `penetrators.csv`.

Physics: `calcConductance2025.m` (Poiseuille conductance), `solvehemodyn2025.m`
(Kirchhoff solve, used for both the haemodynamic and the electrical network),
`calcWSS2025.m` (wall shear stress), `rdotfunCoupling.m` (the adaptation ODE),
`vout.m`, `sens_build.m`, `sens_ss.m`, `sens_cycle_rank.m`, and the electrical model in
`ElectricModel/` (`ElectricEM.m`, `ElectricConductance.m`, `MakeElectricCircuit.m`).

## Governing equations

Hydraulic conductance `G = pi r^4 / (8 mu L)`. Wall shear stress `tau = r (P1 - P2) / (2 L)`.
Shear-gated potassium conductance `gK = ksens |tau| + gK0`, membrane potential
`Em = (gNa ENa + gLeak ELeak + gK EK) / gm`. Adaptation
`dr/dt = r kReg (Vref - gR (r - r*) - Vm)`. Flow deficit
`fd = max((Qreg - Q)/Qreg, 0)`, capillary conductance `gKQ = kSensQ fd`. Steady state is
declared when `max|dr/dt| < 1e-12`.

## redundant/

Superseded scripts, earlier experiments and draft notes were moved to `redundant/` and
are not needed to reproduce any result in the thesis.
