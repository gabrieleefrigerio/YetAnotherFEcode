# Archive of test_04_Dummy_model

Code no longer in use, kept for reference. Everything is also in git
(commit `160b989`, "Snapshot pre-refactor"), so nothing has been lost.

## Why the folder is called `+archive`

The `+` makes it a **MATLAB package**, and packages are excluded from
`genpath`. That is deliberate: `V2_snapshot/` contains old copies of
`RomMC.m`, `RomMCB.m`, `RomMN.m`, `RomRubin.m` and `RomCB.m`, which have the
same names as the live files in `Src/`. With a plain folder, an
`addpath(genpath(...))` would put both versions on the path and which one
wins would depend on the path order, so the July versions could silently be
used instead of the current ones. The `+` makes that impossible.

For the same reason, do **not** rename this folder to `archive` or
`_archive` without first resolving the duplication.

## Contents

| folder         | what it holds |
|----------------|---------------|
| `scripts/`     | superseded mains, tests and plotting scripts (including `test_04_main_V2.m` and `test_04_postProcessing_V2.m`, now merged into the unified mains) |
| `Src/`         | superseded classes: `AbaqusStructure_V2`, `TransientSolverOde_V2`, the Newmark/Leapfrog solvers, the `residual_*` functions, scratch scripts |
| `meshes/`      | Abaqus meshes no longer used (`.inp` V1, V2, V3). The current model is `Src/DummyStructureAbaqus_V4.inp` |
| `Old/`         | pre-existing `Old/` folder, left as it was |
| `V2_snapshot/` | former `Thesis/V2/`: a copy frozen between 2 and 12 July 2026 of files that later evolved in `Src/`. It used to sit on the path alongside the originals |

## Restoring a file

```bash
git mv "Thesis/test_04_Dummy_model/+archive/scripts/name.m" Thesis/test_04_Dummy_model/
```

If it goes back into `Src/`, first check that its name does not clash with a
file already there.
