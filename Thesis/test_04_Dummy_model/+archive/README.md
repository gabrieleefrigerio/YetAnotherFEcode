# Archivio di test_04_Dummy_model

Codice non piu' in uso, conservato per riferimento. Tutto e' anche in git
(commit `160b989`, "Snapshot pre-refactor"), quindi nulla e' andato perso.

## Perche' la cartella si chiama `+archive`

Il `+` la rende un **package MATLAB**, e i package sono esclusi da `genpath`.
E' una scelta deliberata: dentro `V2_snapshot/` ci sono copie vecchie di
`RomMC.m`, `RomMCB.m`, `RomMN.m`, `RomRubin.m` e `RomCB.m`, che hanno lo
stesso nome dei file vivi in `Src/`. Con una cartella normale, un
`addpath(genpath(...))` avrebbe messo entrambe le versioni sul path e quale
delle due vince sarebbe dipeso dall'ordine del path: si sarebbe potuta usare
per sbaglio la versione di luglio invece di quella corrente. Con il `+`
questo non puo' succedere.

Per lo stesso motivo, **non** rinominare questa cartella in `archive` o
`_archive` senza prima aver risolto la duplicazione.

## Contenuto

| cartella       | cosa contiene |
|----------------|---------------|
| `scripts/`     | main, test e script di plot superati (incluso `test_04_main_V2.m` e `test_04_postProcessing_V2.m`, confluiti nei main unificati) |
| `Src/`         | classi superate: `AbaqusStructure_V2`, `TransientSolverOde_V2`, i solutori Newmark/Leapfrog, i `residual_*`, script scratch |
| `meshes/`      | mesh Abaqus non piu' usate (`.inp` V1, V2, V3). Il modello corrente e' `Src/DummyStructureAbaqus_V4.inp` |
| `Old/`         | cartella `Old/` preesistente, lasciata com'era |
| `V2_snapshot/` | ex `Thesis/V2/`: copia congelata al 2-12 luglio 2026 di file poi evoluti in `Src/`. Era sul path insieme agli originali |

## Come recuperare un file

```bash
git mv "Thesis/test_04_Dummy_model/+archive/scripts/nome.m" Thesis/test_04_Dummy_model/
```

Se lo si rimette in `Src/`, controllare prima che non abbia lo stesso nome di
un file gia' presente.
