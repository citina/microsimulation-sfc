# microsimulation-sfc

MATLAB code for the HIV microsimulation model used in:

> Liang C, Suen SC, Nguyen A, Moucheraud C, Hsu L, Holloway IW, Charlebois ED, Steward WT.
> **Impact of COVID-19 Response on the HIV Epidemic in Men Who Have Sex With Men in San Francisco County: The Importance of Rapid Return to Normalcy.**
> *Journal of Acquired Immune Deficiency Syndromes*. 2023;92(5):370–377.
> doi: [10.1097/QAI.0000000000003156](https://doi.org/10.1097/QAI.0000000000003156) · Free full text: [PMC9988211](https://pmc.ncbi.nlm.nih.gov/articles/PMC9988211/)

## What the model does

The model follows men who have sex with men (MSM) in San Francisco County and tracks HIV progression, diagnosis, PrEP and treatment. In the paper, it compares scenarios where HIV services disrupted by COVID-19 (testing, care engagement, PrEP uptake and retention) return to pre-COVID levels by the end of 2022 or by 2025, against a counterfactual in which the disruptions never happened. It also compares prioritizing new patients against retaining existing ones from 2023 to 2025.

## What's in this repository

- `create_init_files/`: builds the starting population of MSM in San Francisco County in 2012.
- `core_simulation/`: the simulation. `wrapper.m` sets up and runs the COVID-19 scenarios, which are defined in its header comments.
- `global_fun/`: helper functions shared by the other scripts.

## Inputs not included

The scripts read an input workbook (`Inputs_*.xlsx`), an initial population file (`init_pop_SF_V5.csv`) and transition tables (`input/transitions/*.csv`). These files are not in this repository, so the code needs them supplied before it can run. They are available on request.
