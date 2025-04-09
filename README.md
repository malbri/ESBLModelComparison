# 🧪 ESBLModelComparison
# Model Comparison: Compartment vs. Source Attribution

This repository provides model codes for the paper:

Brinch, M.L.; Duarte, A.S.R.; Apenteng, O.O.; Hald, T. Modeling the Transmission of ESBL and AmpC-Producing Escherichia coli in Denmark: A Compartmental and Source Attribution Approach. Zoonotic Dis. 2025, 5, 7. https://doi.org/10.3390/zoonoticdis5010007

- **SIR Compartment Model**
- **Source Attribution Model Colonised Individuals**
- **Source Attribution Model Infected Individuals**

These models are used to explore infection dynamics and antimicrobial resistance transmission, with a focus on evaluating intervention strategies.

---

## 📁 Repository Structure
### `/SIRModel`
Contains all files related to the deterministic SIR (Susceptible-Infected-Recovered) compartment model.

**Files:**
- `ESBL_opti_det.R` — Main model script, optimising paramerters
- `Deterministic.RData` — Output from deterministic simulations
- `Optim.RData` — Results from model optimization
- `Scenarios.R` — Scripts for running intervention scenarios
- `Scenarios.RData` — Saved results from scenario runs
- `Sensitivity_analysis.R` — R script to perform sensitivity analysis

### `/ColonisedModel` (OC Model)
Source attribution model estimating the share of colonised individuals attributed to different sources (e.g., food and pets).

**Files:**
- `ESBL_OC.R` — Main R script for the colonised model
- `ESBL_OC` — Model text for Rjags
- `gene_list.xlsx` — Gene data used as model input
- `model_output.xlsx` — Final model output summary
- `obs_pred_counts_OC.xlsx` — Observed vs. predicted counts
- `overall_prev.xlsx` — Overall prevalence data (input)
---

### `/InfectedModel` (CL Model) 
Source attribution model estimating the share of infected individuals attributed to different sources (e.g., food and pets).

**Files:**
- `ESBL_cli.R` — Main R script for the colonised model
- `ESBL_cli` — Model text for Rjags
- `gene_list.xlsx` — Gene data used as model input
- `model_output.xlsx` — Final model output summary
- `obs_pred_counts_CL.xlsx` — Observed vs. predicted counts
- `overall_prev.xlsx` — Overall prevalence data (input)
---

## 📊 Data
Data for the source attribution models is available from the supplementary in the published paper. Following format should be used:

![image](https://github.com/user-attachments/assets/911074d9-f587-4a26-846d-272632df7deb)

