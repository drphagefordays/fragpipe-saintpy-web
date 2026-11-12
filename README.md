# FragPipe → saintPy Studio

**100% Client-Side AP-MS Interactome Scoring & Publication Visualization Suite**

🔗 **Live Web Application (GitHub Pages):** [https://drphagefordays.github.io/fragpipe-saintpy-web/](https://drphagefordays.github.io/fragpipe-saintpy-web/)

---

## Overview

**FragPipe → saintPy Studio** brings the statistical models of [FossatiLab/saintPy](https://github.com/FossatiLab/saintPy) (Significance Analysis of INTeractome) directly into the browser with zero server uploads.

### Key Capabilities

- **Direct FragPipe & SAINT Ingestion:**
  - Drag-and-drop FragPipe `combined_protein.tsv`, `protein.tsv`, `reprint.spc.tsv`, or `reprint.int.tsv`, plus optional `experiment_annotation.tsv` / `fp-manifest`.
  - Also supports classic 3-file (`bait.csv`, `inter.csv`, `prey.csv`) SAINT / saintPy tables.
- **Interactive Experimental Cohort & Prey Hygiene Designer:**
  - Automatically infers Test Baits vs. Negative Controls and verifies bait self-pulldown abundance.
  - Filter out contaminants (`contam_`), reverse decoys (`rev_`), low-probability proteins, and single-peptide hits interactively.
- **Exact In-Browser Statistical Engines:**
  - **SPC Model:** Generalized Poisson Markov Random Field (GP-MRF) with ICM optimization.
  - **MS1 Intensity Model:** Log-Intensity Gaussian MRF with missing-value background estimation.
  - Supports both native Float64 JavaScript execution (instant, offline-capable) and in-browser **Pyodide WebAssembly** (`Python 3.12 + NumPy/SciPy`).
- **Interactive Publication Visualization Suite (HiDPI Canvas + Vector SVG / PNG Export):**
  1. **ProHits-Style Bait × Prey DotPlot Matrix** (color = `SaintScore`, size = relative abundance, ring = empirical `BFDR` tier)
  2. **Volcano & Significance Scatter Plot** (with non-overlapping gene callout pills & enrichment zone)
  3. **Force-Directed Bait–Prey Interactome Network** (hexagonal bait hubs, double-ringed multi-bait complex hub preys, interactive drag)
  4. **Replicate-Level Breakdown & SAINT Posterior Fit Inspector** ($\eta_i$, $\eta_i + d_i$, replicate $Z_{ijr}$ states)
  5. **Global Score Bimodal Distribution & ICM Log-Likelihood Convergence Trace**
