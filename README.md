# rhabdoids-genomic-evolution
This repository contains the necessary code to reproduce the analysis of "Variable latency between the founder genetic event and rhabdoid tumor expansion" Sánchez-Guixé M. et al 2026 (doi: https://doi.org/10.64898/2026.06.30.735306).

## How to run the code from RAW files:

1. Prepare variant calling data folders:
   - Request access to the CRAM files:
      - SJD and PMC rhabdoid samples: EGA repository EGADXXXXXX
      - TARGET rhabdoid samples: dbGaP Study Accession phs000470.v19.p8
      - StJude rhabdoid samples: StJude Cloud
   - Generate variant calling folders following the instructions on the README in `.variant_calling/`
     - Point the output of [sarek pipeline](https://github.com/nf-core/sarek) in the following folders:
       - `./variant_calling/sjd_cohort/sarek_results/sjd_variant_calling/`
       - `./variant_calling/target_cohort/sarek_results/sjd_variant_calling/`
       - `./variant_calling/stjude_cohort/sarek_results/sjd_variant_calling/`
       - `./variant_calling/pmc_case/sarek_results/sjd_variant_calling/`
     - Point the output of [oncoanalyser pipeline](https://github.com/nf-core/oncoanalyser) in the following folders:
       - `./variant_calling/sjd_cohort/oncoanalyser_results/outpout/`
       - `./variant_calling/target_cohort/oncoanalyser_results/output/`
       - `./variant_calling/stjude_cohort/oncoanalyser_results/output/`
       - `./variant_calling/pmc_case/oncoanalyser_results/output/`
   -  Generate the processed mafs folders following the instructions on the README in `./code/qmap_files/`
2. Run figure notebooks:
   - Run the figure notebooks in order (01, 02, 03, ..., 14)
   - The necessary files are stored at the folders in this repository (`./ccf_thresholds`, `./clinical_data`, etc.)
   - The manuscript figure plots will be stored at `./plots/`.
