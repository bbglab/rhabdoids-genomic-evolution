# rhabdoids-genomic-evolution
This repository contains the necessary code to reproduce the analysis of "Variable latency between the founder genetic event and rhabdoid tumor expansion" Sánchez-Guixé M. et al 2026 (doi: https://doi.org/10.64898/2026.06.30.735306).

## How to run the code:

1. Prepare variant calling data folders:
   - Request access to the following datasets in EGA: EGADXXXXXX
   - Download the variant calling folders:
     - From [sarek pipeline](https://github.com/nf-core/sarek):
       - `sjd_variant_calling.tar.gz` and extract folder at `./variant_calling/sjd_cohort/sarek_results/sjd_variant_calling/`
       - `target_variant_calling.tar.gz` and extract folder at `./variant_calling/target_cohort/sarek_results/sjd_variant_calling/`
       - `stjude_variant_calling.tar.gz` and extract folder at `./variant_calling/stjude_cohort/sarek_results/sjd_variant_calling/`
       - `pmc_variant_calling.tar.gz` and extract folder at `./variant_calling/pmc_case/sarek_results/sjd_variant_calling/`
     - From [oncoanalyser pipeline](https://github.com/nf-core/oncoanalyser):
       - `sjd_oncoanalyser_results.tar.gz` and extract folder at `./variant_calling/sjd_cohort/oncoanalyser_results/`
       - `target_oncoanalyser_results.tar.gz` and extract folder at `./variant_calling/target_cohort/oncoanalyser_results/`
       - `stjude_oncoanalyser_results.tar.gz` and extract folder at `./variant_calling/stjude_cohort/oncoanalyser_results/`
       - `pmc_oncoanalyser_results.tar.gz` and extract folder at `./variant_calling/pmc_case/oncoanalyser_results/`
   -  Download the processed mafs folders:
      - `mafs_sjd.tar.gz` and extract folder at `./mafs/mafs_sjd/`
      - `mafs_target.tar.gz` and extract folder at `./mafs/mafs_target/`
      - `mafs_stjude.tar.gz` and extract folder at `./mafs/mafs_stjude/`
      - `mafs_pmc.tar.gz` and extract folder at `./mafs/mafs_pmc/`
2. Run figure notebooks:
   - Run the figure notebooks in order (01, 02, 03, ..., 14)
   - The necessary files are stored at the folders in this repository (`./ccf_thresholds`, `./clinical_data`, etc.)
   - The manuscript figure plots will be stored at `./plots/`.
