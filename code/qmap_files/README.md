# QMAP files to create the processed maf files from the variant calling files form sarek and oncoanalyser ppipelines

The notebook `qmap_process_vcf_files_from_sarek_oncoanalyser.ipynb` has the necessary code to build the 14 qmap files that generate the processed maf files in
`/mafs/maf_sjd/`, `/mafs/maf_target/`, `/mafs/maf_stjude/` and `/mafs/maf_pmc/`. The complete set of files can be also downloaded from EGA repository (see main README in this repo).

QMAP files are a type of job parallelisation files to process multiple runs in a cluster environment. See https://github.com/bbglab/qmap.git for more information.
