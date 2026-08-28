# How to run the variant calling pipelines

You can run directly the variant calling steps with the CRAM files downloaded from EGAXXX (SJD and PMC samples).
TARGET program and StJude Cloud samples need to be converted into FASTQ files and run the alignment with the same reference genome.

## Extract FASTQ files from CRAM files:

Convert to bam if needed:

```
samtools view -b -T reference.fa -o output.bam input.cram 
```
Sort bams before extracting FASTQ files:
```
samtools sort -n path/to/file.bam -o path/to/file.sorted.bam' 
```
Extract FASTQ files:
```
samtools fastq -@ 32 path/to/file.sorted.bam -1 /path/to/file.1.fastq -2 /path/to/file.2.fastq -0 /dev/null -s /dev/null -n'
```


## Run alignment with sarek:

```
nextflow run nf-core/sarek -c path/to/nextflow.conf --input <input_file_name>.csv --outdir . -profile singularity --fasta /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa --fasta_fai /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa.fai --dict /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict --germline_resource /path/to/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz --germline_resource_tbi /path/to/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz.tbi --dbsnp /path/to/Homo_sapiens_assembly38.no_alt.dbsnp138.vcf.gz --dbsnp_tbi /data/bbg/datasets/genomes/GRCh38/Homo_sapiens_assembly38.no_alt.dbsnp138.vcf.gz.tbi --known_indels /path/to/{Homo_sapiens_assembly38.no_alt.known_indels,Mills_and_1000G_gold_standard.indels.hg38.no_alt}.vcf.gz --known_indels_tbi /path/to/{Homo_sapiens_assembly38.no_alt.known_indels,Mills_and_1000G_gold_standard.indels.hg38.no_alt}.vcf.gz.tbi --igenomes_ignore
``` 

## Run variant calling with sarek:

```
nextflow run nf-core/sarek -c path/to/nextflow.conf --input csv/recalibrated.csv --outdir . -profile singularity --fasta /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa --fasta_fai /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa.fai --dict /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict --germline_resource /data/bbg/datasets/genomes/iGenomes/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz --germline_resource_tbi /path/to/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz.tbi --dbsnp /path/to/Homo_sapiens_assembly38.no_alt.dbsnp138.vcf.gz --dbsnp_tbi /path/to/Homo_sapiens_assembly38.no_alt.dbsnp138.vcf.gz.tbi --known_indels /path/to/{Homo_sapiens_assembly38.no_alt.known_indels,Mills_and_1000G_gold_standard.indels.hg38.no_alt}.vcf.gz --known_indels_tbi /path/to/{Homo_sapiens_assembly38.no_alt.known_indels,Mills_and_1000G_gold_standard.indels.hg38.no_alt}.vcf.gz.tbi --save_output_as_bam --step variant_calling --tools 'haplotypecaller,ascat,mutect2,strelka,manta'
``` 

## Run variant calling with oncoanalyser:

``` 
nextflow run nf-core/oncoanalyser -profile singularity -revision 0.2.0 -c path/to/nextflow.conf --mode wgts --genome GRCh38_hmf --input path/to/input.csv --outdir output/ --publish_dir_mode copy
``` 
