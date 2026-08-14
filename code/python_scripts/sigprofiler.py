from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
import click
import os
import json
import pandas as pd

@click.command

@click.option('--project_name',
				'-p',
				required=True,
				help='Project name, e.g. "allmuts"')
@click.option('--genome',
				'-g',
				required=False,
				default='GRCh38',
				help='Reference genome, e.g. GRCh38, GRCh37')
@click.option('--path_to_vcf_files',
				'-vcf',
				required=True,
				help='Path to folder containing vcf files')
@click.option('--scaled',
				'-s',
				is_flag=True,
				required=False,
				default=False,
				help='Whether the generated matrix need to be scaled before running sigprofiler')
@click.option('--cosmic_counts',
				'-cc',
				required=False,
				default=None,
				help='Path to the csv file containing the triplet context counts of the genome used at COSMIC')
@click.option('--target_genome_counts',
				'-tc',
				required=False,
				default=None,
				help='Path to the json file containing the triplet context counts of the target genome used to obtain the mutations')
@click.option('--round_error_correction',
				'-re',
				required=False,
				default=0.0,
				type=float,
				help='Number to substract to the total number of counts in each sample to correct de rounding to 0 error distance. ')
@click.option('--input_type',
	          '-it',
	          required=False,
	          default='matrix',
	          help='type of input data, e.g. matrix')
@click.option('--output',
				'-o',
				required=False,
				default='results',
				help='output folder (if not present it will create a new one)')
@click.option('--min_sigs',
				'-min',
				required=False,
				default=2,
				type=int,
				help='Number of minimum signatures to calculate')
@click.option('--max_sigs',
	          '-max',
	          required=False,
	          default=10,
	          type=int,
	          help='Number of maximum signatures to calculate')
@click.option('--nmf_replicates',
	          '-nmf',
	          required=False,
	          default=100,
	          type=int,
	          help='Number of NMF replicates to use')
@click.option('--cpu',
	          '-c',
	          required=False,
	          default=10,
	          type=int,
	          help='Number of CPUs to use')

def run_sigprofiler(project_name,genome,path_to_vcf_files,scaled,cosmic_counts,target_genome_counts,round_error_correction,input_type,output,min_sigs,max_sigs,nmf_replicates,cpu):
	'''
	Function to run sigProfiler, adapted to python.
	Some values are harcoded.
	'''
	if os.path.exists(path_to_vcf_files):
			matrices = matGen.SigProfilerMatrixGeneratorFunc(project=project_name, reference_genome=genome, 
			path_to_input_files=path_to_vcf_files, exome=False, bed_file=None, chrom_based=False, plot=False, 
			tsb_stat=False, seqInfo=False)

			input_data = 'output/SBS/'+project_name+'.SBS96.all'

			if scaled:
				print('Scaling the matrix to the cosmic trinucleotide contexts...')
				#make dictionaries with the context counts
				target_counts = json.load(open(target_genome_counts,'rb'))
				cosmic_counts_df = pd.read_csv(cosmic_counts)
				cosmic_counts_df.index = cosmic_counts_df[' ']
				cosmic_counts_df.drop([' '],inplace=True,axis=1)
				cosmic_counts_df['total'] = cosmic_counts_df.sum(axis=1)
				cosmic_counts = dict(zip(cosmic_counts_df.index,cosmic_counts_df['total']))

				#read the obtained matrix and create a new matrix scaled to the cosmic context counts
				matrix_df = pd.read_csv(input_data,sep='\t')
				matrix_df['context'] = matrix_df['MutationType'].apply(lambda x: x[0]+x[2]+x[6])
				matrix_df['cosmic_counts'] = matrix_df['context'].map(cosmic_counts)
				matrix_df['target_counts'] = matrix_df['context'].map(target_counts)
				cols = matrix_df.columns
				sample_cols = [col for col in cols if col not in ['MutationType','context','cosmic_counts','target_counts']]
				matrix_scaled_df = matrix_df[['MutationType']]
				matrix_scaled2_df = matrix_df[['MutationType']]
				for col in sample_cols:
				 	#transform counts to cosmic proportions
					matrix_scaled_df[col] = matrix_df[col] * matrix_df['cosmic_counts'] / matrix_df['target_counts']
					total_counts_sample_before = matrix_df[col].sum()
					total_counts_sample_after = matrix_scaled_df[col].sum()
					if total_counts_sample_before == 0:
						matrix_scaled2_df[col] = matrix_df[col]
					else:
						#scale to get the same total number of counts in the sameple (substract to correct for rounding error)
						matrix_scaled2_df[col] = round((matrix_scaled_df[col] * (total_counts_sample_before-float(round_error_correction)) / total_counts_sample_after),0).astype(int)


				input_data = input_data.replace('.all','.all.scaled')
				matrix_scaled2_df.to_csv(input_data,index=None,sep='\t')

			sig.sigProfilerExtractor(input_type=input_type, output=output, input_data=input_data,reference_genome=genome,opportunity_genome=genome,minimum_signatures=min_sigs,maximum_signatures=max_sigs,nmf_replicates=nmf_replicates, cpu=cpu)
	else:
		print('Please put a valid path for --path_to_vcf_files.')
if __name__ == "__main__":
	run_sigprofiler()