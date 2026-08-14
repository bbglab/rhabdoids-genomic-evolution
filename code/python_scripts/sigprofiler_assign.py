from SigProfilerAssignment import Analyzer as Analyze
import click
import os
import json
import pandas as pd

@click.command

@click.option('--samples',
				'-s',
				required=True,
				help='Path to the input somatic mutations file (matrix/vcfs)')
@click.option('--output',
				'-o',
				required=False,
				default='results',
				help='output folder (if not present it will create a new one)')
@click.option('--input_type',
	          '-it',
	          required=False,
	          default='matrix',
	          help='type of input data, e.g. matrix')
@click.option('--context_type',
	          '-ct',
	          required=False,
	          default='96',
	          help='type of context, e.g. 96, 288, etc.')
@click.option('--cosmic_version',
	          '-cv',
	          required=False,
	          default=3.4,
	          type=float,
	          help='Cosmic version (default is 3.4)')
@click.option('--genome',
				'-g',
				required=False,
				default='GRCh38',
				help='Reference genome, e.g. GRCh38, GRCh37')
@click.option('--signature_database',
	          '-sig_db',
	          required=False,
	          default=None,
	          type=str,
	          help='Path to signatures database file')
@click.option('--exclude_signature_subgroups',
	          '-exclude_sigs',
	          required=False,
	          default=None,
	          type=str,
	          help='Subgroup of signatures to remove e.g. MMR_deficiency_signatures')
@click.option('--export_probabilities',
	          '-export_prob',
	          required=False,
	          default=True,
	          type=bool,
	          help='Probability matrix')
@click.option('--export_probabilities_per_mutation',
	          '-export_prob_per_mut',
	          required=False,
	          default=True,
	          type=bool,
	          help='Probability matrices per mutation')
@click.option('--make_plots',
	          '-plots',
	          required=False,
	          default=True,
	          type=bool,
	          help='Generate plots')
@click.option('--sample_reconstruction_plots',
	          '-sr_plots',
	          required=False,
	          default=True,
	          type=bool,
	          help='Select the output format for sample reconstruction plots. Valid inputs are {pdf, png, both, none}.')
@click.option('--verbose',
	          '-v',
	          required=False,
	          default=False,
	          type=bool,
	          help='Prints detailed statements')

def run_sigprofiler_assignment(samples,output,input_type,context_type,
                               cosmic_version,signature_database,genome,
                    exclude_signature_subgroups,export_probabilities,export_probabilities_per_mutation,
                    make_plots,sample_reconstruction_plots, verbose):
	'''
	Function to run sigProfiler assignment adapted to run command line.
	'''
			
	Analyze.cosmic_fit(samples, output, input_type=input_type, context_type=context_type,
				collapse_to_SBS96=True, cosmic_version=cosmic_version, exome=False,
				genome_build=genome, signature_database=signature_database,
				exclude_signature_subgroups=exclude_signature_subgroups, export_probabilities=export_probabilities,
				export_probabilities_per_mutation=export_probabilities_per_mutation, make_plots=make_plots,
				sample_reconstruction_plots=sample_reconstruction_plots, verbose=verbose)	

if __name__ == "__main__":
	run_sigprofiler_assignment()