import pandas as pd
import pybedtools
from io import StringIO
import gzip
import tqdm as tqdm
from multiprocessing import Pool
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
import pickle
import allel

def overlap_tables (t1_pass_df,t2_pass_df,indicator,left,right):
    '''function to get overlapped table between tumors'''
    t1_pass_df['Mutations'] = t1_pass_df['#CHROM'].astype(str)+'_' + t1_pass_df['POS'].astype(str)+ '_' + t1_pass_df['REF'].astype(str)+ '_' + t1_pass_df['ALT'].astype(str)
    t1_pass_df = t1_pass_df[['Mutations']]
    
    t2_pass_df['Mutations'] = t2_pass_df['#CHROM'].astype(str)+'_' + t2_pass_df['POS'].astype(str)+ '_' + t2_pass_df['REF'].astype(str)+ '_' + t2_pass_df['ALT'].astype(str)
    t2_pass_df = t2_pass_df[['Mutations']]
    
    df = pd.merge(t1_pass_df, t2_pass_df, how='outer', indicator=indicator)
    df = df.replace({'left_only':left,'right_only':right})
    
    df = df.groupby(indicator,as_index=False).count()
    
    total = sum(df['Mutations'].tolist())
    
    df['Percent'] = round(df['Mutations']/total*100,2)
    
    return(df)

def venn_tumor_plots (t1_pass_df,t2_pass_df,title):
    '''function to get overlapped table between tumors'''
    t1_pass_df['Mutations'] = t1_pass_df['#CHROM'].astype(str)+'_' + t1_pass_df['POS'].astype(str)+ '_' + t1_pass_df['REF'].astype(str)+ '_' + t1_pass_df['ALT'].astype(str)
    t1_pass_list = t1_pass_df['Mutations'].tolist()
    
    t2_pass_df['Mutations'] = t2_pass_df['#CHROM'].astype(str)+'_' + t2_pass_df['POS'].astype(str)+ '_' + t2_pass_df['REF'].astype(str)+ '_' + t2_pass_df['ALT'].astype(str)
    t2_pass_list = t2_pass_df['Mutations'].tolist()
    
    set1 = set(t1_pass_list)
    set2 = set(t2_pass_list)
    
    plt.title(title, size=20)

    venn2([set1, set2], ('Tumor 1', 'Tumor 2'))
    plt.show()

mutations = ['missense_variant','stop_gained','frameshift_variant']

def boxplot_AF(df1,df2,tumor,sample,title):
    df = pd.concat([df1,df2],ignore_index=True)
    df = df[['Consequence','Feature']][(df['gnomADg_AF']==0)].groupby('Consequence',as_index=False).count()
    s = df['Consequence'].tolist()
    if tumor == 't1':
        data = df1[df1['gnomADg_AF']==0]
    elif tumor == 't2':
        data = df2[df2['gnomADg_AF']==0]
    y = 'Consequence'
    x = sample
    ax = sns.boxplot(x=x,y=y,data=data,order=s)
    ax.set_xlim(-0.1,1)
    plt.axvline(0, 1,0,color='grey',ls='--')
    plt.axvline(0.5, 1,0,color='grey',ls='--')
    ax.set_title(title,size=20,pad=20)
    
    
def cnv_snvs_intersect (cnv_df,snvs_df,caller_cn):
    '''Adds CN information in snvs/indels mafs
    input: cn dataframe from ASCAT output (.cnvs.txt) and snvs/indels dataframe from maf files
    output: dataframe with integrated CN info to vc mafs'''
    
    if caller_cn == 'ascat':
        cnv_bt_df = cnv_df.copy(deep=True)
        cnv_bt_df['CN'] = cnv_bt_df['nMajor']+cnv_bt_df['nMinor']
        cnv_bt_df['chr'] = 'chr' + cnv_bt_df['chr']
        cnv_bt_df = cnv_bt_df.rename(columns={'chr':'chrom','startpos':'start','endpos':'end'})
        cnv_bt_df = cnv_bt_df[['chrom','start','end','CN']]
    elif caller_cn == 'purple':
        cnv_bt_df = cnv_df.copy(deep=True)
        cnv_bt_df = cnv_bt_df.rename(columns={'chromosome':'chrom','copyNumber':'CN'})
        cnv_bt_df = cnv_bt_df[['chrom','start','end','CN']]
    #Prepare snvs_df
    snvs_bt_df = snvs_df.copy(deep=True)
    snvs_bt_df['POS2'] = snvs_bt_df['POS'] + 1
    snvs_bt_df = snvs_bt_df[['#CHROM','POS','POS2','REF','ALT']]

    #pybedtools cnv + snvs
    cnv_bt = pybedtools.BedTool.from_dataframe(cnv_bt_df)
    snvs_bt = pybedtools.BedTool.from_dataframe(snvs_bt_df)

    cnv_and_snvs = cnv_bt.intersect(snvs_bt)

    #Pass intersection to dataframe
    snvs_cnv_df = cnv_and_snvs.to_dataframe(header=None)
    columns= ['#CHROM','POS','POS2','CN'] 
    snvs_cnv_df = snvs_cnv_df.rename(columns={'chrom':'#CHROM','start':'POS','end':'POS2',
                                              'name':'CN'})
    snvs_cnv_df = snvs_cnv_df.drop(columns=['POS2'])

    #Merge intersect with snvs_df
    snvs_cnv_df = pd.merge(snvs_df,snvs_cnv_df,how='left')
    
    return snvs_cnv_df

def calculate_VAF(reads_alt,reads_ref):
    '''Calculate VAF from number of reads
    input: number of reads alternate allele, number of reads reference allele
    output: variant allele frequency'''
    if reads_alt == reads_ref == 0:
        return 0
    else:
        vaf = round (reads_alt/(reads_alt+reads_ref),3)
        return vaf

def calculate_ccf (vaf,cnv,purity):
    '''Calculate CCF from VAF, CN and purity
    input: variant allele frequency, copy number, purity(.purityploidy.txt file from ASCAT)
    output: cancer cell fraction'''
    ccf = vaf*(purity*cnv+ (1-purity)*2)/purity
    return ccf

DAMAGING_CONSEQ= ['transcript_ablation',
'splice_acceptor_variant',
'splice_donor_variant',
'stop_gained',
'frameshift_variant',
'stop_lost',
'start_lost',
'transcript_amplification',
'inframe_insertion',
'inframe_deletion',
'missense_variant',
'protein_altering_variant',
'splice_region_variant',
'incomplete_terminal_codon_variant',
'start_retained_variant',
'stop_retained_variant',
'synonymous_variant',
'coding_sequence_variant',
'mature_miRNA_variant',
'5_prime_UTR_variant',
'3_prime_UTR_variant',
'non_coding_transcript_exon_variant',
'intron_variant',
'NMD_transcript_variant',
'non_coding_transcript_variant',
'upstream_gene_variant',
'downstream_gene_variant',
'TFBS_ablation',
'TFBS_amplification',
'TF_binding_site_variant',
'regulatory_region_ablation',
'regulatory_region_amplification',
'feature_elongation',
'regulatory_region_variant',
'feature_truncation',
'intergenic_variant'
]

most_damaging = ['transcript_ablation',
'splice_acceptor_variant',
'splice_donor_variant',
'stop_gained',
'frameshift_variant',
'stop_lost',
'start_lost',
'transcript_amplification',
'inframe_insertion',
'inframe_deletion',
'missense_variant',
'protein_altering_variant',
'splice_region_variant',
'incomplete_terminal_codon_variant',
'start_retained_variant',
'stop_retained_variant']

def damaging_col(x):
    x = str(x)
    x_list = x.split(',')
    y = list(set(x_list) & set(most_damaging))
    if y != []:       
        return True
    else:
        return False
    
