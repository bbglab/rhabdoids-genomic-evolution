import pandas as pd
from io import StringIO
import gzip
from bgreference import hg19
from aux_data_in_pyvar import CHANNELS

def fix_indels (rw):    
    ref = rw['ref_original']
    alt = rw['alt_original']
    pos = rw['pos_original']
    
    if len(ref)>len(alt): #deletion
        while alt != '-':

            if ref[0] == alt[0]:
                ref = ref[1:]
                alt = alt[1:]
                pos = pos + 1
                if alt == '':
                    alt = '-'
            elif ref[len(ref)-1] == alt[len(alt)-1]:
                ref = ref[0:len(ref)-1]
                alt = alt[0:len(alt)-1]
                pos = pos
                if alt == '':
                    alt = '-'
            elif alt == '*':
                ref = ref
                alt = '-'
                pos = pos
            else:
                print('Another type of deletion?')
                print(ref,alt)
        rw['REF'] = ref
        rw['ALT'] = alt
        rw['POS'] = pos
                
    elif len(ref)<len(alt): #insertion
        while ref != '-':

            if ref[0] == alt[0]:
                ref = ref[1:]
                alt = alt[1:]
                pos = pos + 1
                if ref == '':
                    ref = '-'
            elif ref[len(ref)-1] == alt[len(alt)-1]:
                ref = ref[0:len(ref)-1]
                alt = alt[0:len(alt)-1]
                pos = pos
                if ref == '':
                    ref = '-'
            else:
                print('Another type of insertion?')
                print(ref,alt)
        rw['REF'] = ref
        rw['ALT'] = alt
        rw['POS'] = pos
            
    elif len(ref) == len(alt):
        rw['REF'] = ref
        rw['ALT'] = alt
        rw['POS'] = pos

    else:
        print('Another type of variant?')
        print(ref,alt,pos)
    
    return rw

def calculate_AF (alt,ref):
    if alt+ref == 0:
        af = 0
    else:
        af = alt / (alt+ref)
    return af


def read_vcf(filename, comment='##', sep='\t'):
    """
    #     VCF has a long header starting with ##. The function reads the VCF and converts it into a Pandas dataframe
    #     :param filename: input VCF file
    #     :param comment: ## is the default but it can be any header preceding data in dataframe form
    #     :param sep: tabulated space is the default but it can be anything
    #     :return: Pandas DataFrame
    #     """
    if filename.endswith('vcf.gz') or filename.endswith('txt.gz'):
        lines = ''.join([line for line in gzip.open(filename, 'rt') if not line.startswith(comment)])
    else:
        lines = ''.join([line for line in open(filename) if not line.startswith(comment)])
    try:
        return pd.read_csv(StringIO(lines), sep=sep)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

def get_reads_snvs(rw, status, var=None):
    """
    The name of the function defines the purpose of this function
    :param rw: row of the VCF with the variant (SNVs only)
    :param status: TUMOR or NORMAL columns with the read information
    :param var: ALT or REF depending on the info that one wants to retrieve
    :return: the same row pandas Series object but with added columns
    """
    '''
    rw: row and therefore variant of the dataframe
    '''
    if status == 'TUMOR':
        rw['DP_tumor'] = rw[status].split(':')[0]
    elif status == 'NORMAL':
        rw['DP_normal'] = rw[status].split(':')[0]
    else:
        print('Not a valid status, only TUMOR or NORMAL')

    if rw[var] == 'A':
        rw['reads'] = rw[status].split(':')[-4:-3][0].split(',')[0]
    elif rw[var] == 'C':
        rw['reads'] = rw[status].split(':')[-3:-2][0].split(',')[0]
    elif rw[var] == 'G':
        rw['reads'] = rw[status].split(':')[-2:-1][0].split(',')[0]
    elif rw[var] == 'T':
        rw['reads'] = rw[status].split(':')[-1:][0].split(',')[0]
    else:
        print('You should indicate whether is REF or ALT that you are referring at')
    return rw


def genotype_strelka(rw):
    """
    Gets the genotype estimates of each variant by Strelka. It is different in SNV than InDels
    :param rw: row of the VCF variant
    :return: row with the columns with the genotype information
    """
    if rw['mut_type'] == 'snv':
        info = rw['INFO'].replace('SOMATIC;', '')
        info = info.replace(';OVERLAP', '')
        info_things = dict(item.split("=") for item in info.split(";"))
        rw['GT'] = info_things['SGT']

        #rw['GT'] = rw['INFO'].split(';')[6].split('=')[1]
        if rw['GT'] == str(rw['REF']+rw['REF']+'->'+rw['REF']+rw['ALT']):
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '0/1'
        elif rw['GT'] == str(rw['REF']+rw['REF']+'->'+rw['ALT']+rw['REF']):
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '1/0'
        elif rw['GT'] == str(rw['REF']+rw['REF']+'->'+rw['ALT']+rw['ALT']):
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '1/1'
        else:
            print(rw['GT'])
    else:

        info = rw['INFO'].replace('SOMATIC;', '')
        info = info.replace(';OVERLAP', '')
        info_things = dict(item.split("=") for item in info.split(";"))
        rw['GT'] = info_things['SGT']

        #rw['GT'] = rw['INFO'].split(';')[6].split('=')[1]
        if rw['GT'] == 'ref->het':
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '0/1'
        elif rw['GT'] == 'ref->hom':
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '1/1'
        elif rw['GT'] == 'het->het':
            rw['GT_normal'] = '0/1'
            rw['GT_tumor'] = '1/0'
        elif rw['GT'] == 'ref->ref':
            rw['GT_normal'] = '0/1'
            rw['GT_tumor'] = '0/1'
        elif rw['GT'] == 'hom->hom':
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '1/1'
        else:
            print(rw['GT'])
    rw = rw.drop('GT')
    return rw


def snvs_processing(rw):
    """
    process SNV lines from VCF
    :param rw: row of the tab delimited file (VCF)
    :return: same rw with more columns
    """
    rw['mut_type'] = 'snv'

    # Retrieve read info
    rw = get_reads_snvs(rw, 'TUMOR', 'ALT')
    rw.rename({'reads':'t_alt_reads'}, inplace=True)

    rw = get_reads_snvs(rw, 'TUMOR', 'REF')
    rw.rename({'reads':'t_ref_reads'}, inplace=True)

    rw = get_reads_snvs(rw, 'NORMAL', 'ALT')
    rw.rename({'reads':'n_alt_reads'}, inplace=True)

    rw = get_reads_snvs(rw, 'NORMAL', 'REF')
    rw.rename({'reads':'n_ref_reads'}, inplace=True)

    # Retrieve estimated genotype
    rw = genotype_strelka(rw)

    return rw


def reads_strelka_indels(rw, status):
    """
    The name of the function defines the purpose of this function.
    :param rw: row of the VCF with the variant (InDels only)
    :param status: TUMOR or NORMAL columns with the read information
    :return: the same row pandas Series object but with added columns
    """
    if status == 'TUMOR':
        rw['t_ref_reads'] = rw[status].split(':')[2:3][0].split(',')[0] # [2:3] corresponds to TAR
        rw['t_alt_reads'] = rw[status].split(':')[3:4][0].split(',')[0] # [3:4] corresponds to TIR
        rw['DP_tumor'] = rw[status].split(':')[0]
    elif status == 'NORMAL':
        rw['n_ref_reads'] = rw[status].split(':')[2:3][0].split(',')[0] # [2:3] corresponds to TAR
        rw['n_alt_reads'] = rw[status].split(':')[3:4][0].split(',')[0] # [3:4] corresponds to TIR
        rw['DP_normal'] = rw[status].split(':')[0]
    else:
        print('Not a valid status, only TUMOR or NORMAL')
    return rw


def indels_processing(rw):
    """
    process InDels lines from VCF
    :param rw: row of the tab delimited file (VCF)
    :return: same rw with more columns
    """
    rw['mut_type'] = "indel"

    # Retrieve read info
    rw = reads_strelka_indels(rw, 'TUMOR')
    rw = reads_strelka_indels(rw, 'NORMAL')

    # Retrieve estimated genotype
    rw = genotype_strelka(rw)

    return rw


# GET TRINUCLEOTIDE CONTEXT OF SNVS IN PYRIMIDE REF

def get_context_rev(rw):
    equival_nt = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    pos = rw['POS']
    left, ref, right = hg19(rw['#CHROM'], pos - 1, size=3)
    print(left)
    print(ref)
    print(right)
    alt = rw['ALT']
    rw['TRIPLE'] = left + '[' + ref + '>' + alt + ']' + right
    rw['TRIPLE_COM'] = equival_nt[left] + '[' + equival_nt[ref] + '>' + equival_nt[alt] + ']' + equival_nt[right]
    rw['TRIPLE_COM_REV'] = equival_nt[right] + '[' + equival_nt[ref] + '>' + equival_nt[alt] + ']' + equival_nt[left]

    return rw


def add_pyrimidine_type(df):

    dff = pd.DataFrame()

    for i, rw in df.iterrows():
        if rw['TRIPLE'] in CHANNELS:
            rw['VARIANT_TYPE'] = rw['TRIPLE']
            df_rw_T = rw.to_frame().T
            dff = dff.append(df_rw_T, ignore_index=True)
        else:
            rw['VARIANT_TYPE'] = rw['TRIPLE_COM_REV']
            df_rw_T = rw.to_frame().T
            dff = dff.append(df_rw_T, ignore_index=True)
    dff.drop(columns=['TRIPLE','TRIPLE_COM', 'TRIPLE_COM_REV'], inplace=True)
    return dff


def count_variant_type(df):

    count = {ctxt: 0 for ctxt in CHANNELS}

    for i, rw in df.iterrows():
        count[rw['VARIANT_TYPE']] = count[rw['VARIANT_TYPE']] + 1

    return count


def df_to_dict(df):
    df.rename(columns={'CHROMOSOME': '#CHROM', 'POSITION': 'POS'}, inplace=True)
    df = df[['#CHROM', 'POS', 'REF', 'ALT']]
    df.drop_duplicates(keep='first', inplace=True)

    # make dictionary with counts centered in pyrimidines
    df = df.apply(lambda x: get_context_rev(x), axis=1)
    df = add_pyrimidine_type(df)
    count = count_variant_type(df)

    # from dict to pandas df
    df_96 = pd.DataFrame.from_dict(count, orient='index')
    df_96.reset_index(inplace=True)
    df_96.columns = ['change', 'count']

    # change to dictionary with special key annotation
    df_96[['count']] = df_96[['count']].astype(int)
    #     total = df_96['count'].sum()
    #     df_96['relative_count'] = df_96['count'].apply(lambda x: x/total)
    df_96 = df_96[['count', 'change']]
    df_96.set_index('change', inplace=True)
    df_96 = df_96.sort_index()
    dictionary = df_96.to_dict()['count']

    # fill missing context with 0 counts
    for k in CHANNELS:
        if k not in dictionary.keys():
            dictionary[k] = 0
    return dictionary

