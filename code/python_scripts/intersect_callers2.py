import pandas as pd
import click


def check_callers(df):
    list_callers = list()
    if df.loc[0, 'SAGE'] == True:
        list_callers.extend(['SAGE'])
    if df.loc[0, 'STRELKA'] == True:
        list_callers.extend(['STRELKA'])
    if df.loc[0, 'MUTECT'] == True:
        list_callers.extend(['MUTECT'])
    df['Callers_intersection'] = ",".join(list_callers)
    return (df)


@click.command()

@click.option('--maf_sage',
              '-sa',
              required = True,
              help="Input must be the path to the processed MAF file from SAGE calling")
@click.option('--maf_mutect',
              '-mu',
              required = True,
              help="Input must be the path to the processed MAF file from MuTect2 calling")
@click.option('--maf_strelka',
              '-st',
              required = True,
              help="Input must be the path to the processed MAF file from Strelka calling")
@click.option('--output_dire',
              '-o',
              required = True,
              help="Output MAF with variants in 2 out of the 3 callers")
@click.option('--sample_name',
              '-sn',
              required = True,
              help="Sample name to put in SAMPLE column. This will be as well the name of the file")
@click.option('--by_chrom/--single_file',
              '-c/-s',
              required = False,
              default = False,
              show_default=True,
              help="Whether the output should be in different files by chromosome. Default is False")


def cli(maf_sage, maf_mutect, maf_strelka, output_dire, sample_name, by_chrom):


    # READ RESULTS
    df_maf_sage = pd.read_csv(maf_sage, sep='\t')
    df_maf_mutect = pd.read_csv(maf_mutect, sep='\t')
    df_maf_strelka = pd.read_csv(maf_strelka, sep='\t')

    #Merge by variant and mark variant origin
    cols = ['#CHROM','POS','ID','REF','ALT']

    df = pd.merge(df_maf_sage[cols],df_maf_strelka[cols],on=cols,how='outer',indicator='SAGE')
    df['STRELKA'] = df['SAGE']
    df['SAGE'] = df['SAGE'].replace('both',True)
    df['SAGE'] = df['SAGE'].replace('left_only',True)
    df['SAGE'] = df['SAGE'].replace('right_only',False)
    df['STRELKA'] = df['STRELKA'].replace('both',True)
    df['STRELKA'] = df['STRELKA'].replace('left_only',False)
    df['STRELKA'] = df['STRELKA'].replace('right_only',True)    

    df = pd.merge(df,df_maf_mutect[cols],on=cols,how='outer',indicator='MUTECT')
    df['MUTECT'] = df['MUTECT'].replace('both',True)
    df['MUTECT'] = df['MUTECT'].replace('left_only',False)
    df['MUTECT'] = df['MUTECT'].replace('right_only',True)

    df[['SAGE','STRELKA','MUTECT']] = df[['SAGE','STRELKA','MUTECT']].fillna(False)

    df['num_callers'] = df.apply(lambda row: row['SAGE'] + row['STRELKA'] + row['MUTECT'],axis=1)

    #Select variants with at least 2 callers
    df = df[df['num_callers']>=2]

    df = df.drop('num_callers',axis=1)

    #Get columns from SAGE > STRELKA > MUTECT
    df_sage = pd.merge(df[df['SAGE']==True],df_maf_sage,how='left')
    df_strelka = pd.merge(df[(df['SAGE']==False)&(df['STRELKA']==True)],df_maf_strelka,how='left')
    dff = pd.concat([df_sage,df_strelka],ignore_index=True)

    def callers_intersection_col (row):
        if row['SAGE'] == True and row['STRELKA'] == True and row['MUTECT'] == True:
            return 'SAGE,STRELKA,MUTECT'
        elif row['SAGE'] == True and row['STRELKA'] == True and row['MUTECT'] == False:
            return 'SAGE,STRELKA'
        elif row['SAGE'] == True and row['STRELKA'] == False and row['MUTECT'] == True:
            return 'SAGE,MUTECT'
        elif row['SAGE'] == False and row['STRELKA'] == True and row['MUTECT'] == True:
            return 'STRELKA,MUTECT'

    dff['Callers_intersection'] = dff.apply(lambda row: callers_intersection_col(row),axis=1)

    dff.drop(columns=['SAGE', 'STRELKA', 'MUTECT'], inplace=True)

    
    # ADD SAMPLE COLUMN
    dff['SAMPLE'] = sample_name
    
    if by_chrom == True:
    
        # SAVE TABLES DIVIDED BY CHROMOSOME
        chrom_groups = dff.groupby("#CHROM")

        for chr_ in chrom_groups.groups:
            df_chr = chrom_groups.get_group(chr_)
            df_chr = df_chr.sort_values(by='POS',ascending=True)            
            df_chr.to_csv(output_dire+sample_name+"_"+chr_+'.maf.gz', sep='\t', index=False,compression='gzip')
            
    elif by_chrom == False:

        # SAVE ONE TABLE
        dff.to_csv(output_dire+sample_name+'.maf.gz', sep='\t', index=False,compression='gzip')
        
    else:
        print('by_chrom should be a True or False statement')


if __name__ == '__main__':
    cli()
