import matplotlib.pyplot as plt
import pandas as pd

# path to PAC output from HipSci iPSCs and neurons 
# (PAC outputs standard mapping and PAC results)
ipsc_stand_df = pd.read_csv('', delimiter='\t')
ipsc_df = pd.read_csv('/', delimiter='\t')
neur_stand_df = pd.read_csv('', delimiter='\t')
neur_df = pd.read_csv('', delimiter='\t')
# specify path
output_path = ""
gtf_path = ""


def analysis(ipsc_stand_df, ipsc_df, neur_stand_df, neur_df, output_path, gtf_path):
    """
    Run the analysis in my thesis (linked in relevant thesis section).
    Generates plots comparing the reference allele ratios 
    obtained from PAC and standard mapping. 
    Also generates plots comparing heterozygous sites under ASE detected 
    with binomial test (P<0.05) and Bonferroni-corrected that are unique 
    to cell types and those shared between cell types.
    Outputs genes overlapping ASE sites detected 
    with binomial test (P<0.05) and Bonferroni-corrected.

    Parameters
    ----------
    ipsc_stand_df : df
        PAC output for healthy iPSC HipSci donor from standard mapping.
    ipsc_df : df
        PAC output for healthy iPSC HipSci donor from PAC mapping.
    neur_stand_df : df
        PAC output for healthy neuronal HipSci donor from standard mapping.
    neur_df : df
        PAC output for healthy neuronal HipSci donor from PAC mapping.
    output_path : str
        String of path where the outputs figures and gene annotations will be saved to.
    gtf_path : str
        String of path for GTF file.     

    Returns
    -------
    figure
        6 figures saved in the path specified
        1st figure: reference allele ratio in standard mappign and PAC in iPSCs
        2nd figure: reference allele ratio in standard mappign and PAC in neurons
        3rd figure: cell type specific ASEs detected with binomial test (P<0.05)
        4th figure: cell type shared ASEs detected with binomial test (P<0.05)
        5th figure: cell type specific ASEs detected with binomial test (Bonferroni-corrected)
        6th figure: cell type shared ASEs detected with binomial test (Bonferroni-corrected)
        
    text file
        2 text files with gene annotations for ASE sites
        1st text file is gene annotations for heteroxygous sites for ASE sites 
        (binomial P<0.05) in iPSCs from PAC
        2nd text file is gene annotations for heteroxygous sites for ASE sites 
        (P-val Bonferroni-corrected) in iPSCs from PAC
            

    """
    def columns_to_numeric_and_filter_coverage(df):
        """
        Changes relevant columns to numeric types.
        """
        df['Chr'] = df['Chr'].astype(int)
        df['Pos'] = df['Pos'].astype(int)
        df['MapCov'] = df['MapCov'].astype(float)
        df['TotRead'] = df['TotRead'].astype(int)
        df = df.loc[df.TotRead >= 10]
        
        return df

    def ASE_selection(df_common, df):
        """
        Filters dataframe for significant ASE heterozygous sites. Returns a dataframe 
        with significant ASEv heterozygous sites according to binomial test P<0.05 (sign_df)
        and Bonferroni corrected (sign_mtc_df).
        """
        sign_df = df_common.loc[df_common['p'] < 0.05]
        mtc = 0.05/(len(df))
        sign_mtc_df = df_common.loc[df_common['p'] < mtc]
        return sign_df, sign_mtc_df

    def ASE_shared_sites(df, other_df):
        """
        Filters for heterozygous sites that are shared between two dataframes.
        """
        df_shared = df[df.Chr_Pos.isin(other_df.Chr_Pos)]
        return df_shared

    def ASE_unique_sites(df, other_df):
        """
        Filters for heterozygous sites that are not shared between two dataframes.
        """
        df_unique = df[~df.Chr_Pos.isin(other_df.Chr_Pos)]
        return df_unique
    
    def open_gencode_file(gtf_path):
        """
        Opens GTF file and assigns correct column names
        """
        # open gtf file
        gtf_df = pd.read_csv(gtf_path, sep='\t', comment='#',
                          names=['Chr', 'ignore', 'Feature', 'Start_Pos', 'End_Pos', 
                                  'ignore1', 'ignore2', 'ignore3', 'add_info'])
        
        # add_info column contains a list of columns
        # break up add_info column into separate columns
        gtf_df[['gene_id', 'transcript_id', 'gene_type', 'gene_status', 'gene_name', 
                  'transcript_type', 'transcript_status', 'transcript_name', 'level', 
                  'havana_gene']] = gtf_df['add_info'].str.split('; ', 9, expand=True)
        
        return gtf_df
    

    
    def gtf_file_filtering(gtf_df):
        """
        Filter GTF file for autosomal only genes, and clean up column values.
        """

        gtf_df = open_gencode_file(gtf_path)

        # select only gene in features
        gtf_df =  gtf_df[gtf_df['Feature'].isin(['gene'])]

        # leave only autosomal chromosomes
        gtf_df = gtf_df[~gtf_df['Chr'].isin(['chrM'])]
        gtf_df = gtf_df[~gtf_df['Chr'].isin(['chrY'])]
        gtf_df = gtf_df[~gtf_df['Chr'].isin(['chrX'])]

        # keep only relevant columns
        gtf_df = gtf_df[['Chr', 'Start_Pos', 'End_Pos', 'gene_id', 'gene_name']]

        # strip irrelevant text from column values
        gtf_df['gene_id']=gtf_df.gene_id.str.strip('gene_id')
        gtf_df['gene_id'] = gtf_df['gene_id'].replace('"', '', regex=True)
        gtf_df['gene_id'] = gtf_df['gene_id'].str.split('.').str[0]
        gtf_df['gene_name']=gtf_df.gene_name.str.strip('gene_name')
        gtf_df['gene_name'] = gtf_df['gene_name'].replace('"', '', regex=True)
        gtf_df['Chr']=gtf_df.Chr.str.strip('chr')

        # change to integers
        gtf_df['Chr'] = gtf_df['Chr'].astype(int)
        gtf_df['Start_Pos'] = gtf_df['Start_Pos'].astype(int)
        gtf_df['End_Pos'] = gtf_df['End_Pos'].astype(int)

        return gtf_df

    def gene_annotation_for_df(df, gtf_df):
        """
        Provides gene annotations for the heterozygous sites.
        Merges the PAC output dataframe with the cleaned GTF file,
        outputs a gene that overlaps the heterozygous site.
        """

        gtf_df = gtf_file_filtering(gtf_df)

        # merge dataframe and gtf on Chr
        merged_df = pd.merge(df, gtf_df, how='left', on=['Chr'])
        
        # select rows where position between start and end position columns
        merged_df = merged_df[(merged_df['Pos'].le(merged_df['End_Pos']))]
        merged_df = merged_df[(merged_df['Pos'].ge(merged_df['Start_Pos']))]

        # add columm with removed numbers after .
        merged_df['gene_name_new'] = merged_df['gene_name'].str.split('.').str[0]

        # remove duplicates
        merged_df = merged_df.drop_duplicates(subset=['gene_name_new'], keep='first')

        annotated_df = merged_df[['gene_id']]

        return annotated_df

    # chnage column types to numeric
    ipsc_df = columns_to_numeric_and_filter_coverage(ipsc_df)
    ipsc_stand_df = columns_to_numeric_and_filter_coverage(ipsc_stand_df)
    neur_df = columns_to_numeric_and_filter_coverage(neur_df)
    neur_stand_df = columns_to_numeric_and_filter_coverage(neur_stand_df)

    # dataframe with het sites expressed in both standard alignment and PAC
    ipsc_stand_common_df = (ipsc_stand_df[ipsc_stand_df.Chr_Pos.isin(ipsc_df.Chr_Pos)])
    ipsc_pac_common_df = (ipsc_df[ipsc_df.Chr_Pos.isin(ipsc_stand_df.Chr_Pos)])
    neur_stand_common_df = (neur_stand_df[neur_stand_df.Chr_Pos.isin(neur_df.Chr_Pos)])
    neur_pac_common_df = (neur_df[neur_df.Chr_Pos.isin(neur_stand_df.Chr_Pos)])

    # dataframe with het sites expressed in both tissues in PAC
    ipsc_common_df = (ipsc_df[ipsc_df.Chr_Pos.isin(neur_df.Chr_Pos)])
    neur_common_df = (neur_df[neur_df.Chr_Pos.isin(ipsc_df.Chr_Pos)])

    # filters for significant ASE sites (binomial P<0.05 and Bonferroni-corrected)
    ipsc_sign_df, ipsc_sign_mtc_df = ASE_selection(ipsc_common_df, ipsc_df)
    neur_sign_df, neur_sign_mtc_df = ASE_selection(neur_common_df, neur_df)

    # filters for significant ASE sites shared between cell types
    common_ipsc_neur_sign_df = ASE_shared_sites(ipsc_sign_df, neur_sign_df)
    common_ipsc_neur_sign_mtc_df = ASE_shared_sites(ipsc_sign_mtc_df, neur_sign_mtc_df)

    # filters for significant ASE sites from PAC that are unique for cell type
    uniq_ipsc_sign_df = ASE_unique_sites(ipsc_sign_df, common_ipsc_neur_sign_df)
    uniq_neur_sign_df = ASE_unique_sites(neur_sign_df, common_ipsc_neur_sign_df)
    uniq_ipsc_sign_mtc_df = ASE_unique_sites(ipsc_sign_mtc_df, common_ipsc_neur_sign_mtc_df)
    uniq_neur_sign_mtc_df = ASE_unique_sites(neur_sign_mtc_df, common_ipsc_neur_sign_mtc_df)

    # sort dataframe for the pltos
    uniq_ipsc_sign_df = uniq_ipsc_sign_df.sort_values(by='MapCov')
    uniq_neur_sign_df = uniq_neur_sign_df.sort_values(by='MapCov')
    common_ipsc_neur_sign_df = common_ipsc_neur_sign_df.sort_values(by='MapCov')
    uniq_ipsc_sign_mtc_df = uniq_ipsc_sign_mtc_df.sort_values(by='MapCov')
    uniq_neur_sign_mtc_df = uniq_neur_sign_mtc_df.sort_values(by='MapCov')
    common_ipsc_neur_sign_mtc_df = common_ipsc_neur_sign_mtc_df.sort_values(by='MapCov')

    # open and clean GTF file
    gtf_df = open_gencode_file(gtf_path)
    gtf_df = gtf_file_filtering(gtf_df)
    
    # obtain genes overlapping the heterozygous sites for iPSCs
    annotated_uniq_ipsc_sign_df = gene_annotation_for_df(uniq_ipsc_sign_df, gtf_df)
    annotated_uniq_ipsc_sign_mtc_df = gene_annotation_for_df(uniq_ipsc_sign_mtc_df, gtf_df)

    # save the annotations
    annotated_uniq_ipsc_sign_df = pd.DataFrame(annotated_uniq_ipsc_sign_df).to_csv(output_path + "annotated_ipsc_uniq_ase.txt", sep='\t', index=False, header=None)
    annotated_uniq_ipsc_sign_mtc_df = pd.DataFrame(annotated_uniq_ipsc_sign_mtc_df).to_csv(output_path + "annotated_ipsc_uniq_bonferroni_ase.txt", sep='\t', index=False, header=None)
    
    # plot ipsc stand and pac
    plt.figure(1)
    graph_ipsc_stand_and_pac = ipsc_stand_common_df.MapCov.plot.density(color="red")
    graph_ipsc_stand_and_pac = ipsc_pac_common_df.MapCov.plot.density(color="black")
    graph_ipsc_stand_and_pac.set_xlim(0, 1)
    graph_ipsc_stand_and_pac.legend(['Standard approach', 'PAC'])
    plt.axvline(x=0.5)
    graph_ipsc_stand_and_pac.set_xlabel('Ref allele ratio')
    plt.text(0.02, 2.6, 'Mean & Median' + '\nRef: ' +
             str(round(ipsc_stand_common_df['MapCov'].mean(), 2)) + ' & ' +
             str(round(ipsc_stand_common_df['MapCov'].median(), 2)) + '\nPAC: ' +
             str(round(ipsc_pac_common_df['MapCov'].mean(), 2)) + ' & ' +
             str(round(ipsc_pac_common_df['MapCov'].median(), 2)))
    plt.savefig(output_path + "ipsc_allelic_ratios.png")
    plt.close()

    # plot neur stand and pac
    plt.figure(2)
    graph_neur_stand_and_pac = neur_stand_common_df.MapCov.plot.density(color="red")
    graph_neur_stand_and_pac = neur_pac_common_df.MapCov.plot.density(color="black")
    graph_neur_stand_and_pac.set_xlim(0, 1)
    graph_neur_stand_and_pac.legend(['Standard approach', 'PAC'])
    plt.axvline(x=0.5)
    graph_neur_stand_and_pac.set_xlabel('Ref allele ratio')
    plt.text(0.02, 2.6, 'Mean & Median' + '\nRef: ' +
             str(round(neur_stand_common_df['MapCov'].mean(), 2)) + ' & ' +
             str(round(neur_stand_common_df['MapCov'].median(), 2)) + '\nPAC: ' +
             str(round(neur_pac_common_df['MapCov'].mean(), 2)) + ' & ' +
             str(round(neur_pac_common_df['MapCov'].median(), 2)))
    plt.savefig(output_path + "neur_allelic_ratios.png")
    plt.close()

    # plot ipsc and neur unique ASEs
    plt.figure(3)
    graph_ipsc_neur_uniq_sign = uniq_ipsc_sign_df.MapCov.plot.density(color="red")
    graph_ipsc_neur_uniq_sign = uniq_neur_sign_df.MapCov.plot.density(color="black")
    graph_ipsc_neur_uniq_sign.set_xlim(0, 1)
    graph_ipsc_neur_uniq_sign.legend(['ipsc', 'neuron'])
    graph_ipsc_neur_uniq_sign.set_xlabel('Ratio of reads mapping to ref allele')
    plt.title('iPSC and neurons: unique ASEs (detected at p<0.05)')
    plt.savefig(output_path + "ipsc_neur_uniq_ase.png")
    plt.close()

    # plot common ones
    plt.figure(4)
    graph_common_ipsc_neur_pac_sign = common_ipsc_neur_sign_df.MapCov.plot.density(color="black")
    graph_common_ipsc_neur_pac_sign.set_xlim(0, 1)
    graph_common_ipsc_neur_pac_sign.legend(['common'])
    graph_common_ipsc_neur_pac_sign.set_xlabel('Ratio of reads mapping to ref allele')
    plt.title('iPSC and neurons: common ASEs (detected at p<0.05)')
    plt.savefig(output_path + "common_ipsc_neur_ase.png")
    plt.close()

    # plot ipsc and neur unique ASEs mtc
    plt.figure(5)
    graph_ipsc_neur_uniq_sign_mtc = uniq_ipsc_sign_mtc_df.MapCov.plot.density(color="red")
    graph_ipsc_neur_uniq_sign_mtc = uniq_neur_sign_mtc_df.MapCov.plot.density(color="black")
    graph_ipsc_neur_uniq_sign_mtc.set_xlim(0, 1)
    graph_ipsc_neur_uniq_sign_mtc.legend(['ipsc', 'neuron'])
    graph_ipsc_neur_uniq_sign_mtc.set_xlabel('Ratio of reads mapping to ref allele')
    plt.title('iPSC and neurons: unique ASEs (detection: Bonferroni corrected)')
    plt.savefig(output_path + "ipsc_neur_uniq_mtc_ase.png")
    plt.close()

    # plot common ASEs mtc
    plt.figure(6)
    graph_common_ipsc_neur_sign_mtc = common_ipsc_neur_sign_mtc_df.MapCov.plot.density(color="black")
    graph_common_ipsc_neur_sign_mtc.set_xlim(0, 1)
    graph_common_ipsc_neur_sign_mtc.set_xlabel('Ratio of reads mapping to ref allele')
    plt.title('iPSC and neurons: common ASEs (detection: Bonferroni corrected)')
    plt.savefig(output_path + "common_ipsc_neur_mtc_ase.png")
    plt.close()

    return annotated_uniq_ipsc_sign_df, annotated_uniq_ipsc_sign_mtc_df, graph_ipsc_stand_and_pac, graph_neur_stand_and_pac, graph_ipsc_neur_uniq_sign, graph_common_ipsc_neur_pac_sign, graph_ipsc_neur_uniq_sign_mtc, graph_common_ipsc_neur_sign_mtc

analysis(ipsc_stand_df, ipsc_df, neur_stand_df, neur_df, output_path, gtf_path)
