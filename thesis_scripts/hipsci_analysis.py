import matplotlib.pyplot as plt
import pandas as pd

# path to PAC output from HipSci iPSCs and neurons 
# (PAC outputs standard mapping and PAC results)
ipsc_stand_df = pd.read_csv('', delimiter='\t')
ipsc_df = pd.read_csv('', delimiter='\t')
neur_stand_df = pd.read_csv('/', delimiter='\t')
neur_df = pd.read_csv('', delimiter='\t')

# specify path
analysis(ipsc_stand_df, ipsc_df, neur_stand_df, neur_df, path="")

def analysis(ipsc_stand_df, ipsc_df, neur_stand_df, neur_df, path):
    """
    Run the analysis in my thesis (linked in relevant thesis section).
    Generates plots comparing the reference allele ratios 
    obtained from PAC and standard mapping. 
    Also generates plots comparing heterozygous sites under ASE detected 
    with binomial test (P<0.05) and Bonferroni-corrected that are unique 
    to cell types and those shared between cell types.

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
    path : str
        String of path where the difures will be saved to.

    Returns
    -------
    figure
        6 figures saved in the path specified
        1st figure: reference allele ratio in standard mappign and PAC in iPSCs
        2nd figure: reference allele ratio in standard mappign and PAC in neurons
        3rd figure: 

    """
    def columns_to_numeric(df):
        """
        Changes relevant columns to numeric types.
        """
        df['Chr'] = df['Chr'].astype(int)
        df['Pos'] = df['Pos'].astype(int)
        df['MapCov'] = df['MapCov'].astype(float)
        df['TotRead'] = df['TotRead'].astype(int)
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

    # chnage column types to numeric
    ipsc_df = columns_to_numeric(ipsc_df)
    ipsc_stand_df = columns_to_numeric(ipsc_stand_df)
    neur_df = columns_to_numeric(neur_df)
    neur_stand_df = columns_to_numeric(neur_stand_df)

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
    plt.savefig(path + "ipsc_allelic_ratios.png")
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
    plt.savefig(path + "neur_allelic_ratios.png")
    plt.close()

    # plot ipsc and neur unique ASEs
    plt.figure(3)
    graph_ipsc_neur_uniq_sign = uniq_ipsc_sign_df.MapCov.plot.density(color="red")
    graph_ipsc_neur_uniq_sign = uniq_neur_sign_df.MapCov.plot.density(color="black")
    graph_ipsc_neur_uniq_sign.set_xlim(0, 1)
    graph_ipsc_neur_uniq_sign.legend(['ipsc', 'neuron'])
    graph_ipsc_neur_uniq_sign.set_xlabel('Ratio of reads mapping to ref allele')
    plt.title('iPSC and neurons: unique ASEs (detected at p<0.05)')
    plt.savefig(path + "ipsc_neur_uniq_ase.png")
    plt.close()

    # plot common ones
    plt.figure(4)
    graph_common_ipsc_neur_pac_sign = common_ipsc_neur_sign_df.MapCov.plot.density(color="black")
    graph_common_ipsc_neur_pac_sign.set_xlim(0, 1)
    graph_common_ipsc_neur_pac_sign.legend(['common'])
    graph_common_ipsc_neur_pac_sign.set_xlabel('Ratio of reads mapping to ref allele')
    plt.title('iPSC and neurons: common ASEs (detected at p<0.05)')
    plt.savefig(path + "common_ipsc_neur_ase.png")
    plt.close()

    # plot ipsc and neur unique ASEs mtc
    plt.figure(5)
    graph_ipsc_neur_uniq_sign_mtc = uniq_ipsc_sign_mtc_df.MapCov.plot.density(color="red")
    graph_ipsc_neur_uniq_sign_mtc = uniq_neur_sign_mtc_df.MapCov.plot.density(color="black")
    graph_ipsc_neur_uniq_sign_mtc.set_xlim(0, 1)
    graph_ipsc_neur_uniq_sign_mtc.legend(['ipsc', 'neuron'])
    graph_ipsc_neur_uniq_sign_mtc.set_xlabel('Ratio of reads mapping to ref allele')
    plt.title('iPSC and neurons: unique ASEs (detection: Bonferroni corrected)')
    plt.savefig(path + "ipsc_neur_uniq_mtc_ase.png")
    plt.close()

    # plot common ASEs mtc
    plt.figure(6)
    graph_common_ipsc_neur_sign_mtc = common_ipsc_neur_sign_mtc_df.MapCov.plot.density(color="black")
    graph_common_ipsc_neur_sign_mtc.set_xlim(0, 1)
    graph_common_ipsc_neur_sign_mtc.set_xlabel('Ratio of reads mapping to ref allele')
    plt.title('iPSC and neurons: common ASEs (detection: Bonferroni corrected)')
    plt.savefig(path + "common_ipsc_neur_mtc_ase.png")
    plt.close()

    return graph_ipsc_stand_and_pac, graph_neur_stand_and_pac, graph_ipsc_neur_uniq_sign, graph_common_ipsc_neur_pac_sign, graph_ipsc_neur_uniq_sign_mtc, graph_common_ipsc_neur_sign_mtc
