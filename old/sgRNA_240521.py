### sgRNA.py ###
# Author: Marc Zepeda
# Date: 2024-05-20

# Import Packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# From 230413_sgRNA_library
''' get_csv: Creates a dictionary containing DataFrames corresponding to subset of csvs in a folder.
        path_prefix = folder
        path_suffix = csv subset
        meta = metadata DataFrame
'''
def get_csvs(path_prefix: str(), path_suffix: str(), meta: pd.DataFrame()):
    
    # Generate file paths, create empty dictionary, and obtain metadata.
    file_paths = glob.glob(os.path.join(path_prefix,path_suffix))
    dic = {}
    metadata = meta.drop_duplicates('fasta_ID')

    # Iterate through file paths.
    for file_path in file_paths:
        
        # Obtain DataFrame, fasta_ID, editor, & pam_group from file path. Determine gene from metadata DataFrame.
        file = pd.read_csv(file_path)
        fasta_ID = file_path.split('/')[-1].split('_')[0]
        editor = file_path.split('/')[-1].split('_')[1]
        pam_group = file_path.split('/')[-1].split('_')[2]
        gene = metadata[metadata['fasta_ID'] == fasta_ID]['Gene'].reset_index(drop=True)[0]
        file['Gene'] = gene
        file['Editor'] = editor
        file['PAM Group'] = pam_group

        # Rearrange columns. Remove columns with incorrect data: 'Domain' and 'Plot_site'.
        old_cols = file.columns.tolist()
        new_cols = [old_cols[2]] + old_cols[-2:] + old_cols[0:2] + old_cols[3:-4]
        file = file[new_cols]

        # Add DataFrame to Dictionary
        dic[gene] = file

    return dic

''' create_csvs: Creates csvs from dictionary containing DataFrames.
        dic = dictionary
        out_dir_prefix = main folder
        out_dir_suffix = subfolder
'''
def create_csv(dic: dict(), out_dir_prefix: str(), out_dir_suffix: str()):
    
    # Check for and create output folder for csvs.
    out_dir = os.path.join(out_dir_prefix, out_dir_suffix)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Obtain gene list from dictionary keys.
    gene_list = list(dic.keys())

    # Iterate through genes.
    for gene in gene_list:
        
        # Isolate DataFrame from dictionary.
        file = dic[gene]
        
        # Create file name fome gene, editor, and PAM group.
        editor = file.iloc[0]['Editor']
        pam_group = file.iloc[0]['PAM Group']
        file_name = gene + '_' + editor + '_' + pam_group + '.csv'

        # Create file.
        file_path = os.path.join(out_dir, file_name)
        file.to_csv(file_path, sep=',', index=False)

''' annotate: Filter sgRNA library for missense, silent, and nonsense mutations. Label specified domains within residue ranges.
        dic = dictionary
        meta = metadata DataFrame
'''
def annotate(dic: dict(), meta: pd.DataFrame()):

    # Obtain gene list from dictionary keys and create empty dictionary
    gene_list = list(dic.keys())
    dic2 = {}

    # Iterate through genes.
    for gene in gene_list:

        # Isolate DataFrame from dictionary
        file = dic[gene]

        # Only retain missense, silent, and nonsense mutations.
        file = file[(file['Mut_type'] == 'Missense') |
                (file['Mut_type'] == 'Silent') |
                (file['Mut_type'] == 'Nonsense')].reset_index(drop=True)

        # Create Domain column with Unspecified label.
        file['Domain'] = 'Unspecified'

        # Obtain domain metadata.
        domain_metadata = meta[meta['Gene'] == gene]

        # Iterate through Target_CDS.
        for i_res in range(file.shape[0]):

                # Obtain sgRNA residue coordinates.
                st_crds_res = int(file.iloc[i_res]['Target_resi'].split(",")[0][1:])
                end_crds_res = int(file.iloc[i_res]['Target_resi'].split(",")[1][1:-1])

                # Annotate Domain with metadata label.
                for i_dom in range(domain_metadata.shape[0]):

                        # Obtain domain name and residue range.
                        dom = domain_metadata.iloc[i_dom]['Domain']
                        dom_res_st = domain_metadata.iloc[i_dom]['Start']
                        dom_res_end = domain_metadata.iloc[i_dom]['End']

                        # Specify domain if sgRNA is within the residue range.
                        if (end_crds_res >= dom_res_st) & (st_crds_res <= dom_res_end):
                                file.loc[i_res, 'Domain'] = dom
                                break
        
        # Rearrange columns.
        old_cols = file.columns.tolist()
        new_cols = [old_cols[0]] + [old_cols[-1]] + old_cols[1:-1]
        file = file[new_cols]

        # Add DataFrame to dictionary.
        dic2[gene] = file

    return dic2

''' filter_domains: Filter sgRNA library for specified domains.
        dic = dictionary
'''
def filter_domains(dic: dict()):
    
    # Obtain gene list from dictionary keys and create empty dictionary
    gene_list = list(dic.keys())
    dic2 = {}

    # Iterate through genes.
    for gene in gene_list:

        # Isolate DataFrame from dictionary.
        file = dic[gene]

        # Remove unspecified domains.
        file = file[file['Domain']!='Unspecified'].reset_index(drop=True)

        # Add DataFrame to dictionary.
        dic2[gene] = file

    return dic2

''' domain_counts: Obtain DataFrame containing domain counts for sgRNA library.
        dic = dictionary
'''
def domain_counts(dic: dict()):

    # Obtain gene list from dictionary keys and create empty dictionary
    gene_list = list(dic.keys())
    df = pd.DataFrame(columns=['Gene','Domain','sgRNAs','Editor','PAM Group','Editor_PAM'])

    # Iterate through genes.
    for gene in gene_list:

        # Isolate DataFrame from dictionary
        file = dic[gene]
        
        # Obtain domains and associated sgRNA counts; obtain editor and PAM group; append to DataFrame
        doms = file['Domain'].value_counts().index.tolist()
        dom_cts = file['Domain'].value_counts().tolist()
        editor = file.iloc[0]['Editor']
        pam_group = file.iloc[0]['PAM Group']
        editor_PAM = editor + "_" + pam_group
        dic2 = {'Gene': gene,'Domain': doms,'sgRNAs': dom_cts,'Editor': editor,'PAM Group': pam_group,'Editor_PAM': editor_PAM}
        temp = pd.DataFrame(dic2)
        df = pd.concat([df,temp])

    return df.reset_index(drop=True)

# From 230415_sgRNA_library
''' filter_missense: Filter sgRNA library for missense mutations.
        dic = dictionary
'''
def filter_missense(dic: dict()):

    # Obtain gene list from dictionary keys and create empty dictionary
    gene_list = list(dic.keys())
    dic2 = {}

    # Iterate through genes.
    for gene in gene_list:

        # Isolate DataFrame from dictionary.
        file = dic[gene]

        # Retain missense mutation sgRNAs.
        file = file[file['Mut_type'] == 'Missense'].reset_index(drop=True)

        # Add DataFrame to dictionary.
        dic2[gene] = file

    return dic2

''' summary_df: Create DataFrame containing summary information for sgRNA library: CBE/ABE sgRNAs & CBE/ABE Coverage.
        dic = dictionary CBE
        dic = dictionary ABE
        meta = metadata DataFrame
'''
def summary_df(dic_CBE: dict(), dic_ABE: dict(), meta: pd.DataFrame()):

    # Obtain gene list from dictionary keys and create empty dictionary
    gene_list = list(dic_CBE.keys())
    df = pd.DataFrame()

    # Iterate through genes.
    for gene in gene_list:

        # Isolate DataFrame from dictionary.
        file_CBE = dic_CBE[gene]
        file_ABE = dic_ABE[gene]

        # Obtain domains, residue ranges, and associated sgRNA counts for CBE and ABE
        doms = list(meta[meta['Gene'] == gene]['Domain'])
        start = list(
            meta[meta['Gene'] == gene]['Start'])
        end = list(
            meta[meta['Gene'] == gene]['End'])

        if (file_CBE['Domain'].value_counts().shape[0] == len(start)) & (file_ABE['Domain'].value_counts().shape[0] == len(start)):
            # All domains have sgRNAs for CBE and ABE.
            sgRNAs_CBE = file_CBE['Domain'].value_counts().reindex(doms).tolist()
            sgRNAs_ABE = file_ABE['Domain'].value_counts().reindex(doms).tolist()
        elif file_CBE['Domain'].value_counts().shape[0] == len(start):
            # All domains have sgRNAs for CBE, but not ABE.
            sgRNAs_CBE = file_CBE['Domain'].value_counts().reindex(doms).tolist()
            sgRNAs_ABE = file_ABE['Domain'].value_counts().reindex(doms).tolist()
            sgRNAs_ABE.extend([0] * (len(start) - len(sgRNAs_ABE)))
        elif file_ABE['Domain'].value_counts().shape[0] == len(start):
            # All domains have sgRNAs for ABE, but not CBE.
            sgRNAs_CBE = file_CBE['Domain'].value_counts().reindex(doms).tolist()
            sgRNAs_ABE = file_ABE['Domain'].value_counts().reindex(doms).tolist()
            sgRNAs_CBE.extend([0] * (len(start) - len(sgRNAs_CBE)))
        else:
            # Catch all error message.
            sgRNAs_CBE = file_CBE['Domain'].value_counts().reindex(doms).tolist()
            sgRNAs_ABE = file_ABE['Domain'].value_counts().reindex(doms).tolist()
            print("Error: ", gene, doms, start, end, sgRNAs_CBE, sgRNAs_ABE)

        # Create empty coverage CBE and ABE lists.
        cov_CBE = []
        cov_ABE = []

        # Interate through domains.
        for i_dom, dom in enumerate(doms):
            
            # Isolate a domain
            file_CBE_dom = file_CBE[file_CBE['Domain'] == dom]
            file_ABE_dom = file_ABE[file_ABE['Domain'] == dom]

            # Find unique edit sites for CBE and ABE
            edit_sites_CBE = set()
            for edit_site in file_CBE_dom['Edit_site']:
                if edit_site % 1 != 0:
                    edit_sites_CBE.add(int(round(edit_site+.1)))
                    edit_sites_CBE.add(int(round(edit_site+.1))-1)
                else:
                    edit_sites_CBE.add(int(edit_site))

            edit_sites_ABE = set()
            for edit_site in file_ABE_dom['Edit_site']:
                if edit_site % 1 != 0:
                    edit_sites_ABE.add(int(round(edit_site+.1)))
                    edit_sites_ABE.add(int(round(edit_site+.1))-1)
                else:
                    edit_sites_ABE.add(int(edit_site))

            # Calculate CBE and ABE coverage for domain given unique edit sites and residue range
            cov_CBE_dom = len(edit_sites_CBE)/(end[i_dom]-start[i_dom]+1)
            cov_ABE_dom = len(edit_sites_ABE)/(end[i_dom]-start[i_dom]+1)

            # Add CBE and ABE coverage for domain to CBE and ABE coverage lists.
            cov_CBE.append(cov_CBE_dom)
            cov_ABE.append(cov_ABE_dom)

        # Catch 0 sgRNAs error.
        i = 0
        while i == 0:
            if len(sgRNAs_CBE) < len(cov_CBE):
                sgRNAs_CBE.append(0)
            elif len(sgRNAs_ABE) < len(cov_ABE):
                sgRNAs_ABE.append(0)
            else:
                i = 1

        # Create list of gene repeated by number of domains.
        gene_num_doms = [gene] * len(doms)

        # Create dictionary, convert to temp DataFrame, and append it to output DataFrame.
        dic = {'Gene': gene_num_doms, 'Domain': doms,
            'Start': start, 'End': end,
            'CBE sgRNAs': sgRNAs_CBE, 'ABE sgRNAs': sgRNAs_ABE,
            'CBE Coverage': cov_CBE, 'ABE Coverage': cov_ABE}
        temp = pd.DataFrame(dic)
        df = pd.concat([df, temp])

    return df.sort_values(by=['Gene','Start']).reset_index(drop=True)

''' annotate: Filter sgRNA library for missense, silent, and nonsense mutations. Label specified domains within residue ranges.
        dic = dictionary
        meta = metadata DataFrame
'''
def annotate(dic: dict(), meta: pd.DataFrame()):

    # Obtain gene list from dictionary keys and create empty dictionary
    gene_list = list(dic.keys())
    dic2 = {}

    # Iterate through genes.
    for gene in gene_list:

        # Isolate DataFrame from dictionary
        file = dic[gene]

        # Only retain missense, silent, and nonsense mutations.
        file = file[(file['Mut_type'] == 'Missense') |
                    (file['Mut_type'] == 'Silent') |
                    (file['Mut_type'] == 'Nonsense')].reset_index(drop=True)

        # Create Domain column with Unspecified label.
        file['Domain'] = 'Unspecified'

        # Obtain domain metadata.
        domain_metadata = meta[meta['Gene'] == gene]
        
        # Iterate through Edit_site.
        for i_res in range(file.shape[0]):

            # Obtain sgRNA residue coordinates.
            edit_site = file.iloc[i_res]['Edit_site']
            if edit_site % 1 != 0:
                st_crds_res = int(round(edit_site+.1)-1)
                end_crds_res = int(round(edit_site+.1))
            else:
                st_crds_res = end_crds_res = int(edit_site)
            
            # Annotate Domain with metadata label.
            for i_dom in range(domain_metadata.shape[0]):

                # Obtain domain name and residue range.
                dom = domain_metadata.iloc[i_dom]['Domain']
                dom_res_st = domain_metadata.iloc[i_dom]['Start']
                dom_res_end = domain_metadata.iloc[i_dom]['End']

                # Specify domain if sgRNA is within the residue range.
                if (end_crds_res >= dom_res_st) & (st_crds_res <= dom_res_end):
                    file.loc[i_res, 'Domain'] = dom
                    break

        # Rearrange columns.
        old_cols = file.columns.tolist()
        new_cols = [old_cols[0]] + [old_cols[-1]] + old_cols[1:-1]
        file = file[new_cols]

        # Add DataFrame to dictionary.
        dic2[gene] = file

    return dic2

# From 230417_sgRNA_library
''' get_edit_sites: Creates a list of all edit sites from DataFrame with Edit_site column.
        df = Annotated CRISPOR DataFrame
'''
def get_edit_sites(df: pd.DataFrame()):
    
    # Create empty edit_sites list 
    edit_sites = list()

    # Fill edit_sites list
    for edit_site in df['Edit_site']:
        if edit_site % 1 != 0:
            # sgRNA target two amino acids.
            edit_sites.extend(
                [int(round(edit_site+.1))-1, int(round(edit_site+.1))])
        else:
            # sgRNA targets one amino acid
            edit_sites.extend([int(edit_site)])
    
    return edit_sites

''' histplots: Generates histplots for sgRNA tiling libraries with multiple genes.
        dic_CBE,dic_ABE: dictionary of genes, sgRNA libaries data frames
        dir: specify output directory to save
'''
def histplots(dic_CBE: dict(), dic_ABE: dict(),dir='',title_suffix=' sgRNA library'):
    
    # Check for and create output folder for barplots.
    if dir!='' and not os.path.exists(dir): os.mkdir(dir)

    # Iterate through genes.
    for gene in list(dic_CBE.keys()):

        # Isolate DataFrame from dictionary.
        file_CBE = dic_CBE[gene]
        file_ABE = dic_ABE[gene]

        # Obtain all edit sites for CBE and ABE sgRNAs
        edit_sites_CBE = get_edit_sites(file_CBE)
        edit_sites_ABE = get_edit_sites(file_ABE)

        # Generate CBE and ABE subplots.
        fig, (ax_CBE, ax_ABE) = plt.subplots(nrows=2, ncols=1)
        ax_CBE.hist(edit_sites_CBE, bins=np.arange(
            min(edit_sites_CBE), max(edit_sites_CBE) + 1, 1))
        ax_ABE.hist(edit_sites_ABE, bins=np.arange(
            min(edit_sites_ABE), max(edit_sites_ABE) + 1, 1))
        
        # Add title, x labels, and y labels to subplots.
        ax_CBE.set_title(gene + title_suffix)
        ax_ABE.set_xlabel('Amino Acid')
        ax_CBE.set_ylabel('CBE sgRNAs')
        ax_ABE.set_ylabel('ABE sgRNAs')
        fig.tight_layout()
        fig.show()
        
        # Save barplot.
        plot_name = gene + '.png'
        plot_path = os.path.join(dir, plot_name)
        if dir != '': fig.savefig(plot_path, dpi=600, pad_inches=1)