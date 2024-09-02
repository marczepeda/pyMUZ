### pe.py ###
# Author: Marc Zepeda
# Date: 2024-08-31

# Import packages
import pandas as pd
import numpy as np
import os
from ..pe import pegLIT as pegLIT
import pyMUZ.gen.io as io
import pyMUZ.gen.plot as p

# PrimeDesign Methods
''' PrimeDesign: Run PrimeDesign using Docker
        file (str): Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)
        pbs_length_list (list): List of primer binding site (PBS) lengths for the pegRNA extension. Example: 12 13 14 15
        rtt_length_list (list): List of reverse transcription (RT) template lengths for the pegRNA extension. Example: 10 15 20
        nicking_distance_minimum (int): Minimum nicking distance for designing ngRNAs. (Default: 0 bp)
        nicking_distance_maximum (int): Maximum nicking distance for designing ngRNAs. (Default: 100 bp)
        filter_c1_extension (bool): Filter against pegRNA extensions that start with a C base. (Default: False)
        silent_mutation (bool): Introduce silent mutation into PAM assuming sequence is in-frame. Currently only available with SpCas9. (Default: False)
        genome_wide_design (bool): Whether this is a genome-wide pooled design. This option designs a set of pegRNAs per input without ranging PBS and RTT parameters.
        saturation_mutagenesis (str): Saturation mutagenesis design with prime editing (Options: 'aa', 'base').
        number_of_pegrnas (int): Maximum number of pegRNAs to design for each input sequence. The pegRNAs are ranked by 1) PAM disrupted > PAM intact then 2) distance to edit. (Default: 3)
        number_of_ngrnas (int): Maximum number of ngRNAs to design for each input sequence. The ngRNAs are ranked by 1) PE3b-seed > PE3b-nonseed > PE3 then 2) deviation from nicking_distance_pooled. (Default: 3)
        nicking_distance_pooled (int): The nicking distance between pegRNAs and ngRNAs for pooled designs. PE3b annotation is priority (PE3b seed -> PE3b non-seed), followed by nicking distance closest to this parameter. (Default: 75 bp)
        homology_downstream (int): For pooled designs (genome_wide or saturation_mutagenesis needs to be indicated), this parameter determines the RT extension length downstream of an edit for pegRNA designs. (Default: 10)
        pbs_length_pooled (int): The PBS length to design pegRNAs for pooled design applications. (Default: 14 nt)
        rtt_max_length_pooled (int): Maximum RTT length to design pegRNAs for pooled design applications. (Default: 50 nt)
        out_dir (str): Name of output directory. (Default: ./DATETIMESTAMP_PrimeDesign)
    Dependencies: os,numpy,https://github.com/pinellolab/PrimeDesign
'''
def PrimeDesign(file: str,pbs_length_list: list = [],rtt_length_list: list = [],nicking_distance_minimum: int = 0,
                nicking_distance_maximum: int = 100,filter_c1_extension: bool = False,silent_mutation: bool = False,
                genome_wide_design: bool = False,saturation_mutagenesis: str = None,number_of_pegrnas: int = 3,number_of_ngrnas: int = 3,
                nicking_distance_pooled: int = 75,homology_downstream: int = 10,pbs_length_pooled: int = 14,rtt_max_length_pooled: int = 50,
                out_dir: str = './DATETIMESTAMP_PrimeDesign'):
    
    # Write PrimeDesign Command Line
    cmd = 'docker run -v ${PWD}/:/DATA -w /DATA pinellolab/primedesign primedesign_cli' # prefix
    cmd += f' -f {file}' # Append required parameters
    if pbs_length_list: cmd += f' -pbs {str(np.array(pbs_length_list))[1:-1]}' # Append optional parameters
    if rtt_length_list: cmd += f' -rtt {str(np.array(rtt_length_list))[1:-1]}'
    if nicking_distance_minimum!=0: cmd += f' -nick_dist_min {str(nicking_distance_minimum)}' 
    if nicking_distance_maximum!=100: cmd += f' -nick_dist_max {str(nicking_distance_maximum)}'
    if filter_c1_extension: cmd += f' -filter_c1 {str(filter_c1_extension)}'
    if silent_mutation: cmd += f' - silent_mut'
    if genome_wide_design: cmd += f' -genome_wide'
    if saturation_mutagenesis: cmd += f' -sat_mut {saturation_mutagenesis}'
    if number_of_pegrnas!=3: cmd += f' -n_pegrnas {number_of_pegrnas}'
    if number_of_ngrnas!=3: cmd += f' -n_ngrnas {number_of_ngrnas}'
    if nicking_distance_pooled!=75: cmd += f' -nick_dist_pooled {nicking_distance_pooled}'
    if homology_downstream!=10: cmd += f' -homology_downstream {homology_downstream}'
    if pbs_length_pooled!=14: cmd += f' -pbs_pooled {pbs_length_pooled}'
    if rtt_max_length_pooled!=50: cmd += f' -rtt_pooled {rtt_max_length_pooled}'
    if out_dir!='./DATETIMESTAMP_PrimeDesign': cmd+= f' -out ./DATETIMESTAMP_PrimeDesign'
    print(cmd)

    os.system(cmd) # Execute PrimeDesign Command Line

''' PrimeDesignOutput: Splits peg/ngRNAs from PrimeDesign output & finishes annotations
        pt: path to primeDesign output
        aa_index: 1st amino acid in target sequence index (Optional, Default: start codon = 1)
        scaffold_sequence: sgRNA scaffold sequence (Optional, Default: SpCas9)
    Dependencies: io,numpy
'''
def PrimeDesignOutput(pt: str, aa_index: int=1, 
                      scaffold_sequence: str='GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC'):
    
    # Get PrimeDesign output & seperate pegRNAs and ngRNAs
    primeDesign_output = io.get(pt)
    pegRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='pegRNA']
    ngRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='ngRNA']

    # Generate epegRNAs
    pegRNAs['Edit']=[str(target_name.split('_')[2].split('to')[0]) + # AA Before
                    str(int(target_name.split('_')[1]) + aa_index-1) + # AA Index
                    str(target_name.split('_')[2].split('to')[1]) # AA After
                    for target_name in pegRNAs['Target_name']]
    pegRNAs['Scaffold_sequence']=[scaffold_sequence]*len(pegRNAs)
    pegRNAs['RTT_sequence']=[pegRNAs.iloc[i]['Extension_sequence'][0:int(pegRNAs.iloc[i]['RTT_length'])] for i in range(len(pegRNAs))]
    pegRNAs['PBS_sequence']=[pegRNAs.iloc[i]['Extension_sequence'][int(pegRNAs.iloc[i]['RTT_length']):]  for i in range(len(pegRNAs))]
    pegRNAs=pegRNAs[['pegRNA_number','gRNA_type','Strand','Edit', # Important metadata
                     'Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence', # Sequence information
                     'Target_name','Target_sequence','Spacer_GC_content','PAM_sequence','Extension_sequence','Annotation','pegRNA-to-edit_distance','Nick_index','ngRNA-to-pegRNA_distance','PBS_length','PBS_GC_content','RTT_length','RTT_GC_content','First_extension_nucleotide']] # Less important metadata

    # Generate ngRNAs
    ngRNAs['Edit']=[str(target_name.split('_')[2].split('to')[0]) + # AA Before
                    str(int(target_name.split('_')[1]) + aa_index-1) + # AA Index
                    str(target_name.split('_')[2].split('to')[1]) # AA After
                    for target_name in ngRNAs['Target_name']]
    ngRNAs['Scaffold_sequence']=[scaffold_sequence]*len(ngRNAs)
    ngRNAs['ngRNA_number']=list(np.arange(len(ngRNAs)))
    ngRNAs=ngRNAs[['pegRNA_number','ngRNA_number','gRNA_type','Strand','Edit', # Important metadata
                   'Spacer_sequence','Scaffold_sequence', # Sequence information
                   'Target_name','Target_sequence','Spacer_GC_content','PAM_sequence','Extension_sequence','Annotation','pegRNA-to-edit_distance','Nick_index','ngRNA-to-pegRNA_distance','PBS_length','PBS_GC_content','RTT_length','RTT_GC_content','First_extension_nucleotide']] # Less important metadata

    return pegRNAs,ngRNAs

''' epegRNA_linkers: Generate epegRNA linkers between PBS and 3' hairpin motif & finish annotations
        pegRNAs: pegRNAs DataFrame
        epegRNA_motif_sequence: epegRNA motif sequence (Optional, Default: tevopreQ1)
        checkpoint_dir: Checkpoint directory (Optional)
        checkpoint_file: Checkpoint file name (Optional)
        checkpoint_pt: Previous checkpoint path (Optional)
    Dependencies: pandas,pegLIT,io
'''
def epegRNA_linkers(pegRNAs: pd.DataFrame, epegRNA_motif_sequence: str='CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA',
                    checkpoint_dir: str=None, checkpoint_file=None, checkpoint_pt: str=None):
    
    # Get or make checkpoint DataFrame
    if checkpoint_pt: checkpoint = pd.DataFrame(columns=['pegRNA_number','Linker_sequence'])
    else: checkpoint = io.get(pt=checkpoint_pt)

    # Generate epegRNA linkers between PBS and 3' hairpin motif
    linkers = []
    for i in range(len(pegRNAs)):
        if i>=len(checkpoint):
            linkers.extend(pegLIT.pegLIT(seq_spacer=pegRNAs.iloc[i]['Spacer_sequence'],seq_scaffold=pegRNAs.iloc[i]['Scaffold_sequence'],
                                        seq_template=pegRNAs.iloc[i]['RTT_sequence'],seq_pbs=pegRNAs.iloc[i]['PBS_sequence'],
                                        seq_motif=epegRNA_motif_sequence))
            if checkpoint_dir is not None & checkpoint_file is not None: # Save checkpoints
                checkpoint = pd.concat([checkpoint,pd.DataFrame({'pegRNA_number': [i], 'Linker_sequence': [linkers[i]]})])
                io.save(dir=checkpoint_dir,file=checkpoint_file,obj=checkpoint)
            print(f'Status: {i} out of {len(pegRNAs)}')
    
    # Generate epegRNAs
    pegRNAs['Linker_sequence'] = linkers
    pegRNAs['Motif_sequence'] = [epegRNA_motif_sequence]*len(pegRNAs)
    epegRNAs = pegRNAs[['pegRNA_number','gRNA_type','Strand','Edit', # Important metadata
                        'Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence','Linker_sequence','Motif_sequence', # Sequence information
                        'Target_name','Target_sequence','Spacer_GC_content','PAM_sequence','Extension_sequence','Annotation','pegRNA-to-edit_distance','Nick_index','ngRNA-to-pegRNA_distance','PBS_length','PBS_GC_content','RTT_length','RTT_GC_content','First_extension_nucleotide']] # Less important metadata
    return epegRNAs

''' shared_sequences: Reduce PE library into shared spacers and PBS sequences.
        pegRNAs: pegRNAs DataFrame
        hist_plot: display histogram of reduced PE library (Optional, False)
        hist_dir: directory to save histogram
        hist_file: file name to save histogram
    Dependencies: pandas,plot
'''
def shared_sequences(pegRNAs: pd.DataFrame, hist_plot:bool=True, hist_dir: str=None, hist_file=None, **kwargs):
    
    # Reduce PE library to the set shared of spacers and PBS motifs
    shared = {(pegRNAs.iloc[i]['Spacer_sequence'],pegRNAs.iloc[i]['PBS_sequence']) for i in range(len(pegRNAs))}
    shared_lib = pd.DataFrame(columns=['pegRNA_numbers','Strand','Edits','Spacer_sequence','PBS_sequence'])
    for (spacer,pbs) in shared:
        shared_pegRNAs = pegRNAs[(pegRNAs['Spacer_sequence']==spacer)&(pegRNAs['PBS_sequence']==pbs)]
        shared_lib = pd.concat([shared_lib,pd.DataFrame({'pegRNA_numbers': [shared_pegRNAs['pegRNA_number'].to_list()],
                                                        'Strand': [shared_pegRNAs.iloc[0]['Strand']],
                                                        'Edits': [shared_pegRNAs['Edit'].to_list()],
                                                        'Spacer_sequence': [spacer],
                                                        'PBS_sequence': [pbs]})]).reset_index(drop=True)
    
    # Find shared AAs within the reduced PE library
    aa_numbers_ls=[]
    aa_numbers_min_ls=[]
    aa_numbers_max_ls=[]
    continous_ls=[]
    for edits in shared_lib['Edits']:
        aa_numbers = {int(edit[1:-1]) for edit in edits}
        aa_numbers_min = min(aa_numbers)
        aa_numbers_max = max(aa_numbers)
        if aa_numbers == set(range(aa_numbers_min,aa_numbers_max+1)): continous=True
        else: continous=False
        aa_numbers_ls.append(sorted(aa_numbers))
        aa_numbers_min_ls.append(aa_numbers_min)
        aa_numbers_max_ls.append(aa_numbers_max)
        continous_ls.append(continous)
    shared_lib['AA_numbers']=aa_numbers_ls
    shared_lib['AA_numbers_min']=aa_numbers_min_ls
    shared_lib['AA_numbers_max']=aa_numbers_max_ls
    shared_lib['AA_numbers_continuous']=continous_ls
    
    if hist_plot: # Generate histogram
        shared_hist = pd.DataFrame()
        for i,aa_numbers in enumerate(shared_lib['AA_numbers']):
            shared_hist = pd.concat([shared_hist,pd.DataFrame({'Group_Spacer_PBS': [f'{str(i)}_{shared_lib.iloc[i]["Spacer_sequence"]}_{shared_lib.iloc[i]["PBS_sequence"]}']*len(aa_numbers),
                                                               'AA_number': aa_numbers})]).reset_index(drop=True)
        p.dist(typ='hist',df=shared_hist,x='AA_number',cols='Group_Spacer_PBS',x_axis='AA number',title='Shared Spacers & PBS Sequences in the PE Library',
               x_axis_dims=(min(shared_hist['AA_number']),max(shared_hist['AA_number'])),y_axis_dims=(0,5),
               legend_loc='upper center',legend_bbox_to_anchor=(0.5, -0.1),dir=hist_dir,file=hist_file,legend_ncol=2,**kwargs)

    return shared_lib