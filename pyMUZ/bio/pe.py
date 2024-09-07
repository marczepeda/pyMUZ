### pe.py ###
# Author: Marc Zepeda
# Date: 2024-08-31

# Import packages
import pandas as pd
import numpy as np
import os
import re
from Bio.Seq import Seq
from ..bio import pegLIT as pegLIT
from ..gen import io as io
from ..gen import plot as p

# Biological Dictionaries
dna_aa_codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

aa_dna_codon_table = {
    "F": ["TTT", "TTC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "Y": ["TAT", "TAC"],
    "*": ["TAA", "TAG", "TGA"],  # Stop codons
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],  # Start codon
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "G": ["GGT", "GGC", "GGA", "GGG"]
}

# PrimeDesign Methods
''' PrimeDesignInput: Creates PrimeDesign input file.
        target_name: name of target
        target_sequence: in-frame nucleotide sequence with (saturation mutagensis region)
        dir: name of the output directory 
        file: name of the output file
    Dependencies: pandas,io
    Reference: https://github.com/pinellolab/PrimeDesign/tree/master/PrimeDesign
'''
def PrimeDesignInput(target_name: str, target_sequence: str, dir: str='.', file: str='PrimeDesignInput.csv'):
    io.save(dir=dir,file=file,obj=pd.DataFrame({'target_name': [target_name],'target_sequence': [target_sequence]}))

''' PrimeDesign: Run PrimeDesign using Docker (NEED TO BE RUNNING DESKTOP APP)
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
    if pbs_length_list: cmd += f' -pbs {" ".join([str(val) for val in pbs_length_list])}' # Append optional parameters
    if rtt_length_list: cmd += f' -rtt {" ".join([str(val) for val in rtt_length_list])}'
    if nicking_distance_minimum!=0: cmd += f' -nick_dist_min {str(nicking_distance_minimum)}' 
    if nicking_distance_maximum!=100: cmd += f' -nick_dist_max {str(nicking_distance_maximum)}'
    if filter_c1_extension: cmd += f' -filter_c1 {str(filter_c1_extension)}'
    if silent_mutation: cmd += f' -silent_mut'
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
    pegRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='pegRNA'].reset_index(drop=True)
    ngRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='ngRNA'].reset_index(drop=True)

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

# pegRNA Methods
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
    shared_pegRNAs_lib = pd.DataFrame(columns=['pegRNA_numbers','Strand','Edits','Spacer_sequence','PBS_sequence'])
    for (spacer,pbs) in shared:
        shared_pegRNAs = pegRNAs[(pegRNAs['Spacer_sequence']==spacer)&(pegRNAs['PBS_sequence']==pbs)]
        shared_pegRNAs_lib = pd.concat([shared_pegRNAs_lib,pd.DataFrame({'pegRNA_numbers': [shared_pegRNAs['pegRNA_number'].to_list()],
                                                        'Strand': [shared_pegRNAs.iloc[0]['Strand']],
                                                        'Edits': [shared_pegRNAs['Edit'].to_list()],
                                                        'Spacer_sequence': [spacer],
                                                        'PBS_sequence': [pbs]})]).reset_index(drop=True)
    
    # Find shared AAs within the reduced PE library
    aa_numbers_ls=[]
    aa_numbers_min_ls=[]
    aa_numbers_max_ls=[]
    continous_ls=[]
    for edits in shared_pegRNAs_lib['Edits']:
        aa_numbers = {int(edit[1:-1]) for edit in edits}
        aa_numbers_min = min(aa_numbers)
        aa_numbers_max = max(aa_numbers)
        if aa_numbers == set(range(aa_numbers_min,aa_numbers_max+1)): continous=True
        else: continous=False
        aa_numbers_ls.append(sorted(aa_numbers))
        aa_numbers_min_ls.append(aa_numbers_min)
        aa_numbers_max_ls.append(aa_numbers_max)
        continous_ls.append(continous)
    shared_pegRNAs_lib['AA_numbers']=aa_numbers_ls
    shared_pegRNAs_lib['AA_numbers_min']=aa_numbers_min_ls
    shared_pegRNAs_lib['AA_numbers_max']=aa_numbers_max_ls
    shared_pegRNAs_lib['AA_numbers_continuous']=continous_ls
    
    if hist_plot: # Generate histogram
        shared_hist = pd.DataFrame()
        for i,aa_numbers in enumerate(shared_pegRNAs_lib['AA_numbers']):
            shared_hist = pd.concat([shared_hist,pd.DataFrame({'Group_Spacer_PBS': [f'{str(i)}_{shared_pegRNAs_lib.iloc[i]["Spacer_sequence"]}_{shared_pegRNAs_lib.iloc[i]["PBS_sequence"]}']*len(aa_numbers),
                                                               'AA_number': aa_numbers})]).reset_index(drop=True)
        p.dist(typ='hist',df=shared_hist,x='AA_number',cols='Group_Spacer_PBS',x_axis='AA number',title='Shared Spacers & PBS Sequences in the PE Library',
               x_axis_dims=(min(shared_hist['AA_number']),max(shared_hist['AA_number'])),y_axis_dims=(0,5),
               legend_loc='upper center',legend_bbox_to_anchor=(0.5, -0.1),dir=hist_dir,file=hist_file,legend_ncol=2,**kwargs)

    return shared_pegRNAs_lib


''' get_codons: Returns all codons within a specified frame for a nucleotide sequence
        sequence: nucletide sequence
        frame: codon frame (0, 1, or 2)
'''
def get_codons(sequence,frame:int=0):
    return [sequence[i:i+3] for i in range(frame, len(sequence) - 2, 3)]

''' get_codon_frames: Returns all codon frames for a nucleotide sequence
        seqeuence: nucleotide sequence
'''
def get_codon_frames(sequence):
    return [get_codons(sequence,frame) for frame in range(3)]

''' is_sublist_in_order: Returns if each element in the sub list appears in the correct order in the main list
        main_list: search for it here
        sub_list: find this list
'''
def is_sublist_in_order(main_list, sub_list):
    it = iter(main_list) # Initialize an iterator for the sub_list
    return all(item in it for item in sub_list) # Check if each element in sub_list appears in the correct order in main_list

''' RTT_designer: Design all possible RTT for given spacer & PBS (WT, single insertions, & single deletions)
        pegRNAs: pegRNAs DataFrame
        file (str): Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)
        aa_index: 1st amino acid in target sequence index (Optional, Default: start codon = 1)
        RTT_length: Reverse transcriptase template length (bp)
    Dependencies: io,Bio.Seq.Seq,shared_sequences(),get_codons(),get_codon_frames(),is_sublist_in_order(),aa_dna_codon_table
'''
def RTT_designer(pegRNAs: pd.DataFrame, file: str, rtt_length: int=21, aa_index: int=1):
    
    # Get reference sequence & codons (+ reverse complement)
    target_sequence = io.get(file).iloc[0]['target_sequence']
    seq = Seq(target_sequence.split('(')[1].split(')')[0]) # Break apart target sequences
    flank5 = Seq(target_sequence.split('(')[0])
    flank3 = Seq(target_sequence.split(')')[1])
    seq_nuc = flank5 + seq + flank3  # Join full nucleotide reference sequence
    rc_seq = Seq.reverse_complement(seq) # In-frame reverse complement sequence
    rc_seq_nuc = Seq.reverse_complement(seq_nuc) # Full nucleotide reference reverse complement sequence
    seq_prot = Seq.translate(seq) # In-frame amino acid sequence
    aa_indexes = list(np.arange(aa_index,aa_index+len(seq_prot))) # In-frame amino acid indexes
    codons = get_codons(seq) # Codons
    codons_flank5 = get_codons(flank5,len(flank5)%3) # Codons in-frame flank 5
    codons_flank3 = get_codons(flank3) # Codons in-frame flank 3
    extended_codons = codons_flank5 + codons + codons_flank3 # Codons including flank 5 and flank 3
    extended_codons_nuc = Seq('').join(extended_codons) # Join codons into full in-frame nucleotide sequence
    extended_codons_prot = Seq.translate(extended_codons_nuc) # Translate to full in-frame protein sequence
    extended_codons_aa_indexes = list(np.arange(aa_index-len(codons_flank5),aa_index-len(codons_flank5)+len(extended_codons_prot))) # Obtain full in-frame amino acid indexes


    print(f'FWD Ref: {seq_nuc}')
    print(f'REV Ref: {rc_seq_nuc}')
    print(f'Nucleotides: {seq}')
    print(f'Amino Acids: {seq_prot}\n')

    # Obtain shared spacer and PBS sequences 
    shared_pegRNAs_lib = shared_sequences(pegRNAs=pegRNAs,hist_plot=False)

    # Obtain WT RTT, single insertions, and single deletions
    wildtypes = pd.DataFrame()
    insertions = pd.DataFrame()
    deletions = pd.DataFrame()
    for j,pbs in enumerate(shared_pegRNAs_lib['PBS_sequence']): # Iterate through primer binding sites

        if shared_pegRNAs_lib.iloc[j]['Strand']=='+': # Spacer: + strand; PBS & RTT: - strand
            
            # Obtain WT RTT from - strand
            pbs_j = rc_seq_nuc.find(pbs)
            rtt_wt = rc_seq_nuc[pbs_j-rtt_length:pbs_j]
            wildtypes = pd.concat([wildtypes,
                                   pd.DataFrame({'pegRNA_number': [j],
                                                 'gRNA_type': ['pegRNA'],
                                                 'Strand': [shared_pegRNAs_lib.iloc[j]['Strand']],
                                                 'Edit': [None],
                                                 'Spacer_sequence': [shared_pegRNAs_lib.iloc[j]['Spacer_sequence']],
                                                 'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']],
                                                 'RTT_sequence': [str(rtt_wt)],
                                                 'PBS_sequence': [pbs],
                                                 'Target_name': [None],
                                                 'Target_sequence': [None],
                                                 'Spacer_GC_content': [None], 
                                                 'PAM_sequence': [None],
                                                 'Extension_sequence': [''.join([str(rtt_wt),pbs])], 
                                                 'Annotation': ['wildtype'], 
                                                 'pegRNA-to-edit_distance': [None],
                                                 'Nick_index': [None],
                                                 'ngRNA-to-pegRNA_distance': [None], 
                                                 'PBS_length': [len(pbs)],
                                                 'PBS_GC_content': [None], 
                                                 'RTT_length': [rtt_length], 
                                                 'RTT_GC_content': [None],
                                                 'First_extension_nucleotide': [rtt_wt[0]]})]).reset_index(drop=True)
            
            # Obtain reverse complement WT RTT in-frame from + strand
            rc_rtt_wt = Seq.reverse_complement(rtt_wt) # reverse complement of rtt (+ strand)
            rc_rtt_codon_frames = get_codon_frames(rc_rtt_wt) # codons
            for i,rc_rtt_codon_frame in enumerate(rc_rtt_codon_frames): # Search for in-frame nucleotide sequence
                if is_sublist_in_order(codons,rc_rtt_codon_frame): # Codon frame from reverse complement of rtt matches codons of in-frame nucleotide sequence
                    rc_rtt_wt_inframe_nuc_codons_flank5 = rc_rtt_wt[:i] # Save codon frame flank 5'
                    rc_rtt_wt_inframe_nuc_codons = rc_rtt_codon_frame # Save codon frame
                    rc_rtt_wt_inframe_nuc_codons_flank3 = rc_rtt_wt[i+3*len(rc_rtt_codon_frame):] # Save codon frame flank 3'
                    rc_rtt_wt_inframe_nuc = Seq('').join(rc_rtt_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rc_rtt_wt_inframe_prot = Seq.translate(rc_rtt_wt_inframe_nuc) # Translate to in-frame protein sequence
                    rc_rtt_wt_inframe_prot_indexes = aa_indexes[seq_prot.find(rc_rtt_wt_inframe_prot):seq_prot.find(rc_rtt_wt_inframe_prot)+len(rc_rtt_wt_inframe_prot)] # Obtain correponding aa indexes
                    print('Used codons')
                    break
                elif is_sublist_in_order(extended_codons,rc_rtt_codon_frame): # Codon frame from reverse complement of rtt matches extended codons of in-frame nucleotide sequence
                    rc_rtt_wt_inframe_nuc_codons_flank5 = rc_rtt_wt[:i] # Save codon frame flank 5'
                    rc_rtt_wt_inframe_nuc_codons = rc_rtt_codon_frame # Save codon frame
                    rc_rtt_wt_inframe_nuc_codons_flank3 = rc_rtt_wt[i+3*len(rc_rtt_codon_frame):] # Save codon frame flank 3'
                    rc_rtt_wt_inframe_nuc = Seq('').join(rc_rtt_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rc_rtt_wt_inframe_prot = Seq.translate(rc_rtt_wt_inframe_nuc) # Translate to in-frame protein sequence
                    rc_rtt_wt_inframe_prot_indexes = extended_codons_aa_indexes[extended_codons_prot.find(rc_rtt_wt_inframe_prot):extended_codons_prot.find(rc_rtt_wt_inframe_prot)+len(rc_rtt_wt_inframe_prot)] # Obtain correponding aa indexes
                    print('Used extended codons')
                    break
            
            print(f'Strand: {shared_pegRNAs_lib.iloc[j]["Strand"]}')
            print(f'Nucleotides: {rc_rtt_wt}')
            print(f'Nucleotides 5\' of Codons: {rc_rtt_wt_inframe_nuc_codons_flank5}')
            print(f'Nucleotides Codons: {rc_rtt_wt_inframe_nuc_codons}')
            print(f'Nucleotides 3\' of Codons: {rc_rtt_wt_inframe_nuc_codons_flank3}')
            print(f'Nucleotides In-Frame: {rc_rtt_wt_inframe_nuc}')
            print(f'Amino Acids In-Frame: {rc_rtt_wt_inframe_prot}')
            print(f'Amino Acid #s In-Frame: {rc_rtt_wt_inframe_prot_indexes}\n')

            # Obtain single insertion RTTs from - strand
            edits_in = []
            rtts_in = []
            for i in range(len(rc_rtt_wt_inframe_nuc_codons)): # Iterate through all in-frame codon positions
                for codon_table_aa,codon_table_dna in aa_dna_codon_table.items(): # Obtain all possible codon insertions
                    if codon_table_aa!='*': # Remove stop codons
                        edits_in.append(f'{rc_rtt_wt_inframe_prot[i]}{rc_rtt_wt_inframe_prot_indexes[i]}{rc_rtt_wt_inframe_prot[i]}{codon_table_aa}')
                        rtts_in.append(Seq.reverse_complement(Seq('').join([rc_rtt_wt_inframe_nuc_codons_flank5, # Codon frame flank 5'
                                                              Seq('').join(rc_rtt_wt_inframe_nuc_codons[:i+1]), # Codons before insertion
                                                              Seq(codon_table_dna[0]).lower(), # Insertion codon
                                                              Seq('').join(rc_rtt_wt_inframe_nuc_codons[i+1:]), # Codons after insertion
                                                              rc_rtt_wt_inframe_nuc_codons_flank3]))) # Codon frame flank 3'
            
            print(f'Insertions: {edits_in}')
            print(f'Insertion RTTs: {rtts_in}\n')

            insertions = pd.concat([insertions,
                                    pd.DataFrame({'pegRNA_number': [j]*len(edits_in),
                                                  'gRNA_type': ['pegRNA']*len(edits_in),
                                                  'Strand': [shared_pegRNAs_lib.iloc[j]['Strand']]*len(edits_in),
                                                  'Edit': edits_in,
                                                  'Spacer_sequence': [shared_pegRNAs_lib.iloc[j]['Spacer_sequence']]*len(edits_in),
                                                  'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']]*len(edits_in),
                                                  'RTT_sequence': [str(rtt_in) for rtt_in in rtts_in],
                                                  'PBS_sequence': [pbs]*len(edits_in),
                                                  'Target_name': [None]*len(edits_in),
                                                  'Target_sequence': [None]*len(edits_in),
                                                  'Spacer_GC_content': [None]*len(edits_in),
                                                  'PAM_sequence': [None]*len(edits_in),
                                                  'Extension_sequence': [''.join([str(rtt_in),pbs]) for rtt_in in rtts_in], 
                                                  'Annotation': ['insertion']*len(edits_in),
                                                  'pegRNA-to-edit_distance': [None]*len(edits_in),
                                                  'Nick_index': [None]*len(edits_in),
                                                  'ngRNA-to-pegRNA_distance': [None]*len(edits_in),
                                                  'PBS_length': [len(pbs)]*len(edits_in),
                                                  'PBS_GC_content': [None]*len(edits_in),
                                                  'RTT_length': [len(rtt_in) for rtt_in in rtts_in], 
                                                  'RTT_GC_content': [None]*len(edits_in),
                                                  'First_extension_nucleotide': [rtt_in[0] for rtt_in in rtts_in]})]).reset_index(drop=True)

            # Obtain single deletion RTTs from - strand
            edits_del = [f'{str(aa)}{rc_rtt_wt_inframe_prot_indexes[i]}del' for i,aa in enumerate(rc_rtt_wt_inframe_prot) if i<len(rc_rtt_wt_inframe_prot)-1] # Don't want last AA
            rtts_del = [Seq.reverse_complement(Seq('').join([rc_rtt_wt_inframe_nuc_codons_flank5, # Codon frame flank 5'
                                               Seq('').join(rc_rtt_wt_inframe_nuc_codons[:i]), # Codons before deletion
                                               Seq('').join(rc_rtt_wt_inframe_nuc_codons[i+1:]), # Codons after deletion
                                               rc_rtt_wt_inframe_nuc_codons_flank3])) # Codon frame flank 3'
                                               for i in range(len(rc_rtt_wt_inframe_nuc_codons)) if i<len(rc_rtt_wt_inframe_nuc_codons)-1] # Don't want last AA
            
            print(f'Deletions: {edits_del}')
            print(f'Deletion RTTs: {rtts_del}\n\n')

            deletions = pd.concat([deletions,
                                   pd.DataFrame({'pegRNA_number': [j]*len(edits_del),
                                                 'gRNA_type': ['pegRNA']*len(edits_del),
                                                 'Strand': [shared_pegRNAs_lib.iloc[j]['Strand']]*len(edits_del),
                                                 'Edit': edits_del,
                                                 'Spacer_sequence': [shared_pegRNAs_lib.iloc[j]['Spacer_sequence']]*len(edits_del),
                                                 'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']]*len(edits_del),
                                                 'RTT_sequence': [str(rtt_del) for rtt_del in rtts_del],
                                                 'PBS_sequence': [pbs]*len(edits_del),
                                                 'Target_name': [None]*len(edits_del),
                                                 'Target_sequence': [None]*len(edits_del),
                                                 'Spacer_GC_content': [None]*len(edits_del), 
                                                 'PAM_sequence': [None]*len(edits_del),
                                                 'Extension_sequence': [''.join([str(rtt_del),pbs]) for rtt_del in rtts_del], 
                                                 'Annotation': ['deletion']*len(edits_del), 
                                                 'pegRNA-to-edit_distance': [None]*len(edits_del),
                                                 'Nick_index': [None]*len(edits_del),
                                                 'ngRNA-to-pegRNA_distance': [None]*len(edits_del), 
                                                 'PBS_length': [len(pbs)]*len(edits_del),
                                                 'PBS_GC_content': [None]*len(edits_del),
                                                 'RTT_length': [len(rtt_del) for rtt_del in rtts_del], 
                                                 'RTT_GC_content': [None]*len(edits_del),
                                                 'First_extension_nucleotide': [rtt_del[0] for rtt_del in rtts_del]})]).reset_index(drop=True)
            
        elif shared_pegRNAs_lib.iloc[j]['Strand']=='-': # Spacer: - strand; PBS & RTT: + strand
            
            # Obtain WT RTT from + strand
            pbs_j = seq_nuc.find(pbs)
            rtt_wt = seq_nuc[pbs_j-rtt_length:pbs_j]
            wildtypes = pd.concat([wildtypes,
                                   pd.DataFrame({'pegRNA_number': [j],
                                                 'gRNA_type': ['pegRNA'],
                                                 'Strand': [shared_pegRNAs_lib.iloc[j]['Strand']],
                                                 'Edit': [None],
                                                 'Spacer_sequence': [shared_pegRNAs_lib.iloc[j]['Spacer_sequence']],
                                                 'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']],
                                                 'RTT_sequence': [str(rtt_wt)],
                                                 'PBS_sequence': [pbs],
                                                 'Target_name': [None],
                                                 'Target_sequence': [None],
                                                 'Spacer_GC_content': [None], 
                                                 'PAM_sequence': [None],
                                                 'Extension_sequence': [''.join([str(rtt_wt),pbs])], 
                                                 'Annotation': ['wildtype'], 
                                                 'pegRNA-to-edit_distance': [None],
                                                 'Nick_index': [None],
                                                 'ngRNA-to-pegRNA_distance': [None], 
                                                 'PBS_length': [len(pbs)],
                                                 'PBS_GC_content': [None], 
                                                 'RTT_length': [rtt_length], 
                                                 'RTT_GC_content': [None],
                                                 'First_extension_nucleotide': [rtt_wt[0]]})]).reset_index(drop=True)
            
            # Obtain WT RTT in-frame from + strand
            rtt_codon_frames = get_codon_frames(rtt_wt) # codons
            for i,rtt_codon_frame in enumerate(rtt_codon_frames): # Search for in-frame nucleotide sequence
                if is_sublist_in_order(codons,rtt_codon_frame): # Codon frame from rtt matches codons of in-frame nucleotide sequence
                    rtt_wt_inframe_nuc_codons_flank5 = rtt_wt[:i] # Save codon frame flank 5'
                    rtt_wt_inframe_nuc_codons = rtt_codon_frame # Save codon frame
                    rtt_wt_inframe_nuc_codons_flank3 = rtt_wt[i+3*len(rtt_codon_frame):] # Save codon frame flank 3'
                    rtt_wt_inframe_nuc = Seq('').join(rtt_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rtt_wt_inframe_prot = Seq.translate(rtt_wt_inframe_nuc) # Translate to in-frame protein sequence
                    rtt_wt_inframe_prot_indexes = aa_indexes[seq_prot.find(rtt_wt_inframe_prot):seq_prot.find(rtt_wt_inframe_prot)+len(rtt_wt_inframe_prot)] # Obtain correponding aa indexes
                    print('Used codons')
                    break
                elif is_sublist_in_order(extended_codons,rtt_codon_frame): # Codon frame from rtt matches extended codons of in-frame nucleotide sequence
                    rtt_wt_inframe_nuc_codons_flank5 = rtt_wt[:i] # Save codon frame flank 5'
                    rtt_wt_inframe_nuc_codons = rtt_codon_frame # Save codon frame
                    rtt_wt_inframe_nuc_codons_flank3 = rtt_wt[i+3*len(rtt_codon_frame):] # Save codon frame flank 3'
                    rtt_wt_inframe_nuc = Seq('').join(rtt_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rtt_wt_inframe_prot = Seq.translate(rtt_wt_inframe_nuc) # Translate to in-frame protein sequence
                    rtt_wt_inframe_prot_indexes = extended_codons_aa_indexes[extended_codons_prot.find(rtt_wt_inframe_prot):extended_codons_prot.find(rtt_wt_inframe_prot)+len(rtt_wt_inframe_prot)] # Obtain correponding aa indexes
                    print('Used extended codons')
            
            print(f'Strand: {shared_pegRNAs_lib.iloc[j]["Strand"]}')
            print(f'Nucleotides: {rtt_wt}')
            print(f'Nucleotides 5\' of Codons: {rtt_wt_inframe_nuc_codons_flank5}')
            print(f'Nucleotides Codons: {rtt_wt_inframe_nuc_codons}')
            print(f'Nucleotides 3\' of Codons: {rtt_wt_inframe_nuc_codons_flank3}')
            print(f'Nucleotides In-Frame: {rtt_wt_inframe_nuc}')
            print(f'Amino Acids In-Frame: {rtt_wt_inframe_prot}')
            print(f'Amino Acid #s In-Frame: {rtt_wt_inframe_prot_indexes}\n')

            # Obtain single insertion RTTs from + strand
            edits_in = []
            rtts_in = []
            for i in range(len(rtt_wt_inframe_nuc_codons)): # Iterate through all in-frame codon positions
                for codon_table_aa,codon_table_dna in aa_dna_codon_table.items(): # Obtain all possible codon insertions
                    if codon_table_aa!='*': # Remove stop codons
                        edits_in.append(f'{rtt_wt_inframe_prot[i]}{rtt_wt_inframe_prot_indexes[i]}{rtt_wt_inframe_prot[i]}{codon_table_aa}')
                        rtts_in.append(Seq('').join([rtt_wt_inframe_nuc_codons_flank5, # Codon frame flank 5'
                                       Seq('').join(rtt_wt_inframe_nuc_codons[:i+1]), # Codons before insertion
                                       Seq(codon_table_dna[0]).lower(), # Insertion codon
                                       Seq('').join(rtt_wt_inframe_nuc_codons[i+1:]), # Codons after insertion
                                       rtt_wt_inframe_nuc_codons_flank3])) # Codon frame flank 3'

            print(f'Insertions: {edits_in}')
            print(f'Insertion RTTs: {rtts_in}\n')

            insertions = pd.concat([insertions,
                                    pd.DataFrame({'pegRNA_number': [j]*len(edits_in),
                                                  'gRNA_type': ['pegRNA']*len(edits_in),
                                                  'Strand': [shared_pegRNAs_lib.iloc[j]['Strand']]*len(edits_in),
                                                  'Edit': edits_in,
                                                  'Spacer_sequence': [shared_pegRNAs_lib.iloc[j]['Spacer_sequence']]*len(edits_in),
                                                  'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']]*len(edits_in),
                                                  'RTT_sequence': [str(rtt_in) for rtt_in in rtts_in],
                                                  'PBS_sequence': [pbs]*len(edits_in),
                                                  'Target_name': [None]*len(edits_in),
                                                  'Target_sequence': [None]*len(edits_in),
                                                  'Spacer_GC_content': [None]*len(edits_in),
                                                  'PAM_sequence': [None]*len(edits_in),
                                                  'Extension_sequence': [''.join([str(rtt_in),pbs]) for rtt_in in rtts_in], 
                                                  'Annotation': ['insertion']*len(edits_in),
                                                  'pegRNA-to-edit_distance': [None]*len(edits_in),
                                                  'Nick_index': [None]*len(edits_in),
                                                  'ngRNA-to-pegRNA_distance': [None]*len(edits_in),
                                                  'PBS_length': [len(pbs)]*len(edits_in),
                                                  'PBS_GC_content': [None]*len(edits_in),
                                                  'RTT_length': [len(rtt_in) for rtt_in in rtts_in], 
                                                  'RTT_GC_content': [None]*len(edits_in),
                                                  'First_extension_nucleotide': [rtt_in[0] for rtt_in in rtts_in]})]).reset_index(drop=True)

            # Obtain single deletion RTTs from + strand
            edits_del = [f'{str(aa)}{rtt_wt_inframe_prot_indexes[i]}del' for i,aa in enumerate(rtt_wt_inframe_prot) if i<len(rtt_wt_inframe_prot)-1] # Don't want last AA
            rtts_del = [Seq('').join([rtt_wt_inframe_nuc_codons_flank5, # Codon frame flank 5'
                        Seq('').join(rtt_wt_inframe_nuc_codons[:i]), # Codons before deletion
                        Seq('').join(rtt_wt_inframe_nuc_codons[i+1:]), # Codons after deletion
                        rtt_wt_inframe_nuc_codons_flank3]) # Codon frame flank 3'
                        for i in range(len(rtt_wt_inframe_nuc_codons)) if i<len(rtt_wt_inframe_nuc_codons)-1] # Don't want last AA
            
            print(f'Deletions: {edits_del}')
            print(f'Deletion RTTs: {rtts_del}\n\n')

            deletions = pd.concat([deletions,
                                   pd.DataFrame({'pegRNA_number': [j]*len(edits_del),
                                                 'gRNA_type': ['pegRNA']*len(edits_del),
                                                 'Strand': [shared_pegRNAs_lib.iloc[j]['Strand']]*len(edits_del),
                                                 'Edit': edits_del,
                                                 'Spacer_sequence': [shared_pegRNAs_lib.iloc[j]['Spacer_sequence']]*len(edits_del),
                                                 'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']]*len(edits_del),
                                                 'RTT_sequence': [str(rtt_del) for rtt_del in rtts_del],
                                                 'PBS_sequence': [pbs]*len(edits_del),
                                                 'Target_name': [None]*len(edits_del),
                                                 'Target_sequence': [None]*len(edits_del),
                                                 'Spacer_GC_content': [None]*len(edits_del), 
                                                 'PAM_sequence': [None]*len(edits_del),
                                                 'Extension_sequence': [''.join([str(rtt_del),pbs]) for rtt_del in rtts_del], 
                                                 'Annotation': ['deletion']*len(edits_del), 
                                                 'pegRNA-to-edit_distance': [None]*len(edits_del),
                                                 'Nick_index': [None]*len(edits_del),
                                                 'ngRNA-to-pegRNA_distance': [None]*len(edits_del), 
                                                 'PBS_length': [len(pbs)]*len(edits_del),
                                                 'PBS_GC_content': [None]*len(edits_del),
                                                 'RTT_length': [len(rtt_del) for rtt_del in rtts_del], 
                                                 'RTT_GC_content': [None]*len(edits_del),
                                                 'First_extension_nucleotide': [rtt_del[0] for rtt_del in rtts_del]})]).reset_index(drop=True)

        else: 
            print('Error: Strand column can only have "+" and "-".')
            return

    # Combine wildtype, substitution, insertion, deletion libraries
    return pd.concat([pegRNAs,wildtypes,insertions,deletions]).reset_index(drop=True)

''' compare_RTTs: Compares RTT outcome to WT RTT outcome and returns observed it.
        rtt_prot: RTT outcome, which is AA sequence
        rtt_prot_indexes: RTT outcome indexes, which is AA sequence indexes
        rtt_wt_prot_indexes: RTT WT outcome, which is WT AA sequence
        rtt_wt_prot_indexes: RTT WT outcome indexes, which is WT AA sequence indexes
        annotation: PrimeDesign output, insertion, deletion, or wildtype
        strand: Prime Design output (i.e., "+" or "-")
'''
def compare_RTTs(rtt_prot,rtt_prot_indexes: list,rtt_wt_prot, rtt_wt_prot_indexes: list,annotation: str,strand: str):
    
    # Determine edit category
    if (rtt_prot==rtt_wt_prot)&(rtt_prot_indexes==rtt_wt_prot_indexes): # Wildtype
        print('WT')
        return None
    
    elif annotation=='insertion': # Insertion
        if strand=='+':
            for i in range(len(rtt_wt_prot)): # Find difference starting from 5' end (N-term)
                if rtt_prot[i] != rtt_wt_prot[i]:
                    i_diff = i
                    break
            edit = f'{rtt_wt_prot[i_diff-1]}{rtt_wt_prot_indexes[i_diff-1]}{"".join(rtt_prot[i_diff-1:i_diff+1])}'
            print(f'Insertion Edit: {edit}\n')
            return edit
        elif strand=='-':
            for i in range(len(rtt_wt_prot)-1,-1,-1): # Find difference starting from 3' end (C-term)
                if rtt_prot[i] != rtt_wt_prot[i]:
                    i_diff = i
                    break
            edit = f'{rtt_wt_prot[i_diff]}{rtt_wt_prot_indexes[i_diff]}{"".join(rtt_prot[i_diff-1:i_diff+1])}'
            print(f'Insertion Edit: {edit}\n')
            return edit
    
    elif annotation=='deletion': # Deletion
        if strand=='+':
            for i in range(len(rtt_wt_prot)): # Find difference starting from 5' end (N-term)
                if rtt_prot[i] != rtt_wt_prot[i]:
                    i_diff = i
                    break
            edit = f'{rtt_wt_prot[i_diff]}{rtt_wt_prot_indexes[i_diff]}del'
            print(f'Deletion Edit: {edit}\n')
            return edit
        elif strand=='-':
            for i in range(len(rtt_wt_prot)-1,-1,-1): # Find difference starting from 3' end (C-term)
                if rtt_prot[i] != rtt_wt_prot[i]:
                    i_diff = i
                    break
            edit = f'{rtt_wt_prot[i_diff]}{rtt_wt_prot_indexes[i_diff]}del'
            print(f'Deletion Edit: {edit}\n')
            return edit
        else: print('Error: Strand must be "+" or "-"')

    else: # Substitution
        for i in range(len(rtt_prot)): # Find difference
            if rtt_prot[i] != rtt_wt_prot[i]:
                i_diff = i
                break
        edit = f'{rtt_wt_prot[i_diff]}{rtt_wt_prot_indexes[i_diff]}{rtt_prot[i_diff]}'
        print(f'Substitution Edit: {edit}\n')
        return edit
        

''' pegRNAs_tester: Confirm that pegRNAs should create the predicted edit
        pegRNAs: pegRNAs DataFrame
        file (str): Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)
        aa_index: 1st amino acid in target sequence index (Optional, Default: start codon = 1)
    Dependencies: Bio.Seq.Seq,numpy,pandas,is_sublist_in_order(),get_codons(),get_codon_frames(),compare_RTTs()
    Common errors for calling edits...
        deletions:
            (1) single deletion of repeated AA results in wrong AA #
            (2) N-terminus AA deletion results in wildtype
            (3) AA deletion outside of target sequence boundary results in wrong AA #
'''
def pegRNAs_tester(pegRNAs: pd.DataFrame, file: str, aa_index: int=1):
    
    # Catch all stop codons that are written as "X" instead of "*"
    pegRNAs['Edit'] = pegRNAs['Edit'].replace('X', '*', regex=True)

    # Get reference sequence & codons (+ reverse complement)
    target_sequence = io.get(file).iloc[0]['target_sequence']
    seq = Seq(target_sequence.split('(')[1].split(')')[0]) # Break apart target sequences
    flank5 = Seq(target_sequence.split('(')[0])
    flank3 = Seq(target_sequence.split(')')[1])
    seq_nuc = flank5 + seq + flank3  # Join full nucleotide reference sequence
    rc_seq_nuc = Seq.reverse_complement(seq_nuc) # Full nucleotide reference reverse complement sequence
    seq_prot = Seq.translate(seq) # In-frame amino acid sequence
    aa_indexes = list(np.arange(aa_index,aa_index+len(seq_prot))) # In-frame amino acid indexes
    codons = get_codons(seq) # Codons
    codons_flank5 = get_codons(flank5,len(flank5)%3) # Codons in-frame flank 5
    codons_flank3 = get_codons(flank3) # Codons in-frame flank 3
    extended_codons = codons_flank5 + codons + codons_flank3 # Codons including flank 5 and flank 3
    extended_codons_nuc = Seq('').join(extended_codons) # Join codons into full in-frame nucleotide sequence
    extended_codons_prot = Seq.translate(extended_codons_nuc) # Translate to full in-frame protein sequence
    extended_codons_aa_indexes = list(np.arange(aa_index-len(codons_flank5),aa_index-len(codons_flank5)+len(extended_codons_prot))) # Obtain full in-frame amino acid indexes

    print(f'FWD Ref: {seq_nuc}')
    print(f'REV Ref: {rc_seq_nuc}')
    print(f'Nucleotides: {seq}')
    print(f'Amino Acids: {seq_prot}\n')

    # Get expected edits in the pegRNA library 
    edits_substitutions=[]
    edits_insertions=[]
    edits_deletions=[]
    for i,aa in enumerate(seq_prot):
        edits_substitutions.extend([f'{aa}{str(aa_indexes[i])}{aa2}' for aa2 in aa_dna_codon_table if (aa2!='*')&(aa2!=aa)])
        edits_insertions.extend([f'{aa}{str(aa_indexes[i])}{aa}{aa2}' for aa2 in aa_dna_codon_table if aa2!='*'])
        edits_deletions.append(f'{aa}{str(aa_indexes[i])}del')
    edits_substitutions_set = set(edits_substitutions)
    edits_insertions_set = set(edits_insertions)
    edits_deletions_set = set(edits_deletions)
    
    print(f'Expected Edits in pegRNA library...\nSubstitutions: {edits_substitutions}\nInsertions: {edits_insertions}\nDeletions: {edits_deletions}')
    print(f'All substitions present: {edits_substitutions_set.issubset(set(pegRNAs["Edit"]))}; missing: {edits_substitutions_set-set(pegRNAs["Edit"])}')
    print(f'All insertions present: {edits_insertions_set.issubset(set(pegRNAs["Edit"]))}; missing: {edits_insertions_set-set(pegRNAs["Edit"])}')
    print(f'All deletions present: {edits_deletions_set.issubset(set(pegRNAs["Edit"]))}; missing: {edits_deletions_set-set(pegRNAs["Edit"])}\n')

    edits = [] # Determine all edits in the pegRNA library
    edit_matches = []
    prot_befores = []
    prot_afters = []
    edit_matches_explain = []
    for j in range(len(pegRNAs)): # Iterate through primer binding sites
        
        print(f'Loop: {j}')
        if pegRNAs.iloc[j]['Strand']=='+': # Spacer: + strand; PBS & RTT: - strand
            
            # Obtain reverse complement PBS in-frame from + strand
            rc_pbs = Seq.reverse_complement(Seq(pegRNAs.iloc[j]['PBS_sequence'])) # reverse complement of pbs (+ strand)
            rc_pbs_codon_frames = get_codon_frames(rc_pbs) # codons
            for i,rc_pbs_codon_frame in enumerate(rc_pbs_codon_frames): # Search for in-frame nucleotide sequence
                if is_sublist_in_order(codons,rc_pbs_codon_frame): # Codon frame from reverse complement of pbs matches codons of in-frame nucleotide sequence
                    rc_pbs_inframe_nuc_codons_flank5 = rc_pbs[:i] # Save codon frame flank 5'
                    rc_pbs_inframe_nuc_codons = rc_pbs_codon_frame # Save codon frame
                    rc_pbs_inframe_nuc_codons_flank3 = rc_pbs[i+3*len(rc_pbs_codon_frame):] # Save codon frame flank 3'
                    rc_pbs_inframe_nuc = Seq('').join(rc_pbs_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rc_pbs_inframe_prot = Seq.translate(rc_pbs_inframe_nuc) # Translate to in-frame protein sequence
                    rc_pbs_inframe_prot_indexes = aa_indexes[seq_prot.find(rc_pbs_inframe_prot):seq_prot.find(rc_pbs_inframe_prot)+len(rc_pbs_inframe_prot)] # Obtain correponding aa indexes
                    print('Used codons')
                    break
                elif is_sublist_in_order(extended_codons,rc_pbs_codon_frame): # Codon frame from reverse complement of pbs matches extended codons of in-frame nucleotide sequence
                    rc_pbs_inframe_nuc_codons_flank5 = rc_pbs[:i] # Save codon frame flank 5'
                    rc_pbs_inframe_nuc_codons = rc_pbs_codon_frame # Save codon frame
                    rc_pbs_inframe_nuc_codons_flank3 = rc_pbs[i+3*len(rc_pbs_codon_frame):] # Save codon frame flank 3'
                    rc_pbs_inframe_nuc = Seq('').join(rc_pbs_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rc_pbs_inframe_prot = Seq.translate(rc_pbs_inframe_nuc) # Translate to in-frame protein sequence
                    rc_pbs_inframe_prot_indexes = extended_codons_aa_indexes[extended_codons_prot.find(rc_pbs_inframe_prot):extended_codons_prot.find(rc_pbs_inframe_prot)+len(rc_pbs_inframe_prot)] # Obtain correponding aa indexes
                    print('Used extended codons')

            print(f'Strand: {pegRNAs.iloc[j]["Strand"]}')
            print(f'Annotation: {pegRNAs.iloc[j]["Annotation"]}')
            print(f'PBS Nucleotides: {rc_pbs}')
            print(f'PBS Nucleotides 5\' of Codons: {rc_pbs_inframe_nuc_codons_flank5}')
            print(f'PBS Nucleotides Codons: {rc_pbs_inframe_nuc_codons}')
            print(f'PBS Nucleotides 3\' of Codons: {rc_pbs_inframe_nuc_codons_flank3}')
            print(f'PBS Nucleotides In-Frame: {rc_pbs_inframe_nuc}')
            print(f'PBS Amino Acids In-Frame: {rc_pbs_inframe_prot}')
            print(f'PBS Amino Acid #s In-Frame: {rc_pbs_inframe_prot_indexes}\n')

            # Obtain RTT amino sequence
            rc_rtt = Seq.reverse_complement(Seq(pegRNAs.iloc[j]['RTT_sequence'])) # reverse complement of rtt (+ strand)
            rc_pbs_inframe_flank3_rtt = rc_pbs_inframe_nuc_codons_flank3+rc_rtt
            rc_pbs_inframe_flank3_rtt_codons = get_codons(rc_pbs_inframe_flank3_rtt)
            rc_pbs_inframe_flank3_rtt_prot = Seq.translate(Seq('').join(rc_pbs_inframe_flank3_rtt))
            rc_pbs_inframe_flank3_rtt_prot_indexes = [rc_pbs_inframe_prot_indexes[-1]+1+i for i in range(len(rc_pbs_inframe_flank3_rtt_prot))]

            print(f'RTT (RC) Nucleotides: {rc_rtt}')
            print(f'PBS Flank 3\' + RTT (RC) Nucleotides: {rc_pbs_inframe_flank3_rtt}')
            print(f'PBS Flank 3\' + RTT (RC) NucleotidesCodons: {rc_pbs_inframe_flank3_rtt_codons}')
            print(f'PBS Flank 3\' + RTT (RC) Amino Acids In-Frame: {rc_pbs_inframe_flank3_rtt_prot}')
            print(f'PBS Flank 3\' + RTT (RC) Amino Acids #s In-Frame: {rc_pbs_inframe_flank3_rtt_prot_indexes}\n')

            # Obtain WT RTT from - strand and reverse complement WT RTT in-frame from + strand
            rc_pbs_j = seq_nuc.find(rc_pbs)
            rc_rtt_wt = seq_nuc[rc_pbs_j+len(rc_pbs):rc_pbs_j+len(rc_pbs)+len(rc_rtt)] # reverse complement of rtt (+ strand)
            rc_pbs_inframe_flank3_rtt_wt = rc_pbs_inframe_nuc_codons_flank3+rc_rtt_wt
            rc_pbs_inframe_flank3_rtt_wt_codon_frames = get_codon_frames(rc_pbs_inframe_flank3_rtt_wt) # codons
            for i,rc_pbs_inframe_flank3_rtt_wt_codon_frame in enumerate(rc_pbs_inframe_flank3_rtt_wt_codon_frames): # Search for in-frame nucleotide sequence
                if is_sublist_in_order(codons,rc_pbs_inframe_flank3_rtt_wt_codon_frame): # Codon frame from reverse complement of rtt matches codons of in-frame nucleotide sequence
                    rc_rtt_wt_inframe_nuc_codons_flank5 = rc_rtt_wt[:i] # Save codon frame flank 5'
                    rc_rtt_wt_inframe_nuc_codons = rc_pbs_inframe_flank3_rtt_wt_codon_frame # Save codon frame
                    rc_rtt_wt_inframe_nuc_codons_flank3 = rc_rtt_wt[i+3*len(rc_pbs_inframe_flank3_rtt_wt_codon_frame):] # Save codon frame flank 3'
                    rc_rtt_wt_inframe_nuc = Seq('').join(rc_pbs_inframe_flank3_rtt_wt_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rc_rtt_wt_inframe_prot = Seq.translate(rc_rtt_wt_inframe_nuc) # Translate to in-frame protein sequence
                    rc_rtt_wt_inframe_prot_indexes = aa_indexes[seq_prot.find(rc_rtt_wt_inframe_prot):seq_prot.find(rc_rtt_wt_inframe_prot)+len(rc_rtt_wt_inframe_prot)] # Obtain correponding aa indexes
                    print('Used codons')
                    break
                elif is_sublist_in_order(extended_codons,rc_pbs_inframe_flank3_rtt_wt_codon_frame): # Codon frame from reverse complement of rtt matches extended codons of in-frame nucleotide sequence
                    rc_rtt_wt_inframe_nuc_codons_flank5 = rc_rtt_wt[:i] # Save codon frame flank 5'
                    rc_rtt_wt_inframe_nuc_codons = rc_pbs_inframe_flank3_rtt_wt_codon_frame # Save codon frame
                    rc_rtt_wt_inframe_nuc_codons_flank3 = rc_rtt_wt[i+3*len(rc_pbs_inframe_flank3_rtt_wt_codon_frame):] # Save codon frame flank 3'
                    rc_rtt_wt_inframe_nuc = Seq('').join(rc_pbs_inframe_flank3_rtt_wt_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rc_rtt_wt_inframe_prot = Seq.translate(rc_rtt_wt_inframe_nuc) # Translate to in-frame protein sequence
                    rc_rtt_wt_inframe_prot_indexes = extended_codons_aa_indexes[extended_codons_prot.find(rc_rtt_wt_inframe_prot):extended_codons_prot.find(rc_rtt_wt_inframe_prot)+len(rc_rtt_wt_inframe_prot)] # Obtain correponding aa indexes
                    print('Used extended codons')
            
            print(f'RTT WT (RC) Nucleotides: {rc_rtt_wt}')
            print(f'PBS Flank 3\' + RTT WT (RC) Nucleotides: {rc_pbs_inframe_flank3_rtt_wt}')
            print(f'PBS Flank 3\' + RTT WT (RC) Nucleotides 5\' of Codons: {rc_rtt_wt_inframe_nuc_codons_flank5}')
            print(f'PBS Flank 3\' + RTT WT (RC) Nucleotides Codons: {rc_rtt_wt_inframe_nuc_codons}')
            print(f'PBS Flank 3\' + RTT WT (RC) Nucleotides 3\' of Codons: {rc_rtt_wt_inframe_nuc_codons_flank3}')
            print(f'PBS Flank 3\' + RTT WT (RC) Nucleotides In-Frame: {rc_rtt_wt_inframe_nuc}')
            print(f'PBS Flank 3\' + RTT WT (RC) Amino Acids In-Frame: {rc_rtt_wt_inframe_prot}')
            print(f'PBS Flank 3\' + RTT WT (RC) Amino Acid #s In-Frame: {rc_rtt_wt_inframe_prot_indexes}\n')

            # Determine edit by comparing RTT to WT RTT
            prot_befores.append(''.join([str(aa) for aa in rc_rtt_wt_inframe_prot]))
            prot_afters.append(''.join([str(aa) for aa in rc_pbs_inframe_flank3_rtt_prot]))
            print(f'Expected edit: {pegRNAs.iloc[j]["Edit"]}')
            edit = compare_RTTs(rtt_prot=rc_pbs_inframe_flank3_rtt_prot,
                                rtt_prot_indexes=rc_pbs_inframe_flank3_rtt_prot_indexes,
                                rtt_wt_prot=rc_rtt_wt_inframe_prot,
                                rtt_wt_prot_indexes=rc_rtt_wt_inframe_prot_indexes,
                                annotation=pegRNAs.iloc[j]['Annotation'],
                                strand=pegRNAs.iloc[j]['Strand'])
            edits.append(edit)

            # Compare expected edit to actual edit
            if edit==pegRNAs.iloc[j]['Edit']: 
                edit_matches.append(True)
                edit_matches_explain.append(None) # No need to explain matching edits
            else: 
                edit_matches.append(False)
                if edit is not None:
                    position = [int(num) for num in re.findall(r'\d+', edit)][0] # Find actual edit position
                    if (position==aa_index-1)&(pegRNAs.iloc[j]['Annotation']=='deletion'): edit_matches_explain.append('Common pegRNAs_tester() error: AA deletion outside of target sequence boundary results in wrong AA #')
                    elif pegRNAs.iloc[j]['Annotation']=='deletion': edit_matches_explain.append('Investigate deletion: common error for pegRNAs_tester() is single deletion of repeated AA results in wrong AA #')
                    else: edit_matches_explain.append('Investigate: unknown error')
                elif pegRNAs.iloc[j]['Annotation']=='deletion': edit_matches_explain.append('Common pegRNAs_tester() error: N-terminus AA deletion results in wildtype') 
                else: edit_matches_explain.append('Investigate: unknown error')

        
        elif pegRNAs.iloc[j]['Strand']=='-': # Spacer: - strand; PBS & RTT: + strand
            
            # Obtain PBS in-frame from + strand
            pbs = Seq(pegRNAs.iloc[j]['PBS_sequence']) # pbs (+ strand)
            pbs_codon_frames = get_codon_frames(pbs) # codons
            for i,pbs_codon_frame in enumerate(pbs_codon_frames): # Search for in-frame nucleotide sequence
                if is_sublist_in_order(codons,pbs_codon_frame): # Codon frame from pbs matches codons of in-frame nucleotide sequence
                    pbs_inframe_nuc_codons_flank5 = pbs[:i] # Save codon frame flank 5'
                    pbs_inframe_nuc_codons = pbs_codon_frame # Save codon frame
                    pbs_inframe_nuc_codons_flank3 = pbs[i+3*len(pbs_codon_frame):] # Save codon frame flank 3'
                    pbs_inframe_nuc = Seq('').join(pbs_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    pbs_inframe_prot = Seq.translate(pbs_inframe_nuc) # Translate to in-frame protein sequence
                    pbs_inframe_prot_indexes = aa_indexes[seq_prot.find(pbs_inframe_prot):seq_prot.find(pbs_inframe_prot)+len(pbs_inframe_prot)] # Obtain correponding aa indexes
                    print('Used codons')
                    break
                elif is_sublist_in_order(extended_codons,pbs_codon_frame): # Codon frame from pbs matches extended codons of in-frame nucleotide sequence
                    pbs_inframe_nuc_codons_flank5 = pbs[:i] # Save codon frame flank 5'
                    pbs_inframe_nuc_codons = pbs_codon_frame # Save codon frame
                    pbs_inframe_nuc_codons_flank3 = pbs[i+3*len(pbs_codon_frame):] # Save codon frame flank 3'
                    pbs_inframe_nuc = Seq('').join(pbs_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    pbs_inframe_prot = Seq.translate(pbs_inframe_nuc) # Translate to in-frame protein sequence
                    pbs_inframe_prot_indexes = extended_codons_aa_indexes[extended_codons_prot.find(pbs_inframe_prot):extended_codons_prot.find(pbs_inframe_prot)+len(pbs_inframe_prot)] # Obtain correponding aa indexes
                    print('Used extended codons')
            
            print(f'Strand: {pegRNAs.iloc[j]["Strand"]}')
            print(f'Annotation: {pegRNAs.iloc[j]["Annotation"]}')
            print(f'PBS Nucleotides: {pbs}')
            print(f'PBS Nucleotides 5\' of Codons: {pbs_inframe_nuc_codons_flank5}')
            print(f'PBS Nucleotides Codons: {pbs_inframe_nuc_codons}')
            print(f'PBS Nucleotides 3\' of Codons: {pbs_inframe_nuc_codons_flank3}')
            print(f'PBS Nucleotides In-Frame: {pbs_inframe_nuc}')
            print(f'PBS Amino Acids In-Frame: {pbs_inframe_prot}')
            print(f'PBS Amino Acid #s In-Frame: {pbs_inframe_prot_indexes}\n')

            # Obtain RTT amino sequence
            rtt = Seq(pegRNAs.iloc[j]['RTT_sequence']) # rtt (+ strand)
            rtt_pbs_inframe_flank5 = rtt+pbs_inframe_nuc_codons_flank5
            rtt_pbs_inframe_flank5_codons = get_codons(rtt_pbs_inframe_flank5,len(rtt_pbs_inframe_flank5)%3)
            rtt_pbs_inframe_flank5_prot = Seq.translate(Seq('').join(rtt_pbs_inframe_flank5_codons))
            rtt_pbs_inframe_flank5_prot_indexes = [pbs_inframe_prot_indexes[0]-len(rtt_pbs_inframe_flank5_prot)+i for i in range(len(rtt_pbs_inframe_flank5_prot))]

            print(f'RTT Nucleotides: {rtt}')
            print(f'RTT + PBS Flank 5\' Nucleotides: {rtt_pbs_inframe_flank5}')
            print(f'RTT + PBS Flank 5\'  Nucleotides Codons: {rtt_pbs_inframe_flank5_codons}')
            print(f'RTT + PBS Flank 5\'  Amino Acids In-Frame: {rtt_pbs_inframe_flank5_prot}')
            print(f'RTT + PBS Flank 5\'  Amino Acids #s In-Frame: {rtt_pbs_inframe_flank5_prot_indexes}\n')

            # Obtain WT RTT from + strand and WT RTT in-frame from + strand
            pbs_j = seq_nuc.find(pbs)
            rtt_wt = seq_nuc[pbs_j-len(rtt):pbs_j]
            rtt_wt_pbs_inframe_flank5 = rtt_wt+pbs_inframe_nuc_codons_flank5
            rtt_wt_pbs_inframe_flank5_codon_frames = get_codon_frames(rtt_wt_pbs_inframe_flank5) # codons
            for i,rtt_wt_pbs_inframe_flank5_codon_frame in enumerate(rtt_wt_pbs_inframe_flank5_codon_frames): # Search for in-frame nucleotide sequence
                if is_sublist_in_order(codons,rtt_wt_pbs_inframe_flank5_codon_frame): # Codon frame from rtt matches codons of in-frame nucleotide sequence
                    rtt_wt_inframe_nuc_codons_flank5 = rtt_wt[:i] # Save codon frame flank 5'
                    rtt_wt_inframe_nuc_codons = rtt_wt_pbs_inframe_flank5_codon_frame # Save codon frame
                    rtt_wt_inframe_nuc_codons_flank3 = rtt_wt[i+3*len(rtt_wt_pbs_inframe_flank5_codon_frame):] # Save codon frame flank 3'
                    rtt_wt_inframe_nuc = Seq('').join(rtt_wt_pbs_inframe_flank5_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rtt_wt_inframe_prot = Seq.translate(rtt_wt_inframe_nuc) # Translate to in-frame protein sequence
                    rtt_wt_inframe_prot_indexes = aa_indexes[seq_prot.find(rtt_wt_inframe_prot):seq_prot.find(rtt_wt_inframe_prot)+len(rtt_wt_inframe_prot)] # Obtain correponding aa indexes
                    print('Used codons')
                    break
                elif is_sublist_in_order(extended_codons,rtt_wt_pbs_inframe_flank5_codon_frame): # Codon frame from rtt matches extended codons of in-frame nucleotide sequence
                    rtt_wt_inframe_nuc_codons_flank5 = rtt_wt[:i] # Save codon frame flank 5'
                    rtt_wt_inframe_nuc_codons = rtt_wt_pbs_inframe_flank5_codon_frame # Save codon frame
                    rtt_wt_inframe_nuc_codons_flank3 = rtt_wt[i+3*len(rtt_wt_pbs_inframe_flank5_codon_frame):] # Save codon frame flank 3'
                    rtt_wt_inframe_nuc = Seq('').join(rtt_wt_pbs_inframe_flank5_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rtt_wt_inframe_prot = Seq.translate(rtt_wt_inframe_nuc) # Translate to in-frame protein sequence
                    rtt_wt_inframe_prot_indexes = extended_codons_aa_indexes[extended_codons_prot.find(rtt_wt_inframe_prot):extended_codons_prot.find(rtt_wt_inframe_prot)+len(rtt_wt_inframe_prot)] # Obtain correponding aa indexes
                    print('Used extended codons')         
            
            print(f'RTT WT Nucleotides: {rtt_wt}')
            print(f'RTT WT + PBS Flank 5\' Nucleotides: {rtt_wt_pbs_inframe_flank5}')
            print(f'RTT WT + PBS Flank 5\' Nucleotides 5\' of Codons: {rtt_wt_inframe_nuc_codons_flank5}')
            print(f'RTT WT + PBS Flank 5\' Nucleotides Codons: {rtt_wt_inframe_nuc_codons}')
            print(f'RTT WT + PBS Flank 5\' Nucleotides 3\' of Codons: {rtt_wt_inframe_nuc_codons_flank3}')
            print(f'RTT WT + PBS Flank 5\' Nucleotides In-Frame: {rtt_wt_inframe_nuc}')
            print(f'RTT WT + PBS Flank 5\' Amino Acids In-Frame: {rtt_wt_inframe_prot}')
            print(f'RTT WT + PBS Flank 5\' Amino Acid #s In-Frame: {rtt_wt_inframe_prot_indexes}\n')

            # Determine edit by comparing RTT to WT RTT
            prot_befores.append(''.join([str(aa) for aa in rtt_wt_inframe_prot]))
            prot_afters.append(''.join([str(aa) for aa in rtt_pbs_inframe_flank5_prot]))
            print(f'Expected edit: {pegRNAs.iloc[j]["Edit"]}')
            edit = compare_RTTs(rtt_prot=rtt_pbs_inframe_flank5_prot,
                                rtt_prot_indexes=rtt_pbs_inframe_flank5_prot_indexes,
                                rtt_wt_prot=rtt_wt_inframe_prot,
                                rtt_wt_prot_indexes=rtt_wt_inframe_prot_indexes,
                                annotation=pegRNAs.iloc[j]['Annotation'],
                                strand=pegRNAs.iloc[j]['Strand'])
            edits.append(edit)
            
            # Compare expected edit to actual edit
            if edit==pegRNAs.iloc[j]['Edit']:
                edit_matches.append(True)
                edit_matches_explain.append(None) # No need to explain matching edits
            else: 
                edit_matches.append(False)
                if edit is not None:
                    position = [int(num) for num in re.findall(r'\d+', edit)][0] # Find actual edit position
                    if (position==aa_index-1)&(pegRNAs.iloc[j]['Annotation']=='deletion'): edit_matches_explain.append('Common pegRNAs_tester() error: AA deletion outside of target sequence boundary results in wrong AA #')
                    elif pegRNAs.iloc[j]['Annotation']=='deletion': edit_matches_explain.append('Investigate deletion: common error for pegRNAs_tester() is single deletion of repeated AA results in wrong AA #')
                    else: edit_matches_explain.append('Investigate: unknown error')
                elif pegRNAs.iloc[j]['Annotation']=='deletion': edit_matches_explain.append('Common pegRNAs_tester() error: N-terminus AA deletion results in wildtype') 
                else: edit_matches_explain.append('Investigate: unknown error')

        else: 
            print('Error: Strand column can only have "+" and "-".')
            return
    
    pegRNAs['Edit_check']=edits
    pegRNAs['Edit_matches']=edit_matches
    pegRNAs['AAs_before']=prot_befores
    pegRNAs['AAs_after']=prot_afters
    pegRNAs['Edit_matches_explain']=edit_matches_explain
    return pegRNAs