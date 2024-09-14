### clone.py ###
# Author: Marc Zepeda
# Date: 2024-05-29

# Import packages
import pandas as pd
from Bio.Seq import Seq

# Supporting Methods
''' ord_form: Sigma Alrich ordering formatter
        df: Dataframe
        id: id column name
        spacer: spacer column name
        suf: suffix for oligonucleotide category
    Dependencies: pandas
'''
def ord_form(df:pd.DataFrame,id:str,seq:str,suf:str,pre:str):
    ord = df[[(pre+id+suf),(pre+seq+suf)]]
    ord = ord.rename(columns={(pre+id+suf):'Oligo Name',(pre+seq+suf):'Sequence'})
    scale = []
    bp = []
    for s in ord['Sequence']:
        if len(s)<60: scale.append(0.025)
        else: scale.append(0.05)
        bp.append(len(s))
    ord['Scale (µmol)']=scale
    ord['bp']=bp
    return ord

''' tb: Designs top & bottom oligonucleotides
        df: Dataframe with sequences
        id: id column
        seq: sequence column
        t5: top oligonucleotide 5' overhang
        t3: top oligonucleotide 3' overhang
        b5: bottom oligonucleotide 5' overhang
        b3: bottom oligonucleotide 3' overhang
        tG: add 5' G to spacer if needed
        pre: prefix for ids and id column
    Dependencies: pandas, Bio.Seq
'''
def tb(df:pd.DataFrame,id:str,seq:str,t5:str,t3:str,b5:str,b3:str,tG:bool,pre:str):
    top_ids=[]
    bot_ids=[]
    top_seqs=[]
    bot_seqs=[]
    for i,s in enumerate(df[seq]):
        top_ids.append(pre+str(df.iloc[i][id])+'_t')
        bot_ids.append(pre+str(df.iloc[i][id])+'_b')
        if (tG==True)&(s[0]!='G'): top_seqs.append(t5+'G'+s+t3)
        else: top_seqs.append(t5+s+t3)
        bot_seqs.append(b5+str(Seq(s).reverse_complement()+b3))
    df[pre+id+'_t']=top_ids
    df[pre+id+'_b']=bot_ids
    df[pre+seq+'_t']=top_seqs
    df[pre+seq+'_b']=bot_seqs
    
    return df

# GG cloning
''' sgRNAs: Design GG cloning oligonucleotides for cutting and base editing sgRNAs
        df: Dataframe with sequence information for epegRNAs
        id: id column
        spacer: spacer column (Default: Spacer_sequence)
        t5: top oligonucleotide 5' overhang
        t3: top oligonucleotide 3' overhang
        b5: bottom oligonucleotide 5' overhang (revcom)
        b3: bottom oligonucleotide 3' overhang (revcom)
        tG: add 5' G to spacer if needed (Default: True)
        order: order format
    Dependencies: pandas, top_bot(), ord_form()
'''
def sgRNAs(df:pd.DataFrame,id:str,spacer='Spacer_sequence',
          t5='CACC',t3='',b5='AAAC',b3='',
          tG=True,order=True):
    df=tb(df=df,id=id,seq=spacer,t5=t5,t3=t3,b5=b5,b3=b3,tG=tG,pre='o') # Make top and bottom oligos for spacer inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='o'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='o')]).reset_index(drop=True)
    else: return df # Original dataframe with top and bottom oligos

''' epegRNAs: Design GG cloning oligonucleotides for prime editing epegRNAs
        df: Dataframe with sequence information for epegRNAs
        id: id column
        tG: add 5' G to spacer if needed (Default: True)
        order: order format (Default: True)
        make_extension: concatenate RTT, PBS, and linker to make extension sequence (Default: True)
        }
        spacer: epegRNA spacer column (Default: Spacer_sequence)
        extension: epegRNA extension (Default: Extension_sequence)
            RTT: epegRNA reverse transcripase template column (Default: RTT_sequence)
            PBS: epegRNA primer binding site column (Default: PBS_sequence)
            linker: epegRNA linker column (Default: Linker_sequence)
        }
    Assumptions:
        epegRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
        epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA)
    Dependencies: pandas, top_bot(), ord_form()
'''
def epegRNAs(df: pd.DataFrame,id: str,tG=True, order=True,make_extension=True,
             spacer='Spacer_sequence',spacer_t5='CACC',spacer_t3='GTTTAAGAGC',spacer_b5='',spacer_b3='',
             extension='Extension_sequence',RTT='RTT_sequence',PBS='PBS_sequence',linker='Linker_sequence',extension_t5='',extension_t3='',extension_b5='CGCG',extension_b3='GCACCGACTC'):
    if make_extension==True: df[extension] = df[RTT]+df[PBS]+df[linker] # Make extension by concatenating RTT, PBS, and linker
    else: print(f'Warning: Did not make extension sequence!\nMake sure "{extension}" column includes RTT+PBS+linker for epegRNAs.')
    df=tb(df=df,id=id,seq=spacer,t5=spacer_t5,t3=spacer_t3,b5=spacer_b5,b3=spacer_b3,tG=tG,pre='ps_') # Make top and bottom oligos for spacer inserts
    df=tb(df=df,id=id,seq=extension,t5=extension_t5,t3=extension_t3,b5=extension_b5,b3=extension_b3,tG=False,pre='pe_') # Make top and bottom oligos for extension inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='ps_'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='ps_'),
                                      ord_form(df=df,id=id,seq=extension,suf='_t',pre='pe_'),
                                      ord_form(df=df,id=id,seq=extension,suf='_b',pre='pe_')]).reset_index(drop=True)
    else: return df # Original dataframe with top and bottom oligos

''' ngRNAs: Design GG cloning oligonucleotides for prime editing ngRNAs
        df: Dataframe with spacers
        id: id column
        tG: add 5' G to spacer if needed
        order: order format
        }
        spacer: ngRNA spacer column
            t5: top oligonucleotide 5' overhang
            t3: top oligonucleotide 3' overhang
            b5: bottom oligonucleotide 5' overhang
            b3: bottom oligonucleotide 3' overhang}
        }
    Assumptions:
        ngRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    Dependencies: pandas, top_bot(), ord_form()
'''
def ngRNAs(df: pd.DataFrame,id: str,tG=True, order=True,
           spacer='Spacer_sequence',ngRNA_sp_t5='CACC',ngRNA_sp_t3='GTTTAAGAGC',ngRNA_sp_b5='',ngRNA_sp_b3=''):
    df=tb(df=df,id=id,seq=spacer,t5=ngRNA_sp_t5,t3=ngRNA_sp_t3,b5=ngRNA_sp_b5,b3=ngRNA_sp_b3,tG=tG,pre='ns_') # Make top and bottom oligos for spacer inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='ns_'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='ns_')]).reset_index(drop=True)
    
    else: return df # Original dataframe with top and bottom oligos