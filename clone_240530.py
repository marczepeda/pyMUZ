### clone.py ###
# Author: Marc Zepeda
# Date: 2024-05-29

# Import methods
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
def ord_form(df:pd.DataFrame(),id:str(),seq:str(),suf:str(),pre:str()):
    ord = df[[(pre+id+suf),(pre+seq+suf)]]
    ord = ord.rename(columns={(pre+id+suf):'Oligo Name',(pre+seq+suf):'Sequence'})
    scale = []
    for s in ord['Sequence']:
        if len(s)<60: scale.append(0.025)
        else: scale.append(0.05)
    ord['Scale (µmol)']=scale
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
def tb(df:pd.DataFrame(),id:str(),seq:str(),t5:str(),t3:str(),b5:str(),b3:str(),tG:bool(),pre:str()):
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
''' be: Design GG cloning oligonucleotides for base editing
        df: Dataframe with spacers
        id: id column
        spacer: spacer column
        t5: top oligonucleotide 5' overhang
        t3: top oligonucleotide 3' overhang
        b5: bottom oligonucleotide 5' overhang
        b3: bottom oligonucleotide 3' overhang
        tG: add 5' G to spacer if needed
        order: order format
    Dependencies: pandas, top_bot(), ord_form()
'''
def be(df:pd.DataFrame(),id='id',spacer='spacer',
       t5='CACC',t3='',b5='AAAC',b3='',
       tG=True,order=True):
    df=tb(df=df,id=id,seq=spacer,t5=t5,t3=t3,b5=b5,b3=b3,tG=tG,pre='o')
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='o'),
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='o')])
    else: return df

''' pe: Design GG cloning oligonucleotides for prime editing
        df: Dataframe with spacers
        id: id column
        tG: add 5' G to spacer if needed
        order: order format
        }
        epeg_sp: epegRNA spacer column
        epeg_RTT: epegRNA reverse transcripase template column
        epeg_PBS: epegRNA primer binding site column
        epeg_link: epegRNA linker column
        (epeg_ex: epegRNA extension)
        ng_sp: ngRNA spacer column
            t5: top oligonucleotide 5' overhang
            t3: top oligonucleotide 3' overhang
            b5: bottom oligonucleotide 5' overhang
            b3: bottom oligonucleotide 3' overhang}
        }
    Assumptions:
        epeg/ngRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
        epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA)
    Dependencies: pandas, top_bot(), ord_form()
'''
def pe(df: pd.DataFrame(),id='id',tG=True, order=True,
       epeg_sp='epegRNA_spacer',epeg_sp_t5='CACC',epeg_sp_t3='GTTTAAGAGC',epeg_sp_b5='',epeg_sp_b3='',
       epeg_ex='epegRNA_extension',epeg_ex_t5='',epeg_ex_t3='',epeg_ex_b5='CGCG',epeg_ex_b3='GCACCGACTC',
       ng_sp='ngRNA_spacer',ngRNA_sp_t5='CACC',ngRNA_sp_t3='GTTTAAGAGC',ngRNA_sp_b5='',ngRNA_sp_b3=''):
    df=tb(df=df,id=id,seq=epeg_sp,t5=epeg_sp_t5,t3=epeg_sp_t3,b5=epeg_sp_b5,b3=epeg_sp_b3,tG=tG,pre='ps_')
    df=tb(df=df,id=id,seq=epeg_ex,t5=epeg_ex_t5,t3=epeg_ex_t3,b5=epeg_ex_b5,b3=epeg_ex_b3,tG=False,pre='pe_')
    df=tb(df=df,id=id,seq=ng_sp,t5=ngRNA_sp_t5,t3=ngRNA_sp_t3,b5=ngRNA_sp_b5,b3=ngRNA_sp_b3,tG=tG,pre='ns_')
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=epeg_sp,suf='_t',pre='ps_'),
                                      ord_form(df=df,id=id,seq=epeg_sp,suf='_b',pre='ps_'),
                                      ord_form(df=df,id=id,seq=epeg_ex,suf='_t',pre='pe_'),
                                      ord_form(df=df,id=id,seq=epeg_ex,suf='_b',pre='pe_'),
                                      ord_form(df=df,id=id,seq=ng_sp,suf='_t',pre='ns_'),
                                      ord_form(df=df,id=id,seq=ng_sp,suf='_b',pre='ns_')])
    else: return df