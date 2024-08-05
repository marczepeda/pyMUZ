### tidy.py ###
# Author: Marc Zepeda
# Date: 2024-08-03

# Import methods
import pandas as pd
import re

# Methods for dictionary containing dataframes.
''' split_by: Splits elements of list, set, or series by specified seperator
        series: list, set or series
        by: seperator
'''
def split_by(series, by=', '):
    split_elements = []
    for element in series: 
        if isinstance(element, str): split_elements.extend(element.split(by))
    return split_elements

''' isolate: Isolate rows in dataframes based specified value(s)
        dc: dictionary
        col: df column name
        get: value, set, list, dictionary of dataframes
        want: do you want the value
    Dependencies: re, pandas, split_by()
'''
def isolate(dc: dict(), col: str(), get, get_col='', get_col_split_by='', want=True, exact=True):
    if want==True: 
        if get is None: return {key:df[df[col].isnull()==True].reset_index(drop=True) for key,df in dc.items()}
        elif type(get)==set or type(get)==list or type(get)==pd.Series: 
            if exact==True: 
                if get_col_split_by=='': return {key:df[df[col].isin(set(get))==True].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].isin(set(split_by(get,by=get_col_split_by)))==True].reset_index(drop=True) for key,df in dc.items()}
            else: 
                if get_col_split_by=='': return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in get))==True].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in split_by(get,by=get_col_split_by)))==True].reset_index(drop=True) for key,df in dc.items()}
        elif type(get)==dict: 
            if exact==True: 
                if get_col_split_by=='': return {key:df[df[col].isin(set(get[key][get_col]))==True].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].isin(set(split_by(get[key][get_col],by=get_col_split_by)))==True].reset_index(drop=True) for key,df in dc.items()}
            else: 
                if get_col_split_by=='': return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in get[key][get_col]),case=False, na=False)==True].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in split_by(get[key][get_col],by=get_col_split_by)),case=False, na=False)==True].reset_index(drop=True) for key,df in dc.items()}
        else: return {key:df[df[col]==get].reset_index(drop=True) for key,df in dc.items()}
    else: 
        if get is None: return {key:df[df[col].isnull()==False].reset_index(drop=True) for key,df in dc.items()}
        elif type(get)==set or type(get)==list or type(get)==pd.Series: 
            if exact==True: 
                if get_col_split_by=='': return {key:df[df[col].isin(set(get))==False].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].isin(set(split_by(get,by=get_col_split_by)))==False].reset_index(drop=True) for key,df in dc.items()}
            else: 
                if get_col_split_by=='': return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in get))==False].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in split_by(get,by=get_col_split_by)))==False].reset_index(drop=True) for key,df in dc.items()}
        elif type(get)==dict: 
            if exact==True: 
                if get_col_split_by=='': return {key:df[df[col].isin(set(get[key][get_col]))==True].reset_index(drop=False) for key,df in dc.items()}
                else: return {key:df[df[col].isin(set(split_by(get[key][get_col],by=get_col_split_by)))==True].reset_index(drop=False) for key,df in dc.items()}
            else: 
                if get_col_split_by=='': return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in get[key][get_col]),case=False, na=False)==False].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in split_by(get[key][get_col],by=get_col_split_by)),case=False, na=False)==False].reset_index(drop=True) for key,df in dc.items()}
        else: return {key:df[df[col]!=get].reset_index(drop=True) for key,df in dc.items()}

''' modify: Returns dictionary containing dataframes new or updated column with specified value(s) or function
        dc: dictionary
        col: new/old column name
        val: column value, list, or function.
            e.g., new_val=lambda df: df['AA Mutation'].split('.')[1]
        axis: function is applied to column (1) or row (0)
    Dependencies: pandas
'''
def modify(dc: dict(), col: str(), val, axis=1, **kwargs):
    dc2=dict()
    for key,df in dc.items():
        if callable(val): df2 = df.assign(**{col: df.apply(val, axis=axis, **kwargs)})
        else: df2 = df.assign(**{col: val})
        dc2[key]=df2
    return dc2

''' melt: Returns dictionary containing tidy dataframes
        dc: dictionary of dataframes
        id_vars: metadata columns
    Dependencies: pandas
'''
def melt(dc: dict(),id_vars,**kwargs):
    dc2=dict()
    for key,df in dc.items(): dc2[key]=pd.melt(frame=df,id_vars=id_vars,**kwargs)
    return dc2

''' join: Returns a single dataframe from a dictionary of dataframes
        dc: dictionary of dataframes
        col: name for keys column
    Dependencies: pandas
'''
def join(dc: dict(), col='key'):
    df = pd.DataFrame()
    for key,val in dc.items():
        val[col]=key
        df=pd.concat([df,val])
    return df

''' split: Returns from a dictionary of dataframes from a single dataframe
        df: dataframe
        key: column for spliting dataframe
    Dependencies: pandas
'''
def split(df: pd.DataFrame(), key: str()):
    return {k:df[df[key]==k] for k in list(df[key].value_counts().keys())} 

''' merge: Adds metadata columns to data dataframe using metadata dataframe
        data: data dataframe
        meta: metadata dataframe
        id: id(s) column name(s) [str: both, list: data & meta]
        cols: list of column names in metadata dataframe
'''
def merge(data: pd.DataFrame(), meta: pd.DataFrame(), id, cols):
    if type(id)==str:
        for c in cols: 
            id_c = dict(zip(meta[id],meta[c]))
            data[c] = [id_c[i] for i in data[id]]
    elif (type(id)==list)&(len(id)==2):
        for c in cols: 
            id_c = dict(zip(meta[id[1]],meta[c]))
            data[c] = [id_c[i] for i in data[id[0]]]
    else: print("Error: id needs to be string or list of 2 strings")
    return data