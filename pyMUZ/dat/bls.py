### bls.py ###
# Author: Marc Zepeda
# Date: 2024-08-20

# Import packages
import pandas as pd
import json
import requests

# Series ID methods
''' series_ids: Returns dataframe containing series ids and corresponding metadata
        ls: list from dc_to_ls()
        dropdown_ids: column names
        sep: seperator from sep (optional, default: '.')
        exclude: trailing value from dc_to_ls (optional, default: '.None')
        pre: series id prefix (optional)
    Dependencies: pandas
'''
def series_ids(ls: list(), dropdown_ids: list(), sep='.',exclude='.None', pre=''):
    
    # Clean list by removing trailing value (exclude) and truncated series
    ls = [l.split(exclude)[0] for l in ls]
    print(f'Removed truncated series: {[l for l in ls if len(l.split(sep))!=len(dropdown_ids)]}')
    ls = [l for l in ls if len(l.split(sep))==len(dropdown_ids)] 
    
    # Initialize ids and metas lists
    ids = []
    metas = []
    for l in ls:
        meta = []
        id = pre # Adds prefix if provided
        for s in l.split(sep):
            id += s.split(' ')[0] # Combine ids into a single string
            meta.append(' '.join(s.split(' ')[1:])) # Retain metadata as a list of strings
        ids.append(id)
        metas.append(meta)

    # Generate dataframe using list of dicts from metas and set index using ids
    return pd.DataFrame([dict(zip(dropdown_ids, row)) for row in metas],index=pd.Index(ids,name='series_id')).reset_index()

''' cols_with_subs: Returns list of dataframe columns with substring
        df: Dataframe
        sub: Substring
    Dependencies: pandas
'''
def cols_with_subs(df: pd.DataFrame(), sub: str()):
    return [col for col in df.columns if sub in col]
    
''' series_options: Return dataframe with all categories, options, order, & codes for a set of series IDs
        dc: Dictionary containing dataframes with categories and codes.
        option_col_subs: Substrings to look for options column
        code_col_subs: Substrings to look for code column
    Dependencies: pandas,cols_with_subs()
'''
def series_options(dc: dict(), option_col_subs=['text','name'],code_col_subs=['code']):

    # Get categories as well as their corresponding orders, options, & codes
    orders = [str(k.split('_')[0]) for k in dc.keys()]
    categories = ['_'.join(k.split('_')[1:]) for k in dc.keys()]
    option_cols = [] 
    for v in dc.values():
        for s in option_col_subs: 
            if len(cols_with_subs(v,s))>0: 
                option_cols.extend(cols_with_subs(v,s))
                break
    code_cols = [] 
    for v in dc.values():
        for s in code_col_subs: 
            if len(cols_with_subs(v,s))>0: 
                code_cols.extend(cols_with_subs(v,s))
                break
    
    # Return dataframe with categories, options, order, & codes 
    df = pd.DataFrame() 
    for (cat,order,option_col,code_col) in zip(categories,orders,option_cols,code_cols):
        options = list(dc['_'.join([order,cat])][option_col])
        codes = list(dc['_'.join([order,cat])][code_col])
        df = pd.concat([df,
                        pd.DataFrame({'Category':[cat]*len(options),
                                      'Option':options,
                                      'Order':[order]*len(options),
                                      'Code':codes})]).reset_index(drop=True)
    return df

# API methods
''' api_v1: Returns dataframe containing series_ids data using BLS API V1
        series_ids: list of series_ids
        start_year: start year for series data
        end_year: end year for series data
    Dependencies: requests,json,pandas
'''
def api_v1(series_ids: list(),start_year: str(),end_year: str()):
    
    # Obtain JSON files represented by series_ids from BLS API V1
    json_data = json.loads(requests.post('https://api.bls.gov/publicAPI/v1/timeseries/data/', 
                           data=json.dumps({"seriesid": series_ids,"startyear": start_year, "endyear": end_year}), 
                           headers={'Content-type': 'application/json'}).text)
    print(f'JSON Header:\nStatus: {json_data["status"]}\nResponse Time: {json_data["responseTime"]}\nMessage: {json_data["message"]}')

    if 'series' in json_data['Results']: # Normalize the 'series' list inside 'Results'
        normalized_data = pd.json_normalize(json_data['Results']['series'], 
                                            record_path=['data'],  
                                            meta=['seriesID'],
                                            errors='ignore')
        return normalized_data
    else:
        print("No 'series' data found in 'Results'")
        return pd.DataFrame()  # Return an empty DataFrame if no data is found

''' api_v2: Returns dataframe containing series_ids data using BLS API V2
        series_ids: list of series_ids
        start_year: start year for series data
        end_year: end year for series data
        catalog: true|false
        calculations: true|false
        annual_average: true|false
        aspects: true|false
        registration_key: API registration key
    Dependencies: requests,json,pandas
'''
def api_v2(series_ids: list(),start_year,end_year,catalog=True,calculations=True,annual_average=True,aspects=True,registration_key='b623916dd99845bc8f430711d72c9f38'):
    
    # Obtain JSON files represented by series_ids from BLS API V2
    json_data = json.loads(requests.post('https://api.bls.gov/publicAPI/v2/timeseries/data/', 
                           data=json.dumps({"seriesid":series_ids,"startyear":str(start_year),"endyear":str(end_year),
                                            "catalog":catalog,"calculations":calculations, "annualaverage":annual_average,
                                            "aspects":aspects,"registrationkey":registration_key}), 
                           headers={'Content-type': 'application/json'}).text)
    print(f'JSON Header:\nStatus: {json_data["status"]}\nResponse Time: {json_data["responseTime"]}\nMessage: {json_data["message"]}')

    if 'series' in json_data['Results']: # Normalize the 'series' list inside 'Results'
        normalized_data = pd.json_normalize(json_data['Results']['series'], 
                                            record_path=['data'], 
                                            meta=[
                                                'seriesID', 
                                                ['catalog', 'series_title'],
                                                ['catalog', 'seasonality'],
                                                ['catalog', 'survey_name'],
                                                ['catalog', 'survey_abbreviation'],
                                                ['catalog', 'measure_data_type'],
                                                ['catalog', 'commerce_industry'],
                                                ['catalog', 'commerce_sector'],
                                                ['catalog', 'area']
                                            ], 
                                            errors='ignore')
        return normalized_data
    else:
        print("No 'series' data found in 'Results'")
        return pd.DataFrame()  # Return an empty DataFrame if no data is found