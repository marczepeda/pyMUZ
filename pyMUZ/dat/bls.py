'''
Module: bls.py
Author: Marc Zepeda
Created: 2024-08-20
Description: Bureau of Labor Statistics

Usage:
[Series ID methods]
- series_ids(): returns dataframe containing series ids and corresponding metadata
- cols_with_subs(): returns list of dataframe columns with substring
- series_options(): return dataframe with all categories, options, order, & codes for a set of series IDs

[API methods]
- api_v1(): returns dataframe containing series_ids data using BLS API V1
- api_v2(): returns dataframe containing series_ids data using BLS API V2
- api_batch(): submits BLS API requests in batches

[Data Wrangling Methods]
- convert_str_output(): converts strings to numerical values within BLS output dataframe
- scat_series_ids(): plots times series from BLS output dataframe
'''

# Import packages
import pandas as pd
import json
import requests
import time
from ..gen import tidy as t
from ..gen import io as io
from ..gen import plot as p

# Series ID methods
def series_ids(ls: list, dropdown_ids: list, sep='.',exclude='.None', pre=''):
    ''' 
    series_ids(): returns dataframe containing series ids and corresponding metadata
    
    ls (list): list from dc_to_ls()
    dropdown_ids (list): column names
    sep (str, optional): seperator from sep (Default: '.')
    exclude (str, optional): trailing value from dc_to_ls (Default: '.None')
    pre (str, optional): series id prefix

    Dependencies: pandas
    '''
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

def cols_with_subs(df: pd.DataFrame, sub: str):
    ''' 
    cols_with_subs(): returns list of dataframe columns with substring
    
    Parameters:
    df (dataframe): dataframe
    sub (str): substring
    
    Dependencies: pandas
    '''
    return [col for col in df.columns if sub in col]

def series_options(dc: dict, option_col_subs=['text','name'],code_col_subs=['code']):
    ''' 
    series_options(): return dataframe with all categories, options, order, & codes for a set of series IDs
    
    Parameters:
    dc (dict): Dictionary containing dataframes with categories and codes.
    option_col_subs (list, optional): Substrings to look for options column
    code_col_subs (list, optional): Substrings to look for code column
    
    Dependencies: pandas & cols_with_subs()
    '''
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
def api_v1(series_ids: list, start_year, end_year):
    ''' 
    api_v1(): returns dataframe containing series_ids data using BLS API V1
    
    Parameters:
    series_ids (list): list of series_ids
    start_year: start year for series data
    end_year: end year for series data
    
    Dependencies: requests, json, & pandas
    '''
    # Obtain JSON files represented by series_ids from BLS API V1
    json_data = json.loads(requests.post('https://api.bls.gov/publicAPI/v1/timeseries/data/', 
                           data=json.dumps({"seriesid": series_ids,"startyear": str(start_year), "endyear": str(end_year)}), 
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

def api_v2(series_ids: list,start_year,end_year,catalog=True,calculations=True,annual_average=True,aspects=True,registration_key='b623916dd99845bc8f430711d72c9f38'):
    '''
    api_v2(): returns dataframe containing series_ids data using BLS API V2

    Parameters:
    series_ids: list of series_ids
    start_year: start year for series data
    end_year: end year for series data
    catalog: true|false
    calculations: true|false
    annual_average: true|false
    aspects: true|false
    registration_key: API registration key

    Dependencies: requests, json, & pandas
    '''
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

def api_batch(series_ids: list, start_year, end_year, api_v='2', **kwargs):
    ''' 
    api_batch(): submits BLS API requests in batches
    
    Parameters:
    series_ids (list): list of series_ids
    start_year: start year for series data (int or str)
    end_year: end year for series data (int or str)
    api_v (str, optional): 1 or 2 (Default: 2)
    
    **kwargs:
    catalog (bool, optional): true|false (api_v='2')
    calculations (bool, optional): true|false (api_v='2')
    annual_average (bool, optional): true|false (api_v='2')
    aspects (bool, optional): true|false (api_v='2')
    registration_key (bool, optional): API registration key (api_v='2')
    
    Dependencies: pandas, api_v1(), & api_v2()
    '''
    # Set API Version Rules & Notify User
    api_v_rules = {
    "2": {
        "Service": "Registered",
        "Daily query limit": 500,
        "Series per query limit": 50,
        "Years per query limit": 20,
        "Request rate limit": "50 requests per 10 seconds",
        "Net/Percent Changes": "Yes",
        "Optional annual averages": "Yes",
        "Series description information (catalog)": "Yes"
    },
    "1": {
        "Service": "Unregistered",
        "Daily query limit": 25,
        "Series per query limit": 25,
        "Years per query limit": 10,
        "Request rate limit": "50 requests per 10 seconds",
        "Net/Percent Changes": "No",
        "Optional annual averages": "No",
        "Series description information (catalog)": "No"
    }}
    print(f'BLS API version: {api_v}\n{pd.Series(api_v_rules[api_v])}\n\nSeries ID Batches:')

    # Generate series_ids batches depending on API version
    series_ids_batches = []
    series_ids_batch = []
    for i,series_id in enumerate(series_ids):
        series_ids_batch.append(series_id)
        if (i+1)%api_v_rules[api_v]['Series per query limit']==0:
            print(series_ids_batch)
            series_ids_batches.append(series_ids_batch)
            series_ids_batch=[]
        elif i+1==len(series_ids):
            print(series_ids_batch) 
            series_ids_batches.append(series_ids_batch)

    # Generate years batches depending on API version
    years_batches = []
    if int(end_year)-int(start_year)+1<=api_v_rules[api_v]['Years per query limit']: # Only 1 years_batch
        years_batches.append((start_year,end_year))
    else: # More than 1 years_batch
        for i,year in enumerate(range(int(start_year),int(end_year)+1)):
            if (i+1)%api_v_rules[api_v]['Years per query limit']==0: 
                years_batches.append((str(int(year)-api_v_rules[api_v]['Years per query limit']+1),str(year)))
            elif int(year)==int(end_year):
                years_batches.append((years_batches[len(years_batches)-1][1],str(year)))
    print(f'\nYears Batches: {years_batches}')

    # Check total # of requests before submitting
    r = len(series_ids_batches)*len(years_batches)
    if r>api_v_rules[api_v]['Daily query limit']: # Too many requests
        print(f'\n{r} is above daily query limit for BLS API version {api_v}.\nDid not submit API requests!')
        return
    
    else: # Iterate through years and series ids batches
        print(f'\nProcessing {r} request(s) from BLS API version {api_v}.')
        dc = dict()
        for i,(start_year,end_year) in enumerate(years_batches):
            for j,series_ids_batch in enumerate(series_ids_batches):
                if api_v=='1': dc[f'{str(i)}_{str(j)}'] = api_v1(series_ids=series_ids_batch,start_year=start_year,end_year=end_year,**kwargs)
                elif api_v=='2': dc[f'{str(i)}_{str(j)}'] = api_v2(series_ids=series_ids_batch,start_year=start_year,end_year=end_year,**kwargs)
                time.sleep(0.2) # 50 requests per 10 seconds limit

        try: return t.join(dc,'batch_years-series_ids') # Returns concatenated dataframe
        except Exception: return dc # Returns dictionary of dataframes

# Data Wrangling Methods
def convert_str_output(data: pd.DataFrame, dir: str=None, file: str=None):
    ''' 
    convert_str_output(): converts strings to numerical values within BLS output dataframe
    
    Parameters:
    data (dataframe): BLS output dataframe
    dir (str, optional): save directory
    file (str, optional): save filename
    
    Depedencies: io & pandas
    '''
    data['year']=[int(year) for year in data['year']]
    data['value']=[float(value) for value in data['value']]
    data['year_period']=[int(year)+(float(period[1:])-1)/13 for year,period in zip(data['year'],data['period'])]
    if dir is not None and file is not None: io.save(dir,file,data)
    return data

def scat_series_ids(data: pd.DataFrame, dir: str=None):
    ''' 
    scat_series_ids(): plots times series from BLS output dataframe
    
    Parameters:
    data (pandas): BLS output dataframe
    dir (str, optional): save directory
    
    Depedencies: plot & pandas
    ''' 
    seriesIDs = list(data['seriesID'].value_counts().keys())
    value_max = max(data['value'])

    for seriesID in seriesIDs:
        data_seriesID = data[data['seriesID']==seriesID]
        value_max = max(data_seriesID['value'])
        graph_max = int((value_max+10)/10)*10

        p.scat(typ='line',df=data_seriesID,x='year_period',y='value',
            y_axis_dims=(0,graph_max),
            title=f"{data_seriesID.iloc[0]['catalog.commerce_industry']}:\n{data_seriesID.iloc[0]['catalog.area']}",
            y_axis=f"{data_seriesID.iloc[0]['catalog.measure_data_type']}",
            x_axis=f"Year",
            figsize=(10,5),
            dir=dir,
            file=f"{data_seriesID.iloc[0]['catalog.commerce_industry']}_{data_seriesID.iloc[0]['catalog.area']}.png".replace('/','-'))