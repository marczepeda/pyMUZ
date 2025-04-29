### get_JSON.py ###
# Author: Marc Zepeda
# Date: 2024-01-03

# API extraction packages
from urllib import request
from time import sleep
import json
import pandas as pd
from ..pyMUZ.gen import io

# General packages
import os
from collections import OrderedDict, Counter, defaultdict

### General API ###
# Supporting Methods

def get_json(url: str=None, pt: str=None, file:str=None, dir:str=None):
    ''' 
    get_json(): obtains json file from the internet
    
    Parameters:
    url (str, option 1): website url for json file
    pt (str, option 2): path for json file
    file (str, optional): saved file name
    dir (str, optional): saved file directory
    
    Dependencies: os & requests
    '''
    if url is not None: # url specified
        try: 
            # Send a GET request to the URL
            req = request.Request(url)
            response = request.urlopen(req)
            encoded_response = response.read()
            decoded_response = encoded_response.decode()
            payload = json.loads(decoded_response)

            # Save html file if file and dir are specified
            if file is not None and dir is not None:
                io.save(dir=dir,file=file,obj=pd.DataFrame({'JSON': [payload]}))
                print(f"JSON content successfully saved to {file}")

            return payload
        
        except Exception as e: print(f"{url} error:\n{e}")
    
    elif pt is not None: # pt specified
        try: 
            # Load the JSON data from the file
            with open(pt, 'r') as f: 
                payload = json.load(f)
            
            # Save html file if file and dir are specified
            if file is not None and dir is not None:
                io.save(dir=dir,file=file,obj=pd.DataFrame({'JSON': [payload]}))
                print(f"JSON content successfully saved to {file}")

            return payload

        except Exception as e: print(f"{pt} error:\n{e}")
    
    else: KeyError('Specify url or pt')


## work here...
''' get_results: Returns Pandas dataframe containing "get" values from the "results" section of the JSON file.
        - payload: Python dictionary containing data from JSON file
        - get: list of data
'''
def get_results(payload: dict(), get: list()):
    
    # Fill output data frame with "get" values within "results"
    df = pd.DataFrame(columns=get)
    i=0
    
    for result in payload["results"]:
        
        # Put "get" values into temporary Python dictionary
        values = [result[get[j]] for j in range(len(get))]
        item = dict(zip(get, values))
        
        # Adds "get" values to output dataframe
        obj = pd.DataFrame(item, index=[i])
        df = pd.concat([df, obj], ignore_index=True)
        i+=1
    
    return df

''' set_default_values: Fill in missing keys with default values for list of Python dictionaries or Python dictionary from JSON file
        d: Python dictionary
        default_keys: list of missing keys
        default_values: list of default values
'''
def set_default_values(d: dict, default_keys: list, default_values: list):
    defaults = dict(zip(default_keys,default_values))
    return [dict(defaults, **item) for item in d]

# Main Method

''' get_json_normalize_record: Returns a dictionary of dataframes using pd.json_normalize and specified record keys.
        d: List of Python dictionaries from JSON file
        record_keys: list of record keys
'''
def get_json_normalize_record(d: list, record_keys: list):
    
    # Fill in missing keys with default values
    d_cleaned = set_default_values(d=d, default_keys=record_keys, default_values=[[] for _ in record_keys])

    # Get record meta 
    keys = set()
    for item in d:
        for k in list(item.keys()):
             keys.add((k,type(item[k])))
    m = list()
    for key in keys:
        if key[1]!=type(list()) and key[1]!=type(dict()):
            m.append(key[0])
    
    # Iterate through record_keys with pd.json_normalize
    dfs = dict()
    dfs['no_record_path'] = pd.json_normalize(data=d_cleaned)
    for rc_k in record_keys:
        dfs[rc_k] = pd.json_normalize(data=d_cleaned,
                                      record_path=rc_k,
                                      record_prefix=rc_k+".",
                                      meta=m)
    return dfs

# Not completed nor working on...
# Supporting Methods
''' get_structure: Recursive method that returns a set for a Python dictionary's structure with ":" & ":[]:" seperators indicating levels.
        obj: dict or list objects from JSON file payload
        ls: list containing Python dictionary structure
        k: string that retains previous keys within the Python dictionary structure
'''
def get_structure(obj, s=set, k=''):

    # Checks for instance of dictionary within the obj, appends keys to a set, & recursive search with values
    if isinstance(obj, dict):
        for key, value in obj.items():
            if k=='': # First dictionary level (don't want ":" seperator)
                s.add(key)
                get_structure(value, s, key)
            else: # Subsequent dictionary levels (do want ":" seperator)
                s.add(k+":"+key)
                get_structure(value, s, k+":"+key)
    
    # Checks for instance of list containing dictionary within the obj using recursive search
    if isinstance(obj, list):
        for ls in obj:
            get_structure(ls, s, k+':[]') # List dictionary levels (want ":[]:" seperator)

    # Returns final set with colon seperator indicating dictionary levels.
    return s

''' remove_elements_containing_substring: Remove an element if its string is a substring in another element of the set
        input_set: Python set
        substring: substring
    Credit: ChatGPT 3.5
'''
def remove_elements_containing_substring(input_set, substring):
    elements_to_remove = {element for element in input_set if substring in element}
    input_set -= elements_to_remove

''' remove_substring_elements: Removes elements from a set that contain a specified substring
        input_set: Python set
    Credit: ChatGPT 3.5
'''
def remove_substring_elements(input_set):
    elements_to_remove = set()

    for element in input_set:
        for other_element in input_set:
            if element != other_element and element in other_element:
                elements_to_remove.add(element)

    input_set -= elements_to_remove

''' get_key_sets: Returns sets containing keys with lists and all keys from a list describing Python dictionary's structure
        structure_ls: list containing endpoints for Python dictionary with lists structure minus endpoints containing "exclude" values
'''
def get_key_sets(structure_ls: str):

    # Create sets containing Python dictionary elements that are keys with lists and all keys
    k_ls_set = set()
    k_set = set()

    for element in structure_ls:
        
        # Split element based on ":[]:" seperator and save keys with lists to a set
        parts = set(element.split(":[]:")[:-1])
        for p in parts: # Remove extra information that not part of the list key
            if p.find(":")!=-1:
                parts.discard(p)
                parts.add(p[p.find(":")+1])
        k_ls_set = k_ls_set.union(parts)

        # Split element based on ":" seperator and save all keys to a set
        keys = element.split(":")
        for k in keys:
            if k!='[]':
                k_set.add(k)
    
    return k_ls_set,k_set

''' get_structure_dict: Returns a python dictionary representing the original Python dictionary's structure
        structure_ls: list containing endpoints for Python dictionary with lists structure minus endpoints containing "exclude" values
'''
def get_structure_dict(structure_ls: str): 
    
    # Obtain list of dictionaries based on the Python dictionary's structure
    structure_dict = dict()
    for element in structure_ls:

        # Split element based on ":" seperator (excluding "[]") into a list
        parts_ls = element.split(":")
        parts_ls = [item for item in parts_ls if item != "[]"]

        # Obtain lists containing Python dictionary structure
        temp_ls = list()
        j = 0 # temp_ls counter
        i = len(parts_ls)-1 # parts_ls counter
        
        temp_ls.extend([{parts_ls[i-1]:parts_ls[i]}]) # Start condition
        while i > 1: # Extend condition
            temp_ls.extend([{parts_ls[i-2]:temp_ls[j]}])
            i-=1
            j+=1
        structure_dict[element]=temp_ls[j] # End condition

    return structure_dict

# Main Methods

''' Not completed and stopped working on....
    extract_data: Convert data from JSON file payload (Python dictionary with lists) to dataframe
        payload: JSON file payload (Python dictionary with lists)
        exclude: list of strings
    Dependencies: get_structure, remove_elements_containing_substring, remove_substring_elements, get_key_sets, get_structure_dict
'''
def extract_data(payload: dict, exclude: list): # Not done. Will work on tomorrow.
    
    # Get list containing endpoints for Python dictionary with lists structure minus endpoints containing "exclude" values
    structure_set = get_structure(payload)
    for ex in exclude:
        remove_elements_containing_substring(structure_set, ex)
    remove_substring_elements(structure_set)
    structure_ls = sorted(list(structure_set))
    structure_dict = get_structure_dict(structure_ls)

    # Make sets containing Python dictionary elements that are keys with lists and all keys.
    k_ls_set,k_set = get_key_sets(structure_ls)

    # Fill output data frame with "get" values within "results"
    df = pd.DataFrame(columns=structure_ls)
    i=0

    return df


### NCI SEER API ###
# Supporting Methods

''' get_NCI_SEER_API_payload: Returns Python dictionary from the JSON file encoded in the url.
        - url_prefix: prefix for database
        - url_suffix: suffix for specific JSON file
        - key: NCI SEER API key
'''
def get_NCI_SEER_API_payload(url_prefix: str, url_suffix: str, key: str):
    url = os.path.join(url_prefix,url_suffix)
    req = request.Request(url)
    req.add_header("X-SEERAPI-Key", key)
    response = request.urlopen(req)
    encoded_response = response.read()
    decoded_response = encoded_response.decode()
    payload = json.loads(decoded_response)
    return payload

''' get_NCI_SEER_API_link: Returns url for the next JSON file.
        - url_suffix: suffix for previous JSON file
        - offset: offset value
        - offset_loc: list index for offset based on sep 
        - sep: link seperator value
'''
def get_NCI_SEER_API_link(url_suffix: str, offset: int, offset_loc=3, sep="&"):
    url_suffix_ls = url_suffix.split(sep) # Split url into list
    url_suffix_ls[offset_loc] = url_suffix_ls[offset_loc][0:url_suffix_ls[offset_loc].find('=')+1] + str(offset) # Modify offset
    return sep.join(url_suffix_ls) # Rejoin new url from updated list

''' get_NCI_SEER_API_link_last: Returns url for the last JSON file.
        - url_suffix: suffix for previous JSON file
        - total: total value
        - offset: offset value
        - offset_loc: list index for offset based on sep
        - count_loc: list index for count based on sep 
        - sep: link seperator value
'''
def get_NCI_SEER_API_link_last(url_suffix: str, total: int, offset: int, offset_loc=3, count_loc=2, sep="&"):
    url_suffix_ls = url_suffix.split(sep) # Split url into list
    url_suffix_ls[offset_loc] = url_suffix_ls[offset_loc][0:url_suffix_ls[offset_loc].find('=')+1] + str(offset) # Modify offset
    url_suffix_ls[count_loc] = url_suffix_ls[count_loc][0:url_suffix_ls[count_loc].find('=')+1] + str(total-offset) # Modify count
    return sep.join(url_suffix_ls) # Rejoin new url from updated list

''' get_NCI_SEER_API_results: Returns Pandas dataframe containing "get" values from the "results" section of the JSON file.
        - payload: Python dictionary containing data from JSON file
        - get: list of data
'''
def get_NCI_SEER_API_results(payload: dict, get: list):
    
    # Fill output data frame with "get" values within "results"
    df = pd.DataFrame(columns=get)
    i=0
    
    for result in payload["results"]:
        
        # Put "get" values into temporary Python dictionary
        values = [result[get[j]] for j in range(len(get))]
        item = dict(zip(get, values))
        
        # Adds "get" values to output dataframe
        obj = pd.DataFrame(item, index=[i])
        df = pd.concat([df, obj], ignore_index=True)
        i+=1
    
    return df

# Main Methods
''' get_NCI_SEER_API_disease: Obtain specified data encoded in JSON files from NCI SEER API
        - url_prefix: prefix for database
        - url_suffix: suffix for specific JSON file
        - key: NCI SEER API key
        - get: list of data
    Dependencies: get_NCI_SEER_API_payload, get_NCI_SEER_API_link, get_NCI_SEER_API_link_last, get_NCI_SEER_API_results
'''
def get_NCI_SEER_API_disease(url_prefix: str, url_suffix: str, key: str, get: list):
    
    # Makes Python dictionary from the JSON file encoded in the url.
    payload = get_NCI_SEER_API_payload(url_prefix=url_prefix, url_suffix=url_suffix, key=key)

    # Make Pandas dataframe containing "get" values from the "results" section of the JSON file.
    df = get_NCI_SEER_API_results(payload=payload, get=get)

    # Iterate through JSON files
    cnt = 24 # API step size
    last = False
    while last is False:

        # Get information for next link
        offset = payload['offset']
        total = payload['total']
        offset += cnt 

        # Create next link
        if offset + cnt >= total:
            last = True
            url_suffix = get_NCI_SEER_API_link_last(url_suffix=url_suffix, total=total, offset=offset)
        else:
            url_suffix = get_NCI_SEER_API_link(url_suffix=url_suffix, offset=offset)
        
        # Obtain next payload & append to Pandas dataframe
        payload = get_NCI_SEER_API_payload(url_prefix=url_prefix, url_suffix=url_suffix, key=key)
        print(url_suffix)
        df = pd.concat([df, get_NCI_SEER_API_results(payload=payload, get=get)],ignore_index=True)
    
    return df