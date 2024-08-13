### io.py ###
# Author: Marc Zepeda
# Date: 2024-05-20

# Import methods
import pandas as pd
import os
import csv as c

# CSV methods
''' get: Returns pandas dataframe from csv file.
        dir: output directory path
'''
def get(pt: str(),**kwargs):
    return pd.read_csv(filepath_or_buffer=pt,**kwargs)

''' get_dir: Returns python dictionary of dataframe from csv files within a directory.
        dir: output directory path
        suf: file type suffix
'''
def get_dir(dir: str(),suf='.csv',**kwargs):
    files = os.listdir(dir)
    csv_files = [file for file in files if file[-4:]==suf]
    csvs = dict()
    for csv_file in csv_files:
        csvs[csv_file[:-len(suf)]] = get(pt=os.path.join(dir,csv_file),**kwargs)
    return csvs

''' save: Save .csv file to a specified output directory from obj
        dir: output directory path
        file: file name
        obj: Dataframe, set, list, series
    Dependencies: pandas, os, csv
'''
def save(dir: str(), file: str(), obj, cols=[], id=False, sort=True, **kwargs):
    if not os.path.exists(dir): # Make output directory if it does not exist
        os.mkdir(dir)

    if type(obj)==pd.DataFrame:
        for col in cols: # Check if each element in the list is a string
            if not isinstance(col, str):
                raise ValueError("All elements in the list must be strings.")
        if cols!=[]: obj = obj[cols]
        obj.to_csv(os.path.join(dir,file), index=id,**kwargs)
    elif type(obj)==set or type(obj)==list or type(obj)==pd.Series:
        if sort==True: obj2 = sorted(list(obj))
        else: obj2=list(obj)
        with open(os.path.join(dir,file), 'w', newline='') as csv_file:
            csv_writer = c.writer(csv_file, dialect='excel') # Create a CSV writer object
            csv_writer.writerow(obj2) # Write each row of the list to the CSV file

''' save_dir: Save .csv files to a specified output directory from dictionary of objs
        dir: output directory path
        file_suffix: file name
        dc: dictionary of objs
    Dependencies: pandas, os, csv, save()
'''
def save_dir(dir: str(), file_suffix: str(), dc, **kwargs):
    for key,val in dc.items(): save(dir=dir,file=key+file_suffix,obj=val,**kwargs)