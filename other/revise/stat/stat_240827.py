### stat.py ###
# Author: Marc Zepeda
# Date: 2024-08-26

import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import skew, kurtosis

''' describe: Calculate descriptive statistics for each numerical column in the DataFrame.df: DataFrame
        cols: list of numerical columns to compute statistics (optional)
        id_vars: id column
    Dependencies: pandas, numpy
'''
def describe(df: pd.DataFrame, cols=[], group=''):

    if group!='': df = df.pivot(columns=group) # Splits tidy dataframe
    if len(cols)>0: df = df[cols] # Isolate specified columns

    descriptive = pd.DataFrame()    
    descriptive['mean'] = df.mean() # Mean
    descriptive['median'] = df.median() # Median
    descriptive['variance'] = df.var() # Variance
    descriptive['std_dev'] = df.std() # Standard Deviation
    descriptive['mad'] = df.apply(lambda x: np.median(np.abs(x - x.median()))) # Median Absolute Deviation
    descriptive['min'] = df.min() # Minimum
    descriptive['max'] = df.max() # Maximum
    descriptive['range'] = df.max() - df.min() # Range
    descriptive['skewness'] = df.apply(lambda x: skew(x, nan_policy='omit')) # Skewness
    descriptive['kurtosis'] = df.apply(lambda x: kurtosis(x, nan_policy='omit')) # Kurtosis
    descriptive['count'] = df.count() # Count (non-missing observations)
    descriptive['sum'] = df.sum() # Sum
    descriptive['25%'] = df.quantile(0.25)  # Quantiles (25%, 50%, 75%)
    descriptive['50%'] = df.quantile(0.50)
    descriptive['75%'] = df.quantile(0.75)

    return descriptive
