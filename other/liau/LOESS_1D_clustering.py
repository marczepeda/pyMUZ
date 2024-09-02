# legacy from Kevin, did not delete because I haven't checked through the dependencies of the code
import os, re, sys, time, warnings
from pathlib import Path


import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# this is mainly so that my plots are generated with words as text/in Arial
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = 'Arial'

# specific to the loess function
import scipy.stats as stats
import scipy.spatial.distance as dist
import scipy.interpolate as interp
import scipy.optimize as sp_opt
import scipy.cluster as sp_cl


import statsmodels.api as sm
import statsmodels.stats.multitest as smm
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from Bio.Data import IUPACData

import statsmodels.stats.multitest as smm
import random

# this function calculates the smoothed arbitrary values for each value in your x_out (see below) for the enrichment profile of your POI
# @param: x_obs - the "x-axis" of your plot, usually some summarized aa number of edited window
# @param: y_obs - the "y-axis" of your plot, usually normalized lfc of gRNA enrichment
# @param: x_out - the values you want to interpolate over to get values for, usually a numpy array of discrete values for the whole protein (i.e. np.arange(1,len(POI)+1))
# @param: span - the span that you want to smooth over, requires a little optimization, a good place to start is ~0.05-0.10
# returns df_loess, a dataframe with a column with your x_out and the corresponding loess values in y_loess
def loess_v3(x_obs, y_obs, span, x_out=None, interp_how='quadratic', it=0, loess_kws=dict(missing='raise',return_sorted=False), interp_kws=dict(fill_value='extrapolate')):
    x_obs = x_obs.astype(float)
    if x_out is None:
        x_out = x_obs
    df_loess = pd.DataFrame()
    df_loess['x_vals'] = x_out.astype(float)
    df_loess = df_loess.sort_values(by='x_vals').reset_index(drop=True)
    if interp_how == 'statsmodels':
        #see statsmodel documentation for LOWESS for more info
        df_loess['y_loess'] = lowess(endog=y_obs, exog=x_obs, xvals=df_loess['x_vals'], frac=span, it=it, **loess_kws)
    else:
        df_interp = pd.DataFrame()
        df_interp['x_obs'] = x_obs
        df_interp['y_obs_loess'] = lowess(endog=y_obs, exog=x_obs, xvals=x_obs, frac=span, it=it, missing='drop')#, **loess_kws)
        df_interp = df_interp.groupby('x_obs',as_index=False).agg({'y_obs_loess':'mean'})
        fx_interp = interp.interp1d(x=df_interp['x_obs'], y=df_interp['y_obs_loess'], kind=interp_how, **interp_kws)
        df_loess['y_loess'] = fx_interp(df_loess['x_vals'])
    df_loess['type'] = np.where(df_loess['x_vals'].isin(x_obs), 'loess', 'interp')
    return df_loess

# this function provides the null distribution to compare your actual loess again; it calculates the smoothed arbitrary values for each value in your x_out (see below) for the shuffled enrichment profile of your POI
# @param: x_obs - the "x-axis" of your plot, usually some summarized aa number of edited window
# @param: y_obs - the "y-axis" of your plot, usually normalized lfc of gRNA enrichment
# @param: x_out - the values you want to interpolate over to get values for, usually a numpy array of discrete values for the whole protein (i.e. np.arange(1,len(POI)+1))
# @param: n_repeats - the number of times you want to shuffle, we usually do 10,000x for final data but that takes a while to run so preliminary sample is recommended for optimizatino (I've seen signal with as little as 4000x)
# @param: span - the span that you want to smooth over, requires a little optimization, a good place to start is ~0.05-0.10
# returns df_loess_rand, a dataframe with a column with n_repeats number of calculated loess scores from randomized y_obs
def randomize(x_obs, y_obs, x_out, n_repeats, span):
    df_loess_rand = pd.DataFrame()
    df_loess_rand = pd.concat([loess_v3(x_obs,random.sample(y_obs,len(y_obs)),x_out=x_out, span=span, interp_how='quadratic')[['y_loess']] for i in range(0,n_repeats)], ignore_index=True,axis=1)
    return(df_loess_rand)

# this function calculates the significance of a "spike" in your loess-ed enrichment profile by calculating if it's > 95% of null distribution
# @param: df_loess - output from loess_v3
# @param: df_loess_rand - output from randomize
# @param: n_repeats - the number of times you want to shuffle (see randomized function)
# returns df_pvals, a dataframe with your x_out, corresponding loess values (y_loess) and empirical p-val (sig) as well as the corrected p-values (corr_pvals)
def calculate_sig(df_loess,df_loess_rand,n_repeats):
    df_pvals = df_loess[['x_vals','y_loess']].copy()
    df_pvals['obs_gt'] = df_loess_rand.gt(df_pvals['y_loess'], axis=0).sum(axis=1) # get # of values greater than obs_val
    df_pvals['1t_pval'] = df_pvals['obs_gt'] / n_repeats # divide "rank" of obs val by N to get empirical p-val

    temp = smm.multipletests(pval_df['1t_pval'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False) # apply benjamini-hochberg FDR correction
    df_pvals['sig'] = temp[0]
    df_pvals['corr_pval'] = temp[1]
    return(df_pvals)