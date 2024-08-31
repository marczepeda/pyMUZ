### stat.py ###
# Author: Marc Zepeda
# Date: 2024-08-26

import pandas as pd
import numpy as np
from scipy.stats import skew, kurtosis, ttest_ind, ttest_rel, f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd

''' describe: Calculate descriptive statistics for each numerical column in the DataFrame.df: DataFrame
        df: DataFrame
        cols: list of numerical columns to compute statistics (optional)
        group: column name to split tidy dataframes (optional)
    Dependencies: pandas, numpy, scipy.stats
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

''' difference: Statistical test flow chart
        df: tidy DataFrame
        data_col: data column name
        compare_col: comparisons column name
        compare: list of comparisons
        group: column name to split tidy dataframes (optional)
        factors: # of factors
        same: same (True) or different (False) subjects
        compare: list of comparisons
        para: parameteric (True) or nonparametric (False)
        alpha: significance level for the test (Default: 0.05)
    Dependencies:
'''
def difference(df: pd.DataFrame(), data_col: str(), compare_col: str(), compare: list(), factors=1, same=False, para=True, alpha=0.05):

    if factors==1: # quantity of factors (condition 2 levels?)

        if not same: # Different samples

            if para==True: # Data that follows a normal distribution

                if len(compare)==2: # Comparing 2 samples

                    # Student's T Test Calculation
                    print('Statistical Test: Student\'s T-test \n - Only use if you have two groups or if you are comparing just two of the n groups and are not concerned about inflating the Type I error rate (e.g., False Positives).')
                    t_stat, p_value = ttest_ind(*(df[df[compare_col]==comp][data_col] for comp in compare))
                    inference = pd.DataFrame({'test':['Student\'s T-test']*2,
                                              'comparison': [','.join(compare)]*2, 
                                              'statistic': ['t_stat','p_value'],
                                              'value':[t_stat,p_value]})
                    
                    return inference
                
                elif len(compare)>2: # Comparing more than 2 samples
                    
                    # 1-way Anova Test
                    print('Statistical Test: 1-way Anova\n - If you want to compare all three groups simultaneously to determine if there is a significant difference in their means.\n - If the ANOVA finds a significant difference, then perform post-hoc tests to determine which specific groups differ from each other.')
                    f_stat, p_value = f_oneway(*(df[df[compare_col]==comp][data_col] for comp in compare))
                    inference = pd.DataFrame({'test':['1-way Anova Test']*2,
                                              'comparison': [','.join(compare)]*2, 
                                              'statistic': ['f_stat','p_value'],
                                              'value':[f_stat,p_value]})
                    
                    # Tukey's Honestly Significant Difference Test
                    print('Statistical Test: Tukey\'s Honestly Significant Difference (HSD) Test\n - It\'s suitable when you have equal or unequal group sizes.\n - Other options include Bonferroni Correction, Scheffé\'s Test, Dunnett\'s Test, and Holm\'s Sequential Bonferroni Procedure.')
                    df2 = pd.concat([df[df[compare_col]==comp] for comp in compare]).reset_index(drop=True) # Only include specified data
                    tukey_result = pairwise_tukeyhsd(df2[data_col], df2[compare_col], alpha=alpha)
                    tukey_df = pd.DataFrame(data=tukey_result._results_table.data[1:], columns=tukey_result._results_table.data[0])
                    inference = pd.concat([inference,
                                           pd.DataFrame({'test':['Tukey\'s HSD Test']*len(compare),
                                                         'comparison': [','.join([tukey_df.iloc[i]['group1'],tukey_df.iloc[i]['group2']]) for i in range(len(tukey_df))], 
                                                         'statistic': ['p-adj']*len(compare),
                                                         'value': tukey_df['p-adj']})]).reset_index(drop=True)

                    return inference
                
                else: print('Error: Invalid compare. List needs to contain 2 or more stings')
            
            else: # Data that does not follow a normal distribution
                print('Statistical Test: Mann Whitney U Test\nWIP...')
                
        
        else: # Same sample

            if para==True:  # Data that follows a normal distribution

                if len(compare)==2: # Comparing 2 related (paired) samples
                    
                    # Student's Paired T-test
                    print('Statistical Test: Paired Student\'s T-test \n  - Only use if you have two groups or if you are comparing just two of the n groups and are not concerned about inflating the Type I error rate (e.g., False Positives).\n - Assumes the two groups are related or paired, meaning that each data point in one group has a corresponding data point in the other group.')
                    t_stat, p_value = ttest_rel(*(df[df[compare_col]==comp][data_col] for comp in compare))
                    inference = pd.DataFrame({'test':['Paired Student\'s T-test']*2,
                                              'comparison': [','.join(compare)]*2, 
                                              'statistic': ['t_stat','p_value'],
                                              'value':[t_stat,p_value]})
                    
                    return inference

                elif len(compare)>2: # Comparing 2 or more related (paired) samples
                    
                    # Repeated Anova
                    print('Statistical Test: Repeated 1-way Anova\nWIP...')
                    df2 = pd.concat([df[df[compare_col]==comp] for comp in compare]).reset_index(drop=True) # Only include specified data
                    '''anova = AnovaRM(data=df2, depvar=data_col, subject=compare_col)
                    anova_result = anova.fit()
                    print(anova_result)'''

                else: print('Error: Invalid compare. List needs to contain 2 or more stings')

            else: 
                print('Statistical Test: Wilcocson Paired Test\nWIP...')
    
    # Deal with later...
    elif factors>1: print('Statistical Test: 2 way anova, General Linear (Mixed) Model, etc.\nNot Included...')
    else: print('Error: Invalid factors. Number needs to be greater or equal to 1.')

''' correlation: 
'''