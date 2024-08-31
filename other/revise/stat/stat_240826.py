### stat.py ###
# Author: Marc Zepeda
# Date: 2024-08-26
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import skew, kurtosis, mode, binom, poisson

## Objectives ##
''' 
Combine these methods into smaller groups of methods
Include reminder for what they are for...
Make it easy to include after or along with plot.py
'''

## BST281 Lectures ##
'''
Lecture notes (readings on page 2):
Test statistics summarize data into a single value that captures one or more relevant aspects.
Tests consider a null hypothesis, usually that there is no effect/association/bias/etc.
	The test statistic should have some known distribution, the null distribution, when the null hypothesis is true.
P-value: probability of observing a test statistic at least as extreme if null hypothesis true (from null distribution).
	Null hypothesis rejected if p-value < critical threshold α, calling the result statistically significant.
For normal distributions with known standard deviation, z-test statistic z = (x - μ ̂) / σ is appropriate.
	When the standard deviation is unknown, t = (x - μ ̂) / σ ̂ follows the Student's t-distribution; the t-test.
One-sided/one-tailed tests consider only values in one direction to be "extreme."
	E.g. asking whether a gene's expression is greater than another.
	Two-sided/two-tailed tests consider both directions; whether a gene's expression is different from another.
Non-parametric tests ignore the shape of the distribution, often using a rank transformation.
	But the cost is reduced sensitivity - they are discarding potentially useful information!
	E.g. Mann-Whitney U test, also known as the Wilcoxon rank-sum test.
Permutation test is used when the null distribution is not nice.
	Generate an empirical null distribution directly by calculating test statistic repeatedly in permuted data.
Performance evaluation: perform a test where true values are known (gold standard).
	Error rates: false positive rate (Type I errors), and false negative rate (Type II errors).
	How well the test calls positives: power, precision, and specificity.
	Precision/recall plots and ROC plots capture how well a test does, independent of its parameters.
Area Under the Curve (AUC) commonly used to assess performance (0.5 is random, 1 is perfect).
Multiple hypothesis testing can lead to many false positives.
	Bonferroni correction: conservatively only call α false positives among all tests.
	False Discovery Rate (FDR) correction: instead control the fraction of false positives among all positives.
'''

'''
•	Descriptive statistics summarize information that is already present in a data collection
o	Plots (box plots, histograms, etc.), measures like average, standard deviation, medians, etc.
•	Inferential statistics use a sample of data to make predictions about larger populations or about unobserved/future trends
o	Measurements in the presence of noise, generalizations from a sample to a population, comparisons between datasets: correlations, regression, etc.
•	A statistic is any single value that summarizes an entire dataset
•	Parametric summary statistics describe “well-behaved” data, often normally-distributed
o	Mean, standard deviation, z-scores, etc.
•	Nonparametric summary statistics can be used to describe data regardless of distribution
o	Median, percentiles, quartiles, interquartile range, etc.
•	Summary statistics of paired data allow comparison between datasets
•	Distances quantify differences between datasets; larger means less similar
o	Euclidean (L2 norm), Manhattan (L1 norm)
•	Correlations quantify the similarity of datasets; larger means more similar
o	Pearson, Cosine (parametric), Spearman (nonparametric)
•	Beware of summary statistics: Understand (and visualize) your data before summarizing them!
o	Anscombe’s quartet: a four manually-crafted pairs of datasets with equal means, standard deviations, and Pearson correlations, yet completely different relationships
•	Basics of probability
o	Thinking in terms of the probability of sets of outcomes (events) occurring as subsets of the set of all possible outcomes (sample space)
•	Kolmogorov axioms: one definition of probability that matches reality
•	Bayes’ theorem: Calculating a conditional probability based on the inverse of its condition
'''

## Suggestions from ChatGPT ##
# write a python package for statistical analysis of tidy data sets
def summarize(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return summary statistics for a tidy DataFrame.
    """
    return df.describe(include='all')

def correlation_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return a correlation matrix for numerical columns in the DataFrame.
    """
    return df.corr()

def t_test(sample1, sample2) -> dict:
    """
    Perform an independent t-test between two samples.
    """
    t_stat, p_value = stats.ttest_ind(sample1, sample2)
    return {'t_stat': t_stat, 'p_value': p_value}

def chi_square_test(contingency_table: pd.DataFrame) -> dict:
    """
    Perform a Chi-square test for independence on a contingency table.
    """
    chi2, p, dof, expected = stats.chi2_contingency(contingency_table)
    return {'chi2': chi2, 'p_value': p, 'dof': dof, 'expected': expected}

def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Perform basic data cleaning like handling missing values.
    """
    return df.dropna()

# write code for 1-way anova, wilconox test, etc.
def one_way_anova(*groups) -> dict:
    """
    Perform a one-way ANOVA test.
    
    Parameters:
    *groups: variable number of group arrays (e.g., group1, group2, group3)
    
    Returns:
    dict: containing F-statistic and p-value
    """
    f_stat, p_value = stats.f_oneway(*groups)
    return {'f_stat': f_stat, 'p_value': p_value}

def wilcoxon_test(sample1, sample2) -> dict:
    """
    Perform a Wilcoxon signed-rank test.
    
    Parameters:
    sample1: array-like, first set of paired observations
    sample2: array-like, second set of paired observations
    
    Returns:
    dict: containing test statistic and p-value
    """
    test_stat, p_value = stats.wilcoxon(sample1, sample2)
    return {'test_stat': test_stat, 'p_value': p_value}

def kruskal_wallis_test(*groups) -> dict:
    """
    Perform a Kruskal-Wallis H test.
    
    Parameters:
    *groups: variable number of group arrays (e.g., group1, group2, group3)
    
    Returns:
    dict: containing H-statistic and p-value
    """
    h_stat, p_value = stats.kruskal(*groups)
    return {'h_stat': h_stat, 'p_value': p_value}

def mann_whitney_u_test(sample1, sample2) -> dict:
    """
    Perform a Mann-Whitney U test.
    
    Parameters:
    sample1: array-like, first group of observations
    sample2: array-like, second group of observations
    
    Returns:
    dict: containing U-statistic and p-value
    """
    u_stat, p_value = stats.mannwhitneyu(sample1, sample2)
    return {'u_stat': u_stat, 'p_value': p_value}

def chi_square_test(contingency_table) -> dict:
    """
    Perform a Chi-square test for independence.
    
    Parameters:
    contingency_table: pandas DataFrame, contingency table of observed frequencies
    
    Returns:
    dict: containing Chi-square statistic, p-value, degrees of freedom, and expected frequencies
    """
    chi2_stat, p_value, dof, expected = stats.chi2_contingency(contingency_table)
    return {'chi2_stat': chi2_stat, 'p_value': p_value, 'dof': dof, 'expected': expected}

def paired_t_test(sample1, sample2) -> dict:
    """
    Perform a paired t-test.
    
    Parameters:
    sample1: array-like, first set of paired observations
    sample2: array-like, second set of paired observations
    
    Returns:
    dict: containing t-statistic and p-value
    """
    t_stat, p_value = stats.ttest_rel(sample1, sample2)
    return {'t_stat': t_stat, 'p_value': p_value}

def paired_t_test(sample1, sample2) -> dict:
    """
    Perform a paired t-test.
    
    Parameters:
    sample1: array-like, first set of paired observations
    sample2: array-like, second set of paired observations
    
    Returns:
    dict: containing t-statistic and p-value
    """
    t_stat, p_value = stats.ttest_rel(sample1, sample2)
    return {'t_stat': t_stat, 'p_value': p_value}

def pearson_correlation(x, y) -> dict:
    """
    Calculate the Pearson correlation coefficient.
    
    Parameters:
    x: array-like, first variable
    y: array-like, second variable
    
    Returns:
    dict: containing correlation coefficient and p-value
    """
    corr_coef, p_value = stats.pearsonr(x, y)
    return {'correlation_coefficient': corr_coef, 'p_value': p_value}

# write methods for permutation tests and multiple hypothesis testing

def permutation_test(sample1, sample2, num_permutations=10000, statistic=np.mean) -> dict:
    """
    Perform a permutation test.
    
    Parameters:
    sample1: array-like, first group of observations
    sample2: array-like, second group of observations
    num_permutations: int, number of permutations to perform (default=10000)
    statistic: callable, function to compute the statistic (default=np.mean)
    
    Returns:
    dict: containing the observed statistic, p-value, and the distribution of permuted statistics
    """
    # Calculate the observed statistic
    observed_stat = statistic(sample1) - statistic(sample2)
    
    # Combine samples
    combined = np.concatenate([sample1, sample2])
    
    # Permute and calculate the statistic
    permuted_stats = []
    for _ in range(num_permutations):
        np.random.shuffle(combined)
        permuted_sample1 = combined[:len(sample1)]
        permuted_sample2 = combined[len(sample1):]
        permuted_stat = statistic(permuted_sample1) - statistic(permuted_sample2)
        permuted_stats.append(permuted_stat)
    
    # Calculate p-value
    p_value = np.mean(np.abs(permuted_stats) >= np.abs(observed_stat))
    
    return {'observed_stat': observed_stat, 'p_value': p_value, 'permuted_stats': permuted_stats}

# tidystats/hypothesis.py

def permutation_test(sample1, sample2, num_permutations=10000, statistic=np.mean) -> dict:
    """
    Perform a permutation test.
    
    Parameters:
    sample1: array-like, first group of observations
    sample2: array-like, second group of observations
    num_permutations: int, number of permutations to perform (default=10000)
    statistic: callable, function to compute the statistic (default=np.mean)
    
    Returns:
    dict: containing the observed statistic, p-value, and the distribution of permuted statistics
    """
    # Calculate the observed statistic
    observed_stat = statistic(sample1) - statistic(sample2)
    
    # Combine samples
    combined = np.concatenate([sample1, sample2])
    
    # Permute and calculate the statistic
    permuted_stats = []
    for _ in range(num_permutations):
        np.random.shuffle(combined)
        permuted_sample1 = combined[:len(sample1)]
        permuted_sample2 = combined[len(sample1):]
        permuted_stat = statistic(permuted_sample1) - statistic(permuted_sample2)
        permuted_stats.append(permuted_stat)
    
    # Calculate p-value
    p_value = np.mean(np.abs(permuted_stats) >= np.abs(observed_stat))
    
    return {'observed_stat': observed_stat, 'p_value': p_value, 'permuted_stats': permuted_stats}

def bonferroni_correction(p_values, alpha=0.05) -> dict:
    """
    Perform Bonferroni correction for multiple hypothesis testing.
    
    Parameters:
    p_values: list or array-like, p-values from individual tests
    alpha: float, significance level (default=0.05)
    
    Returns:
    dict: containing the adjusted significance level and a list of rejected hypotheses
    """
    corrected_alpha = alpha / len(p_values)
    rejected = [p <= corrected_alpha for p in p_values]
    return {'corrected_alpha': corrected_alpha, 'rejected': rejected}

def benjamini_hochberg(p_values, alpha=0.05) -> dict:
    """
    Perform Benjamini-Hochberg procedure for multiple hypothesis testing.
    
    Parameters:
    p_values: list or array-like, p-values from individual tests
    alpha: float, significance level (default=0.05)
    
    Returns:
    dict: containing a list of rejected hypotheses
    """
    p_values = np.array(p_values)
    n = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p_values = p_values[sorted_indices]
    
    # Calculate the threshold for each p-value
    thresholds = np.arange(1, n+1) * (alpha / n)
    
    # Determine which p-values are significant
    rejected = sorted_p_values <= thresholds
    
    # Return to the original order
    original_order_rejected = np.zeros(n, dtype=bool)
    original_order_rejected[sorted_indices] = rejected
    
    return {'rejected': original_order_rejected}

# write a methods for descriptive statistics

def descriptive_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate descriptive statistics for each numerical column in the DataFrame.
    
    Parameters:
    df: pandas DataFrame, the dataset for which descriptive statistics are calculated
    
    Returns:
    pd.DataFrame: a DataFrame containing the descriptive statistics
    """
    desc_stats = pd.DataFrame()

    # Mean
    desc_stats['mean'] = df.mean()

    # Median
    desc_stats['median'] = df.median()

    # Mode (returning the first mode in case of multiple modes)
    desc_stats['mode'] = df.apply(lambda x: mode(x, nan_policy='omit')[0][0])

    # Variance
    desc_stats['variance'] = df.var()

    # Standard Deviation
    desc_stats['std_dev'] = df.std()

    # Minimum
    desc_stats['min'] = df.min()

    # Maximum
    desc_stats['max'] = df.max()

    # Range
    desc_stats['range'] = df.max() - df.min()

    # Skewness
    desc_stats['skewness'] = df.apply(lambda x: skew(x, nan_policy='omit'))

    # Kurtosis
    desc_stats['kurtosis'] = df.apply(lambda x: kurtosis(x, nan_policy='omit'))

    # Count (non-missing observations)
    desc_stats['count'] = df.count()

    # Sum
    desc_stats['sum'] = df.sum()

    # Quantiles (25%, 50%, 75%)
    desc_stats['25%'] = df.quantile(0.25)
    desc_stats['50%'] = df.quantile(0.50)
    desc_stats['75%'] = df.quantile(0.75)

    return desc_stats

# write a method for inferential statistics

def inferential_statistics(data, alpha=0.05) -> dict:
    """
    Perform inferential statistics including confidence interval estimation and hypothesis testing.
    
    Parameters:
    data: array-like, the sample data
    alpha: float, significance level for confidence intervals (default=0.05)
    
    Returns:
    dict: containing the sample mean, standard error, confidence interval, and p-value for a hypothesis test
    """
    # Calculate the sample mean
    sample_mean = np.mean(data)
    
    # Calculate the sample standard deviation
    sample_std = np.std(data, ddof=1)
    
    # Calculate the sample size
    n = len(data)
    
    # Calculate the standard error of the mean
    standard_error = sample_std / np.sqrt(n)
    
    # Confidence Interval
    confidence_interval = stats.t.interval(1-alpha, df=n-1, loc=sample_mean, scale=standard_error)
    
    # Hypothesis Test (Two-tailed t-test against population mean 0)
    t_statistic, p_value = stats.ttest_1samp(data, popmean=0)
    
    return {
        'sample_mean': sample_mean,
        'standard_error': standard_error,
        'confidence_interval': confidence_interval,
        't_statistic': t_statistic,
        'p_value': p_value
    }

# write methods for probability

def pmf_binomial(n, p, k) -> float:
    """
    Calculate the probability mass function (PMF) for a binomial distribution.
    
    Parameters:
    n: int, number of trials
    p: float, probability of success in each trial
    k: int, number of successes
    
    Returns:
    float: probability of getting exactly k successes in n trials
    """
    return binom.pmf(k, n, p)

def pmf_poisson(lambda_, k) -> float:
    """
    Calculate the probability mass function (PMF) for a Poisson distribution.
    
    Parameters:
    lambda_: float, average rate (mean) of occurrence
    k: int, number of occurrences
    
    Returns:
    float: probability of observing exactly k occurrences
    """
    return poisson.pmf(k, lambda_)

def cdf_binomial(n, p, k) -> float:
    """
    Calculate the cumulative distribution function (CDF) for a binomial distribution.
    
    Parameters:
    n: int, number of trials
    p: float, probability of success in each trial
    k: int, number of successes
    
    Returns:
    float: cumulative probability of getting up to and including k successes
    """
    return binom.cdf(k, n, p)

def cdf_poisson(lambda_, k) -> float:
    """
    Calculate the cumulative distribution function (CDF) for a Poisson distribution.
    
    Parameters:
    lambda_: float, average rate (mean) of occurrence
    k: int, number of occurrences
    
    Returns:
    float: cumulative probability of observing up to and including k occurrences
    """
    return poisson.cdf(k, lambda_)

from scipy.stats import norm, expon

def pdf_normal(mu, sigma, x) -> float:
    """
    Calculate the probability density function (PDF) for a normal distribution.
    
    Parameters:
    mu: float, mean of the distribution
    sigma: float, standard deviation of the distribution
    x: float, value to evaluate the PDF at
    
    Returns:
    float: the probability density at x
    """
    return norm.pdf(x, mu, sigma)

def pdf_exponential(lambda_, x) -> float:
    """
    Calculate the probability density function (PDF) for an exponential distribution.
    
    Parameters:
    lambda_: float, rate parameter (1/mean)
    x: float, value to evaluate the PDF at
    
    Returns:
    float: the probability density at x
    """
    return expon.pdf(x, scale=1/lambda_)

def cdf_normal(mu, sigma, x) -> float:
    """
    Calculate the cumulative distribution function (CDF) for a normal distribution.
    
    Parameters:
    mu: float, mean of the distribution
    sigma: float, standard deviation of the distribution
    x: float, value to evaluate the CDF at
    
    Returns:
    float: cumulative probability up to and including x
    """
    return norm.cdf(x, mu, sigma)

def cdf_exponential(lambda_, x) -> float:
    """
    Calculate the cumulative distribution function (CDF) for an exponential distribution.
    
    Parameters:
    lambda_: float, rate parameter (1/mean)
    x: float, value to evaluate the CDF at
    
    Returns:
    float: cumulative probability up to and including x
    """
    return expon.cdf(x, scale=1/lambda_)

def joint_probability(p_a, p_b_given_a) -> float:
    """
    Calculate the joint probability of two events A and B.
    
    Parameters:
    p_a: float, probability of event A
    p_b_given_a: float, probability of event B given that A has occurred
    
    Returns:
    float: joint probability of A and B
    """
    return p_a * p_b_given_a

def marginal_probability(p_joint_ab, p_b) -> float:
    """
    Calculate the marginal probability of event A.
    
    Parameters:
    p_joint_ab: float, joint probability of events A and B
    p_b: float, probability of event B
    
    Returns:
    float: marginal probability of A
    """
    return p_joint_ab / p_b

def conditional_probability(p_joint_ab, p_a) -> float:
    """
    Calculate the conditional probability of B given A.
    
    Parameters:
    p_joint_ab: float, joint probability of events A and B
    p_a: float, probability of event A
    
    Returns:
    float: conditional probability of B given A
    """
    return p_joint_ab / p_a