import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import multipletests


def check_intervals_intersect(first_ci, second_ci):
    L1, R1 = first_ci
    L2, R2 = second_ci
    if R1 >= L2 and R2 >= L1:
        return True # confidence intervals intersect
    else:
        return False

def check_dge_with_ci(first_table, second_table):
    genes = list(first_table.columns[:-1])
    ci_test_results = []
    for gene in genes:
        ci_1 = st.t.interval(alpha=0.95,
              df=len(first_table[gene]) - 1,
              loc=np.mean(first_table[gene]),
              scale=st.sem(first_table[gene]))
        ci_2 = st.t.interval(alpha=0.95,
              df=len(second_table[gene]) - 1,
              loc=np.mean(second_table[gene]),
              scale=st.sem(second_table[gene]))
        ci_test_results.append(check_intervals_intersect(ci_1, ci_2))
    return ci_test_results

def check_dge_with_ztest(first_table, second_table, return_pvalue=False):
    genes = list(first_table.columns[:-1])
    z_test_results = []
    z_test_pvalues = []
    for gene in genes:
        _, pvalue = ztest(first_table[gene], second_table[gene])
        z_test_results.append(pvalue < 0.05)
        z_test_pvalues.append(round(pvalue, 4))
    if return_pvalue == False:
        return z_test_results
    else:
        return z_test_pvalues

def mean_diff(first_table, second_table):
    genes = list(first_table.columns[:-1])
    mean_diff = []
    for gene in genes:
        mean_1 = np.mean(first_table[gene])
        mean_2 = np.mean(second_table[gene])
        mean_diff.append(round(mean_2 - mean_1, 1))
    return mean_diff

# Receiving paths to files:
first_cell_type_expressions_path = input('Path to file with expression data for first cell type: ')
second_cell_type_expressions_path = input('Path to file with expression data for second cell type: ')
save_results_table = input('Path to results: ')

# (Used this to check program, don't pay attention)
# first_cell_type_expressions_path = '/Users/amy/Desktop/BI2022/statistics/BI_statistics_2022/hw_6/data/raw_data/b_cells_expression_data.csv'
# second_cell_type_expressions_path = '/Users/amy/Desktop/BI2022/statistics/BI_statistics_2022/hw_6/data/raw_data/nk_cells_expression_data.csv'
# save_results_table = '/Users/amy/Desktop/BI2022/statistics/BI_statistics_2022/hw_6/data/processed_data/results.csv'

# Reading data from files:
first_cell_type_expression_data = pd.read_csv(first_cell_type_expressions_path, index_col=0)
second_cell_type_expression_data = pd.read_csv(second_cell_type_expressions_path, index_col=0)

# Making list of CI test results:
ci_test_results = check_dge_with_ci(first_cell_type_expression_data, second_cell_type_expression_data)

# Making list of Z test results:
z_test_results = check_dge_with_ztest(first_cell_type_expression_data, 
second_cell_type_expression_data)

# Making list of p-values from Z test results:
z_test_p_values = check_dge_with_ztest(first_cell_type_expression_data, 
second_cell_type_expression_data, return_pvalue=True)

# Making list with differences between mean expressions of each gene from first and second tables
mean_diff = mean_diff(first_cell_type_expression_data, second_cell_type_expression_data)

# Asking user to choose the method of multiple correction
method = input('Insert method of multiple correction: ')

# Creating dictionary {'column name': list_of_values} with results
if method == '':
    results = {
    "ci_test_results": ci_test_results,
    "z_test_results": z_test_results,
    "z_test_p_values": z_test_p_values,
    "mean_diff": mean_diff
    }
else:    # Making multiple correction of p-values
    results = {
    "ci_test_results": ci_test_results,
    "z_test_results": z_test_results,
    "z_test_p_values_corrected": multipletests(z_test_p_values, method=method)[1],
    "mean_diff": mean_diff
    }

# Turning dictionary into dataframe
results = pd.DataFrame(results)
results.head()

# Adding first column with gene names
results.insert(loc=0, column='Genes', value=first_cell_type_expression_data.columns[:-1])
results

# Saving dataframe into .csv file
results.to_csv(save_results_table)