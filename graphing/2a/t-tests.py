import pandas as pd
from scipy.stats import ttest_rel
import itertools

if __name__ == '__main__':
    # load in summary file and group data by model col, rows should be folds, each element is r2 score
    summary = pd.read_csv('./results/NaN_test/summary.csv')
        # Pivot the DataFrame
    pivot_df = summary.pivot_table(index='fold', columns='model', values='nan_r2', aggfunc='first')

    # Reorder columns
    unique_models = summary['model'].unique()
    pivot_df = pivot_df[unique_models]

    # Reset index name
    pivot_df.index.name = 'fold'

    # Print the resulting DataFrame
    pivot_df.to_csv('./graphing/summary_pivot.csv')

    # Assuming you have already created the pivot_df DataFrame

    # Create an empty DataFrame to store t-test results
    model_pairs = list(itertools.combinations(pivot_df.columns, 2))
    ttest_results = pd.DataFrame(index=pivot_df.columns, columns=pivot_df.columns)

    # Calculate t-tests and populate the DataFrame
    for model1, model2 in model_pairs:
        t_stat, p_value = ttest_rel(pivot_df[model1], pivot_df[model2], alternative='less')
        # p_value = min(1, p_value * len(model_pairs))
        ttest_results.loc[model1, model2] = p_value
        t_stat, p_value = ttest_rel(pivot_df[model2], pivot_df[model1], alternative='less')
        # p_value = min(1, p_value * len(model_pairs))
        ttest_results.loc[model2, model1] = p_value

    # Convert p-values to numeric values
    # ttest_results = ttest_results.apply(pd.to_numeric, errors='coerce')
    

    # Print the resulting t-test DataFrame
    ttest_results.to_csv('./graphing/ttest_results.csv')

    # for every pair of models, parse model names to get features used, 
    # then if they have only one different feature, create a new dataframe with model vs. model, the differnt feature, the p-value between them, and mean r2 scores of each model
    # then print to csv
    temp = []
    for model1, model2 in model_pairs:
        model1_features = model1.split('-')[1].split('_')
        model2_features = model2.split('-')[1].split('_')
        different_features = []
        for feature in model1_features:
            if feature not in model2_features:
                different_features.append(feature)
        for feature in model2_features:
            if feature not in model1_features:
                different_features.append(feature)
        if len(different_features) == 1:
            model1_r2 = pivot_df[model1].mean()
            model2_r2 = pivot_df[model2].mean()
            temp.append({'model1': model1, 'model2': model2, 'feature': different_features[0], 'p_value': ttest_results.loc[model1, model2], 'model1_r2': model1_r2, 'model2_r2': model2_r2})
            
    feature_comparison_ttest = pd.DataFrame(temp)
    feature_comparison_ttest.to_csv('./graphing/ttest_results_with_features.csv')

