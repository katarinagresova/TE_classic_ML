import pandas as pd
from train import nan_r2, masked_mse_loss
import os
import numpy as np
import argparse

def performance_of_mean_of_preds(preds_file: str, symbol_to_fold_file: str, model_results_file: str, force: bool):
    preds_df = pd.read_csv(preds_file)
    model_results = pd.read_csv(model_results_file, index_col=0)
    # if force:
    #     model_results = model_results[model_results['bio_source'] != 'mean_across_cell_lines']
    #     preds_df = preds_df.drop(columns=['mean_across_cell_lines_pred', 'mean_across_cell_lines_true', 'fold'])


    assert 'mean_across_cell_lines' not in model_results['bio_source'].unique()
    assert 'mean_te_true' in preds_df.columns
    assert 'mean_across_cell_lines_pred' not in preds_df.columns
    assert 'mean_across_cell_lines_true' not in preds_df.columns

    pred_cols = [col for col in preds_df.columns if '_pred' in col]
    pred_cols.remove('mean_te_pred')
    preds_df['mean_across_cell_lines_true'] = preds_df['mean_te_true']
    preds_df['mean_across_cell_lines_pred'] = preds_df[pred_cols].mean(axis=1)
    symbol_to_fold = pd.read_csv(symbol_to_fold_file, index_col=None, sep='\t')
    preds_df = pd.merge(preds_df, symbol_to_fold, left_on='SYMBOL', right_on='SYMBOL', how='left')

    preds_df.to_csv(preds_file, index=False)

    new_results_df = pd.DataFrame(columns=['bio_source', 'fold', 'nan_r2', 'masked_mse_loss', 'total_time(sec)'])
    max_fold = preds_df['fold'].max()
    for fold in range(max_fold + 1):
        labels = np.array(preds_df[preds_df['fold'] == fold]['mean_across_cell_lines_true'])
        mean_of_preds = np.array(preds_df[preds_df['fold'] == fold]['mean_across_cell_lines_pred'])
        new_results_df.loc[len(new_results_df)] = [
            'mean_across_cell_lines', 
            fold, 
            nan_r2(mean_of_preds, labels),
            masked_mse_loss(mean_of_preds, labels).item(),
            0.0
        ]
    new_results_df.loc[len(new_results_df)] = [
        'mean_across_cell_lines', 
        'mean', 
        new_results_df['nan_r2'].mean(),
        new_results_df['masked_mse_loss'].mean(),
        0.0
    ]
    
    model_results = pd.concat([model_results, new_results_df], axis=0, ignore_index=True)
    model_results.to_csv(model_results_file)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--preds_file', '-p', type=str, required=True)
    parser.add_argument('--symbol_to_fold_file', '-s', type=str, required=True)
    parser.add_argument('--model_results_file', '-m', type=str, required=True)
    parser.add_argument('--force', '-f', action='store_true', default=False)
    args = parser.parse_args()

    # preds_file = 'results/human/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/predictions.csv'
    # symbol_to_fold_file = 'data/symbol_to_fold.tsv'
    # model_results_file = 'results/human/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/model_results.csv'
    performance_of_mean_of_preds(args.preds_file, args.symbol_to_fold_file, args.model_results_file, args.force)