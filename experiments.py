import argparse
import os
from train import run_experiment, RANDOM_SEED, get_model_dir_name, nan_r2, masked_mse_loss, create_symbol_to_fold
import numpy as np
from tqdm import tqdm
import pandas as pd

def human_all(args):
    # TE_path = 'human_TE_cellline_all_plain.csv'
    TE_path = 'human_TE_cellline_all_NA_plain.csv'
    symbol_to_fold_path = 'symbol_to_fold.tsv'

    create_symbol_to_fold("./data/" + TE_path, "./data/" + symbol_to_fold_path)

    model_name = 'lgbm'

    feature_set = [
        'LL', 'P5', 'P3', 'CF', 'AAF',
        '3mer_freq_5', 'Struct'
    ]

    # metrics = [r2_score, pearson_corrcoef, spearman_corrcoef, mean_squared_error]
    metrics = [nan_r2, masked_mse_loss] # nan_pearson_corr

    name = get_model_dir_name(model_name, feature_set)
    run_experiment(
        args=args,
        TE_path=TE_path, 
        symbol_to_fold_path=symbol_to_fold_path, 
        model_name=model_name, 
        features_to_extract=feature_set,
        metrics=metrics,
        name=name,
        results_path='results/human/all_cell_lines/',
        bio_source='all',
        species='human',
        )    
    
def human_all_no_struct(args):
    # TE_path = 'human_TE_cellline_all_plain.csv'
    TE_path = 'human_TE_cellline_all_NA_plain.csv'
    symbol_to_fold_path = 'symbol_to_fold.tsv'

    create_symbol_to_fold("./data/" + TE_path, "./data/" + symbol_to_fold_path)

    model_name = 'lgbm'

    feature_set = [
        'LL', 'P5', 'P3', 'CF', 'AAF',
        '3mer_freq_5'
    ]

    # metrics = [r2_score, pearson_corrcoef, spearman_corrcoef, mean_squared_error]
    metrics = [nan_r2, masked_mse_loss] # nan_pearson_corr

    name = get_model_dir_name(model_name, feature_set)
    run_experiment(
        args=args,
        TE_path=TE_path, 
        symbol_to_fold_path=symbol_to_fold_path, 
        model_name=model_name, 
        features_to_extract=feature_set,
        metrics=metrics,
        name=name,
        results_path='results/human/all_cell_lines/',
        bio_source='all',
        species='human',
        )    

def mouse_all(args):
    # TE_path = 'mouse_TE_cellline_all_plain.csv'
    TE_path = 'mouse_TE_cellline_all_plain_NA.csv'
    symbol_to_fold_path = 'mouse_symbol_to_fold.tsv'

    create_symbol_to_fold("./data/" + TE_path, "./data/" + symbol_to_fold_path)

    model_name = 'lgbm'

    feature_set = [
        'LL', 'P5', 'P3', 'CF', 'AAF',
        '3mer_freq_5'
    ]

    # metrics = [r2_score, pearson_corrcoef, spearman_corrcoef, mean_squared_error]
    metrics = [nan_r2, masked_mse_loss] # nan_pearson_corr

    name = get_model_dir_name(model_name, feature_set)
    run_experiment(
        args=args,
        TE_path=TE_path, 
        symbol_to_fold_path=symbol_to_fold_path, 
        model_name=model_name, 
        features_to_extract=feature_set,
        metrics=metrics,
        name=name,
        results_path='results/mouse/all_cell_lines/',
        bio_source='all',
        species='mouse',
        )    

    
def human_model_comparison(args):
    TE_path = 'human_TE_cellline_all_NA_plain.csv'
    symbol_to_fold_path = 'symbol_to_fold.tsv'

    create_symbol_to_fold("./data/" + TE_path, "./data/" + symbol_to_fold_path)

    model_names = [
        'lgbm', 
        'lasso', 
        'elasticnet', 
        'randomforest'
    ]

    feature_set = [
        'LL', 'P5', 'P3', 'CF', 'AAF',
        '3mer_freq_5'
    ]

    # metrics = [r2_score, pearson_corrcoef, spearman_corrcoef, mean_squared_error]
    metrics = [nan_r2, masked_mse_loss] # nan_pearson_corr

    for model_name in model_names:
        name = get_model_dir_name(model_name, feature_set)
        run_experiment(
            args=args,
            TE_path=TE_path, 
            symbol_to_fold_path=symbol_to_fold_path, 
            model_name=model_name, 
            features_to_extract=feature_set,
            metrics=metrics,
            name=name,
            results_path='results/model_compare/',
            bio_source='mean_te',
            species='human',
            )

            
def human_model_comparison_all_features(args):
    TE_path = 'human_TE_cellline_all_NA_plain.csv'
    symbol_to_fold_path = 'symbol_to_fold.tsv'

    create_symbol_to_fold("./data/" + TE_path, "./data/" + symbol_to_fold_path)

    model_names = [
        # 'lgbm', 
        'lasso', 
        'elasticnet', 
        'randomforest'
    ]

    feature_set = ['LL', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'CF', 'AAF', 'DC', 
                   '1mer_freq_5', '1mer_freq_C', '1mer_freq_3', '2mer_freq_5', '2mer_freq_C', '2mer_freq_3', 
                   '3mer_freq_5', '3mer_freq_C', '3mer_freq_3', '4mer_freq_5', '4mer_freq_C', '4mer_freq_3', 
                   '5mer_freq_5', '5mer_freq_C', '5mer_freq_3', '6mer_freq_5', '6mer_freq_C', '6mer_freq_3', 
                   'Struct']


    # metrics = [r2_score, pearson_corrcoef, spearman_corrcoef, mean_squared_error]
    metrics = [nan_r2, masked_mse_loss] # nan_pearson_corr

    for model_name in model_names:
        name = get_model_dir_name(model_name, feature_set)
        run_experiment(
            args=args,
            TE_path=TE_path, 
            symbol_to_fold_path=symbol_to_fold_path, 
            model_name=model_name, 
            features_to_extract=feature_set,
            metrics=metrics,
            name=name,
            results_path='results/model_comparison_all_features/',
            bio_source='mean_te',
            species='human',
            ) 
        

def vst_coposittional_test(args):
    for TE_path in ['human_TE_cellline_all_plain_vst.csv', 'human_TE_cellline_all_NA_plain_vst.csv']:
        # TE_path = 'human_TE_cellline_all_plain_vst.csv'
        # TE_path = 'human_TE_cellline_all_NA_plain_vst.csv'
        symbol_to_fold_path = 'symbol_to_fold.tsv'

        create_symbol_to_fold("./data/" + TE_path, "./data/" + symbol_to_fold_path)

        model_name = 'lgbm'

        feature_set = [
            'LL', 'P5', 'P3', 'CF', 'AAF',
            '3mer_freq_5'
        ]

        # metrics = [r2_score, pearson_corrcoef, spearman_corrcoef, mean_squared_error]
        metrics = [nan_r2, masked_mse_loss] # nan_pearson_corr

        name = get_model_dir_name(model_name, feature_set)
        name = 'with_NA_' + name if 'NA' in TE_path else 'without_NA_' + name 
        run_experiment(
            args=args,
            TE_path=TE_path, 
            symbol_to_fold_path=symbol_to_fold_path, 
            model_name=model_name, 
            features_to_extract=feature_set,
            metrics=metrics,
            name=name,
            results_path='results/vst_copositional/',
            bio_source='mean_te',
            species='human',
            )    




if __name__ == '__main__':
    np.random.seed(RANDOM_SEED)

    tqdm.pandas()
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--experiment', type=str, required=True)
    parser.add_argument('--permutation_fi', '-p', action='store_true', help='Whether to run permutation feature importance')
    parser.add_argument('--save', '-s', action='store_true', help='Whether to save models')
    args = parser.parse_args()

    if args.experiment == 'human_all':
        human_all(args)
    elif args.experiment == 'mouse_all':
        mouse_all(args)
    elif args.experiment == 'human_all_no_struct':
        human_all_no_struct(args)
    elif args.experiment == 'human_model_compare':
        human_model_comparison(args)
    elif args.experiment == 'human_model_compare_all_features':
        human_model_comparison_all_features(args)
    elif args.experiment == 'vst':
        vst_coposittional_test(args)
    # elif args.experiment == 'mouse_model_compare':
    #     mouse_model_compare(args)
    # elif args.experiment == 'human_feature_tests':
    #     human_feature_tests(args)
    # elif args.experiment == 'mouse_feature_tests':
    #     mouse_feature_tests(args)
    
    else:
        raise ValueError(f'Experiment {args.experiment} not found')
