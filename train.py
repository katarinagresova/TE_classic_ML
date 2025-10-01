import argparse
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LassoCV, ElasticNetCV
from sklearn.ensemble import RandomForestRegressor
import lightgbm as lgbm
import time
from lgbm_feature_extract_from_str import dataframe_feature_extract
from data import get_data
from tqdm import tqdm
import os
import numpy as np
from sklearn.model_selection import train_test_split, KFold
from typing import Union
import pickle

RANDOM_SEED = 0

# From Dinghai's code
def nan_r2(
    predictions: np.ndarray, targets: np.ndarray, axis: int = 0
) -> Union[np.ndarray, float]:
    """Calculate R2 between targets and predictions along axis.

    Args:
        targets (np.ndarray): 1D or 2D, may contain nan. If 2D, axis 0 is the
        sample dimension and axis 2 is the target dimenson.
        predictions (np.ndarray): an array with the same shape of targets, may
        contain nan
        axis (int, optional): the axis along which to calculate R2. Defaults
        to 0.

    Returns:
        Union[np.ndarray, float]: R2 (1 less dimension than targets)
    """
    rss = np.nansum((targets - predictions) ** 2, axis)
    tss = np.nansum((targets - np.nanmean(targets, axis, keepdims=True)) ** 2, axis)
    return 1 - rss / tss

# From Dinghai's code
def nan_pearson_corr(predictions: np.ndarray, targets: np.ndarray, axis: int = 0) -> np.ndarray:
    """Calculate Pearson's correlation coefficients between targets and predictions along 'axis'.

    Args:
        targets (np.ndarray): a 2D array, may contain nan
        predictions (np.ndarray): an array with the same shape of targets, may contain nan
        axis (int, optional): the axis along which to calculate correlation coefficients. Defaults to 0.

    Returns:
        np.ndarray: Pearson's correlation coefficients
    """
    if axis == 0:
        targets = targets.T
        predictions = predictions.T
    return np.array(
        [
            np.corrcoef(t[~np.isnan(t)], p[~np.isnan(t)])[0, 1]
            for (t, p) in zip(targets, predictions)
        ]
    )
# From Dinghai's code
def masked_mse_loss(input, target, masked_value="NA", reduction="mean"):
    if masked_value == "NA":
        # target is numpy.ndarray. get mask where target is nan
        mask = np.isnan(target)
    else:
        mask = target == masked_value
    out = (input[~mask] - target[~mask]) ** 2
    if reduction == "mean":
        return out.mean()
    elif reduction == "None":
        return out

def evaluate(y_preds: np.ndarray, y_true: np.ndarray, metrics):
    y_true = y_true.squeeze()
    results = {}
    for metric in metrics:
        x = metric(y_preds, y_true)
        results[metric.__name__] = x
    return results

def train_model(args, data, fold, model_name, features_to_extract, metrics, bio_source='mean_te'):
    start_time = time.time()    

    train_data = data[data['fold'] != fold]
    test_data = data[data['fold'] == fold]
    train_features, train_labels = None, None
    test_features, test_labels = None, None
    if bio_source == 'multi_task':
        assert model_name == 'lgbm'
        all_features, all_labels = dataframe_feature_extract(data, features_to_extract, multi_task=True)
        train_features, test_features, train_labels, test_labels = train_test_split(all_features, all_labels, test_size=0.2, random_state=RANDOM_SEED)
        # train_features, train_labels = multi_task_dataframe_feature_extract(train_data, features_to_extract)
        # test_features, test_labels = multi_task_dataframe_feature_extract(test_data, features_to_extract)
    else:
        train_features, train_labels = dataframe_feature_extract(train_data, features_to_extract, bio_source)
        test_features, test_labels = dataframe_feature_extract(test_data, features_to_extract, bio_source)

    # print('train_features.shape: ', train_features.shape)
    # print('train_labels.shape: ', train_labels.shape)
    # print('test_features.shape: ', test_features.shape)
    # print('test_labels.shape: ', test_labels.shape)

    train_features = train_features.drop(columns=['transcript_id', 'gene_id'])
    test_features = test_features.drop(columns=['transcript_id', 'gene_id'])

    model = None
    if model_name == 'lgbm':
        model = lgbm.LGBMRegressor(n_jobs=12, importance_type='gain', random_state=RANDOM_SEED)
        if bio_source == 'multi_task':
            model = lgbm.LGBMRegressor(n_jobs=12, importance_type='gain', random_state=RANDOM_SEED)
    elif model_name == 'lasso' or model_name == 'elasticnet':
        train_features.replace(np.inf, 1600, inplace=True)
        train_features.replace(-np.inf, -1600, inplace=True)
        train_features.replace(np.nan, 1600, inplace=True)
        test_features.replace(np.inf, 1600, inplace=True)
        test_features.replace(-np.inf, -1600, inplace=True)
        test_features.replace(np.nan, 1600, inplace=True)
        model_type = LassoCV if model_name == 'lasso' else ElasticNetCV
        model = Pipeline([
            ('scaler', StandardScaler()),
            (model_name, model_type(cv=5, n_jobs=12, max_iter=10000, random_state=RANDOM_SEED))
        ])
    elif model_name == 'randomforest':
        train_features.replace(np.inf, 1600, inplace=True)
        train_features.replace(-np.inf, -1600, inplace=True)
        train_features.replace(np.nan, 1600, inplace=True)
        test_features.replace(np.inf, 1600, inplace=True)
        test_features.replace(-np.inf, -1600, inplace=True)
        test_features.replace(np.nan, 1600, inplace=True)
        model = RandomForestRegressor(n_jobs=12, n_estimators=100, max_leaf_nodes=31, random_state=RANDOM_SEED)
    else:
        raise Exception('Model name not supported')

    if bio_source == 'multi_task':
        model.fit(train_features, train_labels.to_numpy().ravel(), categorical_feature=['bio_source'])
    else:
        model.fit(train_features, train_labels.to_numpy().ravel())

    feature_importance = pd.DataFrame(index=train_features.columns)
    if model_name == 'lgbm' or model_name == 'randomforest':
        feature_importance[f'fi_fold_{fold}_{bio_source}'] = model.feature_importances_

    cv_results = None
    if model_name == 'lasso':
        cv_results = pd.Series({'fold': fold, 'best_alpha': model['lasso'].alpha_, 'best_l1_ratio': None})
        feature_importance[f'fi_fold_{fold}_{bio_source}'] = model['lasso'].coef_
    elif model_name == 'elasticnet':
        cv_results = pd.Series({'fold': fold, 'best_alpha': model['elasticnet'].alpha_, 'best_l1_ratio': model['elasticnet'].l1_ratio_})
        feature_importance[f'fi_fold_{fold}_{bio_source}'] = model['elasticnet'].coef_

    y_preds = model.predict(test_features)
    results = evaluate(y_preds, test_labels.to_numpy(), metrics)
    pred = pd.DataFrame({f'{bio_source}_true': test_labels[bio_source], f'{bio_source}_pred': y_preds}, index=test_labels.index, )
    total_time = time.time() - start_time
    results['total_time(sec)'] = total_time
    print("Fold: ", fold)
    print(f"\tResults: {bio_source}, {model_name}, {features_to_extract}")
    for result in results.keys():
        print(f"\t{result}: {results[result]}")

    if args.save:
        save_path = f'{args.new_results_path}/models/{bio_source}'
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        with open(f'{save_path}/model_fold_{fold}.pkl', 'xb') as save_file:
            pickle.dump(model, save_file)

    # plot lgbm tree
    # if model_name == 'lgbm':
    #     import matplotlib.pyplot as plt
    #     lgbm.plot_tree(model, figsize=(20, 20), show_info=['split_gain', 'internal_value', 'internal_count', 'leaf_count'])
    #     plt.savefig(f'./lgbm_tree_fold_{fold}_{bio_source}.png')
    #     assert False

    permutation_fi = None
    if args.permutation_fi:
        print("Running permutation feature importance...")
        permutation_fi = pd.DataFrame(index=test_features.columns)
        for metric in metrics:
            metric_name = metric.__name__
            permutation_fi[f'fold_{fold}_{bio_source}_{metric_name}'] = np.nan
        for feature in tqdm(test_features.columns):
            original_col = test_features[feature].copy()
            test_features[feature] = np.random.permutation(test_features[feature])
            perm_y_preds = model.predict(test_features)
            perm_results = evaluate(perm_y_preds, test_labels.to_numpy(), metrics)
            for metric in metrics:
                metric_name = metric.__name__
                permutation_fi.loc[feature, f'fold_{fold}_{bio_source}_{metric_name}'] = results[metric_name] - perm_results[metric_name]
            test_features[feature] = original_col


    return results, feature_importance, cv_results, pred, permutation_fi

def create_unique_directory(results_dir, name):
    base_path = os.path.join(results_dir, name)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
        return base_path
    
    index = 1
    base_path = f'{base_path}_{index}'
    while os.path.exists(base_path):
        index += 1
        base_path = f'{base_path}_{index}'
    os.makedirs(base_path)
    return base_path

def run_experiment(args, TE_path, symbol_to_fold_path, model_name, features_to_extract, metrics, name='version', data_dir='./data', results_path = 'results', bio_source='mean_te', species='human'):
    new_results_path = create_unique_directory(results_path, name)
    args.new_results_path = new_results_path

    new_data_path = f'data_with_{TE_path}'
    data = get_data(new_data_path, TE_path, symbol_to_fold_path, data_dir, species=species)

    # remove duplicates in SYMBOL from data keeping the largest tx_size
    print(f"Removing duplicates from data. Original size: {data.shape}")
    data = data.sort_values(by=['SYMBOL', 'tx_size'], ascending=False).drop_duplicates(subset=['SYMBOL'], keep='first')
    print(f"New size: {data.shape}")

    assert type(bio_source) == str or type(bio_source) == list
    if type(bio_source) == str:
        if bio_source == 'all':
            bio_source = data.columns[data.columns.str.contains('bio_source_')].to_list()
            bio_source += ['mean_te']
        elif bio_source == 'multi_task':
            bio_source = [bio_source]
        else:
            bio_source = [bio_source]

    results = pd.DataFrame(columns=['bio_source', 'fold'] + [metric.__name__ for metric in metrics] + ['total_time(sec)'])
    feature_importance = None
    permutation_fi = None
    cv_results = pd.DataFrame(columns=['bio_source', 'fold', 'best_alpha', 'best_l1_ratio'])
    predictions = pd.DataFrame(index=data['SYMBOL'])

    num_folds = 10 if bio_source != ['multi_task'] else 1
    for bs in bio_source:
        for fold in range(num_folds):
            rs, fi, cv, pred, p_fi = train_model(args, data, fold, model_name, features_to_extract, metrics, bs)
            results.loc[len(results.index)] = [bs, fold] + [rs[metric.__name__] for metric in metrics] + [rs['total_time(sec)']]
            if feature_importance is None:
                feature_importance = pd.DataFrame(index=fi.index)
            assert feature_importance.index.equals(fi.index)
            feature_importance = pd.concat([feature_importance, fi], axis=1)
            if args.permutation_fi:
                if permutation_fi is None:
                    permutation_fi = pd.DataFrame(index=p_fi.index)
                assert permutation_fi.index.equals(p_fi.index)
                permutation_fi = pd.concat([permutation_fi, p_fi], axis=1)

            for column in pred.columns:
                if column not in predictions.columns:
                    predictions[column] = pd.NA 
            predictions.update(pred) 
            
            if cv is not None:
                cv_results.loc[len(cv_results.index)] = pd.Series({'bio_source': bs, 'fold': fold, 'best_alpha': cv['best_alpha'], 'best_l1_ratio': cv['best_l1_ratio']})
        results_for_bs = results[results['bio_source'] == bs]
        results.loc[len(results.index)] = [bs, 'mean'] + [results_for_bs[metric.__name__].mean() for metric in metrics] + [results_for_bs['total_time(sec)'].mean()]
        print(results)
    
    results.to_csv(new_results_path + '/model_results.csv')
    
    feature_importance['mean_importance'] = feature_importance.mean(axis=1)
    feature_importance['rank'] = feature_importance['mean_importance'].rank(ascending=False)
    if model_name == 'lasso' or model_name == 'elasticnet':
        feature_importance['abs_mean_importance'] = feature_importance['mean_importance'].abs()
        feature_importance['rank'] = feature_importance['abs_mean_importance'].rank(ascending=False)
        cv_results.to_csv(new_results_path + '/cv_results.csv')

    feature_importance[['rank', 'mean_importance']].to_csv(new_results_path + '/feature_importance.csv')

    if args.permutation_fi:
        # for each metric get the mean of the permutation feature importance across all folds and bio_sources
        for metric in metrics:
            permutation_fi[f'mean_{metric.__name__}'] = permutation_fi.filter(like=metric.__name__, axis='columns').mean(axis=1)
        permutation_fi['rank'] = permutation_fi['mean_r2_score'].rank(ascending=False)
        permutation_fi[['rank'] + [f'mean_{metric.__name__}' for metric in metrics]].to_csv(new_results_path + '/permutation_feature_importance.csv')

    predictions.to_csv(new_results_path + '/predictions.csv')

    with open(new_results_path + '/info.txt', 'w') as f:
        f.write(f'TE_path: {TE_path}\n')
        f.write(f'symbol_to_fold_path: {symbol_to_fold_path}\n')
        f.write(f'data_path: {new_data_path}\n')
        f.write(f'model_name: {model_name}\n')
        f.write(f'bio_source: {bio_source}\n')
        f.write(f'features_to_extract: {features_to_extract}\n')
        f.write(f'metrics: {[metric.__name__ for metric in metrics]}\n')
        f.write(f'results_path: {new_results_path}\n')

def get_model_dir_name(model_name, features_to_extract):
    name = f'{model_name}-{"_".join(features_to_extract)}'
    if '1mer5_1merC_1mer3' in name:
        name = name.replace('1mer5_1merC_1mer3', '1mer')
    if '2mer5_2merC_2mer3' in name:
        name = name.replace('2mer5_2merC_2mer3', '2mer')
    if '3mer5_3merC_3mer3' in name:
        name = name.replace('3mer5_3merC_3mer3', '3mer')
    if '4mer5_4merC_4mer3' in name:
        name = name.replace('4mer5_4merC_4mer3', '4mer')
    if '5mer5_5merC_5mer3' in name:
        name = name.replace('5mer5_5merC_5mer3', '5mer')
    if '6mer5_6merC_6mer3' in name:
        name = name.replace('6mer5_6merC_6mer3', '6mer')
    if '1mer_freq_5_1mer_freq_C_1mer_freq_3' in name:
        name = name.replace('1mer_freq_5_1mer_freq_C_1mer_freq_3', '1mer_freq')
    if '2mer_freq_5_2mer_freq_C_2mer_freq_3' in name:
        name = name.replace('2mer_freq_5_2mer_freq_C_2mer_freq_3', '2mer_freq')
    if '3mer_freq_5_3mer_freq_C_3mer_freq_3' in name:
        name = name.replace('3mer_freq_5_3mer_freq_C_3mer_freq_3', '3mer_freq')
    if '4mer_freq_5_4mer_freq_C_4mer_freq_3' in name:
        name = name.replace('4mer_freq_5_4mer_freq_C_4mer_freq_3', '4mer_freq')
    if '5mer_freq_5_5mer_freq_C_5mer_freq_3' in name:
        name = name.replace('5mer_freq_5_5mer_freq_C_5mer_freq_3', '5mer_freq')
    if '6mer_freq_5_6mer_freq_C_6mer_freq_3' in name:
        name = name.replace('6mer_freq_5_6mer_freq_C_6mer_freq_3', '6mer_freq')

    return name

def create_symbol_to_fold(te_path, symbol_to_fold_path):
    if not os.path.exists(symbol_to_fold_path):
        te_data = pd.read_csv(te_path, index_col=0)
        te_data = te_data.T
        te_data.index.name = 'SYMBOL'
        te_data.index = te_data.index.str.replace('-', '.')
        symbol_to_fold = pd.DataFrame(index=te_data.index, columns=['fold'])
        symbol_to_fold = symbol_to_fold
        symbol_to_fold.index.name = 'SYMBOL'
        kf = KFold(n_splits=10, shuffle=True, random_state=RANDOM_SEED)
        for fold, (_, test_index) in enumerate(kf.split(te_data)):
            symbol_to_fold.iloc[test_index, symbol_to_fold.columns.get_loc('fold')] = fold
        print(symbol_to_fold['fold'].value_counts())
        symbol_to_fold.to_csv(symbol_to_fold_path, sep='\t')


if __name__ == '__main__':
    # seed everything
    np.random.seed(RANDOM_SEED)

    tqdm.pandas()
    parser = argparse.ArgumentParser()
    parser.add_argument('--permutation_fi', '-p', action='store_true', help='Whether to run permutation feature importance')
    parser.add_argument('--save', '-s', action='store_true', help='Whether to save models')
    args = parser.parse_args()

    TE_paths = []
    # TE_paths.append('human_TE_cellline_all_plain.csv')
    TE_paths.append('human_TE_cellline_all_NA_plain.csv')
    # TE_paths.append('human_TE_cellline_all_dedup.csv')
    symbol_to_fold_path = 'symbol_to_fold.tsv'

    model_names = [
        'lgbm', 
        # 'lasso', 
        # 'elasticnet', 
        # 'randomforest'
    ]

    feature_tests = [
        # [ # best so far
        #     'LL', 'P5', 'P3', 'CF', 'AAF',
        #     '3mer_freq_5',
        #     'Struct'
        # ],
        # [ # best so far, no struct
        #     'LL', 'P5', 'P3', 'CF', 'AAF',
        #     '3mer_freq_5',
        # ],
        # ['3mer_freq_5'],
        # ['CF'],
        # ['AAF'],
        # ['CF', 'AAF', '3mer_freq_5'],
        # [ # best so far, biochem
        #     'LL', 'P5', 'P3', 'CF', 'AAF',
        #     '3mer_freq_5',
        #     'Struct', 'Biochem'
        # ],
        # ['CF', '3mer_freq_5'],
        # ['AAF', '3mer_freq_5'],
        # ['LL','CF', 'AAF', '3mer_freq_5'],
        # ['P5', 'P3','CF', 'AAF', '3mer_freq_5'],
        # ['LL', 'P5', 'P3', 'CF','3mer_freq_5',],
        # ['CF','AAF'],
        # ['LL','CF', '3mer_freq_5'],
        # ['LL','CF', '3mer_freq_5', 'Struct'],
        # ['CF', 'Struct'],
        # ['LL','CF', 'AAF', '3mer_freq_5', 'Struct']
        # ['Biochem'],
        # [ # all seq features
        #     'LL', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'CF', 'AAF', 'DC', '1mer_freq_5', 
        #     '1mer_freq_C', '1mer_freq_3', '2mer_freq_5', '2mer_freq_C', '2mer_freq_3',
        #     '3mer_freq_5', '3mer_freq_C', '3mer_freq_3', '4mer_freq_5', '4mer_freq_C', '4mer_freq_3',
        #     '5mer_freq_5', '5mer_freq_C', '5mer_freq_3', '6mer_freq_5', '6mer_freq_C', '6mer_freq_3', 'Struct'
        # ],
        [ # all seq features
            'LL', 'P5', 'P3','CF', 'AAF', 'DCF', '3mer_freq_5',
        ],
    ]

    # feature_tests = [
    #     # [ # all features
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC', '1mer5', 
    #     #     '1merC', '1mer3', '2mer5', '2merC', '2mer3',
    #     #     '3mer5', '3merC', '3mer3', '4mer5', '4merC', '4mer3',
    #     #     '5mer5', '5merC', '5mer3', '6mer5', '6merC', '6mer3', 'Struct'
    #     # ],
    #     # [ # no 6mers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC', '1mer5', 
    #     #     '1merC', '1mer3', '2mer5', '2merC', '2mer3',
    #     #     '3mer5', '3merC', '3mer3', '4mer5', '4merC', '4mer3',
    #     #     '5mer5', '5merC', '5mer3', 'Struct'
    #     # ],
    #     # [ # no 5mers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC', '1mer5', 
    #     #     '1merC', '1mer3', '2mer5', '2merC', '2mer3',
    #     #     '3mer5', '3merC', '3mer3', '4mer5', '4merC', '4mer3',
    #     #     'Struct'
    #     # ],
    #     # [ # only 1mers, 2mers, 3mers, no 4mers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC', '1mer5', 
    #     #     '1merC', '1mer3', '2mer5', '2merC', '2mer3',
    #     #     '3mer5', '3merC', '3mer3',
    #     #     'Struct'
    #     # ],

    #     # [ # only 1mers and 2mers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC', '1mer5', 
    #     #     '1merC', '1mer3', '2mer5', '2merC', '2mer3',
    #     #     'Struct'
    #     # ],
    #     # [ # only 1mers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC', '1mer5', 
    #     #     '1merC', '1mer3',
    #     #     'Struct'
    #     # ],
    #     # [ # only 2mers and 3mers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC', 
    #     #     '2mer5', '2merC', '2mer3', '3mer5', '3merC', '3mer3',
    #     #     'Struct'
    #     # ],
    #     # [ # only 2mers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC', 
    #     #     '2mer5', '2merC', '2mer3',
    #     #     'Struct'
    #     # ],
    #     # [ # only 1mers and 3mers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC', 
    #     #     '1mer5', '1merC', '1mer3', '3mer5', '3merC', '3mer3',
    #     #     'Struct'
    #     # ],
    #     # [ # only 3mers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC', 
    #     #     '3mer5', '3merC', '3mer3',
    #     #     'Struct'
    #     # ],
    #     # [ # no 1mers, or any kmers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 'DC',
    #     #     'Struct'
    #     # ],

    #     # [ # no DC, kmers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'C', 
    #     #     'Struct'
    #     # ],
    #     # [ # no K, DC, kmers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'C', 
    #     #     'Struct'
    #     # ],

    #     # [ # try dropping diff regions of 3mers
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'C', 
    #     #     '3mer5', '3merC', '3mer3',
    #     #     'Struct'
    #     # ],

    #     # [ # 3merC and 3mer3
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'C', 
    #     #     '3merC', '3mer3',
    #     #     'Struct'
    #     # ],
    #     # [ # 3mer5 and 3mer3
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'C', 
    #     #     '3mer5', '3mer3',
    #     #     'Struct'
    #     # ],
    #     # [ # Best so far, 3mer5 and 3merC
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # only 3mer5
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'C', 
    #     #     '3mer5',
    #     #     'Struct'
    #     # ],
    #     # [ # only 3merC
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'C', 
    #     #     '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # only 3mer3
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 'C', 
    #     #     '3mer3',
    #     #     'Struct'
    #     # ],
    #     # [ # no C, C is important
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],

    #     # [ # no WP
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # no P3
    #     #     'L', 'P', 'P5', 'PC', 'WP', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # no PC
    #     #     'L', 'P', 'P5', 'P3', 'WP', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # no P5
    #     #     'L', 'P', 'PC', 'P3', 'WP', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # no P, best so far
    #     #     'L', 'P5', 'PC', 'P3', 'WP', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # P5, PC, P3
    #     #     'L', 'P5', 'PC', 'P3', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # PC, P3
    #     #     'L', 'PC', 'P3', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # P5, P3, best so far <=========================================
    #     #     'L', 'P5', 'P3', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # P5, PC
    #     #     'L', 'P5', 'PC', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # P5
    #     #     'L', 'P5', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # PC
    #     #     'L', 'PC', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # P3
    #     #     'L', 'P3', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # no P*
    #     #     'L', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],

    #     # [ # no L
    #     #     'P5', 'P3', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # no Struct
    #     #     'L', 'P5', 'P3', 'C', 
    #     #     '3mer5', '3merC',
    #     # ],

    #     # [ # only C, 3mer5, 3merC, Struct
    #     #     'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # only L, C, 3mer5, 3merC
    #     #     'L', 'C', 
    #     #     '3mer5', '3merC',
    #     # ],
    #     # [ # only C, 3mer5, 3merC
    #     #     'C', 
    #     #     '3mer5', '3merC',
    #     # ],

    #     # [ # from best so far, no L, already ran above
    #     #     'P5', 'P3', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # from best so far, no P5
    #     #     'L', 'P3', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # from best so far, no P3
    #     #     'L', 'P5', 'C', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # from best so far, no C
    #     #     'L', 'P5', 'P3', 
    #     #     '3mer5', '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # from best so far, no 3mer5
    #     #     'L', 'P5', 'P3', 'C', 
    #     #     '3merC',
    #     #     'Struct'
    #     # ],
    #     # [ # from best so far, no 3merC
    #     #     'L', 'P5', 'P3', 'C', 
    #     #     '3mer5',
    #     #     'Struct'
    #     # ],
    #     # [ # from best so far, no Struct
    #     #     'L', 'P5', 'P3', 'C', 
    #     #     '3mer5', '3merC',
    #     # ],
    #     # [
    #     #     'L', 'P', 'P5', 'PC', 'P3', 'WP',
    #     #     'Struct'
    #     # ]
    # ]

    # metrics = [r2_score, pearson_corrcoef, spearman_corrcoef, mean_squared_error]
    metrics = [nan_r2, masked_mse_loss] # nan_pearson_corr

    already_ran = []

    for TE_path in TE_paths:
        for model_name in model_names:
            for features_to_extract in feature_tests:
                name = get_model_dir_name(model_name, features_to_extract)
                if 'NA' in TE_path:
                    name = 'on_NA_' + name
                if name in already_ran:
                    print(f'{name} already ran, skipping...')
                else:
                    already_ran.append(name)
                    run_experiment(
                        args=args,
                        TE_path=TE_path, 
                        symbol_to_fold_path=symbol_to_fold_path, 
                        model_name=model_name, 
                        features_to_extract=features_to_extract, 
                        metrics=metrics,
                        name=name,
                        results_path='results/feature_set_comparison/',
                        bio_source='mean_te'
                        )    

    print("Done!")



