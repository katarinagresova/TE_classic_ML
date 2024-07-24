import pandas as pd
from data import get_data
from lgbm_feature_extract_from_str import dataframe_feature_extract


if __name__ == '__main__':
    TE_path = 'mouse_TE_cellline_all_plain_NA.csv'
    symbol_to_fold_path = 'mouse_symbol_to_fold.tsv'
    new_data_path = f'data_with_{TE_path}'

    features_to_extract = [
        'LL', 'P5', 'P3', 'CF', 'AAF', 
        '3mer_freq_5',
    ]

    data = get_data(new_data_path, TE_path, symbol_to_fold_path, species='mouse')
    all_features, _ = dataframe_feature_extract(data, features_to_extract)
    # bio_sources = data.columns[data.columns.str.contains('bio_source_')].tolist() + ['mean_te']
    # bio_sources = data[bio_sources]
    # bio_sources.index = data['SYMBOL']
    preds = pd.read_csv('./results/mouse/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/predictions.csv', index_col=0)
    true_vals = preds.drop(preds.columns[preds.columns.str.contains('pred')], axis=1)
    all_features = pd.merge(all_features, true_vals, left_index=True, right_index=True)
    
    features = "_".join(features_to_extract)
    all_features.to_csv(f'./graphing/mouse/D/features-{features}.csv')