# unpickle LGBM model and predict TE
import pickle
import pandas as pd
import numpy as np
import os
from lgbm_feature_extract_from_str import dataframe_feature_extract
import yaml
from tqdm import tqdm
from gene_optimize_replace import get_gene_to_optimize
from typing import List

class LGBM_TE_model:
    def __init__(self, models_dir: str):
        self.models = {}
        self.bio_source = None
        files = os.listdir(os.path.join(models_dir, 'models'))
        for file in files:
            model = pickle.load(open(os.path.join(models_dir, 'models', file), 'rb'))
            model_fold = int(file.split('_')[2])
            self.models[model_fold] = model
            model_bio_source = file.split('_')[5].replace('.pkl', '')
            if self.bio_source is None:
                self.bio_source = model_bio_source
            else:
                assert self.bio_source == model_bio_source
        info_file = os.path.join(models_dir, 'info.txt')
        # open info file as yaml
        with open(info_file) as f:
            self.info = yaml.load(f, Loader=yaml.FullLoader)
        self.features_to_extract = self.info['features_to_extract']


    def predict_TE(self, data: pd.DataFrame, pred_name=None) -> pd.DataFrame:
        if pred_name is None:
            pred_name = f'pred_{self.bio_source}'
        assert pred_name not in data.columns
        data[pred_name] = None
        for model_fold, model in tqdm(self.models.items()):
            fold_data = data[data['fold'] == model_fold]
            extracted_features, _ = dataframe_feature_extract(fold_data, self.features_to_extract, te_source='tx_size')
            preds = model.predict(extracted_features)
            data.loc[fold_data.index, pred_name] = preds
            # print(f'fold {model_fold}')
            # print(data.loc[fold_data.index].head())
        return data
    
    def predict_TE_with_single_fold(self, data: pd.DataFrame, fold: int, pred_name=None) -> pd.DataFrame:
        if pred_name is None:
            pred_name = f'pred_{self.bio_source}'
        assert pred_name not in data.columns
        data[pred_name] = None
        model = self.models[fold]
        extracted_features, _ = dataframe_feature_extract(data, self.features_to_extract, te_source='tx_size')
        preds = model.predict(extracted_features)
        data.loc[data.index, pred_name] = preds
        return data



    

if __name__ == '__main__':
    model = LGBM_TE_model('./results/For_Gene_Optimize/bio_source_HEK293T/on_NA_lgbm-LL_P5_P3_CF_AAF_3mer_freq_5')

    # Get Best UTRs
    # data = pd.read_csv('./gene_optimize/data/eGFP_original_5utr_original_3utr.csv')
    # data = model.predict_TE(data)
    # data.to_csv('./gene_optimize/data/eGFP_original_5utr_original_3utr_pred.csv')

    # data = pd.read_csv('./gene_optimize/data/eGFP_original_5utr_CWF19L2_3utr.csv')
    # data = model.predict_TE(data)
    # data.to_csv('./gene_optimize/data/eGFP_original_5utr_CWF19L2_3utr_pred.csv')
    
    

