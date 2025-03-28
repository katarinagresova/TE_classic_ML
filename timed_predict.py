import pickle
import pandas as pd
import os
from lgbm_feature_extract_from_str import dataframe_feature_extract #, calc_struct_data
import yaml
from tqdm import tqdm
from train import nan_r2
import numpy as np
import argparse

import time

class LGBM_TE_model:
    def __init__(self, models_dir: str):
        print(f'Loading models from {models_dir}')
        self.models = {}
        info_file = os.path.join(models_dir, 'info.txt')
        # open info file as yaml
        with open(info_file) as f:
            self.info = yaml.load(f, Loader=yaml.FullLoader)
        self.features_to_extract = self.info['features_to_extract']

        for bio_source_dir in tqdm(os.listdir(os.path.join(models_dir, 'models'))):
            bio_source = os.path.basename(bio_source_dir)
            self.models[bio_source] = {}
            for file in os.listdir(os.path.join(models_dir, 'models', bio_source_dir)):
                model = pickle.load(open(os.path.join(models_dir, 'models', bio_source_dir, file), 'rb'))
                model_fold = int(file.split('_')[2].removesuffix('.pkl'))
                self.models[bio_source][model_fold] = model


    def predict_TE(self, data: pd.DataFrame) -> pd.DataFrame:
        # if "Struct" in self.features_to_extract:
        #     data = calc_struct_data(data)
        number_of_cell_lines = 0
        for bio_source, models in tqdm(self.models.items()):
            number_of_cell_lines += 1
            extracted_features, _ = dataframe_feature_extract(data, self.features_to_extract, te_source='tx_size')
            if "transcript_id" in extracted_features.columns:
                extracted_features = extracted_features.drop(columns=['transcript_id'])
            if "gene_id" in extracted_features.columns:
                extracted_features = extracted_features.drop(columns=['gene_id'])
            pred_name = f'{bio_source}_pred'
            assert pred_name not in data.columns
            all_preds = pd.DataFrame()
            for model_fold, model in models.items():
                all_preds[f"fold_{model_fold}"] = model.predict(extracted_features)
            data[pred_name] = all_preds.mean(axis=1)

        bio_sources = list(self.models.keys())
        cols = ["SYMBOL"]
        if "transcript_id" in data.columns:
            cols.append("transcript_id")
        if "gene_id" in data.columns:
            cols.append("gene_id")
        cols.extend(["tx_size", "utr5_size", "cds_size", "utr3_size", "tx_sequence"])
        for bio_source in bio_sources:
            if f"{bio_source}_true" in data.columns:
                cols.append(f"{bio_source}_true")
            cols.append(f"{bio_source}_pred")

        if "mean_te" in bio_sources:
            bio_sources.remove("mean_te")

        contains_true = True
        for bio_source in bio_sources:
            contains_true = contains_true and f"{bio_source}_true" in data.columns
        if contains_true:
            cols.append("mean_across_cell_lines_true")
            data["mean_across_cell_lines_true"] = data[[f"{bio_source}_true" for bio_source in bio_sources]].mean(axis=1)
        cols.append("mean_across_cell_lines_pred")
        data["mean_across_cell_lines_pred"] = data[[f"{bio_source}_pred" for bio_source in bio_sources]].mean(axis=1)

        data = data[cols]
        return data, number_of_cell_lines
    
    def predict_TE_single_fold(self, data: pd.DataFrame, fold: int) -> pd.DataFrame:
        # if "Struct" in self.features_to_extract:
        #     data = calc_struct_data(data)
        
        for bio_source, models in tqdm(self.models.items()):
            extracted_features, _ = dataframe_feature_extract(data, self.features_to_extract, te_source='tx_size')
            pred_name = f'{bio_source}_pred'
            assert pred_name not in data.columns
            model = models[fold]
            data[pred_name] = model.predict(extracted_features)

        bio_sources = list(self.models.keys())
        cols = ["SYMBOL", "tx_size", "utr5_size", "cds_size", "utr3_size", "tx_sequence"]
        for bio_source in bio_sources:
            if f"{bio_source}_true" in data.columns:
                cols.append(f"{bio_source}_true")
            cols.append(f"{bio_source}_pred")

        if "mean_te" in bio_sources:
            bio_sources.remove("mean_te")

        contains_true = True
        for bio_source in bio_sources:
            contains_true = contains_true and f"{bio_source}_true" in data.columns
        if contains_true:
            cols.append("mean_across_cell_lines_true")
            data["mean_across_cell_lines_true"] = data[[f"{bio_source}_true" for bio_source in bio_sources]].mean(axis=1)
        cols.append("mean_across_cell_lines_pred")
        data["mean_across_cell_lines_pred"] = data[[f"{bio_source}_pred" for bio_source in bio_sources]].mean(axis=1)


        data = data[cols]
        return data


    

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--model_dir', type=str, required=True, help='Directory containing the models, ex: ./results/human/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5')
    argparser.add_argument('--data_path', type=str, required=True, help='Path to the data to predict on, ex: ./examples/predict_input_example.csv')
    argparser.add_argument('--output_path', type=str, required=True, help='Path to save the output, ex: ./examples/predict_output_example.csv')
    args = argparser.parse_args()
   
    start_time = time.time()
    model = LGBM_TE_model(args.model_dir)
    data = pd.read_csv(args.data_path, sep='\t')
    data, number_of_cell_lines = model.predict_TE(data)
    data.to_csv(args.output_path, sep='\t', index=False)
    end_time = time.time()
    total_time = end_time - start_time
    print("Total time taken (sec): ", total_time)
    print("Average time taken per cell line (sec): ", total_time/number_of_cell_lines)



