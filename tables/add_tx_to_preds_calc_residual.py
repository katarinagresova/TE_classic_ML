import pandas as pd
import tqdm
import os

def txs_and_residuals(tx_path, pred_path):
    tx = pd.read_csv(tx_path, sep='\t')
    tx = tx.sort_values(by=['SYMBOL', 'tx_size'], ascending=False).drop_duplicates(subset=['SYMBOL'], keep='first')

    pred = pd.read_csv(pred_path)

    tx = tx[['SYMBOL', 'gene_id', 'transcript_id', 'tx_size','utr5_size','utr3_size','cds_size','tx_sequence', 'fold']]

    assert len(tx) == len(pred), f'{len(tx)} != {len(pred)}'

    merged = pd.merge(tx, pred, on='SYMBOL', how='right')

    merged['mean_te_residual'] = merged['mean_te_true'] - merged['mean_te_pred']
    merged.drop(columns=['mean_te_true', 'mean_te_pred'], inplace=True)

    return merged
    

if __name__ == '__main__':
    tx_path = './data/data_with_human_TE_cellline_all_NA_plain_vst.csv'
    pred_path = './results/vst_copositional/with_NA_lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/predictions.csv'
    with_NA = txs_and_residuals(tx_path, pred_path)

    tx_path = './data/data_with_human_TE_cellline_all_plain_vst.csv'
    pred_path = './results/vst_copositional/without_NA_lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/predictions.csv'
    without_NA = txs_and_residuals(tx_path, pred_path)

    # write to excel
    out_path = './results/vst_copositional/residuals.xlsx'
    with pd.ExcelWriter(out_path) as writer:
        with_NA.to_excel(writer, sheet_name='with_NA', index=False)
        without_NA.to_excel(writer, sheet_name='without_NA', index=False)