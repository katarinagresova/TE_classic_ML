import pandas as pd
import tqdm

def main(species, model):

    true_te = pd.read_excel('./tables/Table_S1.xlsx', sheet_name=species)
    true_te = true_te.drop(columns=['fold'])
    pred_te = pd.read_excel('./tables/Table_S2_with_LGBM.xlsx', sheet_name=f"{model} ({species})")
    pred_te = pred_te.drop(columns=['gene_id','SYMBOL','tx_size','utr5_size','utr3_size','cds_size','tx_sequence','fold'])

    gene_info_cols = ['SYMBOL','gene_id','tx_id','tx_size','utr5_size','utr3_size','cds_size','tx_sequence']
    cell_lines = true_te.columns.str.contains('TE')
    cell_lines = true_te.columns[cell_lines]
    residuals = pd.DataFrame(columns=gene_info_cols)
    #true - pred
    all = pd.merge(true_te, pred_te, on='tx_id', suffixes=('_true', '_pred'), how='inner')
    #print dupes in true_te
    # assert all.shape[0] == pred_te.shape[0], f"Error: {all.shape} != {pred_te.shape}"
    for c in gene_info_cols:
        residuals[c] = all[c]
    for c in cell_lines:
        residuals[f'residual_{c}'] = all[c] - all[f'predicted_{c}']
    # assert residuals.shape[0] == pred_te.shape[0], f"Error: {residuals.shape[0]} != {pred_te.shape[0]}"

    mean_te_residual = pd.DataFrame(columns=gene_info_cols)
    for c in gene_info_cols:
        mean_te_residual[c] = all[c]
    pred_cols = [c for c in all.columns if 'predicted_' in c]
    all[f'predicted_mean_te'] = all[pred_cols].mean(axis=1)
    mean_te_residual['mean_te_residual'] = all['mean_te'] - all['predicted_mean_te']

    print("Correlation between true and predicted mean TE:", all['mean_te'].corr(all['predicted_mean_te']))
    
    sheets = {'residuals': residuals, 'mean_te_residual': mean_te_residual}
    with pd.ExcelWriter(f'./tables/Residuals_{model}_{species}.xlsx') as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)



if __name__ == '__main__':
    # Human or Mouse
    # Multi-task or LGBM
    main('Human', 'LGBM')