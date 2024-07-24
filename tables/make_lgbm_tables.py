import pandas as pd
import tqdm

if __name__ == '__main__':
    tqdm.tqdm.pandas()

    # Load the data
    tableS2 = pd.read_excel('./tables/Table_S2 copy.xlsx', sheet_name=None)

    # assert len(tableS2.keys()) == 3
    # assert 'Multi-task (Human)' in tableS2.keys()
    # assert 'Multi-task (Mouse)' in tableS2.keys()
    # assert 'Single-task Transfer Learning' in tableS2.keys()
    mthuman = tableS2['Multi-task (Human)']
    mtmouse = tableS2['Multi-task (Mouse)']
    sttl = tableS2['Single-task Transfer Learning']

    # assert mthuman['SYMBOL'].equals(sttl['SYMBOL'])

    lgbm_human_preds = pd.read_csv('./results/human/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/predictions.csv')
    lgbm_human_preds = lgbm_human_preds[lgbm_human_preds.columns.drop(list(lgbm_human_preds.filter(regex='_true')))]
    lgbm_human_preds.columns = lgbm_human_preds.columns.str.replace('bio_source', 'predicted_TE')
    lgbm_human_preds.columns = lgbm_human_preds.columns.str.replace('_pred', '')
    lgbm_human_preds.columns = lgbm_human_preds.columns.str.replace(' ', '_', regex=False)
    print("lgbm_human_preds shape: ", lgbm_human_preds.shape)


    human_data = pd.read_csv('./data/data_with_human_TE_cellline_all_NA_plain.csv', sep='\t')
    human_data = human_data[['SYMBOL', 'gene_id', 'transcript_id', 'tx_size','utr5_size','utr3_size','cds_size','tx_sequence', 'fold']]
    human_data = human_data.sort_values(by=['SYMBOL', 'tx_size'], ascending=False).drop_duplicates(subset=['SYMBOL'], keep='first')
    human_data = human_data.rename(columns={'transcript_id': 'tx_id'})
    print("human_data shape: ", human_data.shape)

    lgbm_mouse_preds = pd.read_csv('./results/mouse/all_cell_lines/lgbm-LL_P5_P3_CF_AAF_3mer_freq_5/predictions.csv')
    lgbm_mouse_preds = lgbm_mouse_preds[lgbm_mouse_preds.columns.drop(list(lgbm_mouse_preds.filter(regex='_true')))]
    lgbm_mouse_preds.columns = lgbm_mouse_preds.columns.str.replace('bio_source', 'predicted_TE')
    lgbm_mouse_preds.columns = lgbm_mouse_preds.columns.str.replace('_pred', '')
    lgbm_mouse_preds.columns = lgbm_mouse_preds.columns.str.strip()
    lgbm_mouse_preds.columns = lgbm_mouse_preds.columns.str.replace(' ', '_', regex=False)
    print("lgbm_mouse_preds shape: ", lgbm_mouse_preds.shape)

    mouse_data = pd.read_csv('./data/data_with_mouse_TE_cellline_all_plain_NA.csv', sep='\t')
    mouse_data = mouse_data[['SYMBOL', 'gene_id', 'transcript_id', 'tx_size','utr5_size','utr3_size','cds_size','tx_sequence', 'fold']]
    mouse_data = mouse_data.sort_values(by=['SYMBOL', 'tx_size'], ascending=False).drop_duplicates(subset=['SYMBOL'], keep='first')
    mouse_data = mouse_data.rename(columns={'transcript_id': 'tx_id'})
    print("mouse_data shape: ", mouse_data.shape)

    # assert len(lgbm_human_preds) == len(human_data), f'{len(lgbm_human_preds)} != {len(human_data)}'
    # assert len(lgbm_mouse_preds) == len(mouse_data), f'{len(lgbm_mouse_preds)} != {len(mouse_data)}'

    # assert lgbm_human_preds['SYMBOL'].isin(human_data['SYMBOL']).all()
    # assert lgbm_mouse_preds['SYMBOL'].isin(mouse_data['SYMBOL']).all()
    # assert mthuman['SYMBOL'].isin(lgbm_human_preds['SYMBOL']).all()
    # assert mthuman['SYMBOL'].isin(human_data['SYMBOL']).all()
    # assert mtmouse['SYMBOL'].isin(lgbm_mouse_preds['SYMBOL']).all()
    # assert mtmouse['SYMBOL'].isin(mouse_data['SYMBOL']).all()

    lgbm_table_human = []
    added_symbols = set()
    for i, row in tqdm.tqdm(mthuman.iterrows()):
        symbol = row['SYMBOL']
        preds = lgbm_human_preds.loc[lgbm_human_preds['SYMBOL'] == symbol]
        data = human_data.loc[human_data['SYMBOL'] == symbol]
        assert len(preds) == 1
        assert len(data) == 1
        new_row = {}
        for col in mthuman.columns:
            if col == 'SYMBOL':
                new_row[col] = symbol
            elif col in preds.columns:
                new_row[col] = preds[col].values[0]
            elif col in data.columns:
                new_row[col] = data[col].values[0]
            else:
                raise ValueError(f'Column {col} not found in preds or data')
        
        added_symbols.add(symbol)
        lgbm_table_human.append(new_row)
    
    unadded_symbols = set(lgbm_human_preds['SYMBOL']) - added_symbols
    for symbol in tqdm.tqdm(unadded_symbols):
        preds = lgbm_human_preds.loc[lgbm_human_preds['SYMBOL'] == symbol]
        data = human_data.loc[human_data['SYMBOL'] == symbol]
        assert len(preds) == 1
        assert len(data) == 1
        new_row = {}
        for col in mthuman.columns:
            if col == 'SYMBOL':
                new_row[col] = symbol
            elif col in preds.columns:
                new_row[col] = preds[col].values[0]
            elif col in data.columns:
                new_row[col] = data[col].values[0]
            else:
                raise ValueError(f'Column {col} not found in preds or data')
        
        lgbm_table_human.append(new_row)
    
    lgbm_table_human = pd.DataFrame(lgbm_table_human)

    lgbm_table_mouse = []
    added_symbols = set()
    for i, row in tqdm.tqdm(mtmouse.iterrows()):
        symbol = row['SYMBOL']
        preds = lgbm_mouse_preds.loc[lgbm_mouse_preds['SYMBOL'] == symbol]
        data = mouse_data.loc[mouse_data['SYMBOL'] == symbol]
        assert len(preds) == 1
        assert len(data) == 1
        new_row = {}
        for col in mtmouse.columns:
            if col == 'SYMBOL':
                new_row[col] = symbol
            elif col in preds.columns:
                new_row[col] = preds[col].values[0]
            elif col in data.columns:
                new_row[col] = data[col].values[0]
            else:
                raise ValueError(f'Column {col} not found in preds or data')
        
        added_symbols.add(symbol)
        lgbm_table_mouse.append(new_row)
    
    unadded_symbols = set(lgbm_mouse_preds['SYMBOL']) - added_symbols
    for symbol in tqdm.tqdm(unadded_symbols):
        preds = lgbm_mouse_preds.loc[lgbm_mouse_preds['SYMBOL'] == symbol]
        data = mouse_data.loc[mouse_data['SYMBOL'] == symbol]
        assert len(preds) == 1
        assert len(data) == 1
        new_row = {}
        for col in mtmouse.columns:
            if col == 'SYMBOL':
                new_row[col] = symbol
            elif col in preds.columns:
                new_row[col] = preds[col].values[0]
            elif col in data.columns:
                new_row[col] = data[col].values[0]
            else:
                raise ValueError(f'Column {col} not found in preds or data')
        
        lgbm_table_mouse.append(new_row)

    lgbm_table_mouse = pd.DataFrame(lgbm_table_mouse)

    
    tableS2 = pd.read_excel('./tables/Table_S2 copy.xlsx', sheet_name=None)
    tableS2['LGBM (Human)'] = lgbm_table_human
    tableS2['LGBM (Mouse)'] = lgbm_table_mouse

    with pd.ExcelWriter('./tables/Table_S2 copy.xlsx') as writer:
        for sheet_name, df in tableS2.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print("Done!")



    