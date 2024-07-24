import pandas as pd
from data import get_data
from lgbm_feature_extract_from_str import dataframe_feature_extract
from tqdm import tqdm
from Bio import SeqIO

HUMAN_FA_FILE_PATH = './data/appris_human_v2_selected.fa'
MOUSE_FA_FILE_PATH = './data/appris_mouse_v2_selected.fa'


def extract_lengths(name_list):
    utr5_size, cds_size, utr3_size = None, None, None

    tx_size = int(name_list[6])

    for section in name_list:
        if section.startswith("UTR5:"):
            utr5_start, utr5_end = map(int, section[len("UTR5:"):].split("-"))
            utr5_size = utr5_end - utr5_start + 1
        elif section.startswith("CDS:"):
            cds_start, cds_end = map(int, section[len("CDS:"):].split("-"))
            cds_size = cds_end - cds_start + 1
        elif section.startswith("UTR3:"):
            utr3_start, utr3_end = map(int, section[len("UTR3:"):].split("-"))
            utr3_size = utr3_end - utr3_start + 1

    if utr5_size is None:
        utr5_size = 0
    if utr3_size is None:
        utr3_size = 0
    assert utr5_size is not None and cds_size is not None and utr3_size is not None, f'name_list: {name_list}, utr5_size: {utr5_size}, cds_size: {cds_size}, utr3_size: {utr3_size}'
    assert utr5_size + cds_size + utr3_size == tx_size

    return tx_size, utr5_size, cds_size, utr3_size


def get_genes_from_appris(species='human') -> pd.DataFrame:
    print('Get genes form appris')
    new_data_list = []

    fa_file_path = HUMAN_FA_FILE_PATH if species == 'human' else MOUSE_FA_FILE_PATH

    # missing_list = []

    with open(fa_file_path) as fa_file:
        for record in tqdm(SeqIO.parse(fa_file, 'fasta')):
            name_list = record.name.split('|')
            symbol = name_list[5]
            symbol = symbol.replace('-', '.')
            transcript_id = name_list[0]
            gene_id = name_list[1]
            if symbol[0].isdigit():
                symbol = 'X' + symbol
            tx_size, utr5_size, cds_size, utr3_size = extract_lengths(name_list)
            new_row = {
                'SYMBOL': symbol, 
                'transcript_id': transcript_id,
                'gene_id': gene_id,
                'tx_size': tx_size,
                'utr5_size': utr5_size, 
                'cds_size': cds_size,
                'utr3_size': utr3_size,
                'tx_sequence': str(record.seq),
            }

            new_data_list.append(new_row)
    #         else:
    #             missing_list.append(symbol)

    # with open('missing_genes_in_appris.tsv', 'w') as f:
    #     for item in missing_list:
    #         f.write("%s\n" % item)
            

    new_data = pd.DataFrame(new_data_list)
    # find duplicates in SYMBOL column
    print('Duplicates: \n', new_data[new_data.duplicated(subset=['SYMBOL'], keep=False)])
    print('NaNs getting seqs: ', new_data.isna().sum().sum())

    print('new_data shape after search: ', new_data.shape)
    # There are some missing genes, also need to handle genes with no UTR5
    # struct_utr5_info = pd.read_csv(STRUCT_UTR5_FILE, sep='\t', index_col=0)
    # struct_cds_info = pd.read_csv(STRUCT_CDS_FILE, sep='\t', index_col=0)
    # struct_info = pd.merge(struct_utr5_info, struct_cds_info, left_index=True, right_index=True, how='outer')
    # new_data = pd.merge(new_data, struct_info, left_on='SYMBOL', right_index=True, how='left')

    return new_data


if __name__ == '__main__':
    # TE_path = 'human_TE_cellline_all_NA_plain.csv'
    # symbol_to_fold_path = 'symbol_to_fold.tsv'
    # new_data_path = f'data_with_{TE_path}'

    features_to_extract = [ # all seq features
        'LL', 'P', 'P5', 'PC', 'P3', 'WP', 'K', 'CF', 'AAF', 'DC', '1mer_freq_5', 
        '1mer_freq_C', '1mer_freq_3', '2mer_freq_5', '2mer_freq_C', '2mer_freq_3',
        '3mer_freq_5', '3mer_freq_C', '3mer_freq_3', #'Struct'

        # '4mer_freq_5', '4mer_freq_C', '4mer_freq_3',

        # '5mer_freq_5', '5mer_freq_C', '5mer_freq_3', 

        # '6mer_freq_5', '6mer_freq_C', '6mer_freq_3', 
    ]

    # data = get_data(new_data_path, TE_path, symbol_to_fold_path)
    # data = pd.read_csv("./protein_dataset/pan_cancer_949_alternate_trx/pan_cancer_949_alternate_trx.csv", sep='\t')
    data = get_genes_from_appris()
    data['mean_te'] = 0
    all_features, _ = dataframe_feature_extract(data, features_to_extract)
   
    features = "all_other_features_except_struct"
    all_features.to_csv(f'./protein_dataset/appris_features/{features}.csv')

    features_to_extract = [ # all seq features
        '4mer_freq_5', '4mer_freq_C', '4mer_freq_3',
    ]
    all_features, _ = dataframe_feature_extract(data, features_to_extract)
    features = "4mers"
    all_features.to_csv(f'./protein_dataset/appris_features/{features}.csv')


    features_to_extract = [ # all seq features
        '5mer_freq_5', '5mer_freq_C', '5mer_freq_3', 
    ]
    all_features, _ = dataframe_feature_extract(data, features_to_extract)
    features = "5mers"
    all_features.to_csv(f'./protein_dataset/appris_features/{features}.csv')


    features_to_extract = [ # all seq features
        '6mer_freq_5', '6mer_freq_C', '6mer_freq_3', 
    ]
    all_features, _ = dataframe_feature_extract(data, features_to_extract)
    features = "6mers"
    all_features.to_csv(f'./protein_dataset/appris_features/{features}.csv')