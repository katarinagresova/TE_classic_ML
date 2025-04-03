import pandas as pd
import os
from tqdm import tqdm
import time
from typing import List


def get_true_data(true_data_path: str):
    return pd.read_csv(true_data_path, sep='\t', index_col=None)

def get_gene_to_optimize(gene_path: str):
    df = pd.read_csv(gene_path, index_col=None)
    series = df.iloc[0]
    return series

def extract_regions(sequence, tx_size, utr5_size, cds_size, utr3_size):
    assert tx_size == utr5_size + cds_size + utr3_size
    assert len(sequence) == tx_size
    utr5 = sequence[:utr5_size]
    cds = sequence[utr5_size:utr5_size+cds_size]
    utr3 = sequence[utr5_size+cds_size:utr5_size+cds_size+utr3_size]
    assert len(utr5) == utr5_size
    assert len(cds) == cds_size
    assert len(utr3) == utr3_size
    assert utr5 + cds + utr3 == sequence
    assert len(utr5) + len(cds) + len(utr3) == tx_size
    return utr5, cds, utr3

def replace_utrs(row: pd.Series, gene_to_optimize: pd.Series, utr5_seq, utr3_seq):
    original_utr5, _, original_utr3 = extract_regions(row['tx_sequence'], row['tx_size'], row['utr5_size'], row['cds_size'], row['utr3_size'])
    if utr5_seq is None:
        utr5_seq = original_utr5
        utr5_name = row['SYMBOL']
    if utr3_seq is None:
        utr3_seq = original_utr3
        utr3_name = row['SYMBOL']
    row['UTR5'] = utr5_name
    row['UTR3'] = utr3_name
    row['tx_sequence'] = utr5_seq + gene_to_optimize['NT_SEQ'] + utr3_seq
    row['tx_size'] = len(row['tx_sequence'])
    row['utr5_size'] = len(utr5_seq)
    row['cds_size'] = len(gene_to_optimize['NT_SEQ'])
    row['utr3_size'] = len(utr3_seq)
    return row


def create_replaced_utr_data(gene_to_optimize: pd.Series, true_data: pd.DataFrame, utr5_to_replace=None, utr3_to_replace=None) -> pd.DataFrame:
    assert (utr5_to_replace is None) or (utr3_to_replace is None)

    utr5_seq = None
    if utr5_to_replace is not None:
        row = true_data[true_data['SYMBOL'] == utr5_to_replace].iloc[0]
        utr5_seq, _, _ = extract_regions(row['tx_sequence'], row['tx_size'], row['utr5_size'], row['cds_size'], row['utr3_size'])

    utr3_seq = None
    if utr3_to_replace is not None:
        row = true_data[true_data['SYMBOL'] == utr3_to_replace].iloc[0]
        _, _, utr3_seq = extract_regions(row['tx_sequence'], row['tx_size'], row['utr5_size'], row['cds_size'], row['utr3_size'])

    true_data['UTR5'] = utr5_to_replace
    true_data['UTR3'] = utr3_to_replace
    true_data = true_data[['SYMBOL', 'UTR5', 'UTR3', 'fold', 'tx_size', 'utr5_size', 'cds_size', 'utr3_size', 'tx_sequence']]
    true_data = true_data.progress_apply(lambda row: replace_utrs(row, gene_to_optimize, utr5_seq, utr3_seq), axis=1)

    return true_data


        
AA_TO_CODON_LIST = {
    'M': ['ATG'], 
    'I': ['ATA', 'ATC', 'ATT'], 
    'L': ['TTA', 'TTG', 'CTA', 'CTC', 'CTG', 'CTT'], 
    'V': ['GTA', 'GTC', 'GTG', 'GTT'], 
    'F': ['TTC', 'TTT'], 
    'C': ['TGC', 'TGT'], 
    'A': ['GCA', 'GCC', 'GCG', 'GCT'], 
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'], 
    'T': ['ACA', 'ACC', 'ACG', 'ACT'], 
    'S': ['TCA', 'TCC', 'TCG', 'TCT', 'AGC', 'AGT'], 
    'Y': ['TAC', 'TAT'], 
    'W': ['TGG'], 
    'Q': ['CAA', 'CAG'], 
    'N': ['AAC', 'AAT'], 
    'H': ['CAC', 'CAT'], 
    'E': ['GAA', 'GAG'], 
    'D': ['GAC', 'GAT'], 
    'K': ['AAA', 'AAG'], 
    'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGA', 'AGG'], 
    'X': ['TAA', 'TAG', 'TGA'], 
}

def replace_codons(cds: str, main_codon: str, codons_to_replace: List) -> str:
    cds_codon_list = [cds[i:i+3] for i in range(0, len(cds), 3)]
    for i in range(len(cds_codon_list)):
        if cds_codon_list[i] in codons_to_replace:
            cds_codon_list[i] = main_codon
    new_cds = ''.join(cds_codon_list)
    return new_cds


def create_replaced_codons_data(gene_to_optimize: pd.Series, gene_with_best_utrs: pd.Series) -> (pd.DataFrame, pd.DataFrame):
    assert type(gene_to_optimize) == pd.Series
    assert type(gene_with_best_utrs) == pd.Series

    best_utr5_name = gene_with_best_utrs['UTR5']
    if pd.isna(best_utr5_name):
        best_utr5_name = gene_with_best_utrs['SYMBOL']
    best_utr3_name = gene_with_best_utrs['UTR3']
    if pd.isna(best_utr3_name):
        best_utr3_name = gene_with_best_utrs['SYMBOL']

    utr5, cds, utr3 = extract_regions(gene_with_best_utrs['tx_sequence'], gene_with_best_utrs['tx_size'], gene_with_best_utrs['utr5_size'], gene_with_best_utrs['cds_size'], gene_with_best_utrs['utr3_size'])
    codon_op_df = []
    codon_op_df.append({
        'SYMBOL': 'original',
        'Amino_Acid': 'original',
        'Codon': 'original',
        'UTR5': best_utr5_name,
        'UTR3': best_utr3_name,
        'fold': gene_with_best_utrs['fold'],
        'tx_size': gene_with_best_utrs['tx_size'],
        'utr5_size': gene_with_best_utrs['utr5_size'],
        'cds_size': gene_with_best_utrs['cds_size'],
        'utr3_size': gene_with_best_utrs['utr3_size'],
        'tx_sequence': utr5 + cds + utr3
    })
    
    for aa, codons in AA_TO_CODON_LIST.items():
        for codon in codons:
            new_cds = replace_codons(cds, codon, [c for c in codons if c != codon])
            codon_op_df.append({
                'SYMBOL': aa + "_to_" + codon,
                'Amino_Acid': aa,
                'Codon': codon,
                'UTR5': best_utr5_name,
                'UTR3': best_utr3_name,
                'fold': gene_with_best_utrs['fold'],
                'tx_size': gene_with_best_utrs['tx_size'],
                'utr5_size': gene_with_best_utrs['utr5_size'],
                'cds_size': gene_with_best_utrs['cds_size'],
                'utr3_size': gene_with_best_utrs['utr3_size'],
                'tx_sequence': utr5 + new_cds + utr3
            })
    
    codon_op_df = pd.DataFrame(codon_op_df)
    
    return codon_op_df

if __name__ == '__main__':
    tqdm.pandas()
    true_data_path = './data/data_with_human_TE_cellline_all_NA_plain.csv'
    gene_path = './gene_optimize/data/eGFP.csv'
    true_data = get_true_data(true_data_path)
    gene_to_optimize = get_gene_to_optimize(gene_path)
    # create_replaced_data(gene_to_optimize, true_data, save_name='./gene_optimize/data/eGFP_original_5utr_original_3utr.csv')
    df = create_replaced_utr_data(gene_to_optimize, true_data, utr3_to_replace='CWF19L2')
    df.to_csv('./gene_optimize/data/eGFP_original_5utr_CWF19L2_3utr.csv', index=False)
    
    data = pd.read_csv('./gene_optimize/data/eGFP_original_5utr_original_3utr_pred.csv', dtype={'SYMBOL': str, 'UTR5': str, 'UTR3': str})
    # codon optimize the top sequence
    codon_op_df = (gene_to_optimize, data.sort_values('pred_HEK293T', ascending=False).iloc[0])
    codon_op_df.to_csv('./gene_optimize/data/eGFP_original_5utr_original_3utr_pred_codon_optimized.csv', index=False)