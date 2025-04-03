import argparse
import pandas as pd
import os
import time
import numpy as np
from tqdm import tqdm
import yaml
import random
from typing import Tuple

from gene_optimize_replace import get_gene_to_optimize, get_true_data, create_replaced_utr_data, create_replaced_codons_data, extract_regions, AA_TO_CODON_LIST
from gene_optimize_model import LGBM_TE_model

def get_best_codons(codon_replaced_df: pd.DataFrame, pred_name: str, only_better_performing: bool=False, strictly_increase: bool=False, target_cell_line: str=None) -> pd.DataFrame:
    original_diff = codon_replaced_df[codon_replaced_df['SYMBOL'] == 'original'][pred_name].values[0]
    if strictly_increase:
        assert target_cell_line is not None
        original_pred = codon_replaced_df[codon_replaced_df['SYMBOL'] == 'original'][target_cell_line].values[0]
        codon_replaced_df = codon_replaced_df[codon_replaced_df[target_cell_line] >= original_pred]
    
    grouped_by_aa = codon_replaced_df.groupby('Amino_Acid')
    best_codon_per_aa = grouped_by_aa.apply(lambda x: x.sort_values(pred_name, ascending=False).iloc[0])
    best_codon_per_aa = best_codon_per_aa.set_index('Amino_Acid')[['Codon', pred_name]]
    if only_better_performing:
        best_codon_per_aa = best_codon_per_aa[best_codon_per_aa[pred_name] > original_diff]
    # print(best_codon_per_aa)
    # assert False
    return best_codon_per_aa

def get_replaced_cds(gene_to_optimize: pd.Series, best_codon_per_aa: pd.DataFrame) -> str:
    old_codons = [gene_to_optimize['NT_SEQ'][i:i+3] for i in range(0, len(gene_to_optimize['NT_SEQ']), 3)]
    assert len(old_codons) == len(gene_to_optimize['AA_SEQ'])
    new_cds = ''
    for aa, old_codon in zip(gene_to_optimize['AA_SEQ'], old_codons):
        if aa not in best_codon_per_aa.index or best_codon_per_aa.loc[aa]['Codon'] == 'original':
            new_cds += old_codon
        else:
            new_cds += best_codon_per_aa.loc[aa]['Codon']

    assert len(new_cds) == len(gene_to_optimize['NT_SEQ'])
    return new_cds

def create_best_seq(config, gene_to_optimize: pd.Series, gene_with_best_utrs:pd.Series, codon_replaced_df: pd.DataFrame, 
                    pred_name: str, target_cell_line: str) -> Tuple[pd.Series, dict]:
    
    best_codons_for_cell_line = get_best_codons(codon_replaced_df, pred_name, config['only_better_performing'], config['only_increase_target_cell_line'], target_cell_line)
    new_cds = get_replaced_cds(gene_to_optimize, best_codons_for_cell_line)
    
    utr5, _, utr3 = extract_regions(gene_with_best_utrs['tx_sequence'], gene_with_best_utrs['tx_size'], gene_with_best_utrs['utr5_size'], gene_with_best_utrs['cds_size'], gene_with_best_utrs['utr3_size'])
    seq = utr5 + new_cds + utr3
    best_codons_str = str(best_codons_for_cell_line['Codon'].to_dict())
    new_seq = pd.Series({
        'SYMBOL': 'best_utrs_best_codons',
        'UTR5': gene_with_best_utrs['UTR5'],
        'UTR3': gene_with_best_utrs['UTR3'],
        'replaced_codons': best_codons_str,
        'fold': gene_with_best_utrs['fold'],
        'tx_size': len(seq),
        'utr5_size': len(utr5),
        'cds_size': len(new_cds),
        'utr3_size': len(utr3),
        'tx_sequence': seq
    })
    return new_seq, best_codons_for_cell_line['Codon'].to_dict()

def codon_dgd(config, gene_to_optimize: pd.Series, gene_with_best_utrs:pd.Series, max_model: LGBM_TE_model, min_model: LGBM_TE_model, best_simple_codons: dict) -> pd.Series:

    utr5, _, utr3 = extract_regions(gene_with_best_utrs['tx_sequence'], gene_with_best_utrs['tx_size'], gene_with_best_utrs['utr5_size'], gene_with_best_utrs['cds_size'], gene_with_best_utrs['utr3_size'])
    search_space = AA_TO_CODON_LIST.copy()
    search_space.pop('M')
    search_space.pop('W')
    search_space = {aa: codons + ['original'] for aa, codons in search_space.items()}
    diff_atempts = []
    for num_init in tqdm(range(config['dgd_num_inits']), desc='dgd inits'):
        random.seed(config['dgd_seed'] + num_init)
        best_codons = None
        if config['init_from_simple_codon_op'] == 'random':
            best_codons = {aa:'original' for aa in search_space.keys()}
            for aa, codon in best_simple_codons.items():
                best_codons[aa] = codon
        else:
            best_codons = {aa: random.choice(codons) for aa, codons in search_space.items()}
        best_diff = -np.inf
        for i in tqdm(range(config['dgd_max_iter']), desc='dgd iter'):
            converge = 0
            amino_acids = list(search_space.keys())
            random.shuffle(amino_acids)
            for aa in amino_acids:
                possible_steps = search_space[aa]
                step_results = []
                for step in possible_steps:
                    step_codons = best_codons.copy()
                    step_codons[aa] = step
                    step_codons_df = pd.DataFrame.from_dict(step_codons, orient='index', columns=['Codon'])
                    new_cds = get_replaced_cds(gene_to_optimize, step_codons_df)
                    new_seq = utr5 + new_cds + utr3
                    step_results.append({
                        'SYMBOL': f"{aa}_to_{step}",
                        'replaced_codons': step_codons,
                        'tx_size': len(new_seq),
                        'utr5_size': len(utr5),
                        'cds_size': len(new_cds),
                        'utr3_size': len(utr3),
                        'tx_sequence': new_seq
                    })
                step_results = pd.DataFrame(step_results)
                step_results = max_model.predict_TE_with_single_fold(step_results, gene_with_best_utrs['fold'], pred_name='pred_max')
                step_results = min_model.predict_TE_with_single_fold(step_results, gene_with_best_utrs['fold'], pred_name='pred_min')
                step_results['max-min'] = step_results['pred_max'] - step_results['pred_min']
                if config['dgd_only_increase']:
                    prev_step_best_pred = step_results[step_results['replaced_codons'] == best_codons]['pred_max'].values[0]
                    step_results = step_results[step_results['pred_max'] >= prev_step_best_pred]
                step_results = step_results.sort_values('max-min', ascending=False)[['replaced_codons', 'pred_max', 'pred_min', 'max-min']]
                # print(step_results)
                if best_codons[aa] == step_results.iloc[0]['replaced_codons'][aa]:
                    converge += 1
                else:
                    converge = 0
                best_codons[aa] = step_results.iloc[0]['replaced_codons'][aa]
                best_diff = step_results.iloc[0]['max-min']

            if converge == len(amino_acids):
                print(f'converged after {i} iterations')
                break
        diff_atempts.append({'diffs': best_diff, 'codons': best_codons})
    

    diff_atempts = pd.DataFrame(diff_atempts).sort_values('diffs', ascending=False)
    diff_atempts.to_csv(os.path.join(config['results_path'], f"{gene_to_optimize['SYMBOL']}_dgd_diffs.csv"), index=False, mode='a')
    best_codons = diff_atempts.iloc[0]['codons']
    # print(best_codons)
    new_cds = get_replaced_cds(gene_to_optimize, pd.DataFrame.from_dict(step_codons, orient='index', columns=['Codon']))
    seq = utr5 + new_cds + utr3
    new_seq = pd.Series({
        'SYMBOL': 'best_utrs_dgd_codons',
        'UTR5': gene_with_best_utrs['UTR5'],
        'UTR3': gene_with_best_utrs['UTR3'],
        'replaced_codons': str(best_codons),
        'fold': gene_with_best_utrs['fold'],
        'tx_size': len(seq),
        'utr5_size': len(utr5),
        'cds_size': len(new_cds),
        'utr3_size': len(utr3),
        'tx_sequence': seq
    })
    return new_seq

def minmax(config, max_cell_line: str, min_cell_line: str, max_model: LGBM_TE_model, min_model: LGBM_TE_model, best_utrs: pd.Series, gene_to_optimize: pd.Series):
    # max cell line 1, min cell line 2
    assert max_cell_line != min_cell_line

    # best_utrs = utr_replaced_df.sort_values(f'pred_{max_cell_line}-pred_{min_cell_line}', ascending=False).iloc[0]
    codon_replaced = create_replaced_codons_data(gene_to_optimize, best_utrs)
    codon_replaced = max_model.predict_TE_with_single_fold(codon_replaced, best_utrs['fold'], pred_name=f'co_pred_{max_cell_line}')
    codon_replaced = min_model.predict_TE_with_single_fold(codon_replaced, best_utrs['fold'], pred_name=f'co_pred_{min_cell_line}')
    codon_replaced[f'co_pred_{max_cell_line}-co_pred_{min_cell_line}'] = codon_replaced[f'co_pred_{max_cell_line}'] - codon_replaced[f'co_pred_{min_cell_line}']

    codon_replaced.to_csv(os.path.join(config['results_path'], f"{gene_to_optimize['SYMBOL']}_codon_replaced_for_max_{max_cell_line}_vs_min_{min_cell_line}.csv"), index=False)

    best_seq, simple_best_codons = create_best_seq(
        config, gene_to_optimize, best_utrs, codon_replaced, 
        f'co_pred_{max_cell_line}-co_pred_{min_cell_line}', target_cell_line=f'co_pred_{max_cell_line}'
    )
    best_seq = pd.DataFrame([best_seq])

    original_seq = codon_replaced[codon_replaced['SYMBOL'] == 'original'].drop(columns=['Amino_Acid', 'Codon'])
    original_seq = original_seq.loc[:, ~original_seq.columns.str.contains('pred')]
    original_seq['replaced_codons'] = 'original'
    final_seqs = [best_seq, original_seq]

    if config['discrete_gradient_descent']:
        dgd_seq = pd.DataFrame([codon_dgd(config, gene_to_optimize, best_utrs, max_model, min_model, simple_best_codons)])
        final_seqs.append(dgd_seq)

    best_seq = pd.concat(final_seqs, axis=0, ignore_index=True)

    best_seq = max_model.predict_TE_with_single_fold(best_seq, best_utrs['fold'], pred_name=f'pred_{max_cell_line}')
    best_seq = min_model.predict_TE_with_single_fold(best_seq, best_utrs['fold'], pred_name=f'pred_{min_cell_line}')
    best_seq[f'pred_{max_cell_line}-pred_{min_cell_line}'] = best_seq[f'pred_{max_cell_line}'] - best_seq[f'pred_{min_cell_line}']

    best_seq.to_csv(os.path.join(config['results_path'], f"{gene_to_optimize['SYMBOL']}_best_seq_for_max_{max_cell_line}_vs_min_{min_cell_line}.csv"), index=False)


def get_best_utrs(utr_replaced_df: pd.DataFrame, cell_line_1: str, cell_line_2: str) -> Tuple[pd.Series, pd.Series]:
    # any pair of sequences can not have a difference in length of more than 20%
    # want to find the best pair of sequences that have the smallest difference in length and lowest sum rank
    # if there are multiple pairs with the same sum rank, choose the pair with the highest te diff
    # only search top k
    k = 10
    diff_threshold = 0.2
    top_k_cell_line_1 = utr_replaced_df.sort_values(f'pred_{cell_line_1}-pred_{cell_line_2}', ascending=False).iloc[:k]
    top_k_cell_line_1 = top_k_cell_line_1.reset_index(drop=True)
    top_k_cell_line_2 = utr_replaced_df.sort_values(f'pred_{cell_line_2}-pred_{cell_line_1}', ascending=False).iloc[:k]
    top_k_cell_line_2 = top_k_cell_line_2.reset_index(drop=True)

    best_pair = None
    best_sum_rank = np.inf
    best_te_diff = -np.inf
    for index_1, row_1 in top_k_cell_line_1.iterrows():
        for index_2, row_2 in top_k_cell_line_2.iterrows():
            if row_1['SYMBOL'] == row_2['SYMBOL']:
                continue
            len_diff = abs(row_1['tx_size'] - row_2['tx_size']) / ((row_1['tx_size'] + row_2['tx_size']) / 2)
            if len_diff > diff_threshold:
                continue
            sum_rank = index_1 + index_2
            te_diff = row_1[f'pred_{cell_line_1}-pred_{cell_line_2}']
            if sum_rank < best_sum_rank or (sum_rank == best_sum_rank and te_diff > best_te_diff):
                best_pair = (row_1, row_2)
                best_sum_rank = sum_rank
                best_te_diff = te_diff

    assert best_pair is not None
    return best_pair


if __name__ == "__main__":
    tqdm.pandas()

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', type=str, default='./gene_optimize_minmax_config.yaml', help='config file')
    args = parser.parse_args()

    config = yaml.load(open(args.config, 'r'), Loader=yaml.FullLoader)
    
    cell_line_1 = config['cell_line_1']
    cell_line_2 = config['cell_line_2']
    assert cell_line_1 != cell_line_2

    cell_line_1_model = LGBM_TE_model(config['cell_line_1_model_path'])
    cell_line_2_model = LGBM_TE_model(config['cell_line_2_model_path'])

    assert cell_line_1_model.bio_source == cell_line_1
    assert cell_line_2_model.bio_source == cell_line_2
    assert cell_line_1_model.features_to_extract == cell_line_2_model.features_to_extract

    if not os.path.exists(config['results_path']):
        os.makedirs(config['results_path'])

    true_data = get_true_data(config['true_data_path'])
    gene_to_optimize = get_gene_to_optimize(config['gene_path'])
    utr_replaced_df = create_replaced_utr_data(gene_to_optimize, true_data)

    utr_replaced_df = utr_replaced_df[utr_replaced_df['tx_size'] <= config['max_tx_size']]

    utr_replaced_df = cell_line_1_model.predict_TE(utr_replaced_df, pred_name=f'pred_{cell_line_1}')
    utr_replaced_df = cell_line_2_model.predict_TE(utr_replaced_df, pred_name=f'pred_{cell_line_2}')

    utr_replaced_df[f'pred_{cell_line_1}-pred_{cell_line_2}'] = utr_replaced_df[f'pred_{cell_line_1}'] - utr_replaced_df[f'pred_{cell_line_2}']
    utr_replaced_df[f'pred_{cell_line_2}-pred_{cell_line_1}'] = utr_replaced_df[f'pred_{cell_line_2}'] - utr_replaced_df[f'pred_{cell_line_1}']

    utr_replaced_df = utr_replaced_df.sort_values(f'pred_{cell_line_1}-pred_{cell_line_2}', ascending=False)

    utr_for_cell_line_1, utr_for_cell_line_2 = get_best_utrs(utr_replaced_df, cell_line_1, cell_line_2)

    utr_replaced_df.to_csv(os.path.join(config['results_path'], f"{gene_to_optimize['SYMBOL']}_utr_replaced.csv"), index=False)

    minmax(config, cell_line_1, cell_line_2, cell_line_1_model, cell_line_2_model, utr_for_cell_line_1, gene_to_optimize)
    minmax(config, cell_line_2, cell_line_1, cell_line_2_model, cell_line_1_model, utr_for_cell_line_2, gene_to_optimize)



    
    
    
