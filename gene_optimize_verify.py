import pandas as pd
import argparse
from tqdm import tqdm

from lgbm_feature_extract_from_str import CODON_TO_AMINO_ACID
from gene_optimize_replace import get_gene_to_optimize, extract_regions

def correct_cds(row: pd.Series, gene_to_optimize: pd.Series):
    _, test_cds, _ = extract_regions(row['tx_sequence'], row['tx_size'], row['utr5_size'], row['cds_size'], row['utr3_size'])
    test_codons = [test_cds[i:i+3] for i in range(0, len(test_cds), 3)]
    test_amino_acids = [CODON_TO_AMINO_ACID[codon] for codon in test_codons]
    test_amino_acids = ''.join(test_amino_acids)
    return test_amino_acids == gene_to_optimize['AA_SEQ']

def verify(true_seq: pd.Series, test_seqs: pd.DataFrame):
    test_seqs['correct'] = test_seqs.progress_apply(lambda row: correct_cds(row, true_seq), axis=1)
    print(test_seqs['correct'].value_counts())
    #if any are false then print the false ones
    if False in test_seqs['correct'].value_counts():
        print("False sequences: ")
        print(test_seqs[test_seqs['correct'] == False])
    else:
        print("All sequences are correct!")


if __name__ == '__main__':
    tqdm.pandas()

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gene_to_optimize', type=str, default='./gene_optimize/data/eGFP.csv')
    parser.add_argument('-v', '--seqs_to_verify', type=str, required=True)
    args = parser.parse_args()

    print(f"Target Seq: {args.gene_to_optimize}")
    print(f"Verifying: {args.seqs_to_verify}")

    gene_to_optimize = get_gene_to_optimize('./gene_optimize/data/eGFP.csv')

    seqs_to_verify = pd.read_csv(args.seqs_to_verify)

    verify(gene_to_optimize, seqs_to_verify)