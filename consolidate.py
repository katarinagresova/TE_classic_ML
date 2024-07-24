# read in a csv file to a dataframe, pivot it, and then save it as a csv file
import pandas as pd
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', type=str, help='input model results csv file', required=True)
    
    args = parser.parse_args()


    data = pd.read_csv(args.input_path)
    data = data.pivot(index='bio_source', columns='fold', values='nan_r2')
    output_path = os.path.join(os.path.dirname(args.input_path), 'pivot_model_results.csv')
    data.to_csv(output_path)