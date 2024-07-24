import argparse
import os
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir_path', '-d', type=str, required=True)

    args = parser.parse_args()
    # get dirs in dir_path
    dirs = [d for d in os.listdir(args.dir_path) if os.path.isdir(os.path.join(args.dir_path, d))]

    # summary = pd.DataFrame()
    dfs = []
    for dir in dirs:
        data = pd.read_csv(os.path.join(args.dir_path, dir, 'model_results.csv'), index_col=0)
        data['model'] = dir
        # make dir column first
        cols = data.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        data = data[cols]
        data = data[data['fold'] != 'mean']
        dfs.append(data)
        # mean_data = data[data['fold'] == 'mean']
        # append_len = len(summary)
        # summary.loc[append_len, 'model'] = dir
        # for col in mean_data.columns:
        #     summary.loc[append_len, col] = mean_data.loc[:, col].values[0]
    summary = pd.concat(dfs, axis=0)
    summary.to_csv(os.path.join(args.dir_path, 'summary.csv'), index=False)

    pivot = summary.pivot(index='fold', columns='model', values='nan_r2')
    pivot.index.name = "fold"
    pivot.to_csv(os.path.join(args.dir_path, 'summary_pivot.csv'))
