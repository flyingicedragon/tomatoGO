#!/usr/bin/env python3
"""Pretty GO annotation for GO."""

import pandas as pd
import os.path


def import_old(path: str) -> pd.DataFrame:
    """Import old GO annotation.

    :param path: path to old GO annotation file
    :type path: str
    :return: a dataframe containing GO annotation
    :rtype: pd.DataFrame
    """
    go_anno = pd.read_table(
        path,
        header=None
    )
    go_anno.iloc[:, 2] = go_anno.iloc[:, 2].str[0:14]
    return go_anno.iloc[:, [2, 3]]


def import_go_basic(path: str) -> pd.DataFrame:
    """Import GO basic annotation.

    :param path: path to GO basic file
    :type path: str
    :return: a dataframe containing GO basic
    :rtype: pd.DataFrame
    """
    go_basic = pd.read_table(path)
    return go_basic.loc[:, ['GO', 'level']]


def merge_go(go: pd.DataFrame, go_basic: pd.DataFrame) -> pd.DataFrame:
    """Merge GO and GO basic.

    :param go: GO
    :type go: pd.DataFrame
    :param go_basic: GO basic
    :type go_basic: pd.DataFrame
    :return: dataframe containing GO
    :rtype: pd.DataFrame
    """
    go.columns = ['gene', 'go']
    go_basic.columns = ['go', 'level']
    go: pd.DataFrame = pd.merge(go, go_basic, on='go', how='left')
    return go


def export_go(go: pd.DataFrame, path: str):
    """Export prettied GO annotation.

    :param go: dataframe containing GO annotation
    :type go: pd.DataFrame
    :param path: origin GO path
    :type path: str
    """
    go.columns = ['gene_id', 'GO', 'level']
    go.to_csv(
        path,
        index=False,
        sep='\t'
    )


def main():
    """Entry point for the script pretty_GO_anno."""
    go_basic = import_go_basic('go-basic.tb')
    for file in [
            'ITAG3.2_gene.tsv',
            'ITAG3.2_transcript.tsv',
            'ITAG4_transcript.tsv'
    ]:
        file_path = os.path.join('data', file)
        go = import_old(file_path)
        go = merge_go(go, go_basic)
        export_go(go, file)


if __name__ == '__main__':
    main()
