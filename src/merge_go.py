#!/usr/bin/env python3
"""Merge all GO annotation to one file."""

from functools import wraps

import pandas as pd


class add_export(object):
    """Add export to function.

    :param path: path to save
    :type path: str
    :param sep: sep, defaults to ','
    :type sep: str
    :param header: header, defaults to False
    :type header: bool
    """

    def __init__(self):
        """Construct method."""
        pass

    def __call__(self, f):
        """Add export to function."""
        @wraps(f)
        def decorated(*args, **kwargs):
            dataframe: pd.DataFrame = f(*args, **kwargs)
            dataframe.to_csv(
                kwargs['export_path'],
                sep=kwargs['export_sep'],
                header=kwargs['export_header'],
                index=False
            )
        return decorated


@add_export()
def merge_go(*args, **kwargs):
    """Merge all GO annotation.

    :param *args: path to GO
    :param **kwargs: args for export
    """
    goes = [pd.read_table(path, names=['gene', 'go', 'type']) for path in args]
    go = pd.concat(goes).drop_duplicates()
    return go


def main():
    """Entry point for the script merge_go."""
    merge_go(
        'ITAG3.2_gene.tsv',
        'ITAG3.2_transcript.tsv',
        'ITAG4.0_gene.tsv',
        'ITAG4_transcript.tsv',
        export_path='tomato_go.tsv',
        export_sep='\t',
        export_header=None
    )


if __name__ == '__main__':
    main()
