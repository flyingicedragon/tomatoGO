#!/usr/bin/env python3

# Python script for programtic access to Revigo. Run it with (last output file name is optional):
# python3 revigo.py example.csv result.csv

import time

import click
import requests


@click.command()
@click.option('--go', help='GO results.')
@click.option(
    '--data/--file', default=False, help='GO result is data, not file path.'
)
@click.option('--cutoff', default='0.7', help='Cutoff.')
@click.option(
    '--value-type',
    default='pvalue',
    help='Value type of data, default is pvalue.',
)
@click.option('--species-taxon', default='4081', help='Species taxon ID.')
@click.option('--out', default='revigo.tsv', help='Result file path. Default for system standard output.')
@click.option('--out-type', default='table', help='Result file type.')
@click.option(
    '--namespace',
    default='0',
    help='The namespace for which you are collecting results. 0 represents all, 1 represents BIOLOGICAL_PROCESS, 2 represents CELLULAR_COMPONENT, 3 represents MOLECULAR_FUNCTION.',
)
def revigo(
    go: str,
    data: bool,
    cutoff: str,
    value_type: str,
    species_taxon: str,
    out: str,
    out_type: str,
    namespace: str,
) -> None:
    if data:
        # Read enrichments file
        with open(go, 'r') as f:
            user_data = f.read()
    else:
        user_data = go

    # Submit job to Revigo
    payload = {
        'cutoff': cutoff,
        'valueType': value_type,
        'speciesTaxon': species_taxon,
        'measure': 'SIMREL',
        'goList': user_data,
    }
    if not namespace:
        r = requests.post('http://revigo.irb.hr/Revigo', data=payload)
        with open(out, 'w') as f:
            f.write(r.text)
        return
    r = requests.post('http://revigo.irb.hr/StartJob', data=payload)

    jobid = r.json()['jobid']

    running = True
    while not running:
        r = requests.get(
            'http://revigo.irb.hr/QueryJob',
            data={'jobid': jobid, 'type': 'jstatus'},
        )
        running = r.json()['running']
        time.sleep(1)

    r = requests.get(
        'http://revigo.irb.hr/QueryJob',
        data={'jobid': jobid, 'namespace': namespace, 'type': out_type},
    )
    # Write results to a file
    if out=='':
        print(r.text)
    else:
        with open(out, 'w') as f:
            f.write(r.text)


if __name__ == '__main__':
    revigo()
