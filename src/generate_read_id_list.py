#!/usr/bin/env python3

from pathlib import Path
import os
from snakemake import io
import logging

summary_file = snakemake.input['seqsum']
idlist_file = snakemake.output['idlist']
log_file = snakemake.log[0]


logging.basicConfig(
    filename=log_file,
    format='%(asctime)s %(levelname)-8s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO)


# make sure the basecall file is present
try:
    os.stat(summary_file)
except FileNotFoundError:
    print(f'ERROR. {summary_file} not found.')
    print('       Run basecall.Snakefile')
    raise FileNotFoundError

guppy_dir = Path(summary_file).parent

# need to handle old version of guppy that don't have pass/fail dirs
if Path(guppy_dir, 'pass').is_dir():
    my_read_path = Path(guppy_dir, 'pass', '{read}.fastq')
else:
    my_read_path = Path(guppy_dir, '{read}.fastq')
logging.info(f'Globbing path {my_read_path}')
my_read_names = io.glob_wildcards(my_read_path).read
logging.info('Found read files')
logging.info(my_read_names)
# make sure the fastq file has reads in it
non_empty_read_names = (
    x for x in my_read_names 
    if os.stat(my_read_path.as_posix().format(read=x)).st_size > 0)
# print the read IDs to idlist_file
with open(idlist_file, 'wt') as f:
    for id in non_empty_read_names:
        f.write(f'{id}\n')
