"""Locations of inputs and outputs. Override with environment variables.

SCCNV_DATA     directory holding the raw inputs (Seurat .rds files, WA71 DNA segment table)      default: ./data
SCCNV_WORK     scratch directory (R exports, reference tables)                                     default: ./work
SCCNV_RESULTS  where scripts write result tables / figures                                         default: ./results/run
"""
import os
from pathlib import Path

DATA = Path(os.environ.get('SCCNV_DATA', 'data'))
WORK = Path(os.environ.get('SCCNV_WORK', 'work'))
RESULTS = Path(os.environ.get('SCCNV_RESULTS', 'results/run'))
EXPORT = WORK / 'export'
REF = WORK / 'ref'

# raw input file names expected in SCCNV_DATA (see README)
FILES = dict(
    celllines='pca.with.annotation.rds',
    wa71='wa71.lnmet.rds',
    wa71_segments='WA_71_segment_data.csv',
    wong='Wong.merge.rds',
    healthy='P16_2018_CellRep_Henry_obj.split.list.rds',
)

def results_path(name):
    RESULTS.mkdir(parents=True, exist_ok=True)
    return str(RESULTS / name)
