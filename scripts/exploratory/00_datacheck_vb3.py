"""Exploratory data check (2026-11-11, cloud workspace): overlap of stage-IV clone-table cell
barcodes with each stage-IV library's barcodes. Superseded by scripts/02_library_map.py.
No paper result. Usage: python 00_datacheck_vb3.py <clone_table.xlsx> <barcodes.tsv.gz>..."""
import sys,gzip,pandas as pd
x=pd.read_excel(sys.argv[1],sheet_name='Clone_stage_IV_Exp#3')
cb=set(x['Cell.Barcode'].astype(str).str.replace('-1',''))
for f in sys.argv[2:]:
    s=set(l.strip().split('-')[0] for l in gzip.open(f,'rt'))
    print(f.split('/')[-1], 'cells', len(s), 'overlap with stage-IV barcode table', len(cb & s))
