"""Exploratory data check (2026-11-10, cloud workspace): print each clone-table sheet's shape,
columns and first rows. Superseded by scripts/01_clone_counts.py; produces no result used in
the paper. Usage: python 00_datacheck_vb1.py <clone_table.xlsx>"""
import sys, pandas as pd
x=pd.read_excel(sys.argv[1],sheet_name=None)
for k,v in x.items():
    print('=====',k,v.shape); print(v.columns.tolist()); print(v.head(5).to_string())
