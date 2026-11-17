"""Exploratory data check (2026-11-10, cloud workspace): compare clone definitions (single
viral barcode vs full barcode set) against the published cross-stage counts. Superseded by
scripts/01_clone_counts.py, which implements the chosen definition. No paper result.
Usage: python 00_datacheck_vb2.py <clone_table.xlsx>"""
import sys, pandas as pd
x=pd.read_excel(sys.argv[1],sheet_name=None)
def load(k):
    d=x[k].copy(); d['vb']=d['Virus.Barcode'].astype(str).str.split(';'); return d
pairs=[('Clone_stage_I_Exp#1','Clone_stage_II_Exp#1'),('Clone_stage_II_Exp#2','Clone_stage_III_Exp#2'),('Clone_stage_III_Exp#3','Clone_stage_IV_Exp#3')]
for k in x:
    d=load(k); print(k,'cells',len(d),'unique cells',d['Cell.Barcode'].nunique(),'multi-vb %',round((d.vb.str.len()>1).mean()*100,1),
      'suffixes',d['Cell.Barcode'].astype(str).str.extract(r'(_\d+)$')[0].value_counts(dropna=False).to_dict())
print()
for a,b in pairs:
    A,B=load(a),load(b)
    # rule 1: each individual viral barcode is a clone tag
    ea,eb=A.explode('vb'),B.explode('vb')
    shared=set(ea.vb)&set(eb.vb)
    cells=ea[ea.vb.isin(shared)]['Cell.Barcode'].nunique()+eb[eb.vb.isin(shared)]['Cell.Barcode'].nunique()
    # rule 2: full barcode-set string as clone id
    sa,sb=set(A['Virus.Barcode']),set(B['Virus.Barcode'])
    sh2=sa&sb
    c2=A['Virus.Barcode'].isin(sh2).sum()+B['Virus.Barcode'].isin(sh2).sum()
    print(a,'->',b,'| per-barcode: clones',len(shared),'cells',cells,'| per-set: clones',len(sh2),'cells',c2)
    # clone size distribution under rule 2
    import collections
    print('   per-set clone sizes', pd.Series(list(A[A['Virus.Barcode'].isin(sh2)].groupby('Virus.Barcode').size()+B[B['Virus.Barcode'].isin(sh2)].groupby('Virus.Barcode').size())).describe().round(1).to_dict())
    print('   clusters later stage', B[B['Virus.Barcode'].isin(sh2)]['Cluster.Id'].value_counts().head(12).to_dict())
