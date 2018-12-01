import os
import re
import pandas as pd

parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM4_CaiT/'

# Change working directory
os.chdir(parent_wd)

# Import experimental design, which stores library name, index, codename & experiment setup.
# exp_design = pd.read_excel('experimental_design_all.xlsx')
exp_design = pd.read_csv('exp_table_aug.xlsx', sep="\t")

# This is for easier typing in R
print('\"'+exp_design.codename.str.cat(sep="\",\"")+'\"')
print('\"'+exp_design.treat.str.cat(sep="\",\"")+'\"')
print('\"'+exp_design.dose.str.cat(sep="\",\"")+'\"')
print('\"'+exp_design.period.str.cat(sep="\",\"")+'\"')

saline = exp_design[exp_design.codename.str.startswith('S')]
print('\"'+saline.codename.str.cat(sep="\",\"")+'\"')

L1 = exp_design[exp_design.codename.str.startswith('L')]
L2 = exp_design[exp_design.codename.str.startswith('l')]
L = pd.concat([L1, L2])
print('\"'+L.codename.str.cat(sep="\",\"")+'\"')
print('\"'+L.exp_tag.str.cat(sep="\",\"")+'\"')

P1 = exp_design[exp_design.codename.str.startswith('P')]
P2 = exp_design[exp_design.codename.str.startswith('p')]
P = pd.concat([P1, P2])
print('\"'+P.codename.str.cat(sep="\",\"")+'\"')
print('\"'+P.exp_tag.str.cat(sep="\",\"")+'\"')

T = exp_design[exp_design.codename.str.startswith('T')]
print('\"'+T.codename.str.cat(sep="\",\"")+'\"')
print('\"'+T.exp_tag.str.cat(sep="\",\"")+'\"')

CNO = exp_design[exp_design.exp_tag.str.startswith('CRH')]
print('\"'+CNO.codename.str.cat(sep="\",\"")+'\"')
print('\"'+CNO.exp_tag.str.cat(sep="\",\"")+'\"')

CRS = exp_design[exp_design.exp_tag.str.startswith('CRS')]
print('\"'+CRS.codename.str.cat(sep="\",\"")+'\"')
print('\"'+CRS.exp_tag.str.cat(sep="\",\"")+'\"')

exp_design_new = exp_design
exp_design_new['treat'] = 'null'
for i, row in exp_design_new.iterrows():
    if 'S' in row.codename:
        exp_design_new.at[i, 'treat'] = 'Saline'
    elif ('L' in row.codename) or ('l' in row.codename):
        exp_design_new.at[i, 'treat']  = 'LPS'
    elif ('P' in row.codename) or ('p' in row.codename):
        exp_design_new.at[i, 'treat'] = 'Polyic'
    elif ('3q' in row.codename) or ('mc' in row.codename):
        exp_design_new.at[i, 'treat'] = 'CNO'
    elif 'T' in row.codename:
        exp_design_new.at[i, 'treat'] = 'TNFalpha'
    elif 'crs' in row.codename:
        exp_design_new.at[i, 'treat'] = 'CRS'

exp_design_new['dose'] = 'null'
for i, row in exp_design_new.iterrows():
    if '500ug' in row.experiment:
        exp_design_new.at[i, 'dose'] = '500ug'
    elif '1mg' in row.experiment:
        exp_design_new.at[i, 'dose'] = '1mg'
    elif '2mg' in row.experiment:
        exp_design_new.at[i, 'dose'] = '2mg'
    elif '10mg' in row.experiment:
        exp_design_new.at[i, 'dose'] = '10mg'
    elif '20mg' in row.experiment:
        exp_design_new.at[i, 'dose'] = '20mg'
    elif '50mg' in row.experiment:
        exp_design_new.at[i, 'dose'] = '50mg'
    else:
        exp_design_new.at[i, 'dose'] = '0g'

exp_design_new['period'] = 'null'
for i, row in exp_design_new.iterrows():
    if '6h' in row.experiment:
        exp_design_new.at[i, 'period'] = '6h'
    elif '3h' in row.experiment:
        exp_design_new.at[i, 'period'] = '3h'
    elif '3w' in row.experiment:
        exp_design_new.at[i, 'period'] = '3w'
    elif '4w' in row.experiment:
        exp_design_new.at[i, 'period'] = '4w'
    elif '5w' in row.experiment:
        exp_design_new.at[i, 'period'] = '5w'
    elif '8w' in row.experiment:
        exp_design_new.at[i, 'period'] = '8w'
    elif '14d' in row.experiment:
        exp_design_new.at[i, 'period'] = '2w'
    elif '3d' in row.experiment:
        exp_design_new.at[i, 'period'] = '3d'
    elif '1d' in row.experiment:
        exp_design_new.at[i, 'period'] = '1d'

exp_design_new.to_csv('exp_table_aug.xlsx', sep="\t",index=False)

