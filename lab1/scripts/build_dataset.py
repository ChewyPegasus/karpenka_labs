# -*- coding: utf-8 -*-
"""Шаг 5. Сборка BCR_ABL_dataset.csv / .xlsx.

Объединяет выгрузку активностей (data/chembl_activities.json) и химические
идентификаторы (data/pubchem.json) в одну таблицу по заданию:

    Compound; PubChem_CID; SMILES; Target; Mutation; Activity_type;
    Activity_value; Units; PDB_ID; Source

Правила:
  * одна строка = одно экспериментальное измерение, значения НЕ усредняются;
  * берём только формы WT (assay_variant_mutation = None) и T315I;
  * цензурированные значения сохраняем со знаком: >10000;
  * PDB_ID проставляем там, где соединение реально закристаллизовано с этой
    формой белка (см. раздел 4 отчёта), иначе поле пустое;
  * оставляем только отобранные публикации, чтобы в наборе не было
    случайных клеточных анализов без привязки к форме мишени.

Запуск:  python build_dataset.py
Результат: BCR_ABL_dataset.csv и BCR_ABL_dataset.xlsx в корне работы
"""
import csv
import json
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
DATA = os.path.join(ROOT, 'data')

# Отобранные публикации: для каждой указан человекочитаемый источник
DOCS = {
    'CHEMBL5715937': "ChEMBL; O'Hare T. et al., Cancer Res. 2005, 65:4500-4505 (PMID 15930265)",
    'CHEMBL5678804': 'ChEMBL; Chan W.W. et al., Cancer Cell 2011, 19:556-568 (PMID 21481795)',
    'CHEMBL1240341': 'ChEMBL; Remsing Rix L.L. et al., Leukemia 2009, 23:477-485 (PMID 19039322)',
    'CHEMBL1177713': 'ChEMBL; Huang W.-S. et al., J. Med. Chem. 2010, 53:4701-4719 (PMID 20513156)',
    'CHEMBL1908390': 'ChEMBL; Davis M.I. et al., Nat. Biotechnol. 2011, 29:1046-1051 (PMID 22037378)',
    'CHEMBL4196169': 'ChEMBL; Schoepfer J. et al., J. Med. Chem. 2018, 61:8120-8135 (PMID 30137981)',
    'CHEMBL2311339': 'ChEMBL; Ren X. et al., J. Med. Chem. 2013, 56:879-894 (PMID 23301703)',
    'CHEMBL1144455': 'ChEMBL; Fabian M.A. et al., Nat. Biotechnol. 2005, 23:329-336 (PMID 15711537)',
    'CHEMBL3886395': 'ChEMBL; patent record CHEMBL3886395 (axitinib vs BCR-ABL1 T315I), 2015',
    'CHEMBL5137090': 'ChEMBL; Zhao Y. et al., J. Med. Chem. 2022 (PMID 35944901)',
    'CHEMBL3603836': 'ChEMBL; Yang X. et al., Bioorg. Med. Chem. Lett. 2015, 25:3821 (PMID 26195136)',
}

# С какой структурой PDB связано соединение в каждой форме белка
PDB = {('Imatinib', 'WT'): '2HYY', ('Nilotinib', 'WT'): '3CS9', ('Dasatinib', 'WT'): '2GQG',
       ('Bosutinib', 'WT'): '3UE4', ('Ponatinib', 'WT'): '3OXZ', ('Axitinib', 'WT'): '4WA9',
       ('Asciminib', 'WT'): '5MO4', ('Axitinib', 'T315I'): '4TWP',
       ('Ponatinib', 'T315I'): '3IK3', ('Asciminib', 'T315I'): '5MO4'}

ORDER = ['Imatinib', 'Nilotinib', 'Dasatinib', 'Bosutinib', 'Ponatinib', 'Asciminib', 'Axitinib']
HEADER = ['Compound', 'PubChem_CID', 'SMILES', 'Target', 'Mutation', 'Activity_type',
          'Activity_value', 'Units', 'PDB_ID', 'Source']

acts = json.load(open(os.path.join(DATA, 'chembl_activities.json'), encoding='utf-8'))
pc = json.load(open(os.path.join(DATA, 'pubchem.json'), encoding='utf-8'))

rows = []
for a in acts:
    if a['document_chembl_id'] not in DOCS:
        continue
    if a['assay_variant_mutation'] not in (None, 'T315I'):
        continue
    if not a['standard_value'] or a['standard_units'] != 'nM':
        continue
    name = a['name']
    mutation = 'T315I' if a['assay_variant_mutation'] == 'T315I' else 'WT'
    target = 'BCR-ABL1' if a['target_chembl_id'] == 'CHEMBL2096618' else 'ABL1'
    sign = a['standard_relation'] if a['standard_relation'] not in ('=', None) else ''
    value = '%s%g' % (sign, float(a['standard_value']))
    prop = pc[name]
    rows.append([name, prop['CID'], prop.get('SMILES') or prop.get('ConnectivitySMILES'),
                 target, mutation, a['standard_type'], value, 'nM',
                 PDB.get((name, mutation), ''),
                 '%s; ChEMBL doc %s, assay %s' % (DOCS[a['document_chembl_id']],
                                                  a['document_chembl_id'], a['assay_chembl_id'])])

rows.sort(key=lambda r: (ORDER.index(r[0]), r[4] != 'WT', r[5], r[3]))

# уборка полных дублей (одно и то же значение из одного анализа)
seen, uniq = set(), []
for r in rows:
    key = tuple(r[:8]) + (r[9],)
    if key in seen:
        continue
    seen.add(key)
    uniq.append(r)

csv_path = os.path.join(ROOT, 'BCR_ABL_dataset.csv')
with open(csv_path, 'w', newline='', encoding='utf-8-sig') as f:   # BOM — чтобы Excel не ломал кириллицу
    w = csv.writer(f, delimiter=';')
    w.writerow(HEADER)
    w.writerows(uniq)
print('записано строк:', len(uniq), '->', csv_path)
print(Counter((r[0], r[4]) for r in uniq))

try:
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = 'dataset'
    ws.append(HEADER)
    for r in uniq:
        ws.append(r)
    for col, width in zip('ABCDEFGHIJ', [12, 12, 60, 10, 9, 12, 14, 7, 9, 70]):
        ws.column_dimensions[col].width = width
    ws.freeze_panes = 'A2'
    wb.save(os.path.join(ROOT, 'BCR_ABL_dataset.xlsx'))
    print('также сохранено в BCR_ABL_dataset.xlsx')
except ImportError:
    print('openpyxl не установлен — xlsx не создан (csv достаточно)')
