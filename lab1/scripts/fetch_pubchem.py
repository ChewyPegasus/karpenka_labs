# -*- coding: utf-8 -*-
"""Шаг 4. Химические идентификаторы ингибиторов из PubChem.

Для каждого соединения по названию получаем CID, брутто-формулу, молекулярную
массу, SMILES (изомерный, со стереохимией) и InChIKey.

Запуск:  python fetch_pubchem.py
Результат: data/pubchem.json + таблица в консоль
"""
import json
import os
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')

NAMES = ['Imatinib', 'Nilotinib', 'Dasatinib', 'Bosutinib', 'Ponatinib', 'Asciminib', 'Axitinib']
PROPS = 'MolecularFormula,MolecularWeight,SMILES,ConnectivitySMILES,IUPACName,InChIKey'
URL = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/' + PROPS + '/JSON')

out = {}
for name in NAMES:
    d = json.load(urllib.request.urlopen(URL % urllib.parse.quote(name), timeout=60))
    d = d['PropertyTable']['Properties'][0]
    out[name] = d
    print('%-10s CID %-9s %-22s %7s  %s'
          % (name, d['CID'], d['MolecularFormula'], d['MolecularWeight'], d['InChIKey']))
    print('    SMILES:', d.get('SMILES') or d.get('ConnectivitySMILES'))

os.makedirs(DATA, exist_ok=True)
json.dump(out, open(os.path.join(DATA, 'pubchem.json'), 'w'), indent=1)
print('\nсохранено в data/pubchem.json')
