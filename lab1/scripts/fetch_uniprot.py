# -*- coding: utf-8 -*-
"""Шаг 1. Выгрузка записи ABL1_HUMAN (P00519) из UniProt.

Печатает всё, что нужно для таблицы раздела 1 отчёта: название, ген, организм,
длину последовательности, домены, сайты связывания ATP и остаток в позиции 315.

Запуск:  python fetch_uniprot.py
Результат: data/uniprot_P00519.json + вывод в консоль
"""
import json
import os
import urllib.request

ACC = 'P00519'
OUT = os.path.join(os.path.dirname(__file__), '..', 'data', 'uniprot_%s.json' % ACC)

url = 'https://rest.uniprot.org/uniprotkb/%s.json' % ACC
print('GET', url)
d = json.load(urllib.request.urlopen(url, timeout=120))

os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(d, open(OUT, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

print('accession :', d['primaryAccession'], '|', d['uniProtkbId'], '|', d['entryType'])
print('name      :', d['proteinDescription']['recommendedName']['fullName']['value'])
print('EC        :', [e['value'] for e in d['proteinDescription']['recommendedName'].get('ecNumbers', [])])
print('gene      :', [g['geneName']['value'] for g in d['genes']])
print('organism  :', d['organism']['scientificName'])
print('length    :', d['sequence']['length'])

seq = d['sequence']['value']
print('остаток 315:', seq[314], '| окружение 306-325:', seq[305:325])

print('\n--- домены, мотивы, сайты ---')
INTEREST = ('Domain', 'Region', 'Motif', 'Active site', 'Binding site', 'Site', 'Lipidation')
for f in d['features']:
    if f['type'] in INTEREST:
        loc = f['location']
        print('%-14s %5s-%-5s %s' % (f['type'], loc['start']['value'], loc['end']['value'],
                                     f.get('description', '') or f.get('ligand', {}).get('name', '')))

print('\n--- function / disease ---')
for c in d['comments']:
    if c['commentType'] == 'FUNCTION':
        print('FUNCTION:', c['texts'][0]['value'][:400], '...')
    if c['commentType'] == 'DISEASE':
        if 'disease' in c:
            print('DISEASE :', c['disease']['diseaseId'], '(MIM %s)' % c['disease']['diseaseCrossReference']['id'])
        elif 'note' in c:
            print('NOTE    :', c['note']['texts'][0]['value'][:250])
