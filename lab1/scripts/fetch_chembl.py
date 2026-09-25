# -*- coding: utf-8 -*-
"""Шаг 3. Выгрузка экспериментальных активностей из ChEMBL.

Для семи выбранных ингибиторов забираем ВСЕ записи IC50 / Ki / Kd против двух
мишеней:
    CHEMBL1862    — ABL1 (изолированная киназа)
    CHEMBL2096618 — BCR-ABL1 (слитный белок)
Форму мишени (WT или мутант) даёт поле assay_variant_mutation.

Запуск:  python fetch_chembl.py
Результат: data/chembl_activities.json + краткая сводка в консоль
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
BASE = 'https://www.ebi.ac.uk/chembl/api/data/'

MOLECULES = {
    'Imatinib':  'CHEMBL941',
    'Nilotinib': 'CHEMBL255863',
    'Dasatinib': 'CHEMBL1421',
    'Bosutinib': 'CHEMBL288441',
    'Ponatinib': 'CHEMBL1171837',
    'Asciminib': 'CHEMBL4208229',
    'Axitinib':  'CHEMBL1289926',
}
TARGETS = ['CHEMBL1862', 'CHEMBL2096618']

FIELDS = ['molecule_chembl_id', 'target_chembl_id', 'target_pref_name', 'standard_type',
          'standard_relation', 'standard_value', 'standard_units', 'assay_chembl_id',
          'assay_description', 'assay_type', 'document_chembl_id', 'document_year',
          'assay_variant_mutation', 'assay_variant_accession', 'activity_id',
          'data_validity_comment']


def get(path, **params):
    """GET к ChEMBL API с парой повторов — сервер иногда отвечает медленно."""
    params['format'] = 'json'
    url = BASE + path + '?' + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            return json.load(urllib.request.urlopen(url, timeout=120))
        except Exception as e:
            print('  повтор (%s)' % e, file=sys.stderr)
            time.sleep(3)
    raise RuntimeError('ChEMBL не ответил: ' + url)


def main():
    out = []
    for name, mol in MOLECULES.items():
        for tgt in TARGETS:
            offset, n = 0, 0
            while True:
                d = get('activity', molecule_chembl_id=mol, target_chembl_id=tgt,
                        standard_type__in='IC50,Ki,Kd', limit=1000, offset=offset)
                for a in d['activities']:
                    rec = {k: a.get(k) for k in FIELDS}
                    rec['name'] = name
                    out.append(rec)
                    n += 1
                if not d['page_meta']['next']:
                    break
                offset += 1000
            print('%-10s %-14s %3d записей' % (name, tgt, n))

    os.makedirs(DATA, exist_ok=True)
    json.dump(out, open(os.path.join(DATA, 'chembl_activities.json'), 'w'), indent=0)
    print('\nвсего выгружено:', len(out))

    # сколько из них относится к интересующим нас формам мишени
    wt = sum(1 for a in out if a['assay_variant_mutation'] is None)
    mut = sum(1 for a in out if a['assay_variant_mutation'] == 'T315I')
    print('без пометки мутации (WT):', wt, '| помечено T315I:', mut)


if __name__ == '__main__':
    main()
