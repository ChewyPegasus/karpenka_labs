# -*- coding: utf-8 -*-
"""Шаг 2. Обзор всех структур PDB, относящихся к ABL1.

Что делает:
  1) через Search API RCSB находит все записи PDB, в которых есть ссылка на
     UniProt P00519 (человеческий ABL1) или P00520 (мышиный ортолог);
  2) через Data API (GraphQL) забирает метаданные: метод, разрешение, лиганды,
     цепи, мутации, организм, последовательность;
  3) определяет форму белка по последовательности вокруг остатка 315:
         ...FYIITEF...  -> WT      (Thr315)
         ...FYIIIEF...  -> T315I   (Ile315)
     Это надёжнее, чем читать поле pdbx_mutation: оно заполнено не везде и
     иногда в нумерации изоформы 1b (T334I вместо T315I).

Запуск:  python survey_pdb.py
Результат: data/abl_entries_raw.json, data/abl_rows.json + статистика в консоль
"""
import json
import os
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')

SEARCH = 'https://search.rcsb.org/rcsbsearch/v2/query?json='
GRAPHQL = 'https://data.rcsb.org/graphql'

# Ионы, криопротекторы и прочие «неинтересные» гетероатомы — не считаем их лигандами
IONS = {'SO4', 'GOL', 'CL', 'NA', 'MG', 'NI', 'ZN', 'EDO', 'PO4', 'ACT', 'DMS', 'PEG',
        'TRS', 'CA', 'K', 'IOD', 'MPD', 'FMT', 'BME', 'CIT', 'NO3', '1PE', 'PG4', 'BR',
        'MN', 'SCN', 'AZI', 'ACY', 'UNX', 'HOH', 'IMD', 'TLA', 'MES', 'EPE', 'FLC', 'CD',
        'CO', 'CU', 'GLC', 'NH4', 'PGE', 'P6G', '2PE', 'XPE', 'OXL', 'SIN', 'NHE', 'B3P',
        'MRD', 'DTT', 'TCE', 'PGO', 'SGM', 'BU3', 'PG0', '12P', '15P'}


def search_by_uniprot(acc):
    """Список PDB ID, у которых есть ссылка на указанную запись UniProt."""
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_name",
            "operator": "exact_match", "value": "UniProt"}}]},
        "return_type": "entry", "request_options": {"return_all_hits": True}}
    r = json.load(urllib.request.urlopen(SEARCH + urllib.parse.quote(json.dumps(q)), timeout=120))
    return [x['identifier'] for x in r['result_set']]


def fetch_meta(ids):
    """Метаданные для списка PDB ID одним GraphQL-запросом."""
    gql = '''{ entries(entry_ids:%s){ rcsb_id exptl{method} rcsb_entry_info{resolution_combined}
      struct{title}
      polymer_entities{ rcsb_polymer_entity{pdbx_description pdbx_mutation}
        entity_poly{pdbx_strand_id pdbx_seq_one_letter_code_can}
        rcsb_entity_source_organism{scientific_name}
        rcsb_polymer_entity_align{reference_database_accession aligned_regions{ref_beg_seq_id length}}}
      nonpolymer_entities{ nonpolymer_comp{chem_comp{id name formula_weight}}}}}''' % json.dumps(ids)
    req = urllib.request.Request(GRAPHQL, data=json.dumps({'query': gql}).encode(),
                                 headers={'Content-Type': 'application/json'})
    return json.load(urllib.request.urlopen(req, timeout=300))['data']['entries']


def form_by_sequence(seq):
    """WT / T315I / None — по мотиву вокруг gatekeeper-остатка."""
    if 'FYIITEF' in seq:
        return 'WT'
    if 'FYIIIEF' in seq:
        return 'T315I'
    return None


def main():
    ids = sorted(set(search_by_uniprot('P00519')) | set(search_by_uniprot('P00520')))
    print('найдено записей со ссылкой на ABL1/Abl1:', len(ids))

    entries = []
    for i in range(0, len(ids), 40):          # запрос разбиваем на порции
        entries += fetch_meta(ids[i:i + 40])

    os.makedirs(DATA, exist_ok=True)
    json.dump(entries, open(os.path.join(DATA, 'abl_entries_raw.json'), 'w'), indent=1)

    rows = []
    for e in entries:
        form, org, muts = None, [], []
        for p in e['polymer_entities']:
            seq = (p['entity_poly']['pdbx_seq_one_letter_code_can'] or '').replace('\n', '')
            f = form_by_sequence(seq)
            if f:
                form = f if form in (None, f) else form
            org += [o['scientific_name'] for o in (p['rcsb_entity_source_organism'] or [])
                    if o.get('scientific_name')]
            if p['rcsb_polymer_entity'].get('pdbx_mutation'):
                muts.append(p['rcsb_polymer_entity']['pdbx_mutation'])
        ligs = [n['nonpolymer_comp']['chem_comp']['id'] for n in (e['nonpolymer_entities'] or [])]
        drug = [n['nonpolymer_comp']['chem_comp']['id'] for n in (e['nonpolymer_entities'] or [])
                if n['nonpolymer_comp']['chem_comp']['id'] not in IONS
                and (n['nonpolymer_comp']['chem_comp']['formula_weight'] or 0) >= 200]
        rows.append({'id': e['rcsb_id'], 'title': e['struct']['title'],
                     'method': e['exptl'][0]['method'],
                     'res': (e['rcsb_entry_info']['resolution_combined'] or [None])[0],
                     'form': form, 'mut': '; '.join(muts), 'org': sorted(set(org)),
                     'ligs': ligs, 'drug': drug})
    json.dump(rows, open(os.path.join(DATA, 'abl_rows.json'), 'w'), indent=1)

    kd = [r for r in rows if r['form']]        # записи, покрывающие остаток 315
    print('покрывают позицию 315 (киназный домен):', len(kd))
    print('  WT   :', sum(1 for r in kd if r['form'] == 'WT'))
    print('  T315I:', sum(1 for r in kd if r['form'] == 'T315I'))
    print('  с лигандом / без:', sum(1 for r in kd if r['drug']), '/', sum(1 for r in kd if not r['drug']))
    print('\nвсе структуры T315I:')
    for r in kd:
        if r['form'] == 'T315I':
            print('  %-5s %-16s %4s A  лиганды: %-12s %s'
                  % (r['id'], ','.join(r['org']), r['res'], ','.join(r['drug']), r['mut']))
    print('\nкандидаты WT (человек, есть лиганд, лучшее разрешение):')
    wt = [r for r in kd if r['form'] == 'WT' and r['drug']
          and any(o.lower() == 'homo sapiens' for o in r['org'])]
    for r in sorted(wt, key=lambda r: r['res'] or 99)[:10]:
        print('  %-5s %4s A  %s' % (r['id'], r['res'], ','.join(r['drug'])))


if __name__ == '__main__':
    main()
