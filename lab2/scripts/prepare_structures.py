# -*- coding: utf-8 -*-
"""
Лабораторная работа №2. Подготовка структур BCR::ABL1 WT и T315I.

Делает то же, что последовательность команд ChimeraX из задания
(open / info / select / save), но программно, через gemmi:

  * читает WT_original.cif (2HYY) и T315I_original.cif (3QRJ);
  * определяет цепь ABL1, диапазон остатков, пропуски (missing residues);
  * проверяет нумерацию относительно UniProt P00519 (результат лабораторной №1);
  * находит сокристаллизованный ингибитор, воду, ионы;
  * считает остатки сайта связывания (в радиусе 4.5 A от лиганда);
  * сохраняет рабочие файлы *_complex.pdb (белок + ингибитор) и *_protein.pdb.

Результаты складываются в lab2/data/analysis.json и печатаются в консоль.

Запуск:  python scripts/prepare_structures.py   (из папки lab2)
"""
import json
import os

import gemmi

HERE = os.path.dirname(os.path.abspath(__file__))
LAB2 = os.path.dirname(HERE)
STRUCT = os.path.join(LAB2, 'structures')
DATA = os.path.join(LAB2, 'data')
UNIPROT_JSON = os.path.join(LAB2, '..', 'lab1', 'data', 'uniprot_P00519.json')

# Структуры, выбранные в лабораторной работе №1
TARGETS = [
    # метка,   PDB ID, цепь, код лиганда-ингибитора
    ('WT',    '2hyy', 'A', 'STI'),
    ('T315I', '3qrj', 'A', '919'),
]

KINASE_DOMAIN = (242, 493)   # UniProt P00519, protein kinase domain
GATEKEEPER = 315             # Thr315 / Ile315
CONTACT_CUTOFF = 4.5         # A, радиус для остатков сайта связывания

# Ионы и типичные добавки кристаллизации — не являются ингибитором
NON_LIGAND_HET = {
    'HOH', 'DOD', 'NA', 'K', 'MG', 'CA', 'MN', 'ZN', 'NI', 'CL', 'BR', 'IOD',
    'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'MPD', 'ACT', 'DMS', 'TRS', 'FMT',
}


def uniprot_sequence():
    """Последовательность P00519 из выгрузки лабораторной работы №1."""
    with open(UNIPROT_JSON, encoding='utf-8') as fh:
        entry = json.load(fh)
    return entry['sequence']['value']


def one_letter(res_name):
    code = gemmi.find_tabulated_residue(res_name)
    if code is None or not code.is_amino_acid():
        return None
    return gemmi.find_tabulated_residue(res_name).one_letter_code.upper()


def modelled_residues(chain):
    """Список (номер, имя) реально разрешённых аминокислотных остатков."""
    out = []
    for res in chain:
        if one_letter(res.name) is not None:
            out.append((res.seqid.num, res.name))
    return out


def find_gaps(residues):
    """Внутренние пропуски нумерации в разрешённой части цепи."""
    gaps = []
    for (n1, _), (n2, _) in zip(residues, residues[1:]):
        if n2 > n1 + 1:
            gaps.append((n1 + 1, n2 - 1))
    return gaps


def numbering_offset(residues, uniprot_seq):
    """
    Подбирает сдвиг между нумерацией PDB и UniProt:
    uniprot_position = pdb_seqid + offset. Возвращает (offset, доля совпадений).
    """
    best = (0, -1.0)
    for offset in range(-40, 41):
        hit = total = 0
        for num, name in residues:
            pos = num + offset
            if 1 <= pos <= len(uniprot_seq):
                total += 1
                if one_letter(name) == uniprot_seq[pos - 1]:
                    hit += 1
        if total and hit / total > best[1]:
            best = (offset, hit / total)
    return best


def ligand_residues(chain, ligand_code):
    return [res for res in chain if res.name == ligand_code]


def het_inventory(model):
    """Что ещё есть в структуре помимо белка: вода, ионы, добавки."""
    inv = {}
    for ch in model:
        for res in ch:
            if one_letter(res.name) is None:
                inv[res.name] = inv.get(res.name, 0) + 1
    return inv


def contact_residues(chain, ligand, cutoff):
    """Остатки белка, у которых есть атом ближе cutoff к любому атому лиганда."""
    lig_pos = [atom.pos for atom in ligand]
    hits = []
    for res in chain:
        if one_letter(res.name) is None:
            continue
        best = None
        for atom in res:
            for lp in lig_pos:
                d = atom.pos.dist(lp)
                if best is None or d < best:
                    best = d
        if best is not None and best <= cutoff:
            hits.append((res.seqid.num, res.name, round(best, 2)))
    return hits


def write_selection(src_struct, chain_name, residues_keep, out_path, title):
    """Собирает новую структуру из выбранных остатков и пишет .pdb."""
    out = gemmi.Structure()
    out.name = title
    out.spacegroup_hm = src_struct.spacegroup_hm
    out.cell = src_struct.cell
    model = gemmi.Model('1')
    chain = gemmi.Chain(chain_name)
    for res in residues_keep:
        chain.add_residue(res)
    model.add_chain(chain)
    out.add_model(model)
    out.setup_entities()
    out.write_pdb(out_path)
    n_atoms = sum(len(r) for r in residues_keep)
    return len(residues_keep), n_atoms


def main():
    os.makedirs(DATA, exist_ok=True)
    useq = uniprot_sequence()
    print('UniProt P00519: %d а.о. (из лабораторной №1)' % len(useq))
    print('Остаток %d в UniProt: %s' % (GATEKEEPER, useq[GATEKEEPER - 1]))

    report = {'uniprot_length': len(useq),
              'uniprot_residue_315': useq[GATEKEEPER - 1],
              'kinase_domain': list(KINASE_DOMAIN),
              'structures': {}}

    for label, pdbid, chain_name, lig_code in TARGETS:
        src = os.path.join(STRUCT, pdbid + '.cif')
        original = os.path.join(STRUCT, '%s_original.cif' % label)
        st = gemmi.read_structure(src)
        st.setup_entities()
        # «WT_original.cif» / «T315I_original.cif» из задания
        st.make_mmcif_document().write_file(original)

        model = st[0]
        chain = model[chain_name]
        res_list = modelled_residues(chain)
        gaps = find_gaps(res_list)
        offset, score = numbering_offset(res_list, useq)

        # остаток в позиции gatekeeper (с учётом найденного сдвига)
        want = GATEKEEPER - offset
        gk = [(n, nm) for n, nm in res_list if n == want]

        ligs = ligand_residues(chain, lig_code)
        inv = het_inventory(model)

        contacts = contact_residues(chain, ligs[0], CONTACT_CUTOFF) if ligs else []

        # пропуски внутри киназного домена и внутри сайта связывания
        kd_lo, kd_hi = [x - offset for x in KINASE_DOMAIN]
        gaps_in_kd = [g for g in gaps if g[1] >= kd_lo and g[0] <= kd_hi]
        contact_nums = {c[0] for c in contacts}
        gaps_in_site = [g for g in gaps
                        if any(g[0] - 1 <= n <= g[1] + 1 for n in contact_nums)]

        # --- сохранение рабочих файлов ---
        protein_res = [r for r in chain if one_letter(r.name) is not None]
        complex_res = protein_res + list(ligs)
        p_complex = os.path.join(STRUCT, '%s_complex.pdb' % label)
        p_protein = os.path.join(STRUCT, '%s_protein.pdb' % label)
        nres_c, nat_c = write_selection(st, chain_name, complex_res, p_complex,
                                        '%s_complex' % label)
        nres_p, nat_p = write_selection(st, chain_name, protein_res, p_protein,
                                        '%s_protein' % label)

        info = {
            'label': label,
            'pdb_id': pdbid.upper(),
            'title': st.name,
            'method': 'X-RAY DIFFRACTION',
            'resolution': st.resolution,
            'chains_in_model': [ch.name for ch in model],
            'chain_used': chain_name,
            'n_residues_modelled': len(res_list),
            'range_modelled': [res_list[0][0], res_list[-1][0]],
            'gaps': [list(g) for g in gaps],
            'numbering_offset_to_uniprot': offset,
            'numbering_match_fraction': round(score, 4),
            'gatekeeper_seqid_in_pdb': want,
            'gatekeeper_residue': gk[0][1] if gk else None,
            'ligand_code': lig_code,
            'ligand_copies_in_chain': len(ligs),
            'ligand_atoms': len(ligs[0]) if ligs else 0,
            'het_inventory': inv,
            'binding_site_residues': contacts,
            'gaps_in_kinase_domain': [list(g) for g in gaps_in_kd],
            'gaps_in_binding_site': [list(g) for g in gaps_in_site],
            'files': {
                'original': os.path.basename(original),
                'complex': os.path.basename(p_complex),
                'protein': os.path.basename(p_protein),
            },
            'counts': {
                'complex_residues': nres_c, 'complex_atoms': nat_c,
                'protein_residues': nres_p, 'protein_atoms': nat_p,
            },
        }
        report['structures'][label] = info

        print('\n' + '=' * 72)
        print('%s = %s  (%s)' % (label, pdbid.upper(), st.name))
        print('  разрешение: %.2f A; цепи в модели: %s; рабочая цепь: %s'
              % (st.resolution, ', '.join(info['chains_in_model']), chain_name))
        print('  разрешено остатков: %d, диапазон %d-%d'
              % (len(res_list), res_list[0][0], res_list[-1][0]))
        print('  пропуски: %s' % (gaps if gaps else 'нет'))
        print('  сдвиг нумерации PDB->UniProt: %+d (совпадение %.1f%%)'
              % (offset, 100 * score))
        print('  остаток %d: %s (в PDB номер %d)'
              % (GATEKEEPER, info['gatekeeper_residue'], want))
        print('  лиганд %s: %d копия(и) в цепи, %d атомов'
              % (lig_code, len(ligs), info['ligand_atoms']))
        print('  прочее в структуре: %s' % inv)
        print('  остатков сайта связывания (<= %.1f A): %d'
              % (CONTACT_CUTOFF, len(contacts)))
        print('  пропуски внутри киназного домена: %s'
              % (gaps_in_kd if gaps_in_kd else 'нет'))
        print('  сохранено: %s (%d ост., %d атомов), %s (%d ост., %d атомов)'
              % (os.path.basename(p_complex), nres_c, nat_c,
                 os.path.basename(p_protein), nres_p, nat_p))

    with open(os.path.join(DATA, 'analysis.json'), 'w', encoding='utf-8') as fh:
        json.dump(report, fh, ensure_ascii=False, indent=2)
    print('\n-> data/analysis.json записан')


if __name__ == '__main__':
    main()
