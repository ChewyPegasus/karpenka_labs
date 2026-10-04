# -*- coding: utf-8 -*-
"""
Лабораторная работа №2. Наложение структур WT (2HYY) и T315I (3QRJ).

Вычислительный эквивалент ChimeraX Matchmaker:
  matchmaker #2/A to #1/A pairing ss

Пары остатков строятся по номерам: в лабораторной работе №1 и в
prepare_structures.py проверено, что обе структуры используют одну и ту же
нумерацию изоформы 1a (сдвиг относительно UniProt P00519 равен 0), поэтому
остатки с одинаковым номером — действительно соответственные.

Matchmaker после первичного совмещения выполняет итеративное отбрасывание
(iterative pruning) пар, расходящихся больше чем на cutoff, и печатает RMSD
до и после. Здесь воспроизведена та же процедура.

Результат: data/superposition.json + structures/T315I_superposed_on_WT.pdb

Запуск:  python scripts/superpose.py   (из папки lab2)
"""
import json
import os

import gemmi
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LAB2 = os.path.dirname(HERE)
STRUCT = os.path.join(LAB2, 'structures')
DATA = os.path.join(LAB2, 'data')

REF = ('WT', '2hyy', 'A')        # модель #1 — неподвижная
MOB = ('T315I', '3qrj', 'A')     # модель #2 — совмещаемая

PRUNE_CUTOFF = 2.0               # A, как в ChimeraX Matchmaker по умолчанию
GATEKEEPER = 315
KINASE_DOMAIN = (242, 493)


def ca_atoms(path, chain_name):
    """{номер остатка: (имя, координаты CA)} для цепи."""
    st = gemmi.read_structure(path)
    st.setup_entities()
    chain = st[0][chain_name]
    out = {}
    for res in chain:
        info = gemmi.find_tabulated_residue(res.name)
        if info is None or not info.is_amino_acid():
            continue
        ca = res.find_atom('CA', '*')
        if ca is not None:
            out[res.seqid.num] = (res.name, np.array([ca.pos.x, ca.pos.y, ca.pos.z]))
    return st, out


def kabsch(P, Q):
    """Поворот+сдвиг, совмещающие P с Q. Возвращает (R, t)."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def rmsd(P, Q):
    return float(np.sqrt(((P - Q) ** 2).sum(axis=1).mean()))


def main():
    os.makedirs(DATA, exist_ok=True)

    ref_label, ref_id, ref_chain = REF
    mob_label, mob_id, mob_chain = MOB
    st_ref, ca_ref = ca_atoms(os.path.join(STRUCT, ref_id + '.cif'), ref_chain)
    st_mob, ca_mob = ca_atoms(os.path.join(STRUCT, mob_id + '.cif'), mob_chain)

    common = sorted(set(ca_ref) & set(ca_mob))
    print('CA в %s/%s: %d, в %s/%s: %d, общих номеров: %d'
          % (ref_id.upper(), ref_chain, len(ca_ref),
             mob_id.upper(), mob_chain, len(ca_mob), len(common)))

    Q = np.array([ca_ref[n][1] for n in common])   # неподвижная
    P = np.array([ca_mob[n][1] for n in common])   # подвижная

    # --- первичное совмещение по всем парам ---
    R, t = kabsch(P, Q)
    P_fit = (R @ P.T).T + t
    rmsd_all = rmsd(P_fit, Q)
    print('RMSD по всем %d парам: %.3f A' % (len(common), rmsd_all))

    # --- итеративное отбрасывание (как в Matchmaker) ---
    keep = np.ones(len(common), dtype=bool)
    for it in range(1, 21):
        R, t = kabsch(P[keep], Q[keep])
        P_fit = (R @ P.T).T + t
        dev = np.linalg.norm(P_fit - Q, axis=1)
        new_keep = dev <= PRUNE_CUTOFF
        if new_keep.sum() < 3 or (new_keep == keep).all():
            keep = new_keep if new_keep.sum() >= 3 else keep
            break
        keep = new_keep
    R, t = kabsch(P[keep], Q[keep])
    P_fit = (R @ P.T).T + t
    dev = np.linalg.norm(P_fit - Q, axis=1)
    rmsd_pruned = rmsd(P_fit[keep], Q[keep])
    print('после итеративного отбрасывания (cutoff %.1f A, %d итер.): '
          'RMSD %.3f A по %d парам' % (PRUNE_CUTOFF, it, rmsd_pruned, int(keep.sum())))

    # --- локальные оценки ---
    kd_idx = [i for i, n in enumerate(common) if KINASE_DOMAIN[0] <= n <= KINASE_DOMAIN[1]]
    rmsd_kd = rmsd(P_fit[kd_idx], Q[kd_idx])

    near = [i for i, n in enumerate(common) if abs(n - GATEKEEPER) <= 5]
    rmsd_gk = rmsd(P_fit[near], Q[near])

    gk_i = common.index(GATEKEEPER)
    gk_dev = float(dev[gk_i])

    # 10 самых расходящихся участков
    worst = sorted(zip(common, dev), key=lambda x: -x[1])[:12]

    print('RMSD внутри киназного домена (%d-%d): %.3f A по %d парам'
          % (KINASE_DOMAIN[0], KINASE_DOMAIN[1], rmsd_kd, len(kd_idx)))
    print('RMSD в окрестности остатка %d (+-5): %.3f A по %d парам'
          % (GATEKEEPER, rmsd_gk, len(near)))
    print('отклонение CA остатка %d: %.3f A (%s -> %s)'
          % (GATEKEEPER, gk_dev, ca_ref[GATEKEEPER][0], ca_mob[GATEKEEPER][0]))
    print('наибольшие расхождения:',
          ', '.join('%d:%.1f' % (n, d) for n, d in worst[:8]))

    # --- сохранить совмещённую подвижную структуру ---
    tr = gemmi.Transform()
    tr.mat.fromlist([list(map(float, row)) for row in R])
    tr.vec.fromlist([float(x) for x in t])
    for model in st_mob:
        for ch in model:
            for res in ch:
                for atom in res:
                    atom.pos = gemmi.Position(*tr.apply(atom.pos).tolist())
    out_pdb = os.path.join(STRUCT, 'T315I_superposed_on_WT.pdb')
    st_mob.setup_entities()
    st_mob.write_pdb(out_pdb)
    print('-> %s' % os.path.basename(out_pdb))

    result = {
        'reference': {'label': ref_label, 'pdb_id': ref_id.upper(), 'chain': ref_chain,
                      'n_ca': len(ca_ref)},
        'mobile': {'label': mob_label, 'pdb_id': mob_id.upper(), 'chain': mob_chain,
                   'n_ca': len(ca_mob)},
        'pairing': 'по номерам остатков (нумерация обеих структур совпадает с UniProt)',
        'n_pairs_initial': len(common),
        'rmsd_all_pairs': round(rmsd_all, 3),
        'prune_cutoff': PRUNE_CUTOFF,
        'iterations': it,
        'n_pairs_pruned': int(keep.sum()),
        'rmsd_pruned': round(rmsd_pruned, 3),
        'rmsd_kinase_domain': round(rmsd_kd, 3),
        'n_pairs_kinase_domain': len(kd_idx),
        'rmsd_around_315': round(rmsd_gk, 3),
        'ca_deviation_315': round(gk_dev, 3),
        'residue_315': {'ref': ca_ref[GATEKEEPER][0], 'mobile': ca_mob[GATEKEEPER][0]},
        'largest_deviations': [[int(n), round(float(d), 2)] for n, d in worst],
        'per_residue_deviation': {int(n): round(float(d), 3)
                                  for n, d in zip(common, dev)},
    }
    with open(os.path.join(DATA, 'superposition.json'), 'w', encoding='utf-8') as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)
    print('-> data/superposition.json записан')


if __name__ == '__main__':
    main()
