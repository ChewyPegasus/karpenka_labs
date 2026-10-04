# -*- coding: utf-8 -*-
"""
Лабораторная работа №2. Иллюстрации по реальным координатам структур.

Строит из файлов 2HYY / 3QRJ то же, что показывают команды ChimeraX
(cartoon, style stick, matchmaker), но средствами matplotlib:

  fig_wt_fold.png        ход цепи киназного домена WT + иматиниб + Thr315
  fig_site_wt.png        сайт связывания 2HYY: иматиниб и Thr315
  fig_site_t315i.png     сайт связывания 3QRJ: ребастиниб и Ile315
  fig_315_overlay.png    Thr315 и Ile315 после наложения, крупно
  fig_superposition.png  наложение C-альфа следов WT и T315I
  fig_deviation.png      отклонение C-альфа по остаткам после наложения

Запуск:  python scripts/make_figures.py   (из папки lab2)
"""
import json
import os

import gemmi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
LAB2 = os.path.dirname(HERE)
STRUCT = os.path.join(LAB2, 'structures')
DATA = os.path.join(LAB2, 'data')
FIG = os.path.join(LAB2, 'figures')

# Категориальные слоты 1 и 2 проверенной палитры
C_WT = '#2a78d6'       # синий   — WT / 2HYY
C_MUT = '#eb6834'      # оранжевый — T315I / 3QRJ
INK = '#0b0b0b'
INK2 = '#52514e'
GRID = '#d8d7d2'

GATEKEEPER = 315
KINASE_DOMAIN = (242, 493)
ACT_LOOP = (381, 405)

# Цвета элементов (CPK), углерод подкрашиваем под структуру
ELEM = {'N': '#2f5fd0', 'O': '#d93025', 'S': '#e0a800', 'F': '#1baf7a',
        'CL': '#1baf7a', 'BR': '#8b4513', 'P': '#eda100'}

plt.rcParams.update({
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'savefig.facecolor': 'white', 'font.size': 9,
    'text.color': INK, 'axes.labelcolor': INK,
    'xtick.color': INK2, 'ytick.color': INK2,
})


def load_chain(cif_name, chain_name):
    st = gemmi.read_structure(os.path.join(STRUCT, cif_name))
    st.setup_entities()
    return st, st[0][chain_name]


def is_aa(res):
    info = gemmi.find_tabulated_residue(res.name)
    return info is not None and info.is_amino_acid()


def ca_trace(chain):
    nums, pts = [], []
    for res in chain:
        if not is_aa(res):
            continue
        ca = res.find_atom('CA', '*')
        if ca is not None:
            nums.append(res.seqid.num)
            pts.append([ca.pos.x, ca.pos.y, ca.pos.z])
    return np.array(nums), np.array(pts)


def atoms_of(res, skip_h=True):
    out = []
    for a in res:
        el = a.element.name.upper()
        if skip_h and el == 'H':
            continue
        out.append((a.name, el, np.array([a.pos.x, a.pos.y, a.pos.z])))
    return out


def draw_sticks(ax, atoms, carbon_color, lw=2.0, max_bond=1.95, markersize=0):
    """Рисует связи как отрезки, окрашивая половинки по элементам."""
    n = len(atoms)
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(atoms[i][2] - atoms[j][2])
            if d > max_bond:
                continue
            mid = (atoms[i][2] + atoms[j][2]) / 2
            for k in (i, j):
                col = ELEM.get(atoms[k][1], carbon_color)
                seg = np.array([atoms[k][2], mid])
                ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color=col,
                        lw=lw, solid_capstyle='round', zorder=3)
    if markersize:
        for name, el, p in atoms:
            col = ELEM.get(el, carbon_color)
            ax.plot([p[0]], [p[1]], [p[2]], marker='o', ms=markersize,
                    color=col, zorder=4)


def tidy3d(ax, pts, pad=1.5, zoom=1.55):
    """Равные масштабы, без рамок и осей — как в молекулярном вьюере."""
    pts = np.asarray(pts)
    c = pts.mean(axis=0)
    r = max(np.ptp(pts, axis=0).max() / 2, 1.0) + pad
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_box_aspect((1, 1, 1), zoom=zoom)
    ax.set_axis_off()
    ax.set_position([0.0, 0.0, 1.0, 1.0])


def fig_fold():
    """Ход цепи киназного домена WT, лиганд и Thr315."""
    st, ch = load_chain('2hyy.cif', 'A')
    nums, pts = ca_trace(ch)
    lig = next(r for r in ch if r.name == 'STI')
    gk = next(r for r in ch if r.seqid.num == GATEKEEPER)

    fig = plt.figure(figsize=(7.0, 5.2))
    ax = fig.add_subplot(111, projection='3d')

    # ход цепи: разрывы нумерации не соединяем
    seg = [0]
    for i in range(1, len(nums)):
        if nums[i] != nums[i - 1] + 1:
            ax.plot(pts[seg, 0], pts[seg, 1], pts[seg, 2], color=C_WT,
                    lw=1.6, alpha=.85, zorder=2)
            seg = []
        seg.append(i)
    ax.plot(pts[seg, 0], pts[seg, 1], pts[seg, 2], color=C_WT, lw=1.6,
            alpha=.85, zorder=2)

    draw_sticks(ax, atoms_of(lig), '#444444', lw=2.4)
    draw_sticks(ax, atoms_of(gk), C_MUT, lw=3.0, markersize=3)

    tidy3d(ax, pts, zoom=1.5)
    ax.view_init(elev=16, azim=-62)
    fig.suptitle('2HYY, цепь A: киназный домен ABL1 (WT)\n'
                 'серым — иматиниб (STI), оранжевым — Thr315',
                 fontsize=10, color=INK, y=0.99)
    fig.savefig(os.path.join(FIG, 'fig_wt_fold.png'), dpi=190,
                bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)


def fig_site(cif, chain_name, lig_code, title, out, carbon, site_nums):
    st, ch = load_chain(cif, chain_name)
    lig = next(r for r in ch if r.name == lig_code)
    gk = next(r for r in ch if r.seqid.num == GATEKEEPER)

    fig = plt.figure(figsize=(7.0, 5.0))
    ax = fig.add_subplot(111, projection='3d')

    allpts = []
    # окружение — тонкими линиями
    for res in ch:
        if not is_aa(res) or res.seqid.num not in site_nums:
            continue
        if res.seqid.num == GATEKEEPER:
            continue
        at = atoms_of(res)
        draw_sticks(ax, at, '#b9b8b3', lw=1.0)
        allpts += [a[2] for a in at]

    la = atoms_of(lig)
    draw_sticks(ax, la, carbon, lw=2.6)
    allpts += [a[2] for a in la]

    ga = atoms_of(gk)
    draw_sticks(ax, ga, C_MUT, lw=3.4, markersize=3.5)
    allpts += [a[2] for a in ga]

    cg = np.mean([a[2] for a in ga], axis=0)
    ax.text(cg[0], cg[1], cg[2] + 1.6, '%s%d' % (gk.name.capitalize(), GATEKEEPER),
            color=C_MUT, fontsize=11, weight='bold', ha='center')

    tidy3d(ax, np.array(allpts), pad=0.6, zoom=1.6)
    ax.view_init(elev=14, azim=-70)
    fig.suptitle(title, fontsize=10, color=INK, y=0.99)
    fig.savefig(os.path.join(FIG, out), dpi=190,
                bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)


def fig_overlay_315(sup):
    """Thr315 и Ile315 крупно, после наложения."""
    st_w, ch_w = load_chain('2hyy.cif', 'A')
    stm = gemmi.read_structure(os.path.join(STRUCT, 'T315I_superposed_on_WT.pdb'))
    stm.setup_entities()
    ch_m = stm[0]['A']

    gw = next(r for r in ch_w if r.seqid.num == GATEKEEPER)
    gm = next(r for r in ch_m if r.seqid.num == GATEKEEPER)
    lw_ = next(r for r in ch_w if r.name == 'STI')

    fig = plt.figure(figsize=(7.0, 4.8))
    ax = fig.add_subplot(111, projection='3d')

    aw, am = atoms_of(gw), atoms_of(gm)
    draw_sticks(ax, aw, C_WT, lw=3.4, markersize=4)
    draw_sticks(ax, am, C_MUT, lw=3.4, markersize=4)

    # ближайшая часть иматиниба — чтобы было видно, куда смотрит боковая цепь
    gc = np.mean([a[2] for a in aw], axis=0)
    near_lig = [a for a in atoms_of(lw_) if np.linalg.norm(a[2] - gc) < 8]
    draw_sticks(ax, near_lig, '#b9b8b3', lw=1.6)

    pts = [a[2] for a in aw + am + near_lig]
    tidy3d(ax, np.array(pts), pad=0.5, zoom=1.5)
    ax.view_init(elev=12, azim=-78)
    fig.legend(handles=[Line2D([], [], color=C_WT, lw=3, label='Thr315 (2HYY, WT)'),
                        Line2D([], [], color=C_MUT, lw=3, label='Ile315 (3QRJ, T315I)'),
                        Line2D([], [], color='#b9b8b3', lw=2,
                               label='иматиниб (ближняя часть)')],
               loc='upper left', frameon=False, fontsize=9,
               bbox_to_anchor=(0.02, 0.97))
    fig.suptitle('Остаток 315 после наложения: отклонение C$\\alpha$ %.2f Å'
                 % sup['ca_deviation_315'], fontsize=10, color=INK, y=0.995)
    fig.savefig(os.path.join(FIG, 'fig_315_overlay.png'), dpi=190,
                bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)


def fig_superposition(sup):
    st_w, ch_w = load_chain('2hyy.cif', 'A')
    stm = gemmi.read_structure(os.path.join(STRUCT, 'T315I_superposed_on_WT.pdb'))
    stm.setup_entities()
    ch_m = stm[0]['A']

    nw, pw = ca_trace(ch_w)
    nm, pm = ca_trace(ch_m)

    fig = plt.figure(figsize=(7.2, 5.4))
    ax = fig.add_subplot(111, projection='3d')
    for nums, pts, col in ((nw, pw, C_WT), (nm, pm, C_MUT)):
        seg = [0]
        for i in range(1, len(nums)):
            if nums[i] != nums[i - 1] + 1:
                ax.plot(pts[seg, 0], pts[seg, 1], pts[seg, 2], color=col, lw=1.7,
                        alpha=.9)
                seg = []
            seg.append(i)
        ax.plot(pts[seg, 0], pts[seg, 1], pts[seg, 2], color=col, lw=1.7, alpha=.9)

    tidy3d(ax, np.vstack([pw, pm]), zoom=1.5)
    ax.view_init(elev=16, azim=-62)
    fig.legend(handles=[Line2D([], [], color=C_WT, lw=2.5, label='2HYY, WT'),
                        Line2D([], [], color=C_MUT, lw=2.5, label='3QRJ, T315I')],
               loc='lower left', frameon=False, fontsize=9,
               bbox_to_anchor=(0.03, 0.03))
    fig.suptitle('Наложение C$\\alpha$-следов: RMSD %.3f Å по %d парам\n'
                 '(после итеративного отбрасывания; %.3f Å по всем %d парам)'
                 % (sup['rmsd_pruned'], sup['n_pairs_pruned'],
                    sup['rmsd_all_pairs'], sup['n_pairs_initial']),
                 fontsize=10, color=INK, y=0.995)
    fig.savefig(os.path.join(FIG, 'fig_superposition.png'), dpi=190,
                bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)


def fig_deviation(sup):
    dev = {int(k): v for k, v in sup['per_residue_deviation'].items()}
    xs = np.array(sorted(dev))
    ys = np.array([dev[x] for x in xs])

    fig, ax = plt.subplots(figsize=(7.4, 3.3))
    ax.axvspan(*ACT_LOOP, color='#f0efe9', zorder=0)
    ax.plot(xs, ys, color=C_WT, lw=1.6, zorder=3)
    ax.axhline(sup['prune_cutoff'], color=INK2, lw=1, ls='--', zorder=2)
    ax.axvline(GATEKEEPER, color=C_MUT, lw=1.4, zorder=2)

    ax.annotate('петля активации\n381–405', xy=(np.mean(ACT_LOOP), max(ys) * .93),
                ha='center', fontsize=8, color=INK2)
    ax.annotate('остаток 315\n%.2f Å' % sup['ca_deviation_315'],
                xy=(GATEKEEPER, sup['ca_deviation_315']),
                xytext=(GATEKEEPER - 58, max(ys) * .55), fontsize=8, color=C_MUT,
                arrowprops=dict(arrowstyle='->', color=C_MUT, lw=1))
    ax.annotate('порог отбрасывания %.1f Å' % sup['prune_cutoff'],
                xy=(xs[2], sup['prune_cutoff']), xytext=(xs[2], sup['prune_cutoff'] + .45),
                fontsize=8, color=INK2)

    ax.set_xlabel('номер остатка (нумерация UniProt P00519)')
    ax.set_ylabel('отклонение C$\\alpha$, Å')
    ax.set_title('Локализация различий между 2HYY (WT) и 3QRJ (T315I)',
                 fontsize=10, color=INK)
    ax.grid(axis='y', color=GRID, lw=.7)
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(GRID)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig_deviation.png'), dpi=190)
    plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    with open(os.path.join(DATA, 'analysis.json'), encoding='utf-8') as fh:
        an = json.load(fh)
    with open(os.path.join(DATA, 'superposition.json'), encoding='utf-8') as fh:
        sup = json.load(fh)

    site_wt = {r[0] for r in an['structures']['WT']['binding_site_residues']}
    site_mt = {r[0] for r in an['structures']['T315I']['binding_site_residues']}

    fig_fold()
    print('fig_wt_fold.png')
    fig_site('2hyy.cif', 'A', 'STI',
             '2HYY: сайт связывания, иматиниб (STI) и Thr315',
             'fig_site_wt.png', '#4a3aa7', site_wt)
    print('fig_site_wt.png')
    fig_site('3qrj.cif', 'A', '919',
             '3QRJ: сайт связывания, ребастиниб (919) и Ile315',
             'fig_site_t315i.png', '#008300', site_mt)
    print('fig_site_t315i.png')
    fig_overlay_315(sup)
    print('fig_315_overlay.png')
    fig_superposition(sup)
    print('fig_superposition.png')
    fig_deviation(sup)
    print('fig_deviation.png')


if __name__ == '__main__':
    main()
