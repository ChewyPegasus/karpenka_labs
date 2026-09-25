# -*- coding: utf-8 -*-
"""Шаг 6. График: активность ингибиторов против BCR::ABL1 WT и T315I.

Берёт собранный BCR_ABL_dataset.csv, для каждого соединения считает медиану
IC50 по всем измерениям отдельно для WT и для T315I и строит столбчатую
диаграмму в логарифмической шкале (иначе разницу в 4 порядка не показать).

Если среди измерений есть цензурированные (>5000, >10000 нМ), над столбцом
ставится знак «>»: истинное значение ещё выше.

Запуск:  python plot_activity.py
Результат: figures/fig_activity.png
"""
import csv
import os
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')

rows = list(csv.DictReader(open(os.path.join(ROOT, 'BCR_ABL_dataset.csv'),
                                encoding='utf-8-sig'), delimiter=';'))
COMPOUNDS = ['Imatinib', 'Nilotinib', 'Dasatinib', 'Bosutinib', 'Ponatinib', 'Axitinib']


def values(compound, mutation):
    """Список (значение, цензурировано ли) для IC50 данного соединения и формы."""
    out = []
    for r in rows:
        if r['Compound'] == compound and r['Mutation'] == mutation and r['Activity_type'] == 'IC50':
            v = r['Activity_value']
            out.append((float(v.lstrip('>=<')), v.startswith('>')))
    return out


fig, ax = plt.subplots(figsize=(9, 5))
width = 0.36
x = np.arange(len(COMPOUNDS))

for i, (mutation, color) in enumerate([('WT', '#3b76af'), ('T315I', '#d1603d')]):
    medians, censored = [], []
    for c in COMPOUNDS:
        v = values(c, mutation)
        medians.append(statistics.median([a for a, _ in v]) if v else np.nan)
        censored.append(any(b for _, b in v))
    bars = ax.bar(x + (i - 0.5) * width, medians, width, label='BCR::ABL1 ' + mutation,
                  color=color, edgecolor='black', linewidth=0.5)
    for rect, mv, cen in zip(bars, medians, censored):
        if np.isnan(mv):
            continue
        label = ('%.0f' % mv) if mv >= 10 else ('%g' % round(mv, 2))
        ax.text(rect.get_x() + rect.get_width() / 2, mv * 1.25,
                ('>' if cen else '') + label, ha='center', fontsize=8)

ax.set_yscale('log')
ax.set_ylabel(u'IC$_{50}$, нМ (медиана по всем измерениям, лог. шкала)')
ax.set_xticks(x)
ax.set_xticklabels(COMPOUNDS, rotation=15)
ax.set_title(u'Активность ингибиторов против BCR::ABL1 WT и T315I\n'
             u'(данные ChEMBL, n=%d измерений)' % len(rows))
ax.legend()
ax.grid(axis='y', ls=':', alpha=0.6)
ax.set_axisbelow(True)
fig.tight_layout()

out = os.path.join(ROOT, 'figures', 'fig_activity.png')
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=200)
print('сохранено:', out)

for c in COMPOUNDS:
    print('%-10s WT: %-28s T315I: %s'
          % (c, [v for v, _ in values(c, 'WT')], [v for v, _ in values(c, 'T315I')]))
