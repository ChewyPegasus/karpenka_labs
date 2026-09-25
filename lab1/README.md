# Лабораторная работа №1 — BCR::ABL1 (WT и T315I)

Шпока Владислав, Пясецкий Кирилл — 4 курс, 3 группа.

## Что где лежит

```
report/report.tex        отчёт по заданию (13 пунктов) -> report/report.pdf
guide/guide.tex          подробный разбор: термины, обоснования, ответы на защите -> guide/guide.pdf
BCR_ABL_dataset.csv      набор данных, 78 строк (одна строка = одно измерение)
BCR_ABL_dataset.xlsx     то же в Excel
figures/fig_activity.png график: активность ингибиторов против WT и T315I
scripts/                 скрипты выгрузки и обработки
data/                    сырые ответы UniProt / RCSB / ChEMBL / PubChem (JSON)
```

## Как пересобрать всё с нуля

```
python scripts/fetch_uniprot.py     # UniProt P00519 -> data/uniprot_P00519.json
python scripts/survey_pdb.py        # RCSB PDB       -> data/abl_rows.json
python scripts/fetch_chembl.py      # ChEMBL         -> data/chembl_activities.json
python scripts/fetch_pubchem.py     # PubChem        -> data/pubchem.json
python scripts/build_dataset.py     # -> BCR_ABL_dataset.csv / .xlsx
python scripts/plot_activity.py     # -> figures/fig_activity.png

cd report && pdflatex report.tex && pdflatex report.tex
cd guide  && pdflatex guide.tex  && pdflatex guide.tex
```

Нужны: Python 3.8+, `matplotlib`, `numpy`, `openpyxl` (для .xlsx), MiKTeX или TeX Live.
pdflatex запускается дважды — со второго прохода подставляются номера разделов и ссылок.

## Ключевые результаты

- Мишень: ABL1 (UniProt **P00519**, 1130 а.о.), киназный домен **242–493**, gatekeeper **Thr315**.
- Структуры для дальнейших работ: **2HYY** (WT, 2.40 Å, иматиниб) и **3QRJ** (T315I, 1.82 Å, ребастиниб).
- В PDB: 68 записей с киназным доменом — 60 WT и 8 T315I.
- Датасет: 78 измерений IC50/Kd, 7 соединений, 11 первоисточников.

## Сдача

Папку с материалами переименовать в `Group_N`, где N — порядковый номер в списке группы,
и положить в общую папку `Lab_1`.
