# Fly Crossing Planner

Desktop tool for Drosophila (fruit fly) crosses: type parental genotypes in plain text, list viable F1 genotypes as fractions, and export a one-page PowerPoint crossing scheme.

## Install

```bash
cd CrossingPlanner
python -m pip install -r requirements.txt
```

Requires Python 3.10+ and tkinter (usually included with Python on Windows).

## Run

```bash
python -m crossing_planner
```

## Genotype syntax

- Separate chromosome pairs with `;` (order: sex chromosomes; chr2; chr3; …)
- Separate homologs with `/`
- Multiple loci on one chromatid: spaces (`a b c / d e f`)
- Wild type: `+`
- Males: write X-linked alleles as `allele / Y`, or hemizygous shorthand `w` (expanded to `w / Y`)

### Examples

Female:

```text
w ; Sp / CyO ; TM3, Sb / MKRS
```

Male:

```text
w / Y ; + / CyO ; + / TM6B, Tb
```

Linked loci with recombination (set RF in the GUI, e.g. `0.1`):

```text
♀  a b / + +
♂  a b / Y
```

## Features

- Independent assortment of chromosome pairs
- Sex chromosomes (XX / XY)
- Simple adjacent-locus recombination (global RF, default `0` = complete linkage)
- Males do not recombine. A whole-chromosome balancer on a female chromosome also blocks recombination there. MRS and MKRS only balance proximal chromosome 3, so other loci still use RF
- Marker catalog (bristles, wings, appendages, eye shape, eye color, body) and balancer recognition (FM, CyO, SM, TM, MRS/MKRS, CyO-TM9)
- Homozygous lethal classes dropped; FM7c/FM7c females stay in the ratio and are marked sterile
- F1 classes colored by distinguishable phenotype. The same visible traits share a color
- Fraction and compact genotype views
- Parallel waves and multi-step crosses in one window
- PowerPoint export of the genotypes selected in the results, using the current view checkboxes

## GUI

1. Add stocks, then add one or more crosses. Crosses in the same wave run in parallel
2. Optionally set each cross's RF and any extra lethal alleles
3. Click **计算方案**
4. Click **导出方案 PPT…** to save the plan
