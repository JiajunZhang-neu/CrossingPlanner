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
- Simple adjacent-locus recombination (global RF, default `0` = complete linkage; male Drosophila: no recombination)
- Homozygous lethal exclusion (built-in balancers such as CyO, TM3, TM6B, FM7, … plus optional extra alleles)
- Dominant marker phenotype hints (genotypes remain the primary output)
- F1 ratios as exact `fractions.Fraction` values, re-normalized after removing lethals
- Fraction and compact genotype views
- Single-slide `.pptx` export

## GUI

1. Enter ♀ and ♂ genotypes
2. Optionally set RF and extra lethal alleles
3. Click **计算 F1**
4. Click **导出 PPT…** to save the crossing scheme
