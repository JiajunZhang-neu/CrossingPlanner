"""Fraction-format rendering of F1 genotypes."""

from __future__ import annotations

from fractions import Fraction

from .cross import CrossResult, F1Individual
from .genotype import ChromosomePair, Genotype, format_genotype_compact


def _pair_compact(pair: ChromosomePair) -> str:
    mat, pat = pair.homologs()
    if mat.loci == pat.loci and not (mat.is_y or pat.is_y):
        # Still show as a/a for clarity in F1, or a/a
        return f"{mat.display()}/{pat.display()}"
    return f"{mat.display()}/{pat.display()}"


def format_individual_compact(ind: F1Individual) -> str:
    parts = [_pair_compact(p) for p in ind.pairs]
    sex_mark = "♀" if ind.sex == "female" else "♂"
    return f"{sex_mark} {'; '.join(parts)}  ({ind.probability})"


def format_pair_fraction(pair: ChromosomePair) -> tuple[str, str, str]:
    """
    Return (top, rule, bottom) strings for one chromosome pair,
    with equal width for stacked fraction.
    """
    mat, pat = pair.homologs()
    top = mat.display()
    bot = pat.display()
    width = max(len(top), len(bot), 3)
    top_c = top.center(width)
    bot_c = bot.center(width)
    rule = "-" * width
    return top_c, rule, bot_c


def format_individual_fraction(ind: F1Individual) -> str:
    """
    Multi-line fraction display for one F1 class:

      w        Sp         TM3, Sb
    ------ ; ------ ; ----------   (1/4) ♀
      Y        CyO          +
    """
    if not ind.pairs:
        return f"(empty)  ({ind.probability})"

    tops, rules, bots = [], [], []
    for pair in ind.pairs:
        t, r, b = format_pair_fraction(pair)
        tops.append(t)
        rules.append(r)
        bots.append(b)

    sep = " ; "
    top_line = sep.join(tops)
    rule_line = sep.join(rules)
    bot_line = sep.join(bots)

    sex_mark = "♀" if ind.sex == "female" else "♂"
    # Annotate ratio on the middle (rule) line
    annot = f"   ({ind.probability}) {sex_mark}"
    rule_line = rule_line + annot

    return f"{top_line}\n{rule_line}\n{bot_line}"


def format_cross_result(
    result: CrossResult,
    *,
    show_fraction: bool = True,
    show_compact: bool = True,
    show_phenotype: bool = True,
) -> str:
    """Full text report for GUI / export."""
    lines: list[str] = []
    lines.append("Crossing Scheme")
    lines.append("=" * 48)
    lines.append(f"♀  {result.female_raw}")
    lines.append(f"♂  {result.male_raw}")
    lines.append(f"RF (adjacent) = {result.rf}")
    lines.append("-" * 48)
    lines.append("F1 (viable)")
    lines.append("")

    if not result.offspring:
        lines.append("No viable offspring (all classes lethal or empty).")
    else:
        for i, ind in enumerate(result.offspring, 1):
            lines.append(f"[{i}]")
            if show_fraction:
                lines.append(format_individual_fraction(ind))
            if show_compact:
                lines.append(f"    compact: {format_individual_compact(ind)}")
            if show_phenotype and ind.phenotype_hints:
                hints = ", ".join(ind.phenotype_hints)
                lines.append(f"    phenotype hints (dominant): {hints}")
            lines.append("")

    lines.append("-" * 48)
    lines.append(f"Viable classes: {result.viable_count}")
    # Sum check
    total = sum((ind.probability for ind in result.offspring), Fraction(0))
    lines.append(f"Probability sum (viable): {total}")
    if result.lethal_discarded_probability > 0:
        lines.append(
            f"Lethal discarded (pre-normalization share): "
            f"{result.lethal_discarded_probability}"
            f" of {result.raw_total_before_lethal}"
        )
        lines.append("Note: homozygous lethal balancer/marker classes were excluded.")

    return "\n".join(lines)


def format_parents_line(female: Genotype, male: Genotype) -> str:
    return (
        f"♀ {format_genotype_compact(female)}  ×  "
        f"♂ {format_genotype_compact(male)}"
    )
