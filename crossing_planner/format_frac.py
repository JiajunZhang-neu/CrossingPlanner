"""Fraction-format rendering of F1 genotypes."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .cross import CrossResult, F1Individual
from .genotype import ChromosomePair, Genotype, format_genotype_compact
from .i18n import localize, tr, translate_note
from .phenotype import Phenotype, assign_colors, phenotype_of, read_haplotype


@dataclass(frozen=True)
class TextSpan:
    """One run inside a report line. Used for a checkbox or plus that can highlight on its own."""

    text: str
    color: str | None = None
    tag: str = ""


@dataclass(frozen=True)
class ReportLine:
    text: str
    color: str | None = None  # RRGGBB, without '#'
    tag: str = ""
    spans: tuple[TextSpan, ...] = ()


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


_CHROM_NAMES = {0: "X", 1: "2", 2: "3", 3: "4"}


@dataclass(frozen=True)
class ChromOption:
    """One viable state of a single chromosome, summed across the other chromosomes."""

    pair: ChromosomePair
    probability: Fraction
    phenotype_text: str
    phenotype_key: tuple[tuple[str, str], ...]
    color_name: str
    color_hex: str
    sex: str
    sterile: bool
    cargo: tuple[str, ...]


def chromosome_options(result: CrossResult) -> list[list[ChromOption]]:
    """Marginal viable outcomes per chromosome. Probabilities on one chromosome sum to 1."""
    if not result.offspring:
        return []
    count = len(result.offspring[0].pairs)
    grouped: list[dict[tuple, list]] = [dict() for _ in range(count)]
    for ind in result.offspring:
        for index, pair in enumerate(ind.pairs):
            mat, pat = pair.homologs()
            sex = ind.sex if index == 0 else ""
            key = (sex, mat.loci, pat.loci)
            bucket = grouped[index].setdefault(key, [Fraction(0), pair, sex])
            bucket[0] += ind.probability

    views: list[list[ChromOption]] = []
    for index, buckets in enumerate(grouped):
        pending: list[tuple[Fraction, ChromosomePair, str, Phenotype]] = []
        for probability, pair, sex in buckets.values():
            pheno_sex = sex or "female"
            phenotype = phenotype_of((pair.homologs(),), pheno_sex)
            pending.append((probability, pair, sex, phenotype))
        colors = assign_colors([item[3].key for item in pending])
        options: list[ChromOption] = []
        for probability, pair, sex, phenotype in pending:
            color_name, color_hex = colors.get(phenotype.key, ("灰", "6B6B6B"))
            cargo: list[str] = []
            for haplotype in pair.homologs():
                for token in sorted(read_haplotype(haplotype).extra_tokens):
                    if token not in cargo:
                        cargo.append(token)
            options.append(
                ChromOption(
                    pair=pair,
                    probability=probability,
                    phenotype_text=phenotype.text(),
                    phenotype_key=phenotype.key,
                    color_name=color_name,
                    color_hex=color_hex,
                    sex=sex,
                    sterile=phenotype.sterile,
                    cargo=tuple(cargo),
                )
            )
        options.sort(key=lambda item: (-item.probability, item.sex, _pair_compact(item.pair)))
        views.append(options)
    return views


def match_offspring(
    result: CrossResult, options: list[list[ChromOption]], choices: list[int]
) -> F1Individual | None:
    """Find the full F1 class composed by one option from each chromosome."""
    if not options or len(choices) != len(options):
        return None
    if any(index < 0 or index >= len(column) for index, column in zip(choices, options)):
        return None
    sex = options[0][choices[0]].sex
    chosen = [column[index].pair for column, index in zip(options, choices)]
    for ind in result.offspring:
        if ind.sex != sex:
            continue
        if all(_same_pair(left, right) for left, right in zip(ind.pairs, chosen)):
            return ind
    return None


def _same_pair(left: ChromosomePair, right: ChromosomePair) -> bool:
    left_h = left.homologs()
    right_h = right.homologs()
    return left_h[0].loci == right_h[0].loci and left_h[1].loci == right_h[1].loci


def format_cross_lines(
    result: CrossResult,
    *,
    show_fraction: bool = True,
    show_compact: bool = True,
    show_phenotype: bool = True,
    title: str = "",
    scope: str = "",
    choices: list[int] | None = None,
    picks: list[str] | None = None,
    interactive: bool = False,
    allow_picks: bool = False,
) -> list[ReportLine]:
    """Per-chromosome F1 report. Clickable tags are set when interactive."""
    lines: list[ReportLine] = [
        ReportLine(title or tr("crossing_scheme")),
        ReportLine("=" * 48),
        ReportLine(f"♀  {result.female_raw}"),
        ReportLine(f"♂  {result.male_raw}"),
        ReportLine(tr("rf_value", rf=result.rf)),
        ReportLine("-" * 48),
        ReportLine(tr("f1_help_interactive") if interactive else tr("f1_help_static")),
        ReportLine(""),
    ]
    options = chromosome_options(result)
    if not options:
        lines.append(ReportLine(tr("no_viable")))
    else:
        selected = _clamp_choices(options, choices)
        for index, column in enumerate(options):
            name = _CHROM_NAMES.get(index, str(index + 1))
            lines.append(ReportLine(tr("chromosome", name=name)))
            for option_index, option in enumerate(column):
                lines.extend(
                    _option_lines(
                        option,
                        selected=option_index == selected[index],
                        show_fraction=show_fraction,
                        show_compact=show_compact,
                        show_phenotype=show_phenotype,
                        tag=f"opt:{scope}:{index}:{option_index}" if interactive else "",
                    )
                )
            lines.append(ReportLine(""))
        chosen = match_offspring(result, options, selected)
        lines.append(ReportLine(tr("current_combo")))
        if chosen is None:
            lines.append(ReportLine(tr("combo_missing_viable")))
        else:
            lines.extend(
                _selected_block(
                    chosen,
                    show_fraction=show_fraction,
                    show_compact=True,
                    show_phenotype=show_phenotype,
                )
            )
        if interactive and allow_picks and scope:
            lines.append(
                ReportLine(
                    "＋" + tr("add_combo"),
                    spans=(
                        TextSpan("＋", tag=f"plus:{scope}"),
                        TextSpan(tr("add_combo"), tag=f"add:{scope}"),
                    ),
                )
            )
        if picks:
            lines.append(ReportLine(tr("saved_picks")))
            for pick_index, label in enumerate(picks):
                remove = f"    {tr('remove')}" if interactive else ""
                tag = f"del:{scope}:{pick_index}" if interactive and scope else ""
                lines.append(ReportLine(f"  {pick_index + 1}. {label}{remove}", tag=tag))
            lines.append(ReportLine(""))

    lines.extend(_summary_lines(result))
    return lines


def _selected_block(
    ind: F1Individual,
    *,
    show_fraction: bool,
    show_compact: bool,
    show_phenotype: bool,
) -> list[ReportLine]:
    color = ind.color_hex or None if show_phenotype else None
    lines: list[ReportLine] = []
    if show_compact or not (show_fraction or show_phenotype):
        lines.append(ReportLine(format_individual_compact(ind), color=color))
    if show_fraction:
        for row in format_individual_fraction(ind).splitlines():
            lines.append(ReportLine("  " + row, color=color))
    if show_phenotype:
        sterile = tr("female_sterile_inline") if ind.sterile else ""
        lines.append(
            ReportLine(
                tr("phenotype_line", text=localize(ind.phenotype_text)) + sterile,
                color=color,
            )
        )
    return lines


def _summary_lines(result: CrossResult) -> list[ReportLine]:
    lines = [
        ReportLine("-" * 48),
        ReportLine(tr("viable_classes", count=result.viable_count)),
        ReportLine(
            tr(
                "prob_sum",
                total=sum((ind.probability for ind in result.offspring), Fraction(0)),
            )
        ),
    ]
    sterile_count = sum(1 for ind in result.offspring if ind.sterile)
    if sterile_count:
        lines.append(ReportLine(tr("sterile_kept", count=sterile_count)))
    if result.lethal_discarded_probability > 0:
        lines.append(
            ReportLine(
                tr(
                    "lethal_discarded",
                    discarded=result.lethal_discarded_probability,
                    raw=result.raw_total_before_lethal,
                )
            )
        )
        lines.append(ReportLine(tr("lethal_excluded")))
    for note in result.warnings:
        lines.append(ReportLine(tr("note_prefix", text=translate_note(note))))
    return lines


def _clamp_choices(options: list[list[ChromOption]], choices: list[int] | None) -> list[int]:
    selected = list(choices or [])
    if len(selected) < len(options):
        selected.extend([0] * (len(options) - len(selected)))
    clamped: list[int] = []
    for index, column in enumerate(options):
        choice = selected[index]
        if choice < 0 or choice >= len(column):
            choice = 0
        clamped.append(choice)
    return clamped


def _option_lines(
    option: ChromOption,
    *,
    selected: bool,
    show_fraction: bool,
    show_compact: bool,
    show_phenotype: bool,
    tag: str,
) -> list[ReportLine]:
    mark = "☑" if selected else "☐"
    sex_mark = "♀ " if option.sex == "female" else "♂ " if option.sex == "male" else ""
    color = option.color_hex or None if show_phenotype else None
    if show_compact:
        rest = f" {sex_mark}{_pair_compact(option.pair)}  ({option.probability})"
    else:
        rest = f" {sex_mark}{tr('option_word')}"
    if option.sterile:
        rest += tr("female_sterile_inline")
    box_tag = "box:" + tag[4:] if tag.startswith("opt:") else ""
    if box_tag:
        line = ReportLine(
            mark + rest,
            color=color,
            tag=tag,
            spans=(TextSpan(mark, tag=box_tag), TextSpan(rest, color=color, tag=tag)),
        )
    else:
        line = ReportLine(mark + rest, color=color, tag=tag)
    lines = [line]
    if show_fraction:
        top, rule, bot = format_pair_fraction(option.pair)
        sex = "♀" if option.sex == "female" else "♂" if option.sex == "male" else ""
        lines.append(ReportLine(f"    {top}", color=color, tag=tag))
        lines.append(ReportLine(f"    {rule}   ({option.probability}) {sex}", color=color, tag=tag))
        lines.append(ReportLine(f"    {bot}", color=color, tag=tag))
    if show_phenotype:
        cargo = tr("cargo", items=", ".join(option.cargo)) if option.cargo else ""
        lines.append(
            ReportLine(
                "    " + tr("phenotype_line", text=localize(option.phenotype_text)) + cargo,
                color=color,
                tag=tag,
            )
        )
    return lines


def format_cross_result(
    result: CrossResult,
    *,
    show_fraction: bool = True,
    show_compact: bool = True,
    show_phenotype: bool = True,
) -> str:
    """Full text report for GUI / export."""
    return "\n".join(
        line.text
        for line in format_cross_lines(
            result,
            show_fraction=show_fraction,
            show_compact=show_compact,
            show_phenotype=show_phenotype,
        )
    )


def format_parents_line(female: Genotype, male: Genotype) -> str:
    return (
        f"♀ {format_genotype_compact(female)}  ×  "
        f"♂ {format_genotype_compact(male)}"
    )
