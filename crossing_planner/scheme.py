"""Parallel and multi-step crossing schemes.

A step's parents are stocks, or a chosen offspring class from an earlier step.
Steps that do not depend on each other share a wave and can be set up together.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .cross import CrossResult, cross
from .genotype import Genotype, format_genotype_compact, parse_genotype
from .i18n import tr


@dataclass(frozen=True)
class ParentRef:
    """Where one parent of a step comes from."""

    kind: str  # "stock" or "offspring"
    stock_id: str = ""
    step_id: str = ""
    offspring_index: int = 0

    @staticmethod
    def stock(stock_id: str) -> ParentRef:
        return ParentRef("stock", stock_id=stock_id)

    @staticmethod
    def offspring(step_id: str, index: int) -> ParentRef:
        return ParentRef("offspring", step_id=step_id, offspring_index=index)


@dataclass(frozen=True)
class Stock:
    id: str
    sex: str
    genotype: str
    label: str = ""

    def display_name(self) -> str:
        return self.label or self.id


@dataclass(frozen=True)
class CrossStep:
    id: str
    female: ParentRef
    male: ParentRef
    rf: Fraction = Fraction(0)
    label: str = ""

    def display_name(self) -> str:
        return self.label or self.id


@dataclass(frozen=True)
class StepResult:
    step: CrossStep
    result: CrossResult
    wave: int
    female_text: str
    male_text: str
    female_genotype: Genotype
    male_genotype: Genotype


@dataclass(frozen=True)
class SchemeResult:
    steps: tuple[StepResult, ...]
    waves: tuple[tuple[str, ...], ...]


class SchemeError(ValueError):
    """The scheme cannot be run."""


def wave_ranks(steps: list[CrossStep]) -> dict[str, int]:
    """Topological wave of each step. Independent steps share wave 0."""
    by_id = {step.id: step for step in steps}
    if len(by_id) != len(steps):
        raise SchemeError(tr("err_dup_step"))

    ranks: dict[str, int] = {}

    def rank_of(step_id: str, stack: set[str]) -> int:
        if step_id in ranks:
            return ranks[step_id]
        if step_id not in by_id:
            raise SchemeError(tr("err_missing_step", step_id=step_id))
        if step_id in stack:
            raise SchemeError(tr("err_cycle"))
        stack.add(step_id)
        step = by_id[step_id]
        rank = 0
        for ref in (step.female, step.male):
            if ref.kind == "offspring":
                rank = max(rank, rank_of(ref.step_id, stack) + 1)
        stack.remove(step_id)
        ranks[step_id] = rank
        return rank

    for step in steps:
        rank_of(step.id, set())
    return ranks


def run_scheme(
    stocks: list[Stock],
    steps: list[CrossStep],
    extra_lethals: frozenset[str] | None = None,
) -> SchemeResult:
    """Run every step in wave order. Offspring refs must point at an earlier step."""
    stock_by_id = {stock.id: stock for stock in stocks}
    if len(stock_by_id) != len(stocks):
        raise SchemeError(tr("err_dup_stock"))
    ranks = wave_ranks(steps)
    order = sorted(steps, key=lambda step: (ranks[step.id], step.id))

    done: dict[str, StepResult] = {}
    for step in order:
        female_gt, female_text = _resolve(step.female, "female", stock_by_id, done)
        male_gt, male_text = _resolve(step.male, "male", stock_by_id, done)
        result = cross(female_gt, male_gt, rf=step.rf, extra_lethals=extra_lethals)
        # Keep the resolved text when the genotype raw string is empty (selected F1).
        if not result.female_raw:
            result = _with_parent_text(result, female_text, result.male_raw or male_text)
        if not result.male_raw:
            result = _with_parent_text(result, result.female_raw or female_text, male_text)
        done[step.id] = StepResult(
            step=step,
            result=result,
            wave=ranks[step.id],
            female_text=female_text,
            male_text=male_text,
            female_genotype=female_gt,
            male_genotype=male_gt,
        )

    wave_ids: dict[int, list[str]] = {}
    for step in order:
        wave_ids.setdefault(ranks[step.id], []).append(step.id)
    waves = tuple(tuple(wave_ids[index]) for index in sorted(wave_ids))
    return SchemeResult(steps=tuple(done[step.id] for step in order), waves=waves)


def relabel_scheme(scheme: SchemeResult, steps: list[CrossStep]) -> SchemeResult:
    """Replace stored step labels with the names currently on the plan."""
    by_id = {step.id: step for step in steps}
    updated: list[StepResult] = []
    for item in scheme.steps:
        live = by_id.get(item.step.id, item.step)
        updated.append(
            StepResult(
                step=live,
                result=item.result,
                wave=item.wave,
                female_text=item.female_text,
                male_text=item.male_text,
                female_genotype=item.female_genotype,
                male_genotype=item.male_genotype,
            )
        )
    return SchemeResult(steps=tuple(updated), waves=scheme.waves)


def _resolve(
    ref: ParentRef,
    expected_sex: str,
    stocks: dict[str, Stock],
    done: dict[str, StepResult],
) -> tuple[Genotype, str]:
    if ref.kind == "stock":
        stock = stocks.get(ref.stock_id)
        if stock is None:
            raise SchemeError(tr("err_missing_stock", stock_id=ref.stock_id))
        if stock.sex != expected_sex:
            role = tr("female_role") if expected_sex == "female" else tr("male_role")
            raise SchemeError(tr("err_stock_sex", name=stock.display_name(), role=role))
        genotype = parse_genotype(stock.genotype, stock.sex)
        return genotype, stock.genotype

    if ref.kind != "offspring":
        raise SchemeError(tr("err_bad_kind", kind=ref.kind))

    previous = done.get(ref.step_id)
    if previous is None:
        raise SchemeError(tr("err_step_no_result", step_id=ref.step_id))
    offspring = previous.result.offspring
    if ref.offspring_index < 0 or ref.offspring_index >= len(offspring):
        raise SchemeError(tr("err_no_class", step_id=ref.step_id, n=ref.offspring_index + 1))
    chosen = offspring[ref.offspring_index]
    if chosen.sex != expected_sex:
        role = tr("female_role") if expected_sex == "female" else tr("male_role")
        raise SchemeError(tr("err_class_sex", step_id=ref.step_id, n=ref.offspring_index + 1, role=role))
    # Parent lines already add one sex mark, so the genotype text stays unmarked.
    return chosen.as_genotype(), format_genotype_compact(chosen.as_genotype())


def _with_parent_text(result: CrossResult, female_raw: str, male_raw: str) -> CrossResult:
    return CrossResult(
        offspring=result.offspring,
        raw_total_before_lethal=result.raw_total_before_lethal,
        lethal_discarded_probability=result.lethal_discarded_probability,
        female_raw=female_raw,
        male_raw=male_raw,
        rf=result.rf,
        lethals_used=result.lethals_used,
        warnings=result.warnings,
    )
