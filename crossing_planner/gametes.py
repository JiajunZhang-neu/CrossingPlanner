"""Gamete generation with simple adjacent-locus recombination."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from typing import Iterable

from .genotype import ChromosomePair, Genotype, Haplotype, SexChromRole
from .phenotype import chromosome_suppressed, translocation_indices


@dataclass(frozen=True)
class Gamete:
    """Full haploid contribution: one haplotype per chromosome index."""

    chromosomes: tuple[Haplotype, ...]
    probability: Fraction
    roles: tuple[SexChromRole, ...]

    def key(self) -> tuple[tuple[str, ...], ...]:
        return tuple(h.loci for h in self.chromosomes)


def _single_crossover_products(
    hap_a: tuple[str, ...], hap_b: tuple[str, ...], break_after: int
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Reciprocal products after a single crossover between locus break_after and break_after+1."""
    n = len(hap_a)
    assert len(hap_b) == n
    assert 0 <= break_after < n - 1
    i = break_after + 1
    return hap_a[:i] + hap_b[i:], hap_b[:i] + hap_a[i:]


def _pad_haplotypes(a: Haplotype, b: Haplotype) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """
    Align haplotypes to the same locus count for recombination.

    If lengths differ, pad the shorter with '+' at the end so RF can still apply
    when both have multiple loci; if either is Y, no recombination.
    """
    la, lb = list(a.loci), list(b.loci)
    if a.is_y or b.is_y:
        return tuple(la), tuple(lb)
    n = max(len(la), len(lb), 1)
    while len(la) < n:
        la.append("+")
    while len(lb) < n:
        lb.append("+")
    return tuple(la), tuple(lb)


def haplotypes_from_pair(
    pair: ChromosomePair, rf: Fraction, sex: str
) -> list[tuple[Haplotype, Fraction]]:
    """
    Possible haplotypes transmitted from one chromosome pair, with probabilities.

    rf: recombination fraction between adjacent loci (0 = complete linkage).
    Simple model: at most one crossover, equally likely in each of (n-1) intervals
    when a crossover occurs; probability of crossover = rf (capped use as overall
    recombination rate between ends approximated by adjacent RF for each interval).

    For RF interpretation: each adjacent interval recombines independently with
    probability `rf`, but we use the simple single-exchange approximation:
    - With probability (1 - rf_eff): parental haplotypes 1:1
    - With probability rf_eff: one of the (n-1) intervals chosen uniformly,
      reciprocal recombinants each get half of that share

    For a single interval, rf_eff = rf.
    For multiple intervals, rf_eff = min(Fraction(1, 2), rf * (n-1)) is NOT used;
    instead each interval has recombination probability rf, and we only model
    exactly 0 or 1 crossover total:
      P(no CO) = (1-rf)^(n-1) approximated as (1 - rf) for simplicity when n=2,
      and for n>2: P(CO in interval i only) proportional — plan says
      "单次交换近似".

    Implementation (simple):
      If n <= 1 or rf == 0 or Y involved:
        parental only, each 1/2 (or 1 for XX homozygous same, still 1/2+1/2 same)
      Else:
        P(parental type set) = 1 - rf
          each parental hap gets (1-rf)/2
        P(recombinant) = rf
          for each of (n-1) breakpoints, probability rf/(n-1)
          each reciprocal product gets half of that = rf/(2*(n-1))
    """
    mat, pat = pair.homologs()

    # Male sex chromosome: X or Y each 1/2, no recombination in male Drosophila (classic)
    if pair.role == SexChromRole.XY and sex == "male":
        results: list[tuple[Haplotype, Fraction]] = []
        for hap in (mat, pat):
            results.append((hap.normalized(), Fraction(1, 2)))
        return _merge_hap_probs(results)

    # Female XX or autosomes: recombination allowed
    if mat.is_y or pat.is_y:
        # Shouldn't happen for female; treat as parental
        return _merge_hap_probs(
            [(mat.normalized(), Fraction(1, 2)), (pat.normalized(), Fraction(1, 2))]
        )

    if rf <= 0:
        return _merge_hap_probs(
            [
                (mat.normalized(), Fraction(1, 2)),
                (pat.normalized(), Fraction(1, 2)),
            ]
        )

    a, b = _pad_haplotypes(mat, pat)
    n = len(a)

    if n <= 1:
        return _merge_hap_probs(
            [
                (mat.normalized(), Fraction(1, 2)),
                (pat.normalized(), Fraction(1, 2)),
            ]
        )

    # Cap rf at 1/2
    if rf > Fraction(1, 2):
        rf = Fraction(1, 2)

    intervals = n - 1
    parental_share = Fraction(1) - rf
    out: list[tuple[Haplotype, Fraction]] = [
        (Haplotype(a).normalized(), parental_share / 2),
        (Haplotype(b).normalized(), parental_share / 2),
    ]

    per_break = rf / intervals
    for br in range(intervals):
        r1, r2 = _single_crossover_products(a, b, br)
        out.append((Haplotype(r1).normalized(), per_break / 2))
        out.append((Haplotype(r2).normalized(), per_break / 2))

    return _merge_hap_probs(out)


def _merge_hap_probs(
    items: Iterable[tuple[Haplotype, Fraction]],
) -> list[tuple[Haplotype, Fraction]]:
    acc: dict[tuple[str, ...], Fraction] = {}
    for hap, p in items:
        key = hap.normalized().loci
        acc[key] = acc.get(key, Fraction(0)) + p
    return [(Haplotype(k), v) for k, v in acc.items()]


def generate_gametes(genotype: Genotype, rf: Fraction) -> list[Gamete]:
    """Cartesian product of per-chromosome haplotype contributions."""
    per_chr: list[list[tuple[Haplotype, Fraction]]] = []
    roles: list[SexChromRole] = []
    suppressed = translocation_indices(genotype.pairs)

    for index, pair in enumerate(genotype.pairs):
        roles.append(pair.role)
        rf_here = rf
        if genotype.sex == "male" or chromosome_suppressed(pair) or index in suppressed:
            rf_here = Fraction(0)
        per_chr.append(haplotypes_from_pair(pair, rf_here, genotype.sex))

    gametes: list[Gamete] = []
    for combo in product(*per_chr):
        haps = tuple(h for h, _ in combo)
        prob = Fraction(1)
        for _, p in combo:
            prob *= p
        gametes.append(
            Gamete(chromosomes=haps, probability=prob, roles=tuple(roles))
        )

    return _merge_gametes(gametes)


def _merge_gametes(gametes: list[Gamete]) -> list[Gamete]:
    acc: dict[tuple[tuple[str, ...], ...], tuple[Fraction, tuple[SexChromRole, ...]]] = {}
    for g in gametes:
        key = g.key()
        if key in acc:
            old_p, roles = acc[key]
            acc[key] = (old_p + g.probability, roles)
        else:
            acc[key] = (g.probability, g.roles)
    return [
        Gamete(
            chromosomes=tuple(Haplotype(loci) for loci in key),
            probability=p,
            roles=roles,
        )
        for key, (p, roles) in acc.items()
    ]
