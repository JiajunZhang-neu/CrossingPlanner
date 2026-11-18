"""Allele-token helpers and compatibility views of the marker catalog."""

from __future__ import annotations

from .catalog import BALANCERS, MARKERS, split_tokens


def parse_extra_lethals(text: str) -> frozenset[str]:
    """Parse comma/semicolon-separated extra lethal allele names."""
    if not text or not text.strip():
        return frozenset()
    parts = []
    for chunk in text.replace(";", ",").split(","):
        name = chunk.strip()
        if name:
            parts.append(name)
    return frozenset(parts)


def allele_tokens(allele: str) -> list[str]:
    """
    Split a haplotype locus token into comparable allele pieces.

    Examples: 'TM3, Sb' -> ['TM3', 'Sb']; 'CyO' -> ['CyO']; 'w[*]' -> ['w[*]']
    """
    return split_tokens(allele)


DEFAULT_LETHAL_ALLELES: frozenset[str] = frozenset(
    balancer.symbol
    for balancer in BALANCERS
    if balancer.homozygous == "lethal"
) | frozenset(marker.symbol for marker in MARKERS if marker.homozygous == "lethal")

DEFAULT_DOMINANT_MARKERS: frozenset[str] = frozenset(
    marker.symbol for marker in MARKERS if marker.dominance == "dominant"
) | frozenset(balancer.symbol for balancer in BALANCERS)


def haplotype_contains_lethal(haplotype_loci: tuple[str, ...], lethals: frozenset[str]) -> set[str]:
    """Return lethal allele names found anywhere in the haplotype loci."""
    found: set[str] = set()
    for locus in haplotype_loci:
        for tok in allele_tokens(locus):
            if tok in lethals:
                found.add(tok)
            for lethal in lethals:
                if tok == lethal or tok.startswith(lethal + ",") or tok.startswith(lethal + " "):
                    found.add(lethal)
    return found


def dominant_visible_in_pair(
    maternal: tuple[str, ...],
    paternal: tuple[str, ...],
    dominants: frozenset[str] = DEFAULT_DOMINANT_MARKERS,
) -> list[str]:
    """List dominant markers present on either homolog."""
    visible: list[str] = []
    seen: set[str] = set()
    for locus in (*maternal, *paternal):
        for tok in allele_tokens(locus):
            if tok in ("+", "Y"):
                continue
            if tok in dominants and tok not in seen:
                seen.add(tok)
                visible.append(tok)
    return visible
