"""Built-in lethal and dominant marker tables for Drosophila genetics."""

from __future__ import annotations

# Balancers / rearrangements typically lethal when homozygous.
DEFAULT_LETHAL_ALLELES: frozenset[str] = frozenset(
    {
        "CyO",
        "SM5",
        "SM6",
        "SM6a",
        "SM6b",
        "TM3",
        "TM6",
        "TM6B",
        "TM6C",
        "FM7",
        "FM7a",
        "FM7c",
        "FM6",
        "Binsn",
        "Ser",
    }
)

# Dominant visible markers (name used for phenotype hints).
# Keys are allele tokens as they commonly appear in genotypes.
DEFAULT_DOMINANT_MARKERS: frozenset[str] = frozenset(
    {
        "Cy",
        "CyO",
        "Sp",
        "Sco",
        "Pm",
        "Bla",
        "L",
        "Gla",
        "SM5",
        "SM6",
        "SM6a",
        "SM6b",
        "Sb",
        "Ser",
        "Hu",
        "D",
        "Dr",
        "Ly",
        "Ubx",
        "TM3",
        "TM6",
        "TM6B",
        "TM6C",
        "Tb",
        "MKRS",
        "Bar",
        "B",
        "Binsn",
        "FM7",
        "FM7a",
        "FM7c",
        "FM6",
        "w[+]",
        "If",
        "KrIf-1",
        "Gl",
    }
)


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
    allele = allele.strip()
    if not allele:
        return []
    # Comma-separated compound balancer labels
    if "," in allele:
        return [p.strip() for p in allele.split(",") if p.strip()]
    return [allele]


def haplotype_contains_lethal(haplotype_loci: tuple[str, ...], lethals: frozenset[str]) -> set[str]:
    """Return lethal allele names found anywhere in the haplotype loci."""
    found: set[str] = set()
    for locus in haplotype_loci:
        for tok in allele_tokens(locus):
            if tok in lethals:
                found.add(tok)
            # Also match bare balancer if token starts with it (e.g. TM6B, Tb[*])
            for lethal in lethals:
                if tok == lethal or tok.startswith(lethal + ",") or tok.startswith(lethal + " "):
                    found.add(lethal)
    return found


def dominant_visible_in_pair(
    maternal: tuple[str, ...],
    paternal: tuple[str, ...],
    dominants: frozenset[str] = DEFAULT_DOMINANT_MARKERS,
) -> list[str]:
    """List dominant markers present on either homolog (phenotype hint)."""
    visible: list[str] = []
    seen: set[str] = set()
    for locus in (*maternal, *paternal):
        for tok in allele_tokens(locus):
            if tok in ("+", "Y"):
                continue
            if tok in dominants and tok not in seen:
                seen.add(tok)
                visible.append(tok)
            else:
                # Partial match for compound labels like 'TM6B, Tb'
                for d in dominants:
                    if (tok == d or d in allele_tokens(tok)) and d not in seen:
                        seen.add(d)
                        visible.append(d)
    return visible
