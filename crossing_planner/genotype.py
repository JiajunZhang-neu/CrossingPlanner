"""Parse and represent Drosophila genotypes from plain text."""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum, auto


class SexChromRole(Enum):
    AUTOSOME = auto()
    XX = auto()  # female X pair (or first chromosome pair treated as sex)
    XY = auto()  # male X/Y


class ParseError(ValueError):
    """Raised when genotype text cannot be parsed."""


@dataclass(frozen=True)
class Haplotype:
    """Ordered loci on one chromatid / chromosome copy."""

    loci: tuple[str, ...]

    @property
    def is_y(self) -> bool:
        return len(self.loci) == 1 and self.loci[0].upper() == "Y"

    def display(self) -> str:
        if not self.loci:
            return "+"
        return " ".join(self.loci)

    def normalized(self) -> Haplotype:
        """Canonical empty -> '+'."""
        if not self.loci:
            return Haplotype(("+",))
        return self


@dataclass(frozen=True)
class ChromosomePair:
    """A homologous pair (or X/Y)."""

    maternal: Haplotype  # first haplotype written (left of /)
    paternal: Haplotype  # second haplotype (right of /); for hemizygous may be synthetic
    role: SexChromRole = SexChromRole.AUTOSOME
    homozygous_shorthand: bool = False  # True if input was a single haplotype (no /)

    def homologs(self) -> tuple[Haplotype, Haplotype]:
        return (self.maternal.normalized(), self.paternal.normalized())


@dataclass(frozen=True)
class Genotype:
    """Full genotype: ordered chromosome pairs."""

    pairs: tuple[ChromosomePair, ...]
    sex: str  # 'female' or 'male'
    raw: str = ""

    def chromosome_count(self) -> int:
        return len(self.pairs)


def _split_chromosomes(text: str) -> list[str]:
    parts = [p.strip() for p in text.split(";")]
    return [p for p in parts if p]


def _parse_haplotype(text: str) -> Haplotype:
    """
    Parse one chromatid.

    Spaces separate loci. Comma-separated labels stay one locus:
    'TM3, Sb' -> ('TM3, Sb',); 'a b' -> ('a', 'b').
    """
    text = text.strip()
    if not text:
        return Haplotype(("+",))

    # "TM3, Sb" / "TM6B, Tb[*]" remain a single locus: collapse comma-space, split, restore
    compact = re.sub(r",\s+", ",", text)
    raw_tokens = [t for t in compact.split() if t]
    tokens = [t.replace(",", ", ") for t in raw_tokens]
    if not tokens:
        return Haplotype(("+",))
    return Haplotype(tuple(tokens))


def _parse_pair(segment: str, is_first: bool, sex: str) -> ChromosomePair:
    segment = segment.strip()
    if not segment:
        raise ParseError("Empty chromosome segment")

    if "/" in segment:
        left, right = segment.split("/", 1)
        mat = _parse_haplotype(left)
        pat = _parse_haplotype(right)
        shorthand = False
    else:
        hap = _parse_haplotype(segment)
        mat = hap
        pat = hap
        shorthand = True

    role = SexChromRole.AUTOSOME
    if is_first:
        if sex == "male":
            if pat.is_y or mat.is_y:
                role = SexChromRole.XY
                if mat.is_y and not pat.is_y:
                    mat, pat = pat, mat
            else:
                role = SexChromRole.XY
                if shorthand:
                    pat = Haplotype(("Y",))
                    shorthand = False
        else:
            role = SexChromRole.XX

    return ChromosomePair(
        maternal=mat,
        paternal=pat,
        role=role,
        homozygous_shorthand=shorthand,
    )


def parse_genotype(text: str, sex: str) -> Genotype:
    """
    Parse plain-text Drosophila genotype.

    sex: 'female' or 'male'
    """
    sex = sex.lower().strip()
    if sex not in ("female", "male"):
        raise ParseError("sex must be 'female' or 'male'")

    cleaned = " ".join(text.replace("\n", " ").split())
    if not cleaned:
        raise ParseError("Genotype text is empty")

    segments = _split_chromosomes(cleaned)
    if not segments:
        raise ParseError("No chromosome segments found")

    pairs: list[ChromosomePair] = []
    for i, seg in enumerate(segments):
        pairs.append(_parse_pair(seg, is_first=(i == 0), sex=sex))

    if sex == "male" and pairs:
        first = pairs[0]
        if first.role == SexChromRole.XY:
            mat, pat = first.homologs()
            if not (mat.is_y or pat.is_y):
                raise ParseError(
                    "Male sex chromosome should include Y (e.g. 'w / Y'). "
                    "Hemizygous 'w' is also accepted and expanded to 'w / Y'."
                )

    return Genotype(pairs=tuple(pairs), sex=sex, raw=cleaned)


def genotypes_compatible(female: Genotype, male: Genotype) -> None:
    """Raise ParseError if chromosome counts differ."""
    if female.sex != "female":
        raise ParseError("First parent must be female")
    if male.sex != "male":
        raise ParseError("Second parent must be male")
    if len(female.pairs) != len(male.pairs):
        raise ParseError(
            f"Chromosome count mismatch: female has {len(female.pairs)}, "
            f"male has {len(male.pairs)}"
        )


def format_genotype_compact(gt: Genotype) -> str:
    """Single-line compact genotype string."""
    parts: list[str] = []
    for pair in gt.pairs:
        mat, pat = pair.homologs()
        if pair.homozygous_shorthand and mat == pat and pair.role != SexChromRole.XY:
            parts.append(mat.display())
        else:
            parts.append(f"{mat.display()} / {pat.display()}")
    return " ; ".join(parts)


def sort_homologs_for_display(
    a: Haplotype, b: Haplotype
) -> tuple[Haplotype, Haplotype]:
    """Stable display order: non-Y first; then lexicographic."""
    if a.is_y and not b.is_y:
        return b, a
    if b.is_y and not a.is_y:
        return a, b
    if a.display() <= b.display():
        return a, b
    return b, a
