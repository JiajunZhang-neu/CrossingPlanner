"""Cross female x male gametes, filter lethals, normalize F1 ratios."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .gametes import Gamete, generate_gametes
from .genotype import (
    ChromosomePair,
    Genotype,
    Haplotype,
    SexChromRole,
    genotypes_compatible,
    sort_homologs_for_display,
)
from .markers import (
    DEFAULT_DOMINANT_MARKERS,
    DEFAULT_LETHAL_ALLELES,
    allele_tokens,
    dominant_visible_in_pair,
    haplotype_contains_lethal,
)


@dataclass(frozen=True)
class F1Individual:
    """One distinct viable F1 genotype class."""

    pairs: tuple[ChromosomePair, ...]
    probability: Fraction  # normalized among viable
    phenotype_hints: tuple[str, ...]
    sex: str  # 'female' or 'male'

    def as_genotype(self) -> Genotype:
        return Genotype(pairs=self.pairs, sex=self.sex, raw="")


@dataclass(frozen=True)
class CrossResult:
    offspring: tuple[F1Individual, ...]
    raw_total_before_lethal: Fraction
    lethal_discarded_probability: Fraction
    female_raw: str
    male_raw: str
    rf: Fraction
    lethals_used: frozenset[str]

    @property
    def viable_count(self) -> int:
        return len(self.offspring)


def _sex_from_sex_chrom(maternal_x: Haplotype, paternal_sex: Haplotype) -> str:
    """F1 sex from sex-chromosome contribution: Y from father -> male."""
    if paternal_sex.is_y:
        return "male"
    return "female"


def _pair_has_lethal_homozygote(
    mat: Haplotype, pat: Haplotype, lethals: frozenset[str]
) -> bool:
    """True if both homologs carry the same lethal allele token."""
    if mat.is_y or pat.is_y:
        # Hemizygous X/Y: balancer on X alone is viable (FM7/Y males exist)
        return False

    mat_lethals = haplotype_contains_lethal(mat.loci, lethals)
    pat_lethals = haplotype_contains_lethal(pat.loci, lethals)
    shared = mat_lethals & pat_lethals
    if not shared:
        # Also check raw locus string equality for lethal tokens
        mat_toks = set()
        pat_toks = set()
        for locus in mat.loci:
            mat_toks.update(allele_tokens(locus))
        for locus in pat.loci:
            pat_toks.update(allele_tokens(locus))
        shared = (mat_toks & pat_toks) & set(lethals)
    return bool(shared)


def _zygote_is_lethal(
    fem_gamete: Gamete, male_gamete: Gamete, lethals: frozenset[str]
) -> bool:
    for fh, mh in zip(fem_gamete.chromosomes, male_gamete.chromosomes):
        if _pair_has_lethal_homozygote(fh, mh, lethals):
            return True
    return False


def _build_pairs(
    fem_gamete: Gamete, male_gamete: Gamete, offspring_sex: str
) -> tuple[ChromosomePair, ...]:
    pairs: list[ChromosomePair] = []
    for i, (fh, mh) in enumerate(
        zip(fem_gamete.chromosomes, male_gamete.chromosomes)
    ):
        role = fem_gamete.roles[i] if i < len(fem_gamete.roles) else SexChromRole.AUTOSOME
        if i == 0:
            role = SexChromRole.XY if offspring_sex == "male" else SexChromRole.XX
            top, bot = sort_homologs_for_display(fh, mh)
            # For males, keep X / Y order
            if offspring_sex == "male":
                if mh.is_y:
                    top, bot = fh.normalized(), mh.normalized()
                elif fh.is_y:
                    top, bot = mh.normalized(), fh.normalized()
            pairs.append(
                ChromosomePair(
                    maternal=top.normalized(),
                    paternal=bot.normalized(),
                    role=role,
                )
            )
        else:
            top, bot = sort_homologs_for_display(fh, mh)
            pairs.append(
                ChromosomePair(
                    maternal=top.normalized(),
                    paternal=bot.normalized(),
                    role=SexChromRole.AUTOSOME,
                )
            )
    return tuple(pairs)


def _phenotype_hints(pairs: tuple[ChromosomePair, ...]) -> tuple[str, ...]:
    hints: list[str] = []
    seen: set[str] = set()
    for pair in pairs:
        mat, pat = pair.homologs()
        for h in dominant_visible_in_pair(mat.loci, pat.loci, DEFAULT_DOMINANT_MARKERS):
            if h not in seen:
                seen.add(h)
                hints.append(h)
    return tuple(hints)


def _genotype_key(pairs: tuple[ChromosomePair, ...], sex: str) -> tuple:
    key_pairs = []
    for idx, pair in enumerate(pairs):
        mat, pat = pair.homologs()
        # ordered for sex chromosome; unordered for autosomes
        if idx == 0 and pair.role in (SexChromRole.XY, SexChromRole.XX):
            key_pairs.append((mat.loci, pat.loci))
        else:
            a, b = sorted((mat.loci, pat.loci))
            key_pairs.append((a, b))
    return (sex, tuple(key_pairs))


def cross(
    female: Genotype,
    male: Genotype,
    rf: Fraction = Fraction(0),
    extra_lethals: frozenset[str] | None = None,
) -> CrossResult:
    """
    Compute viable F1 genotype classes with normalized Fraction probabilities.
    """
    genotypes_compatible(female, male)
    lethals = DEFAULT_LETHAL_ALLELES | (extra_lethals or frozenset())

    fem_gametes = generate_gametes(female, rf)
    male_gametes = generate_gametes(male, rf)

    # Accumulate raw probabilities before normalization
    raw: dict[tuple, tuple[Fraction, tuple[ChromosomePair, ...], str, tuple[str, ...]]] = {}
    total_raw = Fraction(0)
    lethal_prob = Fraction(0)

    for fg in fem_gametes:
        for mg in male_gametes:
            p = fg.probability * mg.probability
            total_raw += p
            if _zygote_is_lethal(fg, mg, lethals):
                lethal_prob += p
                continue

            # Sex from father's sex-chromosome contribution (index 0)
            paternal_sex_chr = mg.chromosomes[0]
            maternal_sex_chr = fg.chromosomes[0]
            offspring_sex = _sex_from_sex_chrom(maternal_sex_chr, paternal_sex_chr)
            pairs = _build_pairs(fg, mg, offspring_sex)
            hints = _phenotype_hints(pairs)
            key = _genotype_key(pairs, offspring_sex)
            if key in raw:
                old_p, old_pairs, old_sex, old_hints = raw[key]
                raw[key] = (old_p + p, old_pairs, old_sex, old_hints)
            else:
                raw[key] = (p, pairs, offspring_sex, hints)

    viable_raw = total_raw - lethal_prob
    offspring_list: list[F1Individual] = []

    if viable_raw > 0:
        for _key, (p, pairs, sex, hints) in raw.items():
            offspring_list.append(
                F1Individual(
                    pairs=pairs,
                    probability=p / viable_raw,
                    phenotype_hints=hints,
                    sex=sex,
                )
            )

    offspring_list.sort(key=lambda x: (-x.probability, x.sex, str(x.pairs)))

    return CrossResult(
        offspring=tuple(offspring_list),
        raw_total_before_lethal=total_raw,
        lethal_discarded_probability=lethal_prob,
        female_raw=female.raw,
        male_raw=male.raw,
        rf=rf,
        lethals_used=lethals,
    )
