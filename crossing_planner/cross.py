"""Cross female x male gametes, filter lethals, normalize F1 ratios."""

from __future__ import annotations

from dataclasses import dataclass, field
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
from .phenotype import (
    Phenotype,
    assign_colors,
    balancer_chromosome_warnings,
    chromosome_suppressed,
    confusable_warnings,
    has_proximal_balancer,
    phenotype_of,
    translocation_indices,
    zygote_status,
)


@dataclass(frozen=True)
class F1Individual:
    """One distinct viable F1 genotype class."""

    pairs: tuple[ChromosomePair, ...]
    probability: Fraction  # normalized among viable
    phenotype_hints: tuple[str, ...]
    sex: str  # 'female' or 'male'
    phenotype_text: str = ""
    phenotype_key: tuple[tuple[str, str], ...] = ()
    color_name: str = ""
    color_hex: str = ""
    sterile: bool = False

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
    warnings: tuple[str, ...] = field(default_factory=tuple)

    @property
    def viable_count(self) -> int:
        return len(self.offspring)


def _sex_from_sex_chrom(maternal_x: Haplotype, paternal_sex: Haplotype) -> str:
    """F1 sex from sex-chromosome contribution: Y from father -> male."""
    if paternal_sex.is_y:
        return "male"
    return "female"


def _build_pairs(
    fem_gamete: Gamete, male_gamete: Gamete, offspring_sex: str
) -> tuple[ChromosomePair, ...]:
    pairs: list[ChromosomePair] = []
    for i, (fh, mh) in enumerate(zip(fem_gamete.chromosomes, male_gamete.chromosomes)):
        role = fem_gamete.roles[i] if i < len(fem_gamete.roles) else SexChromRole.AUTOSOME
        if i == 0:
            role = SexChromRole.XY if offspring_sex == "male" else SexChromRole.XX
            top, bot = sort_homologs_for_display(fh, mh)
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


def _hint_symbols(phenotype: Phenotype) -> tuple[str, ...]:
    symbols: list[str] = []
    for trait in phenotype.traits:
        for symbol in sorted(trait.markers):
            if symbol not in symbols:
                symbols.append(symbol)
    return tuple(symbols)


def _genotype_key(pairs: tuple[ChromosomePair, ...], sex: str) -> tuple:
    key_pairs = []
    for idx, pair in enumerate(pairs):
        mat, pat = pair.homologs()
        if idx == 0 and pair.role in (SexChromRole.XY, SexChromRole.XX):
            key_pairs.append((mat.loci, pat.loci))
        else:
            a, b = sorted((mat.loci, pat.loci))
            key_pairs.append((a, b))
    return (sex, tuple(key_pairs))


def _cross_warnings(female: Genotype, male: Genotype, rf: Fraction) -> list[str]:
    notes: list[str] = []
    if rf > 0:
        notes.append("雄蝇各染色体不重组。")
        female_suppressed = any(chromosome_suppressed(pair) for pair in female.pairs)
        female_suppressed = female_suppressed or bool(translocation_indices(female.pairs))
        if female_suppressed:
            notes.append("雌性亲本里，带全臂平衡子的染色体对按完全连锁计算。")
        if has_proximal_balancer(female.pairs) or has_proximal_balancer(male.pairs):
            notes.append("MRS / MKRS 只抑制 3 号染色体近端交换，其余位点仍按设定的 RF 计算。")
    notes.extend(balancer_chromosome_warnings(female.pairs))
    notes.extend(balancer_chromosome_warnings(male.pairs))
    return notes


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
    lethals = extra_lethals or frozenset()

    fem_gametes = generate_gametes(female, rf)
    male_gametes = generate_gametes(male, rf)

    raw: dict[tuple, tuple[Fraction, tuple[ChromosomePair, ...], str, Phenotype]] = {}
    total_raw = Fraction(0)
    lethal_prob = Fraction(0)

    for fg in fem_gametes:
        for mg in male_gametes:
            p = fg.probability * mg.probability
            total_raw += p
            paternal_sex_chr = mg.chromosomes[0]
            maternal_sex_chr = fg.chromosomes[0]
            offspring_sex = _sex_from_sex_chrom(maternal_sex_chr, paternal_sex_chr)
            status = zygote_status(fg.chromosomes, mg.chromosomes, offspring_sex, lethals)
            if status == "lethal":
                lethal_prob += p
                continue

            pairs = _build_pairs(fg, mg, offspring_sex)
            homologs = tuple(pair.homologs() for pair in pairs)
            phenotype = phenotype_of(homologs, offspring_sex)
            if status == "sterile":
                phenotype = Phenotype(traits=phenotype.traits, sterile=True)
            key = _genotype_key(pairs, offspring_sex)
            if key in raw:
                old_p, old_pairs, old_sex, old_pheno = raw[key]
                raw[key] = (old_p + p, old_pairs, old_sex, old_pheno)
            else:
                raw[key] = (p, pairs, offspring_sex, phenotype)

    viable_raw = total_raw - lethal_prob
    pending: list[tuple[Fraction, tuple[ChromosomePair, ...], str, Phenotype]] = []
    if viable_raw > 0:
        for _key, (p, pairs, sex, phenotype) in raw.items():
            pending.append((p / viable_raw, pairs, sex, phenotype))

    colors = assign_colors([item[3].key for item in pending])
    offspring_list: list[F1Individual] = []
    for probability, pairs, sex, phenotype in pending:
        color_name, color_hex = colors.get(phenotype.key, ("灰", "6B6B6B"))
        offspring_list.append(
            F1Individual(
                pairs=pairs,
                probability=probability,
                phenotype_hints=_hint_symbols(phenotype),
                sex=sex,
                phenotype_text=phenotype.text(),
                phenotype_key=phenotype.key,
                color_name=color_name,
                color_hex=color_hex,
                sterile=phenotype.sterile,
            )
        )

    offspring_list.sort(key=lambda item: (-item.probability, item.sex, str(item.pairs)))

    warnings = _cross_warnings(female, male, rf)
    warnings.extend(confusable_warnings([item[3] for item in pending]))

    return CrossResult(
        offspring=tuple(offspring_list),
        raw_total_before_lethal=total_raw,
        lethal_discarded_probability=lethal_prob,
        female_raw=female.raw,
        male_raw=male.raw,
        rf=rf,
        lethals_used=lethals,
        warnings=tuple(dict.fromkeys(warnings)),
    )
