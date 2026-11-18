"""Recognize balancers and markers, then score visible phenotypes."""

from __future__ import annotations

from dataclasses import dataclass

from .catalog import (
    CHROM_BY_INDEX,
    Balancer,
    Marker,
    is_wild_type_token,
    lookup_balancer,
    lookup_marker,
    markers_confusable,
    split_fused_markers,
    split_tokens,
)
from .genotype import ChromosomePair, Haplotype, SexChromRole


_AXIS_ORDER = ("眼色", "眼形", "翅", "刚毛", "体色", "幼虫", "附肢")

_COLOR_PALETTE: tuple[tuple[str, str], ...] = (
    ("蓝", "1F4E79"),
    ("橙", "C65911"),
    ("绿", "548235"),
    ("紫", "7030A0"),
    ("红", "C00000"),
    ("金", "BF8F00"),
    ("青", "0D7377"),
    ("棕", "833C0C"),
    ("品红", "A1256B"),
    ("藏青", "2F5496"),
    ("橄榄", "7F6000"),
    ("石板", "5B5B5B"),
)


@dataclass(frozen=True)
class HapContent:
    """Markers and balancers present on one homolog."""

    balancers: tuple[Balancer, ...]
    markers: tuple[Marker, ...]
    extra_tokens: frozenset[str]


@dataclass(frozen=True)
class TraitState:
    axis: str
    label: str
    markers: frozenset[str]
    stage: str


@dataclass(frozen=True)
class Phenotype:
    traits: tuple[TraitState, ...]
    sterile: bool

    @property
    def key(self) -> tuple[tuple[str, str], ...]:
        return tuple((trait.axis, trait.label) for trait in self.traits)

    def text(self) -> str:
        if not self.traits:
            return "无明显突变性状"
        return " · ".join(f"{trait.axis}:{trait.label}" for trait in self.traits)


def haplotype_tokens(haplotype: Haplotype) -> list[str]:
    tokens: list[str] = []
    for locus in haplotype.loci:
        tokens.extend(split_tokens(locus))
    return tokens


def read_haplotype(haplotype: Haplotype) -> HapContent:
    balancers: list[Balancer] = []
    markers: list[Marker] = []
    extras: set[str] = set()
    seen_bal: set[str] = set()
    seen_marker: set[str] = set()

    for token in haplotype_tokens(haplotype):
        if is_wild_type_token(token):
            continue
        balancer = lookup_balancer(token)
        if balancer is not None:
            if balancer.symbol not in seen_bal:
                seen_bal.add(balancer.symbol)
                balancers.append(balancer)
            for symbol in (*balancer.dominant_markers, *balancer.recessive_markers):
                marker = lookup_marker(symbol)
                if marker is not None and marker.symbol not in seen_marker:
                    seen_marker.add(marker.symbol)
                    markers.append(marker)
            continue
        marker = lookup_marker(token)
        if marker is not None:
            if marker.symbol not in seen_marker:
                seen_marker.add(marker.symbol)
                markers.append(marker)
            continue
        fused = split_fused_markers(token)
        if fused:
            for symbol in fused:
                fused_marker = lookup_marker(symbol)
                if fused_marker is not None and fused_marker.symbol not in seen_marker:
                    seen_marker.add(fused_marker.symbol)
                    markers.append(fused_marker)
            continue
        extras.add(token)

    return HapContent(
        balancers=tuple(balancers),
        markers=tuple(markers),
        extra_tokens=frozenset(extras),
    )


def _is_sex_pair(pair: ChromosomePair, index: int) -> bool:
    return index == 0 or pair.role in (SexChromRole.XX, SexChromRole.XY)


def pair_status(
    maternal: Haplotype,
    paternal: Haplotype,
    *,
    sex: str,
    index: int,
    extra_lethals: frozenset[str],
) -> str:
    """Return 'lethal', 'sterile', or 'viable' for one chromosome pair."""
    if maternal.is_y and paternal.is_y:
        return "lethal"

    mat = read_haplotype(maternal)
    pat = read_haplotype(paternal)
    hemizygous = _is_sex_pair_hemizygous(maternal, paternal, sex, index)

    status = "viable"
    mat_bal = {item.symbol: item for item in mat.balancers}
    pat_bal = {item.symbol: item for item in pat.balancers}
    for symbol, balancer in mat_bal.items():
        if symbol not in pat_bal:
            continue
        if balancer.homozygous == "lethal":
            return "lethal"
        if balancer.homozygous == "female_sterile" and sex == "female":
            status = "sterile"

    mat_markers = {item.symbol: item for item in mat.markers}
    pat_markers = {item.symbol: item for item in pat.markers}

    if hemizygous:
        present = pat_markers if maternal.is_y else mat_markers
        for marker in present.values():
            if marker.male_lethal and sex == "male":
                return "lethal"
    else:
        for symbol, marker in mat_markers.items():
            if symbol not in pat_markers:
                continue
            if marker.homozygous == "lethal":
                return "lethal"
            if marker.homozygous == "female_sterile" and sex == "female":
                status = "sterile"

    if extra_lethals and not hemizygous:
        mat_extra = _lethal_names(maternal, extra_lethals)
        pat_extra = _lethal_names(paternal, extra_lethals)
        if mat_extra & pat_extra:
            return "lethal"
    return status


def _is_sex_pair_hemizygous(
    maternal: Haplotype, paternal: Haplotype, sex: str, index: int
) -> bool:
    if sex != "male":
        return False
    if index != 0:
        return False
    return maternal.is_y or paternal.is_y


def _lethal_names(haplotype: Haplotype, extra: frozenset[str]) -> set[str]:
    found: set[str] = set()
    for token in haplotype_tokens(haplotype):
        if token in extra:
            found.add(token)
        for name in extra:
            if token == name or token.startswith(name + " ") or token.startswith(name + ","):
                found.add(name)
    return found


def zygote_status(
    maternal_haps: tuple[Haplotype, ...],
    paternal_haps: tuple[Haplotype, ...],
    sex: str,
    extra_lethals: frozenset[str],
) -> str:
    """Worst status across chromosome pairs: lethal before sterile."""
    saw_sterile = False
    for index, (maternal, paternal) in enumerate(zip(maternal_haps, paternal_haps)):
        status = pair_status(
            maternal,
            paternal,
            sex=sex,
            index=index,
            extra_lethals=extra_lethals,
        )
        if status == "lethal":
            return "lethal"
        if status == "sterile":
            saw_sterile = True
    return "sterile" if saw_sterile else "viable"


def _visible_markers(maternal: HapContent, paternal: HapContent, hemizygous: bool) -> list[Marker]:
    if hemizygous:
        if _content_is_y_only(maternal) and not _content_is_y_only(paternal):
            content = paternal
        else:
            content = maternal
        return [marker for marker in content.markers if marker.scorable]

    visible: list[Marker] = []
    seen: set[str] = set()
    mat = {item.symbol: item for item in maternal.markers}
    pat = {item.symbol: item for item in paternal.markers}
    for symbol, marker in mat.items():
        if not marker.scorable or symbol in seen:
            continue
        if marker.dominance == "dominant" or symbol in pat:
            seen.add(symbol)
            visible.append(marker)
    for symbol, marker in pat.items():
        if not marker.scorable or symbol in seen:
            continue
        if marker.dominance == "dominant":
            seen.add(symbol)
            visible.append(marker)
    return visible


def _content_is_y_only(content: HapContent) -> bool:
    return not content.balancers and not content.markers and not content.extra_tokens


def phenotype_of(
    homolog_pairs: tuple[tuple[Haplotype, Haplotype], ...],
    sex: str,
) -> Phenotype:
    """Visible trait states for one viable genotype."""
    by_axis: dict[str, list[tuple[str, str, str]]] = {}
    sterile = False

    for index, (maternal, paternal) in enumerate(homolog_pairs):
        mat = read_haplotype(maternal)
        pat = read_haplotype(paternal)
        hemizygous = _is_sex_pair_hemizygous(maternal, paternal, sex, index)
        for marker in _visible_markers(mat, pat, hemizygous):
            for axis, label in marker.effects:
                by_axis.setdefault(axis, []).append((label, marker.symbol, marker.stage))
        if sex == "female":
            shared = {item.symbol for item in mat.markers} & {item.symbol for item in pat.markers}
            for marker in mat.markers:
                if marker.symbol in shared and marker.homozygous == "female_sterile":
                    sterile = True
            for balancer in mat.balancers:
                if balancer.homozygous == "female_sterile" and any(
                    item.symbol == balancer.symbol for item in pat.balancers
                ):
                    sterile = True

    traits = _traits_from_axes(by_axis)
    return Phenotype(traits=tuple(traits), sterile=sterile)


def _traits_from_axes(by_axis: dict[str, list[tuple[str, str, str]]]) -> list[TraitState]:
    traits: list[TraitState] = []
    for axis in _AXIS_ORDER:
        entries = by_axis.get(axis)
        if not entries:
            continue
        if axis == "眼色":
            label, symbols, stage = _resolve_eye_color(entries)
            traits.append(
                TraitState(axis=axis, label=label, markers=frozenset(symbols), stage=stage)
            )
            continue
        labels = sorted({label for label, _symbol, _stage in entries})
        symbols = {symbol for _label, symbol, _stage in entries}
        stage = "larva" if any(item[2] == "larva" for item in entries) else "adult"
        traits.append(
            TraitState(
                axis=axis,
                label="+".join(labels),
                markers=frozenset(symbols),
                stage=stage,
            )
        )
    return traits


def _resolve_eye_color(
    entries: list[tuple[str, str, str]],
) -> tuple[str, frozenset[str], str]:
    symbols = {symbol for _label, symbol, _stage in entries}
    labels = {symbol: label for label, symbol, _stage in entries}
    if "w" in symbols:
        return "白眼", frozenset({"w"}), "adult"
    if "bw" in symbols and "st" in symbols:
        return "白眼", frozenset({"bw", "st"}), "adult"
    if "bw" in symbols and "cn" in symbols:
        return "白眼", frozenset({"bw", "cn"}), "adult"
    ordered = sorted(labels[symbol] for symbol in symbols)
    # Identical labels (cn and v both 朱红眼) collapse.
    unique = list(dict.fromkeys(ordered))
    return "+".join(unique), frozenset(symbols), "adult"


def assign_colors(
    keys: list[tuple[tuple[str, str], ...]],
) -> dict[tuple[tuple[str, str], ...], tuple[str, str]]:
    """Map each phenotype key to (color name, hex). Identical keys share a color."""
    wild = ("灰", "6B6B6B")
    distinct = sorted({key for key in keys if key})
    mapping: dict[tuple[tuple[str, str], ...], tuple[str, str]] = {(): wild}
    for index, key in enumerate(distinct):
        if index < len(_COLOR_PALETTE):
            mapping[key] = _COLOR_PALETTE[index]
        else:
            mapping[key] = (f"色{index + 1}", _extra_hex(index))
    return mapping


def _extra_hex(index: int) -> str:
    hue = (index * 47) % 360
    saturation = 0.55
    lightness = 0.38
    red, green, blue = _hsl_to_rgb(hue, saturation, lightness)
    return f"{red:02X}{green:02X}{blue:02X}"


def _hsl_to_rgb(hue: float, saturation: float, lightness: float) -> tuple[int, int, int]:
    c = (1 - abs(2 * lightness - 1)) * saturation
    hp = hue / 60
    x = c * (1 - abs(hp % 2 - 1))
    if hp < 1:
        r1, g1, b1 = c, x, 0
    elif hp < 2:
        r1, g1, b1 = x, c, 0
    elif hp < 3:
        r1, g1, b1 = 0, c, x
    elif hp < 4:
        r1, g1, b1 = 0, x, c
    elif hp < 5:
        r1, g1, b1 = x, 0, c
    else:
        r1, g1, b1 = c, 0, x
    m = lightness - c / 2
    return (round((r1 + m) * 255), round((g1 + m) * 255), round((b1 + m) * 255))


def confusable_warnings(phenotypes: list[Phenotype]) -> tuple[str, ...]:
    """Warn when two classes differ by a single hard-to-score trait."""
    unique: dict[tuple[tuple[str, str], ...], Phenotype] = {}
    for phenotype in phenotypes:
        unique.setdefault(phenotype.key, phenotype)
    items = list(unique.values())
    notes: list[str] = []
    seen: set[tuple[str, str]] = set()
    for left_index, left in enumerate(items):
        left_axes = {trait.axis: trait for trait in left.traits}
        for right in items[left_index + 1 :]:
            right_axes = {trait.axis: trait for trait in right.traits}
            if set(left_axes) != set(right_axes):
                continue
            differed: list[tuple[TraitState, TraitState]] = []
            for axis, trait in left_axes.items():
                other = right_axes[axis]
                if trait.label != other.label:
                    differed.append((trait, other))
            if len(differed) != 1:
                continue
            a, b = differed[0]
            if not any(markers_confusable(x, y) for x in a.markers for y in b.markers):
                continue
            token = tuple(sorted((a.label, b.label)))
            if token in seen:
                continue
            seen.add(token)
            notes.append(f"{a.axis}：{token[0]} 与 {token[1]} 容易看混，这两类只差这一项。")
    return tuple(notes)


def chromosome_suppressed(pair: ChromosomePair) -> bool:
    """True when either homolog carries a whole-chromosome balancer."""
    for haplotype in pair.homologs():
        for balancer in read_haplotype(haplotype).balancers:
            if balancer.suppression == "chromosome" and len(balancer.chromosomes) == 1:
                return True
    return False


def translocation_indices(pairs: tuple[ChromosomePair, ...]) -> set[int]:
    """Chromosome indexes whose recombination is suppressed by a multi-chromosome balancer."""
    found = False
    targets: set[str] = set()
    for pair in pairs:
        for haplotype in pair.homologs():
            for balancer in read_haplotype(haplotype).balancers:
                if len(balancer.chromosomes) > 1 and balancer.suppression == "chromosome":
                    found = True
                    targets.update(balancer.chromosomes)
    if not found:
        return set()
    return {index for index, name in CHROM_BY_INDEX.items() if name in targets and index < len(pairs)}


def has_proximal_balancer(pairs: tuple[ChromosomePair, ...]) -> bool:
    for pair in pairs:
        for haplotype in pair.homologs():
            for balancer in read_haplotype(haplotype).balancers:
                if balancer.suppression == "proximal":
                    return True
    return False


def balancer_chromosome_warnings(pairs: tuple[ChromosomePair, ...]) -> tuple[str, ...]:
    notes: list[str] = []
    for index, pair in enumerate(pairs):
        expected = CHROM_BY_INDEX.get(index)
        if expected is None:
            continue
        for haplotype in pair.homologs():
            for balancer in read_haplotype(haplotype).balancers:
                if len(balancer.chromosomes) == 1 and balancer.chromosomes[0] != expected:
                    notes.append(
                        f"{balancer.symbol} 通常在 {balancer.chromosomes[0]} 染色体上，"
                        f"这里写在第 {index + 1} 对。"
                    )
    return tuple(dict.fromkeys(notes))
