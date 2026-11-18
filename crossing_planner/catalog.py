"""Marker and balancer catalog.

Phenotype labels are short original phrases. FlyBase accession numbers identify
balancers; cytology and atlas prose are not copied.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Marker:
    """One scorable (or lethal-only) allele family."""

    symbol: str
    chromosome: str
    dominance: str  # "dominant" or "recessive"
    homozygous: str  # "lethal", "viable", or "female_sterile"
    effects: tuple[tuple[str, str], ...]  # (axis, label)
    aliases: tuple[str, ...] = ()
    stage: str = "adult"  # "adult" or "larva"
    scorable: bool = True
    male_lethal: bool = False


@dataclass(frozen=True)
class Balancer:
    """A rearranged chromosome that can be recognized in genotype text."""

    symbol: str
    chromosomes: tuple[str, ...]
    dominant_markers: tuple[str, ...]
    recessive_markers: tuple[str, ...]
    homozygous: str  # "lethal" or "female_sterile"
    male_viable: bool
    suppression: str  # "chromosome" or "proximal"
    flybase_id: str
    aliases: tuple[str, ...] = ()
    note: str = ""


def _m(
    symbol: str,
    chromosome: str,
    dominance: str,
    homozygous: str,
    effects: tuple[tuple[str, str], ...],
    aliases: tuple[str, ...] = (),
    *,
    stage: str = "adult",
    scorable: bool = True,
    male_lethal: bool = False,
) -> Marker:
    return Marker(
        symbol=symbol,
        chromosome=chromosome,
        dominance=dominance,
        homozygous=homozygous,
        effects=effects,
        aliases=aliases,
        stage=stage,
        scorable=scorable,
        male_lethal=male_lethal,
    )


# Axes: 眼色 眼形 翅 刚毛 体色 幼虫 附肢
MARKERS: tuple[Marker, ...] = (
    # Bristles
    _m("Bl", "2", "dominant", "lethal", (("刚毛", "刚毛钝短"),)),
    _m("dm", "X", "recessive", "viable", (("刚毛", "刚毛细小"),)),
    _m("emc", "3", "recessive", "viable", (("刚毛", "额外刚毛"),)),
    _m("f", "X", "recessive", "viable", (("刚毛", "刚毛分叉"),)),
    _m("h", "3", "recessive", "viable", (("刚毛", "额外小刚毛"),)),
    _m("Hu", "3", "dominant", "lethal", (("刚毛", "肩鬃增多"),)),
    _m("jv", "3", "recessive", "viable", (("刚毛", "刚毛如矛"),)),
    _m("M", "multi", "dominant", "lethal", (("刚毛", "刚毛细短"),), scorable=False),
    _m("Pin", "2", "dominant", "lethal", (("刚毛", "胸部刚毛聚拢"),)),
    _m("sc", "X", "recessive", "viable", (("刚毛", "缺盾片鬃"),)),
    _m("Sco", "2", "dominant", "lethal", (("刚毛", "盾片鬃缺失"),)),
    _m("sn", "X", "recessive", "viable", (("刚毛", "刚毛卷曲"),)),
    _m("Sp", "2", "dominant", "lethal", (("刚毛", "胸侧多刚毛"),)),
    _m("ss", "3", "recessive", "viable", (("刚毛", "刚毛缩短"),)),
    _m("Sb", "3", "dominant", "lethal", (("刚毛", "短刚毛"),)),
    # Wings
    _m("ci", "4", "recessive", "viable", (("翅", "肘脉中断"),)),
    _m("ct", "X", "recessive", "viable", (("翅", "翅端切削"),)),
    _m("cu", "3", "recessive", "viable", (("翅", "翅上卷"),)),
    _m("Cy", "2", "dominant", "lethal", (("翅", "卷翅"),)),
    _m("D", "3", "dominant", "lethal", (("翅", "翅展开"), ("刚毛", "缺背侧鬃"))),
    _m("dp", "2", "recessive", "viable", (("翅", "翅截短"),)),
    _m("Ly", "3", "dominant", "lethal", (("翅", "翅缘削去"),)),
    _m("m", "X", "recessive", "viable", (("翅", "小翅"),)),
    _m("N", "X", "dominant", "lethal", (("翅", "翅缘缺刻"),), male_lethal=True),
    _m("ri", "3", "recessive", "viable", (("翅", "径脉不全"),)),
    _m("ve", "3", "recessive", "viable", (("翅", "翅脉缩短"),)),
    _m("vg", "2", "recessive", "viable", (("翅", "残翅"),)),
    _m("wg", "2", "recessive", "viable", (("翅", "无翅"),)),
    _m("Xa", "2", "dominant", "lethal", (("翅", "翅呈X形"),), ("apXa",)),
    _m("Ser", "3", "dominant", "lethal", (("翅", "翅缘锯齿"),)),
    # Appendages
    _m("al", "2", "recessive", "viable", (("附肢", "无触角芒"),)),
    _m("Antp", "3", "dominant", "lethal", (("附肢", "触角成足"),)),
    _m("Dll", "2", "dominant", "lethal", (("附肢", "远端附肢缩短"),)),
    _m("Scr", "3", "dominant", "lethal", (("附肢", "性梳减少"),)),
    _m("th", "3", "recessive", "viable", (("附肢", "触角芒变细"),)),
    _m("Ubx", "3", "dominant", "lethal", (("附肢", "平衡棒成翅"),)),
    # Eye shape
    _m("B", "X", "dominant", "viable", (("眼形", "棒眼"),)),
    _m("Dfd", "3", "dominant", "lethal", (("眼形", "眼畸形"),)),
    _m("Dr", "3", "dominant", "lethal", (("眼形", "眼缩小"),)),
    _m("eyD", "4", "dominant", "lethal", (("眼形", "无眼"),), ("ey^D",)),
    _m("Gl", "3", "dominant", "lethal", (("眼形", "眼面发粘"),)),
    _m("gl", "3", "recessive", "viable", (("眼形", "玻璃眼"),)),
    _m("Gla", "2", "dominant", "lethal", (("眼形", "釉眼"),), ("wgGla",)),
    _m("If", "2", "dominant", "lethal", (("眼形", "小眼不规则"),), ("KrIf-1",)),
    _m("L", "2", "dominant", "lethal", (("眼形", "叶状眼"),)),
    _m("oc", "X", "recessive", "female_sterile", (("眼形", "无单眼"),)),
    _m("ro", "3", "recessive", "viable", (("眼形", "眼面粗糙"),)),
    _m("Roi", "2", "dominant", "lethal", (("眼形", "眼面粗糙"),)),
    _m("ru", "3", "recessive", "viable", (("眼形", "眼面粗糙"),)),
    # Eye color
    _m("bw", "2", "recessive", "viable", (("眼色", "棕眼"),)),
    _m("ca", "3", "recessive", "viable", (("眼色", "酒红眼"),)),
    _m("cn", "2", "recessive", "viable", (("眼色", "朱红眼"),)),
    _m("g", "X", "recessive", "viable", (("眼色", "石榴红眼"),)),
    _m("Hn", "3", "dominant", "viable", (("眼色", "深色眼"),)),
    _m("kar", "3", "recessive", "viable", (("眼色", "朱红眼"),)),
    _m("lt", "2", "recessive", "viable", (("眼色", "浅色眼"),)),
    _m("Mé", "3", "dominant", "lethal", (("眼色", "云纹眼"),), ("Me",)),
    _m("p", "3", "recessive", "viable", (("眼色", "粉红眼"),)),
    _m("pr", "2", "recessive", "viable", (("眼色", "紫眼"),)),
    _m("ry", "3", "recessive", "viable", (("眼色", "玫瑰眼"),)),
    _m("se", "3", "recessive", "viable", (("眼色", "深棕眼"),)),
    _m("st", "3", "recessive", "viable", (("眼色", "猩红眼"),)),
    _m("v", "X", "recessive", "viable", (("眼色", "朱红眼"),)),
    _m("w", "X", "recessive", "viable", (("眼色", "白眼"),)),
    # Body
    _m("Bc", "2", "dominant", "lethal", (("体色", "黑细胞"),)),
    _m("e", "3", "recessive", "viable", (("体色", "乌木体"),)),
    _m("ptg", "X", "recessive", "viable", (("体色", "胸部五边形斑"),)),
    _m("sp", "2", "recessive", "viable", (("体色", "翅基黑点"),)),
    _m("Tb", "3", "dominant", "lethal", (("幼虫", "幼虫粗短"),), stage="larva"),
    _m("y", "X", "recessive", "viable", (("体色", "黄体"),)),
    # Referenced by balancer genotypes, not a separate atlas plate
    _m("ds", "2", "recessive", "viable", (("翅", "翅短"),)),
)


BALANCERS: tuple[Balancer, ...] = (
    Balancer("FM6", ("X",), ("B",), ("y", "sc", "dm"), "lethal", True, "chromosome", "FBba0000003"),
    Balancer("FM7", ("X",), ("B",), (), "lethal", True, "chromosome", "", ("In(1)FM7",), "未写明亚型时只认显性棒眼"),
    Balancer("FM7a", ("X",), ("B",), ("y", "sc", "w", "v"), "lethal", True, "chromosome", "FBba0000007", note="雄性可育"),
    Balancer("FM7c", ("X",), ("B",), ("y", "sc", "w", "sn", "v", "g"), "female_sterile", True, "chromosome", "FBba0000009", note="纯合雌性可活但不育"),
    Balancer("FM7d", ("X",), ("B",), ("y", "sc"), "lethal", True, "chromosome", "FBba0000216", note="雄性和杂合雌性可育"),
    Balancer("FM7i", ("X",), ("B",), ("y", "sc", "w", "oc", "ptg"), "lethal", True, "chromosome", "FBba0000226", note="雄性可育"),
    Balancer("Binsn", ("X",), ("B",), (), "lethal", True, "chromosome", "", note="经典 X 平衡子，图谱平衡子章未单列"),
    Balancer("CyO", ("2",), ("Cy",), ("dp", "pr", "cn"), "lethal", False, "chromosome", "FBba0000025", ("In(2LR)O",)),
    Balancer("SM1", ("2",), ("Cy",), ("al", "cn", "sp"), "lethal", False, "chromosome", "FBba0000037"),
    Balancer("SM5", ("2",), ("Cy",), ("al", "ds", "lt", "cn", "sp"), "lethal", False, "chromosome", "FBba0000038", note="杂合体生活力低于 SM1"),
    Balancer("SM6", ("2",), ("Cy",), (), "lethal", False, "chromosome", "", note="未写明 a/b 时只认卷翅"),
    Balancer("SM6a", ("2",), ("Cy",), ("al", "dp", "cn", "sp"), "lethal", False, "chromosome", "FBba0000039"),
    Balancer("SM6b", ("2",), ("Cy", "Roi"), ("al", "dp", "cn", "sp"), "lethal", False, "chromosome", "FBba0000040"),
    Balancer(
        "CyO-TM9",
        ("2", "3"),
        ("Cy", "Sb"),
        ("cn", "dp", "e", "pr", "st", "th"),
        "lethal",
        False,
        "chromosome",
        "FBab0023073",
        ("T(2;3)CyO-TM9", "T23CyO-TM9"),
        "同时抑制 2 号和 3 号交换",
    ),
    Balancer("MRS", ("3",), ("Sb", "M"), ("ry",), "lethal", False, "proximal", "FBba0000068", note="只平衡 3 号近端；Minute 几乎看不出"),
    Balancer("MKRS", ("3",), ("Sb", "M"), ("kar", "ry"), "lethal", False, "proximal", "FBba0000066", note="只平衡 3 号近端；Minute 几乎看不出"),
    Balancer("TM1", ("3",), ("Mé",), ("ri",), "lethal", False, "chromosome", "FBba0000042"),
    Balancer("TM2", ("3",), ("Ubx",), ("emc", "e"), "lethal", False, "chromosome", "FBba0000062"),
    Balancer("TM3", ("3",), ("Ser",), ("ri", "e"), "lethal", False, "chromosome", "FBba0000047", note="显性标记是 Ser；Sb 仅在基因型写明时计入"),
    Balancer("TM6", ("3",), ("Hn", "Ubx"), ("e",), "lethal", False, "chromosome", "FBba0000056"),
    Balancer("TM6B", ("3",), ("Hu",), ("e",), "lethal", False, "chromosome", "FBba0000057", note="常写 Tb，偶尔写 D；未写明则不自动加上"),
    Balancer("TM6C", ("3",), (), ("e",), "lethal", False, "chromosome", "FBba0000071", note="常带 Sb 或 Tb，未写明则不自动加上"),
    Balancer("TM8", ("3",), ("Sb",), ("th", "st", "e"), "lethal", False, "chromosome", "FBba0000060"),
    Balancer("TM9", ("3",), ("Sb",), ("th", "st", "e"), "lethal", False, "chromosome", "FBba0000061", note="标记组成同 TM8，断点不同"),
)


# Markers that are hard to tell apart when no second trait separates the classes.
CONFUSABLE_GROUPS: tuple[tuple[str, ...], ...] = (
    ("cn", "st", "v", "kar"),
    ("se", "ca", "ry"),
    ("Sb", "Bl"),
    ("Ser", "N", "Ly", "ct"),
    ("Cy", "cu"),
    ("ro", "ru", "Roi"),
    ("gl", "Gl", "Gla"),
)


def _register(table: dict[str, object], key: str, value: object) -> None:
    if key in table and table[key] is not value:
        raise RuntimeError(f"duplicate catalog token {key}")
    table[key] = value


_MARKERS: dict[str, Marker] = {}
for _marker in MARKERS:
    _register(_MARKERS, _marker.symbol, _marker)
    for _alias in _marker.aliases:
        _register(_MARKERS, _alias, _marker)

_BALANCERS: dict[str, Balancer] = {}
for _balancer in BALANCERS:
    _register(_BALANCERS, _balancer.symbol, _balancer)
    for _alias in _balancer.aliases:
        _register(_BALANCERS, _alias, _balancer)

_CONFUSABLE: dict[str, frozenset[str]] = {}
for _group in CONFUSABLE_GROUPS:
    others = frozenset(_group)
    for _symbol in _group:
        _CONFUSABLE[_symbol] = others


def split_tokens(locus: str) -> list[str]:
    """Split a locus label. Commas keep compound balancer names as pieces."""
    text = locus.strip()
    if not text:
        return []
    if "," in text:
        return [part.strip() for part in text.split(",") if part.strip()]
    return [text]


def is_wild_type_token(token: str) -> bool:
    text = token.strip()
    if text in {"+", "Y"}:
        return True
    return "[+]" in text


def _strip_superscript(token: str) -> str:
    stem = token.rstrip("0123456789")
    return stem if stem and stem != token else token


def lookup_balancer(token: str) -> Balancer | None:
    if is_wild_type_token(token):
        return None
    found = _BALANCERS.get(token)
    if found is not None:
        return found
    return _BALANCERS.get(_strip_superscript(token))


def lookup_marker(token: str) -> Marker | None:
    if is_wild_type_token(token):
        return None
    if lookup_balancer(token) is not None:
        return None
    found = _MARKERS.get(token)
    if found is not None:
        return found
    return _MARKERS.get(_strip_superscript(token))


_MARKER_SYMBOLS_DESC: tuple[str, ...] = tuple(
    sorted(_MARKERS, key=len, reverse=True)
)


def split_fused_markers(token: str) -> tuple[str, ...] | None:
    """Split a compacted allele string such as 'yw' or 'y1w1118'.

    Exact markers and balancers are left untouched. A fusion is accepted only
    when the whole token is a sequence of known symbols, so 'wg' stays wingless
    and 'UAS-SytHA' is not cut apart.
    """
    if is_wild_type_token(token) or lookup_balancer(token) is not None:
        return None
    if token in _MARKERS or _strip_superscript(token) in _MARKERS:
        return None
    if not token or not token[0].isalpha():
        return None
    if any(not (char.isalpha() or char.isdigit()) for char in token):
        return None

    parts: list[str] = []
    index = 0
    while index < len(token):
        hit = ""
        for symbol in _MARKER_SYMBOLS_DESC:
            if token.startswith(symbol, index):
                hit = symbol
                break
        if not hit:
            return None
        index += len(hit)
        while index < len(token) and token[index].isdigit():
            index += 1
        parts.append(hit)
    if len(parts) < 2:
        return None
    return tuple(parts)


def markers_confusable(left: str, right: str) -> bool:
    if left == right:
        return False
    group = _CONFUSABLE.get(left)
    return group is not None and right in group


CHROM_BY_INDEX: dict[int, str] = {0: "X", 1: "2", 2: "3", 3: "4"}
