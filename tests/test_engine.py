"""Engine tests for markers, balancers, colors, and multi-step schemes."""

from __future__ import annotations

import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from crossing_planner.catalog import MARKERS, split_fused_markers
from crossing_planner.cross import cross
from crossing_planner.format_frac import chromosome_options, format_cross_lines, match_offspring
from crossing_planner.gametes import generate_gametes
from crossing_planner.genotype import parse_genotype
from crossing_planner.i18n import PHENOTYPE_EN, set_language
from crossing_planner.pptx_export import export_crossing_pptx, export_scheme_pptx
from crossing_planner.gui import _ref_label
from crossing_planner.project_file import dumps_project, freeze, loads_project
from crossing_planner.scheme import CrossStep, ParentRef, SchemeError, Stock, run_scheme


def _cross(female: str, male: str, rf: str = "0"):
    return cross(
        parse_genotype(female, "female"),
        parse_genotype(male, "male"),
        rf=Fraction(rf),
    )


class SampleCrossTests(unittest.TestCase):
    def test_sample_counts_and_colors(self) -> None:
        result = _cross(
            "w ; Sp / CyO ; TM3, Sb / MKRS",
            "w / Y ; + / CyO ; + / TM6B, Tb",
        )
        self.assertEqual(len(result.offspring), 24)
        self.assertEqual(result.lethal_discarded_probability, Fraction(1, 4))
        self.assertEqual(sum(ind.probability for ind in result.offspring), 1)
        self.assertTrue(all(ind.probability == Fraction(1, 24) for ind in result.offspring))
        self.assertTrue(all("白眼" in ind.phenotype_text for ind in result.offspring))
        compact = [ind.pairs[1].homologs()[0].display() + "/" + ind.pairs[1].homologs()[1].display() for ind in result.offspring]
        self.assertFalse(any("CyO/CyO" in text or "CyO / CyO" in text for text in compact))
        curly = [ind for ind in result.offspring if "卷翅" in ind.phenotype_text]
        self.assertTrue(curly)
        shared = {ind.color_hex for ind in result.offspring if ind.phenotype_key == curly[0].phenotype_key}
        self.assertEqual(len(shared), 1)


class BalancerTests(unittest.TestCase):
    def test_cyo_over_sm6a_is_lethal(self) -> None:
        result = _cross("+ ; CyO / +", "+ / Y ; SM6a / +")
        self.assertEqual(len(result.offspring), 6)
        for ind in result.offspring:
            text = " ".join(hap.display() for hap in ind.pairs[1].homologs())
            self.assertFalse("CyO" in text and "SM6a" in text)

    def test_fm7c_homozygous_females_stay_and_are_sterile(self) -> None:
        result = _cross("FM7c / +", "FM7c / Y")
        self.assertEqual(sum(ind.probability for ind in result.offspring), 1)
        self.assertEqual(result.lethal_discarded_probability, 0)
        sterile = [ind for ind in result.offspring if ind.sterile]
        self.assertEqual(len(sterile), 1)
        self.assertEqual(sterile[0].sex, "female")
        self.assertIn("白眼", sterile[0].phenotype_text)
        self.assertIn("棒眼", sterile[0].phenotype_text)
        hetero = [
            ind
            for ind in result.offspring
            if ind.sex == "female" and not ind.sterile and "棒眼" in ind.phenotype_text
        ]
        self.assertTrue(hetero)
        self.assertNotIn("白眼", hetero[0].phenotype_text)

    def test_notch_males_die(self) -> None:
        result = _cross("N / +", "+ / Y")
        males = [ind for ind in result.offspring if ind.sex == "male"]
        self.assertTrue(males)
        self.assertTrue(all("翅缘缺刻" not in ind.phenotype_text for ind in males))
        self.assertGreater(result.lethal_discarded_probability, 0)

    def test_balancer_blocks_recombination(self) -> None:
        female = parse_genotype("+ ; a b / CyO", "female")
        gametes = generate_gametes(female, Fraction(1, 2))
        chrom2 = {gamete.chromosomes[1].display() for gamete in gametes}
        self.assertEqual(chrom2, {"a b", "CyO"})

    def test_free_chromosome_still_recombines(self) -> None:
        female = parse_genotype("a b / + + ; CyO / +", "female")
        gametes = generate_gametes(female, Fraction(1, 2))
        chrom_x = {gamete.chromosomes[0].display() for gamete in gametes}
        chrom2 = {gamete.chromosomes[1].display() for gamete in gametes}
        self.assertIn("a +", chrom_x)
        self.assertEqual(chrom2, {"CyO", "+"})

    def test_male_autosomes_do_not_recombine(self) -> None:
        male = parse_genotype("+ / Y ; a b / + +", "male")
        gametes = generate_gametes(male, Fraction(1, 2))
        chrom2 = {gamete.chromosomes[1].display() for gamete in gametes}
        self.assertEqual(chrom2, {"a b", "+ +"})

    def test_translocation_symbol_parses_and_suppresses_both(self) -> None:
        female = parse_genotype("+ ; T(2;3)CyO-TM9 / a b ; c d / + +", "female")
        self.assertEqual(len(female.pairs), 3)
        gametes = generate_gametes(female, Fraction(1, 2))
        chrom2 = {gamete.chromosomes[1].display() for gamete in gametes}
        chrom3 = {gamete.chromosomes[2].display() for gamete in gametes}
        self.assertEqual(chrom2, {"T(2;3)CyO-TM9", "a b"})
        self.assertEqual(chrom3, {"c d", "+ +"})


class PhenotypeTests(unittest.TestCase):
    def test_white_masks_cinnabar(self) -> None:
        result = _cross("w ; cn / cn", "w / Y ; cn / cn")
        self.assertTrue(result.offspring)
        self.assertTrue(all(ind.phenotype_text.startswith("眼色:白眼") for ind in result.offspring))
        self.assertTrue(all("朱红眼" not in ind.phenotype_text for ind in result.offspring))

    def test_brown_and_scarlet_are_white(self) -> None:
        result = _cross("+ ; bw / bw ; st / st", "+ / Y ; bw / bw ; st / st")
        self.assertTrue(all("眼色:白眼" in ind.phenotype_text for ind in result.offspring))

    def test_same_visible_trait_shares_a_color(self) -> None:
        result = _cross("+ ; CyO / + ; TM3, Sb / +", "+ / Y ; + / + ; + / +")
        colors: dict[tuple, set[str]] = {}
        for ind in result.offspring:
            colors.setdefault(ind.phenotype_key, set()).add(ind.color_hex)
        self.assertTrue(colors)
        self.assertTrue(all(len(values) == 1 for values in colors.values()))
        tb = _cross("+ ; + / + ; TM6B, Tb / +", "+ / Y ; + / + ; + / +")
        plain = {ind.phenotype_key for ind in result.offspring}
        tubby = {
            ind.phenotype_key
            for ind in tb.offspring
            if any(axis == "幼虫" for axis, _label in ind.phenotype_key)
        }
        self.assertTrue(tubby)
        self.assertTrue(tubby.isdisjoint(plain))

    def test_confusable_eye_colors_warn(self) -> None:
        result = _cross("+ ; cn / + ; st / st", "+ / Y ; cn / + ; st / st")
        self.assertTrue(any("容易看混" in note for note in result.warnings))


class FusedAlleleTests(unittest.TestCase):
    def test_yw_is_yellow_and_white(self) -> None:
        self.assertEqual(split_fused_markers("yw"), ("y", "w"))
        self.assertEqual(split_fused_markers("y1w1118"), ("y", "w"))
        self.assertIsNone(split_fused_markers("wg"))
        self.assertIsNone(split_fused_markers("UAS-SytHA"))
        result = _cross("yw UAS-SytHA / yw", "yw / Y")
        males = [ind for ind in result.offspring if ind.sex == "male"]
        self.assertTrue(males)
        for ind in males:
            self.assertIn("黄体", ind.phenotype_text)
            self.assertIn("白眼", ind.phenotype_text)
            self.assertNotEqual(ind.phenotype_text, "无明显突变性状")
        displays = [
            hap.display()
            for ind in males
            for pair in ind.pairs
            for hap in pair.homologs()
        ]
        self.assertTrue(any("UAS-SytHA" in text for text in displays))

    def test_wg_is_not_split_into_white(self) -> None:
        result = _cross("wg / wg", "wg / Y")
        self.assertTrue(any("无翅" in ind.phenotype_text for ind in result.offspring))
        self.assertTrue(all("白眼" not in ind.phenotype_text for ind in result.offspring))


class ChromosomeViewTests(unittest.TestCase):
    def test_marginals_multiply_back_to_each_class(self) -> None:
        result = _cross(
            "w ; Sp / CyO ; TM3, Sb / MKRS",
            "w / Y ; + / CyO ; + / TM6B, Tb",
        )
        options = chromosome_options(result)
        self.assertEqual(len(options), 3)
        for column in options:
            self.assertEqual(sum(item.probability for item in column), 1)
        self.assertLess(sum(len(column) for column in options), len(result.offspring))
        for index, column in enumerate(options):
            choices = [0] * len(options)
            choices[index] = 0
            # Every combination of the first option of X with each later option exists.
        choices = [0, 0, 0]
        matched = match_offspring(result, options, choices)
        self.assertIsNotNone(matched)
        product = Fraction(1)
        for column, choice in zip(options, choices):
            product *= column[choice].probability
        self.assertEqual(matched.probability, product)


class SchemeTests(unittest.TestCase):
    def test_parallel_waves_then_a_cross(self) -> None:
        stocks = [
            Stock("a", "female", "w ; Sp / CyO"),
            Stock("b", "male", "w / Y ; + / +"),
            Stock("c", "female", "+ ; + / +"),
            Stock("d", "male", "+ / Y ; + / +"),
        ]
        first = run_scheme(
            stocks,
            [
                CrossStep("c1", ParentRef.stock("a"), ParentRef.stock("b")),
                CrossStep("c2", ParentRef.stock("c"), ParentRef.stock("d")),
            ],
        )
        self.assertEqual(first.waves, (("c1", "c2"),))
        female_index = next(i for i, ind in enumerate(first.steps[0].result.offspring) if ind.sex == "female")
        male_index = next(i for i, ind in enumerate(first.steps[1].result.offspring) if ind.sex == "male")
        second = run_scheme(
            stocks,
            [
                CrossStep("c1", ParentRef.stock("a"), ParentRef.stock("b")),
                CrossStep("c2", ParentRef.stock("c"), ParentRef.stock("d")),
                CrossStep(
                    "c3",
                    ParentRef.offspring("c1", female_index),
                    ParentRef.offspring("c2", male_index),
                ),
            ],
        )
        self.assertEqual(second.waves, (("c1", "c2"), ("c3",)))
        later = second.steps[-1]
        self.assertGreater(later.result.viable_count, 0)
        self.assertNotIn("♀", later.female_text)
        self.assertNotIn("♂", later.male_text)
        parent_lines = [
            line.text
            for line in format_cross_lines(later.result)
            if line.text.startswith("♀  ") or line.text.startswith("♂  ")
        ]
        self.assertEqual([line.count("♀") + line.count("♂") for line in parent_lines], [1, 1])

    def test_cycle_is_rejected(self) -> None:
        stocks = [
            Stock("a", "female", "+"),
            Stock("b", "male", "+ / Y"),
        ]
        with self.assertRaises(SchemeError):
            run_scheme(
                stocks,
                [
                    CrossStep("c1", ParentRef.offspring("c2", 0), ParentRef.stock("b")),
                    CrossStep("c2", ParentRef.stock("a"), ParentRef.offspring("c1", 0)),
                ],
            )


class ProjectTests(unittest.TestCase):
    def tearDown(self) -> None:
        set_language("en")

    def test_step_list_uses_the_current_cross_name(self) -> None:
        set_language("zh")
        stocks = [Stock("s1", "female", "w", "雌")]
        steps = [CrossStep("c1", ParentRef.stock("s1"), ParentRef.stock("s2"), label="常染色体")]
        label = _ref_label(ParentRef.offspring("c1", 0), stocks, steps)
        self.assertEqual(label, "常染色体 第1类")
        renamed = [CrossStep("c1", ParentRef.stock("s1"), ParentRef.stock("s2"), label="甲组")]
        self.assertEqual(_ref_label(ParentRef.offspring("c1", 0), stocks, renamed), "甲组 第1类")

    def test_project_file_roundtrip_keeps_fingerprints(self) -> None:
        fingerprint = ("female", ((("yw", "UAS-SytHA"), ("yw",)),))
        text = dumps_project({"picks": {"c1": [{"fingerprint": fingerprint}]}})
        loaded = loads_project(text)
        self.assertEqual(freeze(loaded["picks"]["c1"][0]["fingerprint"]), fingerprint)


class ExportTests(unittest.TestCase):
    def setUp(self) -> None:
        set_language("zh")

    def tearDown(self) -> None:
        set_language("en")

    def test_pptx_files_are_written(self) -> None:
        result = _cross("w ; CyO / +", "w / Y ; + / +")
        stocks = [
            Stock("a", "female", "w ; CyO / +"),
            Stock("b", "male", "w / Y ; + / +"),
        ]
        scheme = run_scheme(stocks, [CrossStep("c1", ParentRef.stock("a"), ParentRef.stock("b"))])
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            one = export_crossing_pptx(result, root / "one.pptx")
            plan = export_scheme_pptx(scheme, root / "plan.pptx")
            self.assertGreater(one.stat().st_size, 1000)
            self.assertGreater(plan.stat().st_size, 1000)

    def test_pptx_uses_the_checked_genotype_and_view_flags(self) -> None:
        result = _cross("B / +", "+ / Y")
        options = chromosome_options(result)
        chosen = match_offspring(result, options, [1])
        other = match_offspring(result, options, [0])
        assert chosen is not None and other is not None
        with tempfile.TemporaryDirectory() as folder:
            path = export_crossing_pptx(
                result,
                Path(folder) / "one.pptx",
                choices=[1],
                show_fraction=False,
                show_compact=True,
                show_phenotype=True,
            )
            text = _pptx_text(path)
            colored = export_crossing_pptx(
                result,
                Path(folder) / "fraction.pptx",
                choices=[1],
                show_fraction=True,
                show_compact=True,
                show_phenotype=False,
            )
            fraction_text = _pptx_text(colored)
        self.assertIn("染色体 X", text)
        self.assertIn("☑", text)
        for option in options[0]:
            mat, pat = option.pair.homologs()
            self.assertIn(f"{mat.display()}/{pat.display()}", text)
        self.assertNotIn("\n    ---", text)
        self.assertIn("表型:", text)
        self.assertIn("染色体 X", fraction_text)
        self.assertIn("\n    ---", fraction_text)
        self.assertNotIn("表型:", fraction_text)
        self.assertIsNotNone(chosen)
        self.assertIsNotNone(other)
        multi = _cross("w ; CyO / +", "w / Y ; + / +")
        with tempfile.TemporaryDirectory() as folder:
            compact_only = export_crossing_pptx(
                multi,
                Path(folder) / "compact.pptx",
                show_fraction=False,
                show_compact=True,
                show_phenotype=False,
            )
            compact_text = _pptx_text(compact_only)
        self.assertIn("染色体 X", compact_text)
        self.assertIn("染色体 2", compact_text)
        self.assertIn("CyO", compact_text)
        self.assertNotIn("\n    ---", compact_text)
        self.assertNotIn("表型:", compact_text)

    def test_english_report_and_pptx(self) -> None:
        result = _cross("yw", "yw / Y")
        set_language("en")
        self.assertIn("eye color", PHENOTYPE_EN["眼色"] and "eye color")
        for marker in MARKERS:
            for axis, label in marker.effects:
                self.assertIn(axis, PHENOTYPE_EN)
                self.assertIn(label, PHENOTYPE_EN)
        with tempfile.TemporaryDirectory() as folder:
            path = export_crossing_pptx(
                result,
                Path(folder) / "en.pptx",
                choices=[0],
                show_fraction=False,
                show_compact=True,
                show_phenotype=True,
            )
            text = _pptx_text(path)
        self.assertIn("Phenotype:", text)
        self.assertIn("yellow body", text)
        self.assertIn("white eyes", text)
        self.assertNotIn("黄体", text)


def _pptx_text(path: Path) -> str:
    from pptx import Presentation

    presentation = Presentation(str(path))
    chunks: list[str] = []
    for slide in presentation.slides:
        for shape in slide.shapes:
            if not shape.has_text_frame:
                continue
            for paragraph in shape.text_frame.paragraphs:
                chunks.append(paragraph.text)
    return "\n".join(chunks)


if __name__ == "__main__":
    unittest.main()
