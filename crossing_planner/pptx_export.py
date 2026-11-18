"""Export crossing schemes to PowerPoint, one slide per step plus an overview."""

from __future__ import annotations

from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

from .cross import CrossResult
from .format_frac import ReportLine, format_cross_lines
from .i18n import tr
from .scheme import SchemeResult


_LINES_PER_SLIDE = 28


def export_crossing_pptx(
    result: CrossResult,
    path: str | Path,
    *,
    choices: list[int] | None = None,
    picks: list[str] | None = None,
    show_fraction: bool = True,
    show_compact: bool = True,
    show_phenotype: bool = True,
) -> Path:
    """Slides matching the result window: every chromosome option, plus the current checks."""
    path = Path(path)
    prs = _blank_presentation()
    lines = format_cross_lines(
        result,
        title=tr("ppt_single_title"),
        choices=choices,
        picks=picks,
        show_fraction=show_fraction,
        show_compact=show_compact,
        show_phenotype=show_phenotype,
    )
    _add_line_slides(prs, lines, tr("ppt_single_title"))
    prs.save(str(path))
    return path


def export_scheme_pptx(
    scheme: SchemeResult,
    path: str | Path,
    *,
    choices: dict[str, list[int]] | None = None,
    picks: dict[str, list[str]] | None = None,
    show_fraction: bool = True,
    show_compact: bool = True,
    show_phenotype: bool = True,
) -> Path:
    """Overview, then each step's chromosome options drawn with the current view flags."""
    path = Path(path)
    prs = _blank_presentation()
    _add_line_slides(prs, _overview_lines(scheme), tr("ppt_plan_title"))
    choices = choices or {}
    picks = picks or {}
    for step in scheme.steps:
        title = tr("wave_step", wave=step.wave + 1, name=step.step.display_name())
        lines = format_cross_lines(
            step.result,
            title=title,
            scope=step.step.id,
            choices=choices.get(step.step.id),
            picks=picks.get(step.step.id),
            show_fraction=show_fraction,
            show_compact=show_compact,
            show_phenotype=show_phenotype,
        )
        _add_line_slides(prs, lines, title)
    prs.save(str(path))
    return path


def _blank_presentation() -> Presentation:
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    return prs


def _overview_lines(scheme: SchemeResult) -> list[ReportLine]:
    lines = [ReportLine(tr("overview_title")), ReportLine(tr("overview_hint")), ReportLine("")]
    by_id = {step.step.id: step for step in scheme.steps}
    for wave_number, wave in enumerate(scheme.waves, 1):
        lines.append(ReportLine(tr("wave_heading", n=wave_number)))
        for step_id in wave:
            step = by_id[step_id]
            lines.append(ReportLine(f"  {step.step.display_name()}"))
            lines.append(ReportLine(f"    ♀  {step.female_text}"))
            lines.append(ReportLine(f"    ♂  {step.male_text}"))
            lines.append(
                ReportLine(tr("overview_parents", rf=step.step.rf, count=step.result.viable_count))
            )
        lines.append(ReportLine(""))
    return lines


def _add_line_slides(prs: Presentation, lines: list[ReportLine], title: str) -> None:
    if not lines:
        lines = [ReportLine(title)]
    chunks = [lines[index : index + _LINES_PER_SLIDE] for index in range(0, len(lines), _LINES_PER_SLIDE)]
    for chunk_index, chunk in enumerate(chunks):
        slide = prs.slides.add_slide(prs.slide_layouts[6])
        heading = title if len(chunks) == 1 else f"{title}  ({chunk_index + 1}/{len(chunks)})"
        _add_textbox(
            slide,
            Inches(0.45),
            Inches(0.22),
            Inches(12.4),
            Inches(0.4),
            [ReportLine(heading)],
            font_size=22,
            bold=True,
            align=PP_ALIGN.LEFT,
        )
        _add_textbox(
            slide,
            Inches(0.45),
            Inches(0.7),
            Inches(12.4),
            Inches(6.5),
            chunk,
            font_size=12,
            align=PP_ALIGN.LEFT,
            font_name="Consolas",
        )


def _add_textbox(
    slide,
    left,
    top,
    width,
    height,
    lines: list[ReportLine],
    *,
    font_size: int = 14,
    bold: bool = False,
    align=PP_ALIGN.LEFT,
    font_name: str = "Calibri",
) -> None:
    box = slide.shapes.add_textbox(left, top, width, height)
    frame = box.text_frame
    frame.word_wrap = True
    for index, line in enumerate(lines):
        paragraph = frame.paragraphs[0] if index == 0 else frame.add_paragraph()
        paragraph.alignment = align
        run = paragraph.add_run()
        run.text = line.text
        run.font.size = Pt(font_size)
        run.font.bold = bold
        run.font.name = font_name
        if line.color:
            run.font.color.rgb = RGBColor.from_string(line.color)
