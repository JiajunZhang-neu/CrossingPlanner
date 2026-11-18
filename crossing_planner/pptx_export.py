"""Export a single-page PowerPoint crossing scheme."""

from __future__ import annotations

from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

from .cross import CrossResult
from .format_frac import format_individual_compact, format_individual_fraction


def export_crossing_pptx(result: CrossResult, path: str | Path) -> Path:
    """
    Write a one-slide PPTX summarizing parents and viable F1 genotypes.
    """
    path = Path(path)
    prs = Presentation()
    # Widescreen 16:9
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    blank = prs.slide_layouts[6]  # blank
    slide = prs.slides.add_slide(blank)

    # Title
    _add_textbox(
        slide,
        Inches(0.5),
        Inches(0.25),
        Inches(12.3),
        Inches(0.5),
        "Crossing Scheme",
        font_size=28,
        bold=True,
        align=PP_ALIGN.CENTER,
    )

    # Parents
    parent_text = f"♀  {result.female_raw}\n×\n♂  {result.male_raw}"
    _add_textbox(
        slide,
        Inches(0.5),
        Inches(0.85),
        Inches(12.3),
        Inches(1.0),
        parent_text,
        font_size=16,
        align=PP_ALIGN.CENTER,
        font_name="Consolas",
    )

    # Arrow / F1 label
    _add_textbox(
        slide,
        Inches(0.5),
        Inches(1.85),
        Inches(12.3),
        Inches(0.35),
        "↓  F1  ↓",
        font_size=18,
        bold=True,
        align=PP_ALIGN.CENTER,
    )

    # F1 body
    body_lines: list[str] = []
    if not result.offspring:
        body_lines.append("No viable offspring.")
    else:
        for i, ind in enumerate(result.offspring, 1):
            body_lines.append(f"[{i}]  {format_individual_compact(ind)}")
            # Include stacked fraction (may be multi-line)
            frac = format_individual_fraction(ind)
            for line in frac.splitlines():
                body_lines.append("     " + line)
            if ind.phenotype_hints:
                body_lines.append(
                    "     phenotype: " + ", ".join(ind.phenotype_hints)
                )
            body_lines.append("")

    body = "\n".join(body_lines).rstrip() + "\n"
    _add_textbox(
        slide,
        Inches(0.6),
        Inches(2.25),
        Inches(12.1),
        Inches(4.5),
        body,
        font_size=11,
        align=PP_ALIGN.LEFT,
        font_name="Consolas",
    )

    # Footer
    footer_parts = [
        f"Viable classes: {result.viable_count}",
        f"RF={result.rf}",
    ]
    if result.lethal_discarded_probability > 0:
        footer_parts.append(
            "Homozygous lethal classes excluded and ratios re-normalized."
        )
    _add_textbox(
        slide,
        Inches(0.5),
        Inches(6.9),
        Inches(12.3),
        Inches(0.4),
        "  |  ".join(footer_parts),
        font_size=10,
        align=PP_ALIGN.CENTER,
        color=(0x55, 0x55, 0x55),
    )

    prs.save(str(path))
    return path


def _style_paragraph(p, *, font_size: int, bold: bool, align, font_name: str, color) -> None:
    p.alignment = align
    for run in p.runs:
        run.font.size = Pt(font_size)
        run.font.bold = bold
        run.font.name = font_name
        if color is not None:
            run.font.color.rgb = RGBColor(*color)


def _add_textbox(
    slide,
    left,
    top,
    width,
    height,
    text: str,
    *,
    font_size: int = 14,
    bold: bool = False,
    align=PP_ALIGN.LEFT,
    font_name: str = "Calibri",
    color: tuple[int, int, int] | None = None,
) -> None:
    box = slide.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame
    tf.word_wrap = True
    lines = text.split("\n")
    first = True
    for line in lines:
        if first:
            p = tf.paragraphs[0]
            first = False
        else:
            p = tf.add_paragraph()
        p.text = line
        _style_paragraph(
            p,
            font_size=font_size,
            bold=bold,
            align=align,
            font_name=font_name,
            color=color,
        )
