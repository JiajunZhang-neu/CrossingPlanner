"""tkinter GUI for Fly Crossing Planner."""

from __future__ import annotations

import json
import tkinter as tk
from dataclasses import dataclass
from fractions import Fraction
from tkinter import filedialog, messagebox, scrolledtext, ttk

from .cross import CrossResult, F1Individual
from .format_frac import (
    ReportLine,
    TextSpan,
    chromosome_options,
    format_cross_lines,
    format_individual_compact,
    match_offspring,
)
from .genotype import ParseError, parse_genotype
from .i18n import get_language, localize, set_language, tr
from .markers import parse_extra_lethals
from .pptx_export import export_scheme_pptx
from .project_file import freeze, read_project, write_project
from .scheme import (
    CrossStep,
    ParentRef,
    SchemeError,
    SchemeResult,
    StepResult,
    Stock,
    relabel_scheme,
    run_scheme,
)


@dataclass(frozen=True)
class SavedPick:
    fingerprint: tuple
    summary: str
    sex: str
    phenotype_text: str


class CrossingPlannerApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.geometry("1040x820")
        self.minsize(860, 640)

        self._stocks: list[Stock] = []
        self._steps: list[CrossStep] = []
        self._picks: dict[str, list[SavedPick]] = {}
        self._draft: dict[str, list[int]] = {}
        self._stock_serial = 0
        self._step_serial = 0
        self._last_scheme: SchemeResult | None = None
        self._i18n_widgets: list[tuple[tk.Widget, str]] = []

        self._build_ui()
        for var in (self.show_fraction, self.show_compact, self.show_pheno):
            var.trace_add("write", self._on_view_toggle)

    def _build_ui(self) -> None:
        bar = ttk.Frame(self)
        bar.pack(fill=tk.X, padx=8, pady=(6, 0))
        self._track(ttk.Label(bar, text=tr("lang_label")), "lang_label").pack(side=tk.LEFT)
        self.lang_var = tk.StringVar(value="en")
        ttk.Radiobutton(
            bar, text="中文", value="zh", variable=self.lang_var, command=self._on_language
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Radiobutton(
            bar, text="English", value="en", variable=self.lang_var, command=self._on_language
        ).pack(side=tk.LEFT, padx=(8, 0))

        scheme = ttk.Frame(self, padding=8)
        scheme.pack(fill=tk.BOTH, expand=True)
        self._build_scheme(scheme)
        self.title(tr("app_title"))

    def _track(self, widget: tk.Widget, key: str) -> tk.Widget:
        self._i18n_widgets.append((widget, key))
        return widget

    def _apply_language(self) -> None:
        self.title(tr("app_title"))
        for widget, key in self._i18n_widgets:
            widget.configure(text=tr(key))

    def _on_language(self) -> None:
        set_language(self.lang_var.get())
        self._apply_language()
        self._refresh_scheme_lists()
        if self._last_scheme is not None:
            self._render_scheme()

    def _build_scheme(self, main: ttk.Frame) -> None:
        pad = {"padx": 6, "pady": 4}

        stocks = self._track(ttk.LabelFrame(main, text=tr("stocks"), padding=8), "stocks")
        stocks.pack(fill=tk.X, **pad)
        self.stock_sex = tk.StringVar(value="female")
        ttk.Radiobutton(stocks, text="♀", value="female", variable=self.stock_sex).grid(row=0, column=0, sticky=tk.W)
        ttk.Radiobutton(stocks, text="♂", value="male", variable=self.stock_sex).grid(row=0, column=1, sticky=tk.W)
        self.stock_label = tk.StringVar()
        ttk.Entry(stocks, textvariable=self.stock_label, width=18).grid(row=0, column=2, padx=4)
        self._track(ttk.Label(stocks, text=tr("name")), "name").grid(row=0, column=3, sticky=tk.W)
        self.stock_genotype = tk.StringVar()
        ttk.Entry(stocks, textvariable=self.stock_genotype).grid(row=1, column=0, columnspan=4, sticky=tk.EW, pady=4)
        self._track(ttk.Button(stocks, text=tr("add_stock"), command=self.on_add_stock), "add_stock").grid(
            row=1, column=4, padx=4
        )
        self._track(ttk.Button(stocks, text=tr("rename"), command=self.on_rename_stock), "rename").grid(
            row=1, column=5, padx=4
        )
        self._track(
            ttk.Button(stocks, text=tr("delete_stock"), command=self.on_remove_stock), "delete_stock"
        ).grid(row=1, column=6)
        self.stock_list = tk.Listbox(stocks, height=4)
        self.stock_list.grid(row=2, column=0, columnspan=7, sticky=tk.EW, pady=(4, 0))
        self.stock_list.bind("<<ListboxSelect>>", self._on_stock_selected)
        stocks.columnconfigure(3, weight=1)

        self._track(ttk.Label(main, text=tr("scheme_hint")), "scheme_hint").pack(fill=tk.X, padx=8)

        steps = self._track(ttk.LabelFrame(main, text=tr("steps"), padding=8), "steps")
        steps.pack(fill=tk.X, **pad)
        self._track(ttk.Label(steps, text=tr("name")), "name").grid(row=0, column=0, sticky=tk.W)
        self.step_name = tk.StringVar()
        ttk.Entry(steps, textvariable=self.step_name, width=16).grid(row=0, column=1, sticky=tk.W, padx=4)
        self._track(ttk.Button(steps, text=tr("rename"), command=self.on_rename_step), "rename").grid(
            row=0, column=2, padx=4
        )
        ttk.Label(steps, text="♀").grid(row=1, column=0, sticky=tk.W, pady=(6, 0))
        self.step_female = ttk.Combobox(steps, state="readonly", width=28)
        self.step_female.grid(row=1, column=1, columnspan=2, sticky=tk.EW, padx=4, pady=(6, 0))
        ttk.Label(steps, text="♂").grid(row=1, column=3, sticky=tk.W, pady=(6, 0))
        self.step_male = ttk.Combobox(steps, state="readonly", width=28)
        self.step_male.grid(row=1, column=4, sticky=tk.EW, padx=4, pady=(6, 0))
        ttk.Label(steps, text="RF").grid(row=1, column=5, sticky=tk.W, pady=(6, 0))
        self.step_rf = tk.StringVar(value="0")
        ttk.Entry(steps, textvariable=self.step_rf, width=8).grid(row=1, column=6, padx=4, pady=(6, 0))
        self._track(ttk.Button(steps, text=tr("add_step"), command=self.on_add_step), "add_step").grid(
            row=1, column=7, padx=4, pady=(6, 0)
        )
        self.step_list = tk.Listbox(steps, height=4)
        self.step_list.grid(row=2, column=0, columnspan=7, sticky=tk.EW, pady=(4, 0))
        self.step_list.bind("<<ListboxSelect>>", self._on_step_selected)
        self._track(
            ttk.Button(steps, text=tr("delete_step"), command=self.on_remove_step), "delete_step"
        ).grid(row=2, column=7, sticky=tk.N)
        steps.columnconfigure(4, weight=1)

        self.show_fraction = tk.BooleanVar(value=True)
        self.show_compact = tk.BooleanVar(value=True)
        self.show_pheno = tk.BooleanVar(value=True)

        lethal = ttk.Frame(main)
        lethal.pack(fill=tk.X, **pad)
        self._track(ttk.Label(lethal, text=tr("extra_lethals")), "extra_lethals").pack(side=tk.LEFT)
        self.lethal_var = tk.StringVar(value="")
        ttk.Entry(lethal, textvariable=self.lethal_var).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(8, 0))

        actions = ttk.Frame(main)
        actions.pack(fill=tk.X, **pad)
        self._track(
            ttk.Button(actions, text=tr("calc_scheme"), command=self.on_scheme_calculate), "calc_scheme"
        ).pack(side=tk.LEFT, padx=(0, 8))
        self._track(
            ttk.Button(actions, text=tr("export_scheme"), command=self.on_scheme_export), "export_scheme"
        ).pack(side=tk.LEFT, padx=(0, 8))
        self._track(
            ttk.Button(actions, text=tr("load_parallel"), command=self.on_scheme_example), "load_parallel"
        ).pack(side=tk.LEFT, padx=(0, 8))
        self._track(
            ttk.Button(actions, text=tr("save_project"), command=self.on_save_project), "save_project"
        ).pack(side=tk.LEFT, padx=(0, 8))
        self._track(
            ttk.Button(actions, text=tr("open_project"), command=self.on_open_project), "open_project"
        ).pack(side=tk.LEFT)
        self._add_view_checks(actions).pack(side=tk.LEFT, padx=(16, 0))

        results = self._track(ttk.LabelFrame(main, text=tr("scheme_results"), padding=8), "scheme_results")
        results.pack(fill=tk.BOTH, expand=True, **pad)
        self.scheme_text = self._result_widget(results, cursor="arrow")
        self.scheme_text.pack(fill=tk.BOTH, expand=True)

    def _add_view_checks(self, parent: tk.Widget) -> ttk.Frame:
        checks = ttk.Frame(parent)
        self._track(
            ttk.Checkbutton(checks, text=tr("view_fraction"), variable=self.show_fraction), "view_fraction"
        ).pack(side=tk.LEFT, padx=(0, 12))
        self._track(
            ttk.Checkbutton(checks, text=tr("view_compact"), variable=self.show_compact), "view_compact"
        ).pack(side=tk.LEFT, padx=(0, 12))
        self._track(
            ttk.Checkbutton(checks, text=tr("view_phenotype"), variable=self.show_pheno), "view_phenotype"
        ).pack(side=tk.LEFT)
        return checks

    def _result_widget(self, parent: tk.Widget, cursor: str = "arrow") -> scrolledtext.ScrolledText:
        widget = scrolledtext.ScrolledText(parent, wrap=tk.NONE, font=("Consolas", 10), cursor=cursor)
        widget.bind("<Button-1>", self._on_text_click)
        widget.bind("<Key>", self._swallow_key)
        widget.bind("<Motion>", self._on_hover)
        widget.bind("<Leave>", self._on_hover_leave)
        widget._hover_tag = ""  # type: ignore[attr-defined]
        return widget

    def _parse_rf(self, raw: str) -> Fraction:
        raw = raw.strip()
        if not raw:
            return Fraction(0)
        try:
            rf = Fraction(raw) if "/" in raw else Fraction(str(float(raw))).limit_denominator(1000)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(tr("invalid_rf", value=raw)) from exc
        if rf < 0 or rf > Fraction(1, 2):
            raise ValueError(tr("rf_range"))
        return rf

    def _view_flags(self) -> dict[str, bool]:
        return {
            "show_fraction": self.show_fraction.get(),
            "show_compact": self.show_compact.get(),
            "show_phenotype": self.show_pheno.get(),
        }

    def on_add_stock(self) -> None:
        text = self.stock_genotype.get().strip()
        sex = self.stock_sex.get()
        try:
            parse_genotype(text, sex)
        except ParseError as exc:
            messagebox.showerror(tr("stock_parse_error"), str(exc))
            return
        self._stock_serial += 1
        stock = Stock(
            id=f"s{self._stock_serial}",
            sex=sex,
            genotype=text,
            label=self.stock_label.get().strip(),
        )
        self._stocks.append(stock)
        self.stock_genotype.set("")
        self.stock_label.set("")
        self._refresh_scheme_lists()

    def on_rename_stock(self) -> None:
        selection = self.stock_list.curselection()
        if not selection:
            messagebox.showerror(tr("no_stock_selected"), tr("select_stock_first"))
            return
        index = selection[0]
        current = self._stocks[index]
        self._stocks[index] = Stock(
            id=current.id,
            sex=current.sex,
            genotype=current.genotype,
            label=self.stock_label.get().strip(),
        )
        self._refresh_scheme_lists()

    def on_remove_stock(self) -> None:
        selection = self.stock_list.curselection()
        if not selection:
            return
        removed = self._stocks.pop(selection[0])
        self._steps = [
            step
            for step in self._steps
            if not _ref_uses_stock(step.female, removed.id) and not _ref_uses_stock(step.male, removed.id)
        ]
        self._refresh_scheme_lists()

    def on_add_step(self) -> None:
        female = self._selected_ref(self.step_female, "female")
        male = self._selected_ref(self.step_male, "male")
        if female is None or male is None:
            messagebox.showerror(tr("cannot_add"), tr("pick_both_parents"))
            return
        try:
            rf = self._parse_rf(self.step_rf.get())
        except ValueError as exc:
            messagebox.showerror("RF", str(exc))
            return
        self._step_serial += 1
        self._steps.append(
            CrossStep(
                id=f"c{self._step_serial}",
                female=female,
                male=male,
                rf=rf,
                label=self.step_name.get().strip(),
            )
        )
        self.step_name.set("")
        self._refresh_scheme_lists()

    def on_rename_step(self) -> None:
        selection = self.step_list.curselection()
        if not selection:
            messagebox.showerror(tr("no_step_selected"), tr("select_step_first"))
            return
        index = selection[0]
        current = self._steps[index]
        self._steps[index] = CrossStep(
            id=current.id,
            female=current.female,
            male=current.male,
            rf=current.rf,
            label=self.step_name.get().strip(),
        )
        self._refresh_scheme_lists()
        if self._last_scheme is not None:
            self._render_scheme()

    def on_remove_step(self) -> None:
        selection = self.step_list.curselection()
        if not selection:
            return
        removed = self._steps.pop(selection[0])
        self._picks.pop(removed.id, None)
        self._draft.pop(removed.id, None)
        self._steps = [
            step
            for step in self._steps
            if step.female.step_id != removed.id and step.male.step_id != removed.id
        ]
        self._refresh_scheme_lists()

    def on_scheme_example(self) -> None:
        self._stocks = [
            Stock("s1", "female", "w ; Sp / CyO ; TM3, Sb / MKRS", tr("ex_double_female")),
            Stock("s2", "male", "w / Y ; + / CyO ; + / TM6B, Tb", tr("ex_tester_male")),
            Stock("s3", "female", "B / +", tr("ex_bar_female")),
            Stock("s4", "male", "+ / Y", tr("ex_wild_male")),
        ]
        self._steps = [
            CrossStep("c1", ParentRef.stock("s1"), ParentRef.stock("s2"), label=tr("ex_autosome")),
            CrossStep("c2", ParentRef.stock("s3"), ParentRef.stock("s4"), label=tr("ex_x")),
        ]
        self._picks = {}
        self._draft = {}
        self._stock_serial = 4
        self._step_serial = 2
        self._last_scheme = None
        self._refresh_scheme_lists()

    def on_scheme_calculate(self) -> None:
        if not self._steps:
            messagebox.showerror(tr("no_steps"), tr("add_stock_and_step"))
            return
        try:
            scheme = run_scheme(
                self._stocks,
                self._steps,
                extra_lethals=parse_extra_lethals(self.lethal_var.get()),
            )
        except (SchemeError, ParseError, ValueError) as exc:
            messagebox.showerror(tr("scheme_failed"), str(exc))
            return
        self._last_scheme = scheme
        for step in scheme.steps:
            self._draft[step.step.id] = _fit_choices(step.result, self._draft.get(step.step.id))
        self._render_scheme()
        self._refresh_scheme_lists()

    def on_scheme_export(self) -> None:
        if self._last_scheme is None:
            self.on_scheme_calculate()
            if self._last_scheme is None:
                return
        path = filedialog.asksaveasfilename(
            title=tr("export_scheme_title"),
            defaultextension=".pptx",
            filetypes=[("PowerPoint", "*.pptx"), ("All files", "*.*")],
            initialfile="crossing_plan.pptx",
        )
        if not path:
            return
        scheme = relabel_scheme(self._last_scheme, self._steps)
        try:
            export_scheme_pptx(
                scheme,
                path,
                choices={step.step.id: list(self._draft.get(step.step.id, [])) for step in scheme.steps},
                picks={step.step.id: self._pick_labels(step) for step in scheme.steps},
                **self._view_flags(),
            )
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror(tr("export_failed"), str(exc))
            return
        messagebox.showinfo(tr("exported"), tr("saved_path", path=path))

    def _on_stock_selected(self, _event: object) -> None:
        selection = self.stock_list.curselection()
        if not selection:
            return
        self.stock_label.set(self._stocks[selection[0]].label)

    def _on_step_selected(self, _event: object) -> None:
        selection = self.step_list.curselection()
        if not selection:
            return
        self.step_name.set(self._steps[selection[0]].label)

    def _on_view_toggle(self, *_args: object) -> None:
        if self._last_scheme is not None:
            self._render_scheme()

    def _render_scheme(self) -> None:
        if self._last_scheme is None:
            return
        lines: list[ReportLine] = [ReportLine(tr("waves_independent")), ReportLine("")]
        for wave_number, wave in enumerate(self._last_scheme.waves, 1):
            names = ", ".join(self._step_title(step_id) for step_id in wave)
            lines.append(ReportLine(tr("wave_line", n=wave_number, names=names)))
        lines.append(ReportLine(""))
        flags = self._view_flags()
        for step in self._last_scheme.steps:
            lines.append(ReportLine("=" * 48))
            picks = self._pick_labels(step)
            live = next((item for item in self._steps if item.id == step.step.id), step.step)
            lines.extend(
                format_cross_lines(
                    step.result,
                    title=tr(
                        "wave_step_parents",
                        wave=step.wave + 1,
                        name=live.display_name(),
                        female=step.female_text,
                        male=step.male_text,
                    ),
                    scope=step.step.id,
                    choices=self._draft.get(step.step.id),
                    picks=picks,
                    interactive=True,
                    allow_picks=True,
                    **flags,
                )
            )
            lines.append(ReportLine(""))
        self._show_lines(self.scheme_text, lines)

    def _pick_labels(self, step: StepResult) -> list[str]:
        labels: list[str] = []
        for pick in self._picks.get(step.step.id, []):
            suffix = "" if self._resolve_pick(step, pick) is not None else tr("not_in_result")
            labels.append(pick.summary + suffix)
        return labels

    def _choose(self, scope: str, chrom: int, option: int) -> None:
        result = self._result_for_scope(scope)
        if result is None:
            return
        choices = _fit_choices(result, self._draft.get(scope))
        if chrom >= len(choices):
            return
        choices[chrom] = option
        self._draft[scope] = choices
        self._render_scheme()

    def _add_pick(self, scope: str) -> None:
        step = self._step_result(scope)
        if step is None:
            return
        choices = _fit_choices(step.result, self._draft.get(scope))
        options = chromosome_options(step.result)
        chosen = match_offspring(step.result, options, choices)
        if chosen is None:
            messagebox.showerror(tr("cannot_pick"), tr("combo_missing_viable"))
            return
        fingerprint = _fingerprint(chosen)
        current = self._picks.setdefault(scope, [])
        if any(pick.fingerprint == fingerprint for pick in current):
            return
        current.append(
            SavedPick(
                fingerprint=fingerprint,
                summary=format_individual_compact(chosen),
                sex=chosen.sex,
                phenotype_text=chosen.phenotype_text,
            )
        )
        self._refresh_scheme_lists()
        self._render_scheme()

    def _remove_pick(self, scope: str, index: int) -> None:
        picks = self._picks.get(scope, [])
        if index < 0 or index >= len(picks):
            return
        del picks[index]
        self._refresh_scheme_lists()
        self._render_scheme()

    def _on_text_click(self, event: tk.Event) -> str | None:
        widget = event.widget
        index = widget.index(f"@{event.x},{event.y}")
        for tag in widget.tag_names(index):
            if tag.startswith("opt:") or tag.startswith("box:"):
                _kind, scope, chrom, option = tag.split(":")
                self._choose(scope, int(chrom), int(option))
                return "break"
            if tag.startswith("add:") or tag.startswith("plus:"):
                self._add_pick(tag.split(":", 1)[1])
                return "break"
            if tag.startswith("del:"):
                _del, scope, number = tag.split(":")
                self._remove_pick(scope, int(number))
                return "break"
        return None

    def _swallow_key(self, event: tk.Event) -> str | None:
        if event.state & 0x4:
            return None
        if event.keysym in {
            "Left",
            "Right",
            "Up",
            "Down",
            "Prior",
            "Next",
            "Home",
            "End",
            "Shift_L",
            "Shift_R",
        }:
            return None
        return "break"

    def _result_for_scope(self, scope: str) -> CrossResult | None:
        step = self._step_result(scope)
        return None if step is None else step.result

    def _step_result(self, step_id: str) -> StepResult | None:
        if self._last_scheme is None:
            return None
        for step in self._last_scheme.steps:
            if step.step.id == step_id:
                return step
        return None

    def _resolve_pick(self, step: StepResult, pick: SavedPick) -> int | None:
        for index, ind in enumerate(step.result.offspring):
            if _fingerprint(ind) == pick.fingerprint:
                return index
        return None

    def _selected_ref(self, combo: ttk.Combobox, sex: str) -> ParentRef | None:
        label = combo.get()
        for text, ref in self._source_choices(sex):
            if text == label:
                return ref
        return None

    def _source_choices(self, sex: str) -> list[tuple[str, ParentRef]]:
        choices: list[tuple[str, ParentRef]] = []
        for stock in self._stocks:
            if stock.sex != sex:
                continue
            mark = "♀" if sex == "female" else "♂"
            choices.append(
                (
                    tr("stock_source", name=stock.display_name(), mark=mark, genotype=stock.genotype),
                    ParentRef.stock(stock.id),
                )
            )
        if self._last_scheme is None:
            return choices
        for step in self._last_scheme.steps:
            for number, pick in enumerate(self._picks.get(step.step.id, []), 1):
                if pick.sex != sex:
                    continue
                index = self._resolve_pick(step, pick)
                if index is None:
                    continue
                mark = "♀" if sex == "female" else "♂"
                choices.append(
                    (
                        tr(
                            "pick_source",
                            name=self._step_title(step.step.id),
                            n=number,
                            mark=mark,
                            phenotype=localize(pick.phenotype_text),
                        ),
                        ParentRef.offspring(step.step.id, index),
                    )
                )
        return choices

    def _refresh_scheme_lists(self) -> None:
        self.stock_list.delete(0, tk.END)
        for stock in self._stocks:
            mark = "♀" if stock.sex == "female" else "♂"
            self.stock_list.insert(tk.END, f"{stock.display_name()}  {mark}  {stock.genotype}")

        self.step_list.delete(0, tk.END)
        for step in self._steps:
            count = len(self._picks.get(step.id, []))
            suffix = tr("picks_suffix", count=count) if count else ""
            self.step_list.insert(
                tk.END,
                f"{step.display_name()}  ♀ {_ref_label(step.female, self._stocks, self._steps)}  ×  "
                f"♂ {_ref_label(step.male, self._stocks, self._steps)}  RF={step.rf}{suffix}",
            )

        female_labels = [text for text, _ref in self._source_choices("female")]
        male_labels = [text for text, _ref in self._source_choices("male")]
        self.step_female["values"] = female_labels
        self.step_male["values"] = male_labels
        if self.step_female.get() not in female_labels:
            self.step_female.set(female_labels[0] if female_labels else "")
        if self.step_male.get() not in male_labels:
            self.step_male.set(male_labels[0] if male_labels else "")

    def _show_lines(self, widget: scrolledtext.ScrolledText, lines: list[ReportLine]) -> None:
        first, _last = widget.yview()
        self._clear_hover(widget)
        widget.configure(state=tk.NORMAL)
        widget.delete("1.0", tk.END)
        for line in lines:
            pieces = line.spans or (TextSpan(line.text, line.color, line.tag),)
            for piece_index, piece in enumerate(pieces):
                tags: list[str] = []
                if piece.color:
                    color_tag = "c" + piece.color
                    if color_tag not in widget.tag_names():
                        widget.tag_configure(color_tag, foreground="#" + piece.color)
                    tags.append(color_tag)
                if piece.tag:
                    tags.append(piece.tag)
                ending = "\n" if piece_index == len(pieces) - 1 else ""
                widget.insert(tk.END, piece.text + ending, tuple(tags))
        widget.yview_moveto(first)

    def _on_hover(self, event: tk.Event) -> None:
        widget = event.widget
        index = widget.index(f"@{event.x},{event.y}")
        target = ""
        for tag in widget.tag_names(index):
            if tag.startswith("box:") or tag.startswith("plus:"):
                target = tag
                break
        previous = getattr(widget, "_hover_tag", "")
        if target == previous:
            return
        self._clear_hover(widget)
        if target:
            widget.tag_configure(target, background="#D6E6F5")
            widget.tag_raise(target)
            widget._hover_tag = target

    def _on_hover_leave(self, event: tk.Event) -> None:
        self._clear_hover(event.widget)

    def _clear_hover(self, widget: tk.Text) -> None:
        previous = getattr(widget, "_hover_tag", "")
        if not previous:
            return
        try:
            widget.tag_configure(previous, background="")
        except tk.TclError:
            pass
        widget._hover_tag = ""

    def on_save_project(self) -> None:
        path = filedialog.asksaveasfilename(
            title=tr("save_title"),
            defaultextension=".json",
            filetypes=[(tr("project_filetype"), "*.json"), (tr("all_files"), "*.*")],
            initialfile="crossing_project.json",
        )
        if not path:
            return
        try:
            write_project(path, self._project_state())
        except OSError as exc:
            messagebox.showerror(tr("save_failed"), str(exc))
            return
        messagebox.showinfo(tr("saved"), path)

    def on_open_project(self) -> None:
        if self._stocks or self._steps:
            if not messagebox.askyesno(tr("open_title"), tr("open_replace")):
                return
        path = filedialog.askopenfilename(
            title=tr("open_title"),
            filetypes=[(tr("project_filetype"), "*.json"), (tr("all_files"), "*.*")],
        )
        if not path:
            return
        try:
            data = read_project(path)
            self._apply_project(data)
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
            messagebox.showerror(tr("cannot_open"), str(exc))

    def _project_state(self) -> dict:
        return {
            "language": self.lang_var.get(),
            "extra_lethals": self.lethal_var.get(),
            "view": {
                "fraction": self.show_fraction.get(),
                "compact": self.show_compact.get(),
                "phenotype": self.show_pheno.get(),
            },
            "serials": {"stock": self._stock_serial, "step": self._step_serial},
            "stocks": [
                {"id": stock.id, "sex": stock.sex, "genotype": stock.genotype, "label": stock.label}
                for stock in self._stocks
            ],
            "steps": [
                {
                    "id": step.id,
                    "label": step.label,
                    "rf": str(step.rf),
                    "female": _ref_state(step.female),
                    "male": _ref_state(step.male),
                }
                for step in self._steps
            ],
            "picks": {
                step_id: [
                    {
                        "fingerprint": pick.fingerprint,
                        "summary": pick.summary,
                        "sex": pick.sex,
                        "phenotype_text": pick.phenotype_text,
                    }
                    for pick in picks
                ]
                for step_id, picks in self._picks.items()
            },
            "draft": self._draft,
        }

    def _apply_project(self, data: dict) -> None:
        set_language(str(data.get("language", "en")))
        self.lang_var.set("en" if get_language() == "en" else "zh")
        self._apply_language()
        single = data.get("single", {})
        view = data.get("view", {})
        self.lethal_var.set(str(data.get("extra_lethals", single.get("extra_lethals", ""))))
        self.show_fraction.set(bool(view.get("fraction", True)))
        self.show_compact.set(bool(view.get("compact", True)))
        self.show_pheno.set(bool(view.get("phenotype", True)))

        self._stocks = [
            Stock(
                id=item["id"],
                sex=item["sex"],
                genotype=item["genotype"],
                label=item.get("label", ""),
            )
            for item in data.get("stocks", [])
        ]
        self._steps = [
            CrossStep(
                id=item["id"],
                female=_ref_from_state(item["female"]),
                male=_ref_from_state(item["male"]),
                rf=Fraction(item.get("rf", "0")),
                label=item.get("label", ""),
            )
            for item in data.get("steps", [])
        ]
        self._picks = {}
        for step_id, picks in data.get("picks", {}).items():
            self._picks[step_id] = [
                SavedPick(
                    fingerprint=freeze(item["fingerprint"]),
                    summary=item.get("summary", ""),
                    sex=item.get("sex", ""),
                    phenotype_text=item.get("phenotype_text", ""),
                )
                for item in picks
            ]
        self._draft = {
            scope: [int(index) for index in choices] for scope, choices in data.get("draft", {}).items()
        }
        serials = data.get("serials", {})
        self._stock_serial = int(serials.get("stock", len(self._stocks)))
        self._step_serial = int(serials.get("step", len(self._steps)))
        self._last_scheme = None
        self.scheme_text.delete("1.0", tk.END)
        self._refresh_scheme_lists()
        if self._steps:
            self.on_scheme_calculate()

    def _step_title(self, step_id: str) -> str:
        for step in self._steps:
            if step.id == step_id:
                return step.display_name()
        return step_id


def _fit_choices(result: CrossResult, previous: list[int] | None) -> list[int]:
    width = len(result.offspring[0].pairs) if result.offspring else 0
    choices = list(previous or [])
    if len(choices) < width:
        choices.extend([0] * (width - len(choices)))
    return choices[:width]


def _fingerprint(ind: F1Individual) -> tuple:
    pairs = []
    for pair in ind.pairs:
        maternal, paternal = pair.homologs()
        pairs.append((maternal.loci, paternal.loci))
    return (ind.sex, tuple(pairs))


def _ref_uses_stock(ref: ParentRef, stock_id: str) -> bool:
    return ref.kind == "stock" and ref.stock_id == stock_id


def _ref_label(ref: ParentRef, stocks: list[Stock], steps: list[CrossStep]) -> str:
    if ref.kind == "stock":
        for stock in stocks:
            if stock.id == ref.stock_id:
                return stock.display_name()
        return ref.stock_id
    name = ref.step_id
    for step in steps:
        if step.id == ref.step_id:
            name = step.display_name()
            break
    return tr("offspring_class", name=name, n=ref.offspring_index + 1)


def _ref_state(ref: ParentRef) -> dict[str, str | int]:
    return {
        "kind": ref.kind,
        "stock_id": ref.stock_id,
        "step_id": ref.step_id,
        "offspring_index": ref.offspring_index,
    }


def _ref_from_state(data: dict) -> ParentRef:
    return ParentRef(
        kind=str(data.get("kind", "stock")),
        stock_id=str(data.get("stock_id", "")),
        step_id=str(data.get("step_id", "")),
        offspring_index=int(data.get("offspring_index", 0)),
    )


def run_app() -> None:
    app = CrossingPlannerApp()
    app.mainloop()
