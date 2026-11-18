"""tkinter GUI for Fly Crossing Planner."""

from __future__ import annotations

import tkinter as tk
from fractions import Fraction
from tkinter import filedialog, messagebox, scrolledtext, ttk

from .cross import CrossResult, cross
from .format_frac import format_cross_result
from .genotype import ParseError, parse_genotype
from .markers import parse_extra_lethals
from .pptx_export import export_crossing_pptx


class CrossingPlannerApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Fly Crossing Planner")
        self.geometry("900x720")
        self.minsize(720, 560)

        self._last_result: CrossResult | None = None

        self._build_ui()
        self._load_examples()

    def _build_ui(self) -> None:
        pad = {"padx": 8, "pady": 4}

        main = ttk.Frame(self, padding=8)
        main.pack(fill=tk.BOTH, expand=True)

        # Parents
        parents = ttk.LabelFrame(main, text="Parents", padding=8)
        parents.pack(fill=tk.X, **pad)

        ttk.Label(parents, text="♀ Female genotype").grid(row=0, column=0, sticky=tk.W)
        self.female_text = scrolledtext.ScrolledText(parents, height=3, wrap=tk.WORD)
        self.female_text.grid(row=1, column=0, columnspan=2, sticky=tk.EW, pady=2)

        ttk.Label(parents, text="♂ Male genotype").grid(row=2, column=0, sticky=tk.W, pady=(8, 0))
        self.male_text = scrolledtext.ScrolledText(parents, height=3, wrap=tk.WORD)
        self.male_text.grid(row=3, column=0, columnspan=2, sticky=tk.EW, pady=2)

        parents.columnconfigure(0, weight=1)

        # Options
        opts = ttk.LabelFrame(main, text="Options", padding=8)
        opts.pack(fill=tk.X, **pad)

        ttk.Label(opts, text="Recombination RF (adjacent loci, 0–0.5):").grid(
            row=0, column=0, sticky=tk.W
        )
        self.rf_var = tk.StringVar(value="0")
        ttk.Entry(opts, textvariable=self.rf_var, width=10).grid(
            row=0, column=1, sticky=tk.W, padx=4
        )

        ttk.Label(opts, text="Extra lethal alleles (comma-separated):").grid(
            row=1, column=0, sticky=tk.W, pady=(6, 0)
        )
        self.lethal_var = tk.StringVar(value="")
        ttk.Entry(opts, textvariable=self.lethal_var).grid(
            row=1, column=1, sticky=tk.EW, padx=4, pady=(6, 0)
        )
        opts.columnconfigure(1, weight=1)

        self.show_fraction = tk.BooleanVar(value=True)
        self.show_compact = tk.BooleanVar(value=True)
        self.show_pheno = tk.BooleanVar(value=True)
        checks = ttk.Frame(opts)
        checks.grid(row=2, column=0, columnspan=2, sticky=tk.W, pady=(8, 0))
        ttk.Checkbutton(checks, text="Fraction view", variable=self.show_fraction).pack(
            side=tk.LEFT, padx=(0, 12)
        )
        ttk.Checkbutton(checks, text="Compact view", variable=self.show_compact).pack(
            side=tk.LEFT, padx=(0, 12)
        )
        ttk.Checkbutton(checks, text="Phenotype hints", variable=self.show_pheno).pack(
            side=tk.LEFT
        )

        # Buttons
        btns = ttk.Frame(main)
        btns.pack(fill=tk.X, **pad)
        ttk.Button(btns, text="计算 F1", command=self.on_calculate).pack(
            side=tk.LEFT, padx=(0, 8)
        )
        ttk.Button(btns, text="导出 PPT…", command=self.on_export).pack(side=tk.LEFT)
        ttk.Button(btns, text="Load example", command=self._load_examples).pack(
            side=tk.RIGHT
        )

        # Results
        results = ttk.LabelFrame(main, text="F1 results", padding=8)
        results.pack(fill=tk.BOTH, expand=True, **pad)
        self.result_text = scrolledtext.ScrolledText(
            results, wrap=tk.NONE, font=("Consolas", 10), state=tk.DISABLED
        )
        self.result_text.pack(fill=tk.BOTH, expand=True)
        xscroll = ttk.Scrollbar(results, orient=tk.HORIZONTAL, command=self.result_text.xview)
        self.result_text.configure(xscrollcommand=xscroll.set)

    def _load_examples(self) -> None:
        female = "w ; Sp / CyO ; TM3, Sb / MKRS"
        male = "w / Y ; + / CyO ; + / TM6B, Tb"
        self.female_text.delete("1.0", tk.END)
        self.female_text.insert(tk.END, female)
        self.male_text.delete("1.0", tk.END)
        self.male_text.insert(tk.END, male)

    def _parse_rf(self) -> Fraction:
        raw = self.rf_var.get().strip()
        if not raw:
            return Fraction(0)
        try:
            if "/" in raw:
                rf = Fraction(raw)
            else:
                rf = Fraction(str(float(raw))).limit_denominator(1000)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"Invalid RF value: {raw}") from exc
        if rf < 0 or rf > Fraction(1, 2):
            raise ValueError("RF must be between 0 and 0.5")
        return rf

    def on_calculate(self) -> None:
        try:
            female = parse_genotype(self.female_text.get("1.0", tk.END), "female")
            male = parse_genotype(self.male_text.get("1.0", tk.END), "male")
            rf = self._parse_rf()
            extra = parse_extra_lethals(self.lethal_var.get())
            result = cross(female, male, rf=rf, extra_lethals=extra)
        except (ParseError, ValueError) as exc:
            messagebox.showerror("Error", str(exc))
            return
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror("Unexpected error", str(exc))
            return

        self._last_result = result
        report = format_cross_result(
            result,
            show_fraction=self.show_fraction.get(),
            show_compact=self.show_compact.get(),
            show_phenotype=self.show_pheno.get(),
        )
        self.result_text.configure(state=tk.NORMAL)
        self.result_text.delete("1.0", tk.END)
        self.result_text.insert(tk.END, report)
        self.result_text.configure(state=tk.DISABLED)

    def on_export(self) -> None:
        if self._last_result is None:
            # Try calculate first
            self.on_calculate()
            if self._last_result is None:
                return

        path = filedialog.asksaveasfilename(
            title="Export PowerPoint",
            defaultextension=".pptx",
            filetypes=[("PowerPoint", "*.pptx"), ("All files", "*.*")],
            initialfile="crossing_scheme.pptx",
        )
        if not path:
            return
        try:
            export_crossing_pptx(self._last_result, path)
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror("Export failed", str(exc))
            return
        messagebox.showinfo("Exported", f"Saved:\n{path}")


def run_app() -> None:
    app = CrossingPlannerApp()
    app.mainloop()
