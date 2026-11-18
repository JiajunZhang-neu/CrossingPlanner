"""Save and reopen a crossing plan as JSON."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .i18n import tr


PROJECT_VERSION = 1


def dumps_project(state: dict[str, Any]) -> str:
    payload = {"version": PROJECT_VERSION, **state}
    return json.dumps(payload, ensure_ascii=False, indent=2)


def loads_project(text: str) -> dict[str, Any]:
    data = json.loads(text)
    if not isinstance(data, dict) or data.get("version") != PROJECT_VERSION:
        raise ValueError(tr("bad_project"))
    return data


def write_project(path: str | Path, state: dict[str, Any]) -> Path:
    path = Path(path)
    path.write_text(dumps_project(state), encoding="utf-8")
    return path


def read_project(path: str | Path) -> dict[str, Any]:
    return loads_project(Path(path).read_text(encoding="utf-8"))


def freeze(value: Any) -> Any:
    """Turn JSON lists back into tuples so saved fingerprints match live ones."""
    if isinstance(value, list):
        return tuple(freeze(item) for item in value)
    return value
