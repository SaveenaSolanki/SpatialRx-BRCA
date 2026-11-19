"""Antibody (protein-feature) name handling.

The Visium CytAssist antibody panel stores features as ``SYMBOL-<n>``
(e.g. ``CD68-1``, ``PTPRC-2``) to disambiguate panels and isoforms, **except**
for symbols that legitimately contain a dash (e.g. ``HLA-DRA``) which appear
unsuffixed.

The original pipeline mapped antibodies to genes with ``name.split("-")[0]``,
which is wrong for ``HLA-DRA`` (yields ``HLA``) and silently dropped that pair
from every matched-marker analysis (ISSUE-001). These helpers replace that
heuristic and are covered by ``tests/test_proteins.py``.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence

__all__ = ["protein_to_symbol", "antibody_index", "isotype"]

# only a trailing "-<digits>" is a panel/isoform suffix
_SUFFIX = re.compile(r"-\d+$")

# isotype controls shipped with the 10x Immuno-Oncology panel
_ISOTYPES = frozenset(
    {"mouse_IgG1k", "mouse_IgG2a", "mouse_IgG2bk", "rat_IgG2a"}
)


def protein_to_symbol(name: str, known_symbols: Iterable[str] | None = None) -> str:
    """Map an antibody feature name to its gene symbol.

    A trailing ``-<digits>`` is treated as a panel/isoform suffix (``CD68-1`` ->
    ``CD68``, ``PTPRC-2`` -> ``PTPRC``), so dashed symbols such as ``HLA-DRA`` are
    preserved. Ambiguous cases — genes whose *own* symbol ends in ``-<digits>``
    (``NKX2-1``, ``KRTAP5-1``) — are resolved by passing the gene-symbol universe
    (``known_symbols``): a name present verbatim in it is never stripped.

    >>> protein_to_symbol("CD68-1")
    'CD68'
    >>> protein_to_symbol("HLA-DRA")
    'HLA-DRA'
    >>> protein_to_symbol("NKX2-1")                       # ambiguous by default
    'NKX2'
    >>> protein_to_symbol("NKX2-1", known_symbols={"NKX2-1"})
    'NKX2-1'
    """
    name = str(name)
    if known_symbols is not None and name in set(known_symbols):
        return name
    return _SUFFIX.sub("", name)


def antibody_index(
    names: Sequence[str], known_symbols: Iterable[str] | None = None
) -> dict[str, int]:
    """Map gene symbol -> column index, keeping the first isoform encountered.

    Isoforms (``PTPRC-1``, ``PTPRC-2``) and isotypes are handled deterministically:
    isoforms collapse to the first column, isotype controls are excluded.
    Pass ``known_symbols`` (e.g. the GEX ``var_names``) to protect gene symbols
    that themselves end in ``-<digits>``.
    """
    ks = set(known_symbols) if known_symbols is not None else None
    idx: dict[str, int] = {}
    for i, name in enumerate(names):
        name = str(name)
        if name in _ISOTYPES:
            continue
        idx.setdefault(protein_to_symbol(name, ks), i)
    return idx


def isotype(name: str) -> bool:
    """True if the antibody feature is an isotype control."""
    return str(name) in _ISOTYPES


def protein_symbols(
    names: Iterable[str], known_symbols: Iterable[str] | None = None
) -> list[str]:
    """Vectorised ``protein_to_symbol`` over a feature-name iterable."""
    ks = known_symbols
    return [protein_to_symbol(n, ks) for n in names]
