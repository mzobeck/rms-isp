"""Pluggable annotator backends for phase 1 (variant -> consequence/HGVSp).

See docs/adr/0003-pluggable-annotators.md for the design rationale.
"""
from __future__ import annotations

from .base import Annotator, Variant, VariantAnnotation, get_annotator

__all__ = ["Annotator", "Variant", "VariantAnnotation", "get_annotator"]
