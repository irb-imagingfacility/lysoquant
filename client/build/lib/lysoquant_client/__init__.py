"""Segment images with a LysoQuant server from Python, without Fiji."""

from .client import (
    DEFAULT_FOLDERS,
    LysoQuantError,
    Segmentation,
    check_server,
    make_rgb,
    segment,
)

__all__ = ["DEFAULT_FOLDERS", "LysoQuantError", "Segmentation", "check_server", "make_rgb", "segment"]
