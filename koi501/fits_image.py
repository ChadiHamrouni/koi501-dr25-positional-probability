"""Reader for single-image FITS files.

The four follow-up frames are primary-HDU images with a linear CD matrix, which
is all a contrast curve needs, so the format is read directly rather than
adding astropy to the requirements.
"""

from __future__ import annotations

import numpy as np

from .models import FitsImageHeader

BLOCK = 2880
CARD = 80
DTYPES = {8: ">u1", 16: ">i2", 32: ">i4", 64: ">i8", -32: ">f4", -64: ">f8"}


def _value(raw: str):
    text = raw.split("/", 1)[0].strip() if not raw.strip().startswith("'") else raw.strip()
    if text.startswith("'"):
        return text[1:text.rfind("'")].rstrip()
    if text in ("T", "F"):
        return text == "T"
    try:
        return int(text)
    except ValueError:
        try:
            return float(text.replace("D", "E"))
        except ValueError:
            return text


def read_image(data: bytes) -> tuple[np.ndarray, dict]:
    """Return the primary image as float64 and its header cards as a dict."""
    cards, offset = {}, 0
    while True:
        block = data[offset:offset + BLOCK].decode("ascii")
        offset += BLOCK
        for i in range(0, BLOCK, CARD):
            card = block[i:i + CARD]
            key = card[:8].strip()
            if key == "END":
                header = FitsImageHeader.model_validate(cards)
                count = header.naxis1 * header.naxis2
                dtype = np.dtype(DTYPES[header.bitpix])
                pixels = np.frombuffer(data, dtype=dtype, count=count, offset=offset)
                image = pixels.astype(np.float64).reshape(header.naxis2, header.naxis1)
                return image * header.bscale + header.bzero, cards
            if card[8:10] == "= " and key not in cards:
                cards[key] = _value(card[10:])
