#!/usr/bin/env python3
"""Decode encoder PBMs with zxing-cpp and require an exact payload match.

Looks for pairs ``<tag>.pbm`` + ``<tag>.txt`` in the working directory.
Every version 1..40 must have both the short sample (``vNN``) and the
near-capacity sample (``vNNf``). Extra pairs (short strings) are checked too.
A missing required file is a failure. The encoder's own read-back is not used.
"""

import glob
import os
import sys

import numpy as np
import zxingcpp
from PIL import Image


def decode_pbm(path):
    toks = open(path).read().split()
    if not toks or toks[0] != "P1":
        raise SystemExit("not a P1 PBM: %s" % path)
    w, h = int(toks[1]), int(toks[2])
    bits = np.array([int(x) for x in toks[3 : 3 + w * h]], dtype=np.uint8).reshape(h, w)
    # PBM 1 is dark. Scale so the detector has a stable module grid.
    img = np.where(bits, 0, 255).astype(np.uint8)
    big = np.kron(img, np.ones((8, 8), dtype=np.uint8))
    res = zxingcpp.read_barcodes(Image.fromarray(big))
    return res[0].text if res else ""


def main():
    required = []
    for v in range(1, 41):
        required.append("v%02d" % v)
        required.append("v%02df" % v)

    ok = 0
    bad = []
    for tag in required + ["short", "short2"]:
        pbm = tag + ".pbm"
        txt = tag + ".txt"
        if not os.path.exists(pbm) or not os.path.exists(txt):
            bad.append(tag + " (missing)")
            print("  %-8s MISSING" % tag)
            continue
        expect = open(txt, encoding="utf-8").read()
        got = decode_pbm(pbm)
        if got == expect:
            ok += 1
            print("  %-8s OK %d chars" % (tag, len(got)))
        else:
            bad.append(tag)
            print("  %-8s FAIL expect %d got %d" % (tag, len(expect), len(got)))

    # Anything else the encoder dropped in this directory is checked as well.
    for pbm in sorted(glob.glob("*.pbm")):
        tag = os.path.splitext(os.path.basename(pbm))[0]
        if tag in required or tag in ("short", "short2"):
            continue
        txt = tag + ".txt"
        if not os.path.exists(txt):
            continue
        expect = open(txt, encoding="utf-8").read()
        got = decode_pbm(pbm)
        if got == expect:
            ok += 1
            print("  %-8s OK %d chars" % (tag, len(got)))
        else:
            bad.append(tag)
            print("  %-8s FAIL expect %d got %d" % (tag, len(expect), len(got)))

    print()
    print("passed %d, failed: %s" % (ok, bad if bad else "none"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
