"""numpy だけで 4 近傍の連結成分ラベリング (scipy.ndimage.label の代わり)"""
import numpy as np


def label(mask):
    lab = np.zeros(mask.shape, int)
    cur = 0
    H, W = mask.shape
    for i0, j0 in zip(*np.nonzero(mask)):
        if lab[i0, j0]:
            continue
        cur += 1
        stack = [(i0, j0)]
        lab[i0, j0] = cur
        while stack:
            i, j = stack.pop()
            for a, b in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
                if 0 <= a < H and 0 <= b < W and mask[a, b] and not lab[a, b]:
                    lab[a, b] = cur
                    stack.append((a, b))
    return lab
