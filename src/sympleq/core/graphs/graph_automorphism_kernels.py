from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, fastmath=True)
def _consistent_numba(
    S_mod: np.ndarray,
    phi: np.ndarray,
    mapped_stack: np.ndarray,
    mapped_len: int,
    i: int,
    y: int,
) -> bool:
    """Check edge-colour consistency against all already-mapped indices."""
    for t in range(mapped_len):
        j = mapped_stack[t]
        yj = phi[j]
        if S_mod[i, j] != S_mod[y, yj]:
            return False
        if S_mod[j, i] != S_mod[yj, y]:
            return False
    return True


def _build_bitrows_binary(S_mod: np.ndarray) -> tuple[np.ndarray, int]:
    """For p=2 only. Pack each row's 0/1 into chunks of 64 bits."""
    n = S_mod.shape[0]
    C = (n + 63) // 64
    bits = np.zeros((n, C), dtype=np.uint64)
    for i in range(n):
        row = bits[i]
        for j in range(n):
            if S_mod[i, j] & 1:
                row[j >> 6] |= (np.uint64(1) << np.uint64(j & 63))
    return bits, C


@njit(cache=True, fastmath=True)
def _consistent_bitset(
    bits: np.ndarray,
    phi: np.ndarray,
    mapped_stack: np.ndarray,
    mapped_len: int,
    i: int,
    y: int,
) -> bool:
    """Same logic as _consistent_numba but reading single bits from packed rows."""
    for t in range(mapped_len):
        j = mapped_stack[t]
        yj = phi[j]
        bi = (bits[i, j >> 6] >> (j & 63)) & 1
        by = (bits[y, yj >> 6] >> (yj & 63)) & 1
        if bi != by:
            return False
        bji = (bits[j, i >> 6] >> (i & 63)) & 1
        byy = (bits[yj, y >> 6] >> (y & 63)) & 1
        if bji != byy:
            return False
    return True


class _ConsistencyChecker:
    """Reusable consistency kernel; chooses bitset or direct variant once."""

    def __init__(self, S_mod: np.ndarray, p2_bitset: bool):
        self.S_mod = S_mod
        if p2_bitset:
            self.bits, _ = _build_bitrows_binary(S_mod)
            self._fn = self._bitset
        else:
            self._fn = self._direct

    def __call__(
        self,
        phi: np.ndarray,
        mapped_stack: np.ndarray,
        mapped_len: int,
        i: int,
        y: int,
    ) -> bool:
        return self._fn(phi, mapped_stack, int(mapped_len), int(i), int(y))

    def _bitset(self, phi: np.ndarray, mapped_stack: np.ndarray, mapped_len: int, i: int, y: int) -> bool:
        return _consistent_bitset(self.bits, phi, mapped_stack, mapped_len, i, y)

    def _direct(self, phi: np.ndarray, mapped_stack: np.ndarray, mapped_len: int, i: int, y: int) -> bool:
        return _consistent_numba(self.S_mod, phi, mapped_stack, mapped_len, i, y)
