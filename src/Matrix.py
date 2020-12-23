from __future__ import annotations
from typing import List

from galois_field import ElementInGFp, ElementInGFpn, GFp


class Matrix:
    """A Matrix in GF(p)
    """

    def __init__(self, mat: List[List[ElementInGFp]]):
        self.__values = mat
        self.__m = len(mat)
        self.__n = len(mat[0])

    def __getitem__(self, i: int) -> List[ElementInGFp]:
        return self.__values[i]

    def __setitem__(self, i: int, value: List[ElementInGFp]) -> None:
        self.__values[i] = value

    @property
    def m(self) -> int:
        return self.__m

    @property
    def n(self) -> int:
        return self.__n

    @property
    def rank(self) -> int:
        return _seek_rank_over_Field(self)

    @staticmethod
    def from_list(li: List[ElementInGFpn]) -> Matrix:
        print("li", li)
        mat = _conv_GFpn_list_to_GFp_matrix(li)
        return Matrix(mat)


def _conv_GFpn_list_to_GFp_matrix(li: List[ElementInGFpn]) -> List[List[ElementInGFp]]:
    Fp = GFp(li[0].p)
    result = []
    max_len = max(map(lambda gi: len(gi.coeffs), li))
    for gi in li:
        result.append([Fp.elm(0)]*(max_len-len(gi.coeffs)) +
                      list(map(lambda g: Fp.elm(int(g)), gi.coeffs)))

    return result


def _seek_rank_over_Field(mat: Matrix) -> int:
    """Find the rank of a matrix over GF(p).

    Args:
        mat (Matrix): M_ij in GF(p).

    Returns:
        int: The rank of a matrix over GF(p).
    """
    # must be m < n
    m = mat.m
    n = mat.n

    mi, ni = 0, 0

    # Convert the matrix to row echelon form.
    while mi < m and ni < n:
        # Get the line number after the mi line where the ni column is not 0.
        mj = _find_none_zero_mi(mat, ni, mi)
        if mj == -1:
            ni += 1
            continue
        # Exchange mi and mj lines
        if mi != mj:
            mat = _swap_rows(mat, mi, mj)

        # Set column ni of line mi to 1
        mat[mi] = list(map(lambda a: a*mat[mi][ni].inverse(), mat[mi]))
        # Set column ni of other rows to 0
        for mk in range(mi + 1, m):
            for nk in range(ni, n)[::-1]:
                mat[mk][nk] -= mat[mi][nk]*mat[mk][ni]

        mi += 1
        ni += 1

    # Find a row where all elements is not 0.
    for mi, row in enumerate(mat[::-1]):
        for a in row:
            if a != 0:
                return m - mi


def _find_none_zero_mi(mat: Matrix, ni: int, start_mi: int) -> int:
    for mi in range(start_mi, mat.m):
        if mat[mi][ni] != 0:
            return mi

    return -1


def _swap_rows(mat: Matrix, mi: int, mj: int) -> Matrix:
    mat[mi], mat[mj] = mat[mj], mat[mi]
    return mat
