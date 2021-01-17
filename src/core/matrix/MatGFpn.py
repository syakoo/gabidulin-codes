from __future__ import annotations
from typing import List, Union

from galois_field import ElementInGFpn, GFp, GFpn

from . import rank
from .types import ElmInGF, MatValues


class MatGFpn:
    def __init__(self, mat: MatValues):
        self.__values = mat
        self.__m = len(mat)
        self.__n = len(mat[0])

    def __getitem__(self, i: int) -> List[ElmInGF]:
        return self.__values[i]

    def __setitem__(self, i: int, value: List[ElmInGF]) -> None:
        self.__values[i] = value

    @property
    def m(self) -> int:
        return self.__m

    @property
    def n(self) -> int:
        return self.__n

    @property
    def rank(self) -> int:
        return rank.seek_rank_over_Field(self.__values)

    def __add__(self, other: MatGFpn) -> MatGFpn:
        if not (self.m == other.m and self.n == other.n):
            raise ValueError("Two matrices have different sizes.")

        result_values = [[self[j][i] + other[j][i]
                          for i in range(self.n)] for j in range(self.m)]
        return MatGFpn(result_values)

    def __sub__(self, other: MatGFpn) -> MatGFpn:
        if not (self.m == other.m and self.n == other.n):
            raise ValueError("Two matrices have different sizes.")

        result_values = [[self[j][i] - other[j][i]
                          for i in range(self.n)] for j in range(self.m)]
        return MatGFpn(result_values)

    def __mul__(self, other: MatGFpn) -> MatGFpn:
        if not (self.n == other.m):
            raise ValueError("Two matrices have different sizes.")

        result_values = []
        for j in range(self.m):
            row = []
            for i in range(other.n):
                val = 0
                for k in range(self.n):
                    val = val + self[j][k]*other[k][i]
                row.append(val)
            result_values.append(row)
        return MatGFpn(result_values)

    def __str__(self) -> str:
        str_rows = map(lambda row: " ".join(map(str, row)), self.__values)
        return "[ " + "\n  ".join(str_rows) + " ]"

    def transpose(self) -> MatGFpn:
        result_values = [[self[j][i]
                          for j in range(self.m)] for i in range(self.n)]
        return MatGFpn(result_values)

    @staticmethod
    def from_int_values(values: List[List[int]], F: Union[GFp, GFpn]) -> MatGFpn:
        if isinstance(F, GFp):
            applyed_values = list(
                map(lambda row: list(map(F.elm, row)), values))
        elif isinstance(F, GFpn):
            applyed_values = list(
                map(lambda row: list(map(lambda val: F.elm([val]), row)), values))
        else:
            raise TypeError("F must be an instance of GFp or GFpn.")

        return MatGFpn(applyed_values)

    @staticmethod
    def from_vect_over_GFpn(vector: List[ElementInGFpn]) -> MatGFpn:
        F = GFp(vector[0].p)
        result = []
        max_len = len(vector[0].mod_poly.coeffs) - 1
        for vi in vector:
            result.append([F.elm(0)]*(max_len - len(vi.coeffs)) +
                          list(map(F.elm, vi.coeffs)))

        return MatGFpn(result)
