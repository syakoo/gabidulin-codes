from __future__ import annotations
from typing import List, Union

from galois_field import ElementInGFp, ElementInGFpn, GFp, GFpn

ElmInGF = Union[ElementInGFp, ElementInGFpn]


class MatGFpn:
    def __init__(self, mat: List[List[ElmInGF]]):
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


if __name__ == "__main__":
    GF = GFp(5)

    Mat1 = MatGFpn.from_int_values([[1, 2, 3], [4, 5, 6]], GF)
    Mat2 = MatGFpn.from_int_values([[1, 3, 5], [2, 4, 6]], GF)
    print(Mat1)
    print(Mat2)
    print(f"add:\n{Mat1 + Mat2}")
    print(f"sub:\n{Mat1 - Mat2}")
    print(f"tMat2:\n{Mat2.transpose()}")
    print(f"mul:\n{Mat1 * Mat2.transpose()}")

    GF2 = GFpn(5, [1, 0, 0, 0, 2])
    values = [[1, 2, 3], [4, 5, 6]]
    Mat3 = MatGFpn.from_int_values(values, GF2)
    Mat4 = MatGFpn.from_vect_over_GFpn(list(map(GF2.elm, values)))
    print(Mat3)
    print(Mat4)
