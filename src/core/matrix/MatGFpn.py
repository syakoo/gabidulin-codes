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

    def __str__(self) -> str:
        str_rows = map(lambda row: " ".join(map(str, row)), self.__values)
        return "[ " + "\n  ".join(str_rows) + " ]"

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


if __name__ == "__main__":
    GF = GFp(5)

    Mat1 = MatGFpn.from_int_values([[1, 2, 3], [4, 5, 6]], GF)
    Mat2 = MatGFpn.from_int_values([[1, 3, 5], [2, 4, 6]], GF)
    print(Mat1)
    print(Mat2)
