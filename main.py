from typing import List

from galois_field import GFp, GFpn, ElementInGFp

from GabidulinCodes import GabidulinCodes
from src.core.matrix.MatGFpn import MatGFpn


def print_mat_Fp(mat: List[List[ElementInGFp]]) -> None:
    print(list(map(lambda row: list(map(lambda a: int(a), row)), mat)))


def main():
    F = GFpn(5, [1, 0, 0, 0, 2])
    gs = [F.elm([1, 2]), F.elm([1, 0, 1]), F.elm([3, 2, 0])]
    Gab = GabidulinCodes(F, gs, 3)
    print(F, Gab.G, sep="\n")

    # gs_as_mat = Matrix.from_list(gs)
    # print(gs_as_mat.rank)
    message = [F.elm([1, 0, 0]), F.elm([0, 2, 0]), F.elm([0, 0, 3])]
    codeword = Gab.encode(message)
    print(codeword)


def foo():
    GF = GFp(5)

    Mat1 = MatGFpn.from_int_values([[1, 2, 3], [4, 5, 6]], GF)
    Mat2 = MatGFpn.from_int_values([[1, 3, 5], [2, 4, 6]], GF)
    print(Mat1)
    print(Mat1.rank)
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
    print("="*20)
    print("inverse:")
    Mat = MatGFpn.from_int_values(([1, 2, 3], [2, 4, 3], [7, 8, 9]), GF)
    inv_Mat = Mat.inverse()
    print(Mat)
    print(inv_Mat)
    print(Mat*inv_Mat)


if __name__ == "__main__":
    main()
    # foo()
