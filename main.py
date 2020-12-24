from typing import List

from galois_field import GFpn, ElementInGFp

from GabidulinCodes import GabidulinCodes
from src.Matrix import Matrix


def print_mat_Fp(mat: List[List[ElementInGFp]]) -> None:
    print(list(map(lambda row: list(map(lambda a: int(a), row)), mat)))


def main():
    F = GFpn(5, [1, 0, 0, 0, 2])
    gs = [F.elm([1, 2]), F.elm([1, 0, 1]), F.elm([3, 2, 0])]
    Gab = GabidulinCodes(F, gs, 3)
    print(F, Gab.G)

    # gs_as_mat = Matrix.from_list(gs)
    # print(gs_as_mat.rank)
    message = [F.elm([1, 0, 0]), F.elm([0, 2, 0]), F.elm([0, 0, 3])]
    codeword = Gab.encode(message)
    print(codeword)


if __name__ == "__main__":
    main()
