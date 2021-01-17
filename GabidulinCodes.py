from typing import List

from galois_field import GFpn, ElementInGFpn

from src import encode
from src.core.matrix import MatGFpn


class GabidulinCodes:
    def __init__(self, F: GFpn, gs: List[ElementInGFpn], k: int):
        self.F = F
        self.gs = gs
        self.__n = len(gs)
        self.__k = k
        self.__d = self.__n - self.__k + 1

        G_values = [gs]
        for i in range(1, k):
            G_row = []
            exp_pi = F.p ** i
            for g in gs:
                # g^[i] = g^{p^i}
                frob_pow = g ** exp_pi
                G_row.append(frob_pow)

            G_values.append(G_row)
        self.__G = MatGFpn(G_values)

    def encode(self, x: List[ElementInGFpn]) -> List[ElementInGFpn]:
        codeword = encode.encode(x, self.__G)
        return codeword

    @property
    def G(self) -> MatGFpn:
        return self.__G

    @property
    def n(self) -> int:
        return self.__n

    @property
    def k(self) -> int:
        return self.__k

    @property
    def d(self) -> int:
        return self.__d
