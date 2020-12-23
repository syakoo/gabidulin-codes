from typing import List

from galois_field import GFpn, ElementInGFpn


class GabidulinCodes:
    def __init__(self, F: GFpn, gs: List[ElementInGFpn], k: int):
        self.F = F
        self.gs = gs
        self.__n = len(gs)
        self.__k = k
        self.__d = self.__n - self.__k + 1

        self.__G = [gs]
        for i in range(1, k):
            G_row = []
            exp_pi = F.p ** i
            for g in gs:
                # g^[i] = g^{p^i}
                frob_pow = g ** exp_pi
                G_row.append(frob_pow)

            self.__G.append(G_row)

    @property
    def G(self) -> List[List[ElementInGFpn]]:
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
