from typing import List

from galois_field import ElementInGFpn


def encode(x: List[ElementInGFpn], G: List[List[ElementInGFpn]]) -> List[ElementInGFpn]:
    result = []

    for ni in range(len(G[0])):
        for mi, xi in enumerate(x):
            if mi == 0:
                result.append(xi*G[mi][ni])
            else:
                result[ni] += xi*G[mi][ni]

    return result
