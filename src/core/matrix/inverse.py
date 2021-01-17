
from .types import MatValues


def seek_inverse(values: MatValues) -> MatValues:
    n = len(values)
    inv_values = [[1 if i == j else 0 for i in range(n)] for j in range(n)]

    for i in range(n):
        # Get the line number after the mi line where the ni column is not 0
        mj = _find_none_zero_mi(values, i, i)
        if mj == -1:
            raise ValueError("Could not find the inverse of the matrix.")

        # Exchange i and mj lines
        if i != mj:
            values = _swap_rows(values, i, mj)
            inv_values = _swap_rows(inv_values, i, mj)

        # Set column ni of line mi to 1
        inv_a = values[i][i].inverse()
        values[i] = list(map(lambda val: val*inv_a, values[i]))
        inv_values[i] = list(map(lambda val: val*inv_a, inv_values[i]))

        # Set column ni of other rows to 0
        for mi in range(n):
            if mi == i:
                continue
            val = values[mi][i]
            for ni in range(n):
                if ni >= i:
                    values[mi][ni] = values[mi][ni] - values[i][ni]*val

                inv_values[mi][ni] = inv_values[mi][ni] - inv_values[i][ni]*val

    return inv_values


def _find_none_zero_mi(values: MatValues, ni: int, start_mi: int) -> int:
    for mi in range(start_mi, len(values)):
        if values[mi][ni] != 0:
            return mi

    return -1


def _swap_rows(values: MatValues, mi: int, mj: int) -> MatValues:
    values[mi], values[mj] = values[mj], values[mi]
    return values
