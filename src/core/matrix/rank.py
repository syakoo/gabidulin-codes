from .types import MatValues


def _seek_rank_over_Field(values: MatValues) -> int:
    """Find the rank of a matrix over GF(p).

    Args:
        values (MatValues): M_ij in GF(p^n).

    Returns:
        int: The rank of a matrix over GF(p).
    """
    # must be m < n
    m = len(values)
    n = len(values[0])

    mi, ni = 0, 0

    # Convert the matrix to row echelon form.
    while mi < m and ni < n:
        # Get the line number after the mi line where the ni column is not 0.
        mj = _find_none_zero_mi(values, ni, mi)
        if mj == -1:
            ni += 1
            continue
        # Exchange mi and mj lines
        if mi != mj:
            values = _swap_rows(values, mi, mj)

        # Set column ni of line mi to 1
        values[mi] = list(
            map(lambda a: a*values[mi][ni].inverse(), values[mi]))
        # Set column ni of other rows to 0
        for mk in range(mi + 1, m):
            for nk in range(ni, n)[::-1]:
                values[mk][nk] -= values[mi][nk]*values[mk][ni]

        mi += 1
        ni += 1

    # Find a row where all elements is not 0.
    for mi, row in enumerate(values[::-1]):
        for a in row:
            if a != 0:
                return m - mi


def _find_none_zero_mi(values: MatValues, ni: int, start_mi: int) -> int:
    for mi in range(start_mi, len(values)):
        if values[mi][ni] != 0:
            return mi

    return -1


def _swap_rows(values: MatValues, mi: int, mj: int) -> MatValues:
    values[mi], values[mj] = values[mj], values[mi]
    return values
