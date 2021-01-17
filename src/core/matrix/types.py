from typing import List, Union

from galois_field import ElementInGFp, ElementInGFpn


ElmInGF = Union[ElementInGFp, ElementInGFpn]
MatValues = List[List[ElmInGF]]
