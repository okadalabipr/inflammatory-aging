from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "TNF",
    "TNFR",
    "TTR",
    "aTTR",
    "IKK",
    "pIKK",
    "iIKK",
    "IkBa",
    "NFkB",
    "IkBa_NFkB",
    "IkBan",
    "NFkBn",
    "IkBa_NFkBn",
    "pIkBa_NFkB",
    "pIkBa",
    "ikbamRNA",
    "a20mRNA",
    "A20",
]

NUM: int = len(NAMES)

Species = make_dataclass(
    cls_name="Species",
    fields=[(name, int) for name in NAMES],
    namespace={"NAMES": NAMES, "NUM": NUM},
    frozen=True,
)

name2idx: Dict[str, int] = {k: v for v, k in enumerate(NAMES)}

V = Species(**name2idx)

del name2idx
