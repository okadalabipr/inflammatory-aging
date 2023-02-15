from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "kf1",
    "kf2",
    "kf3",
    "kf4",
    "kf5",
    "kf6",
    "kf7",
    "kf8",
    "kf9",
    "kr9",
    "kf10",
    "kr10",
    "kf11",
    "kf12",
    "kf13",
    "kf14",
    "kf15",
    "kf16",
    "kf17",
    "kf18",
    "kf19",
    "V20",
    "K20",
    "n20",
    "kf21",
    "kf22",
    "kf23",
    "kf24",
    "V25",
    "K25",
    "n25",
    "kf26",
    "kf27",
    "kf28",
]

NUM: int = len(NAMES)

Parameters = make_dataclass(
    cls_name="Parameters",
    fields=[(name, int) for name in NAMES],
    namespace={"NAMES": NAMES, "NUM": NUM},
    frozen=True,
)

name2idx: Dict[str, int] = {k: v for v, k in enumerate(NAMES)}

C = Parameters(**name2idx)

del name2idx
