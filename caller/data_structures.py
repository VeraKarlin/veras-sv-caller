from dataclasses import dataclass, field
from enum import Enum

@dataclass
class Alignment:
    query_name: str
    is_primary: bool
    chrom: str
    query_start: int
    query_end: int
    ref_start: int
    ref_end: int
    strand: str
    cigar: str
    mapq: float
    phase: int = 0

class SVType(Enum):
    INV = "inversion"
    DEL = "deletion"
    INS = "insertion"
    DUP = "duplication"
    BND = "breakend of translocation"
    TRANSLOCATION_BALANCED = "translocation_balanced"
    TRANSLOCATION_UNBALANCED = "translocation_unbalanced"

@dataclass
class BPInfo:
    chrom: str
    pos: int
    strand: float
    direction: float
    cov_before: float
    cov_after: float

@dataclass
class SVInfo:
    first: BPInfo
    second: BPInfo
    sv_type: SVType
    sequence: str = ""
    support: int = 1
    phase: int = 0
    phase_ratio: float = 1.0
    sv_pipeline: str = ""
    is_inverted: bool = False