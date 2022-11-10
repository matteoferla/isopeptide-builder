import enum
import types
from typing import (Optional, Dict)
import pyrosetta

prc: types.ModuleType = pyrosetta.rosetta.core


# ---------------------------------------------------------------------------------------------------------

class IsopetideType(enum.Enum):
    ASP_LYS_NONE = enum.auto()  # noqa
    ASP_LYS_GLU = enum.auto()  # noqa

def create_weighted_scorefxn(
                             weights: Dict[prc.scoring.ScoreType, float],
                             score_name: str = 'ref2015_cart') -> pyrosetta.ScoreFunction:
    scorefxn = pyrosetta.create_score_function(score_name)
    for weight_name, weight in weights.items():
        scorefxn.set_weight(weight_name, weight)
    return scorefxn
