"""Exact trimmed CAD shape bounds and pairwise distance, in millimetres."""
import numpy as np
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepExtrema import BRepExtrema_DistShapeShape

def bounds(shape):
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, box, False, False)
    return np.asarray(box.Get()).reshape(2, 3)


def gap(first, second):
    judge = BRepExtrema_DistShapeShape(first, second)
    assert judge.IsDone() and judge.NbSolution()
    return {
        "distance_mm": judge.Value(),
        "point_1_mm": list(judge.PointOnShape1(1).Coord()),
        "point_2_mm": list(judge.PointOnShape2(1).Coord()),
        "inner_solution": judge.InnerSolution(),
    }
