"""Read placed STEP assembly occurrences, keeping original labels on the private side."""

from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path
import random

import numpy as np
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GeomAbs import GeomAbs_Cylinder
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPCAFControl import STEPCAFControl_Reader
from OCP.TCollection import TCollection_AsciiString, TCollection_ExtendedString
from OCP.TColStd import TColStd_SequenceOfAsciiString
from OCP.TDataStd import TDataStd_Name
from OCP.TDF import TDF_Label, TDF_LabelSequence, TDF_Tool
from OCP.TDocStd import TDocStd_Document
from OCP.TopAbs import TopAbs_FACE, TopAbs_REVERSED, TopAbs_SOLID
from OCP.TopExp import TopExp
from OCP.TopLoc import TopLoc_Location
from OCP.TopTools import TopTools_IndexedMapOfShape
from OCP.TopoDS import TopoDS
from OCP.XCAFDoc import XCAFDoc_DocumentTool, XCAFDoc_ShapeTool
from OCP.gp import gp_Trsf

from common import sha
from measure_geometry import bounds


def indexed(shape, kind):
    indices = TopTools_IndexedMapOfShape()
    TopExp.MapShapes_s(shape, kind, indices)
    return [indices.FindKey(i) for i in range(1, indices.Extent() + 1)]


def properties(shape, volume=False):
    props = GProp_GProps()
    if volume:
        BRepGProp.VolumeProperties_s(shape, props)
    else:
        BRepGProp.SurfaceProperties_s(shape, props)
    return props.Mass(), list(props.CentreOfMass().Coord())


def placed(shape, orientation, offset):
    tr = gp_Trsf()
    tr.SetValues(*np.column_stack((orientation,offset)).ravel())
    return shape.Moved(TopLoc_Location(tr))


def scan(shape):
    box = bounds(shape)
    faces = indexed(shape,TopAbs_FACE)
    cylinders = []
    for number,raw in enumerate(faces,1):
        face = TopoDS.Face_s(raw)
        surface = BRepAdaptor_Surface(face,True)
        if surface.GetType() != GeomAbs_Cylinder:
            continue
        cylinder = surface.Cylinder()
        axis = np.asarray(cylinder.Axis().Direction().Coord())
        point = np.asarray(cylinder.Location().Coord())
        first,last = surface.FirstVParameter(),surface.LastVParameter()
        angle = abs(surface.LastUParameter()-surface.FirstUParameter())
        if np.isfinite([first,last,angle]).all() and last-first>.025 and angle>.02:
            cylinders.append({"a":(point+first*axis).tolist(),"b":(point+last*axis).tolist(),
                              "r":cylinder.Radius(),"angle":angle,"face":number,
                              "inner":face.Orientation()==TopAbs_REVERSED})
    return {"bbox":box.tolist(),"dims":np.ptp(box,axis=0).tolist(),"volume":abs(properties(shape,True)[0]),
            "area":properties(shape)[0],"solid_count":len(indexed(shape,TopAbs_SOLID)),"faces":len(faces),
            "is_valid":BRepCheck_Analyzer(shape).IsValid(),"cylinders":cylinders}


def label_name(label):
    name = TDataStd_Name()
    return name.Get().ToExtString() if label.FindAttribute(TDataStd_Name.GetID_s(),name) else ""


def label_key(label):
    key = TCollection_AsciiString()
    TDF_Tool.Entry_s(label,key)
    return key.ToCString()


class Model:
    def __init__(self, source):
        self.source = Path(source).resolve(strict=True)
        if self.source.suffix.lower() not in {".step", ".stp"}:
            raise ValueError("The geometry analysis input must be a millimetre STEP assembly")
        self.source_hash = sha(self.source)
        reader = STEPCAFControl_Reader()
        reader.SetNameMode(True)
        reader.SetColorMode(False)
        reader.SetLayerMode(False)
        if reader.ReadFile(str(self.source)) != IFSelect_RetDone:
            raise ValueError("STEP parsing failed")
        units = [TColStd_SequenceOfAsciiString() for _ in range(3)]
        reader.Reader().FileUnits(*units)
        length = {units[0].Value(i).ToCString().lower() for i in range(1,units[0].Length()+1)}
        if not length or not length <= {"millimetre", "millimeter", "mm"}:
            raise ValueError(f"Convert or independently verify STEP length units first: {sorted(length)}")
        self._doc = TDocStd_Document(TCollection_ExtendedString("XmlXCAF"))
        if not reader.Transfer(self._doc):
            raise ValueError("STEP assembly transfer failed")
        self.geometry = {}
        self.rows = []
        parents = {""}

        def descend(label, location, path, names):
            definition = TDF_Label()
            if XCAFDoc_ShapeTool.GetReferredShape_s(label,definition):
                location = location.Multiplied(XCAFDoc_ShapeTool.GetLocation_s(label))
            else:
                definition = label
            caption = label_name(definition) or label_name(label)
            names = [*names,caption]
            children = TDF_LabelSequence()
            XCAFDoc_ShapeTool.GetComponents_s(definition,children)
            if children.Length():
                parents.add(path)
                for n in range(1,children.Length()+1):
                    descend(children.Value(n),location,f"{path}/{n:06d}",names)
                return
            key = label_key(definition)
            if key not in self.geometry:
                self.geometry[key] = XCAFDoc_ShapeTool.GetShape_s(definition)
            matrix = location.Transformation()
            rotation = np.asarray([[matrix.Value(i,j) for j in (1,2,3)] for i in (1,2,3)])
            if not np.allclose(rotation.T @ rotation,np.eye(3),atol=1e-7) or abs(np.linalg.det(rotation)-1)>1e-7:
                raise ValueError(f"Unsupported scaled or mirrored placement for input occurrence {path}")
            offset = matrix.TranslationPart().Coord()
            self.rows.append({"path":path,"definition":key,"origin_mm":offset,"rotation":rotation.T.ravel().tolist(),
                              "original_labels":names})

        roots = TDF_LabelSequence()
        XCAFDoc_DocumentTool.ShapeTool_s(self._doc.Main()).GetFreeShapes(roots)
        for n in range(1,roots.Length()+1):
            descend(roots.Value(n),TopLoc_Location(),f"{n:06d}",[])
        if not self.rows:
            raise ValueError("STEP has no inspectable occurrences")
        self.by_path = {r["path"]:r for r in self.rows}
        keys = sorted(self.by_path)
        random.Random(709313).shuffle(keys)
        self.ids = {key:f"I{i:05d}" for i,key in enumerate(keys,1)}
        self.paths = {v:k for k,v in self.ids.items()}
        branches = sorted(parents)
        random.Random(6977).shuffle(branches)
        self.parent = {key:f"T{i:04d}" for i,key in enumerate(branches,1)}
        self.model_id = "Sample"
        self.features = {}
        print(f"STEP occurrences {len(self.rows)} / unique definitions {len(self.geometry)}",flush=True)

    def prepare_definitions(self):
        for definition,shape in self.geometry.items():
            if shape.IsNull():
                self.features[definition] = {"error":"empty shape"}
                continue
            try:
                self.features[definition] = scan(shape)
            except (RuntimeError, ValueError, AssertionError) as exc:
                self.features[definition] = {"error":type(exc).__name__}
        print("STEP definition reads",Counter(v.get("error","read") for v in self.features.values()),flush=True)

    def shape(self, definition):
        return self.geometry[definition]

    def pose(self, ident):
        row = self.by_path[self.paths[ident]]
        return np.asarray(row["rotation"]).reshape(3,3).T,np.asarray(row["origin_mm"])

    @lru_cache(maxsize=1200)
    def local_shape(self, ident):
        row = self.by_path[self.paths[ident]]
        rotation,shift = self.pose(ident)
        return placed(self.shape(row["definition"]),self.axes@rotation,self.axes@(shift-self.origin))

    def restore_basis(self, base):
        self.origin = np.asarray(base["basis"]["origin"])
        self.axes = np.asarray(base["basis"]["axes"])

    def signatures(self, ident, require_complete=False):
        rotation,offset = self.pose(ident)
        source = self.shape(self.by_path[self.paths[ident]]["definition"])
        body = placed(source,rotation,offset)
        solids = indexed(body,TopAbs_SOLID)
        if require_complete:
            covered = TopTools_IndexedMapOfShape()
            for solid in solids:
                TopExp.MapShapes_s(solid,TopAbs_FACE,covered)
            if covered.Extent() != len(indexed(body,TopAbs_FACE)):
                return []
        signatures = []
        for solid in solids:
            v, xyz = properties(solid,True)
            signatures.append({"tight_box":bounds(solid).ravel().tolist(),"center":xyz,"volume":abs(v),
                               "faces":len(indexed(solid,TopAbs_FACE)),"area":properties(solid)[0]})
        return signatures

    def private_trace(self):
        return {ident:{"occurrence":path,"labels":self.by_path[path]["original_labels"]} for ident,path in self.paths.items()}

    def basis(self, external=None):
        rotation, shift = external if external is not None else (np.eye(3), np.zeros(3))
        candidates, centers, guides, plates = [], [], [], []
        for row in self.rows:
            data = self.features[row["definition"]]
            if "bbox" not in data:
                continue
            orientation, offset = self.pose(self.ids[row["path"]])
            matrix = rotation @ orientation
            offset = rotation@offset+shift
            box = np.asarray(data["bbox"])
            center = matrix @ box.mean(axis=0) + offset
            dims = np.ptp(box, axis=0)
            short, middle, long = sorted(dims)
            major = int(dims.argmax())
            if long >= 350:
                centers.append(center)
            if long >= 2400 and long >= 5*middle:
                candidates.append((float(long),matrix[:,major],matrix))
            if long >= 1000 and middle <= 85:
                if any(not c["inner"] and c["angle"] >= 2.6 and c["r"] >= 7 and
                       np.linalg.norm(np.asarray(c["b"])-c["a"]) >= .7*long for c in data["cylinders"]):
                    guides.append(matrix[:,major])
            if short <= 14 and middle >= 350 and long >= 550 and data["area"] > 50000:
                plates.append((center, matrix[:,int(dims.argmin())], middle*long, short))
        if not candidates or not guides:
            raise ValueError("Insufficient geometry for unsigned box axes")
        q = np.asarray(centers)
        _, vectors = np.linalg.eigh(np.cov(q.T))
        # Snap the major distribution direction to an actual long member.
        _, x, ref = max(candidates, key=lambda p: (round(abs(p[1]@vectors[:,-1]),4),
                                                   round(p[0],5)))
        support = np.abs(np.asarray(guides)@ref).sum(axis=0)
        support[np.argmax(np.abs(ref.T@x))] = -1
        z = ref[:,int(support.argmax())]
        coord = q@z
        middle = (coord.max()+coord.min())*.5
        halves = np.zeros(2)
        for center, normal, area, thick in plates:
            if abs(normal@z)>.98:
                # Broad horizontal deck and lower continuous closure contribute more
                # plate section than the upper skins in these samples.
                halves[int(center@z > middle)] += area*min(thick,12)
        if min(halves) <= 0 or max(halves)/min(halves)<1.13 or support.max() < 4:
            raise ValueError(f"Unresolved top/bottom polarity: masses {halves}, guides {support}")
        if halves.argmax() == 1:
            z = -z
        if np.linalg.norm(x) < .99 or np.linalg.norm(z) < .99 or abs(x@z)>1.e-7:
            raise ValueError("Box axes not verified")
        y = np.cross(z,x)
        axes = np.array([x,y,z])
        assert abs(np.linalg.det(axes)-1)<1.e-7
        origin = q.mean(axis=0)
        return origin, axes, {"long_member_mm":round(max(t[0] for t in candidates),4),
                              "guide_axis_support":round(float(support.max()),4),
                              "upper_lower_plate_proxy":sorted(np.round(halves,3).tolist()),
                              "axes_source":"geometry and occurrence transforms"}

    def make_base(self, datum=None, external=None):
        if datum is None:
            origin, axes, basis_info = self.basis(external)
        else:
            origin = np.asarray(datum["origin"], dtype=float)
            axes = np.asarray(datum["axes"], dtype=float)
            if (origin.shape != (3,) or axes.shape != (3,3) or not np.isfinite(origin).all()
                    or not np.allclose(axes @ axes.T,np.eye(3),atol=1e-7) or abs(np.linalg.det(axes)-1)>1e-7):
                raise ValueError("Datum must be a right-handed orthonormal millimetre reference")
            basis_info = {"axes_source":"operator datum"}
        rotation, shift = external if external is not None else (np.eye(3), np.zeros(3))
        original_origin = rotation.T@(origin-shift)
        original_axes = axes@rotation
        self.origin, self.axes = original_origin, original_axes
        instances = {}
        for path, ident in self.ids.items():
            row = self.by_path[path]
            d = row["definition"]
            data = self.features[d]
            item = {"parent": self.parent[path.rpartition("/")[0]], "availability": "step_visible"}
            if "bbox" not in data:
                item["availability"] = "read_error"
            else:
                ori, off = self.pose(ident)
                matrix, offset = original_axes@ori, original_axes@(off-original_origin)
                # Original topology-vertex envelope is used for screening. Use
                # analytic faces for paired tests, and tight BRep bounds for export.
                local_box = np.asarray(data["bbox"])
                xyz = np.asarray(list(product(*local_box.T)))@matrix.T + offset
                limits = np.array([xyz.min(axis=0),xyz.max(axis=0)])
                if np.ptp(limits,axis=0).max()>170:
                    # A rotated AABB loses e.g. diagonal orientation. Tighten only
                    # rotated candidates with appreciable extents.
                    if np.min(np.max(np.abs(matrix),axis=0)) < .999999:
                        shape = placed(self.shape(d),matrix,offset)
                        limits = bounds(shape)
                size = np.ptp(limits,axis=0)
                a = sorted(data["dims"])
                cyl = data["cylinders"]
                axis_lengths=[np.linalg.norm(np.array(c["b"])-c["a"]) for c in cyl]
                full = [c for c,s in zip(cyl,axis_lengths) if c["angle"]>=2.55 and not c["inner"]
                        and s>=.70*max(a)]
                holes=[]
                for c,s in zip(cyl,axis_lengths):
                    if c["inner"] and c["r"]>=8 and c["angle"]>1.05 and s>.5:
                        direction = matrix@(np.asarray(c["b"])-c["a"])/s
                        holes.append(tuple(np.round(np.abs(direction),3)))
                item.update({"box": np.round(limits,5).tolist(), "shape_mm": a,
                             "direction": np.round(matrix[:,int(np.argmax(data["dims"]))],6).tolist(),
                             "solids": data["solid_count"], "faces": data["faces"],
                             "volume": round(data["volume"],4), "area": round(data["area"],4),
                             "valid": data["is_valid"],
                             "shaft": round(max((c["r"] for c in full),default=0),4),
                             "opening_dirs": sorted(set(holes)),
                             "concentric_candidate": any(c["angle"]>=2.8 for c in cyl)})
            instances[ident] = item
        # Use true physical extents of medium/large geometry. The long reference
        # is based on these surfaces, not on the origin of an assembly file.
        spans = [r["box"] for r in instances.values() if "box" in r and max(r["shape_mm"])>=400]
        frame = np.array([np.asarray(spans)[:,0].min(axis=0),np.asarray(spans)[:,1].max(axis=0)])
        tree = {v:(self.parent[p.rpartition("/")[0]] if p else None) for p,v in self.parent.items()}
        packet = {"model":self.model_id,"configuration":"C01","unit":"mm","frame_box":frame.tolist(),
                  "basis":{"origin":origin.tolist(),"axes":axes.tolist(),**basis_info},
                  "parts":instances,"tree":tree}
        return packet
