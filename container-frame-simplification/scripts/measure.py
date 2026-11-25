"""Label-free measurements of through-bolted stacks and thin-panel interfaces."""

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations

import numpy as np
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCP.BRepTools import BRepTools
from OCP.GeomAbs import GeomAbs_Circle, GeomAbs_Line, GeomAbs_Plane
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_REVERSED, TopAbs_WIRE
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

from intake import indexed, placed, properties
from recognize import coords, geometry
from relations import Graph, intervals

FIT_MM = .12
TOUCH_MM = .03
PARALLEL = .9999


def segment_distance(point, segment):
    vec = segment[1] - segment[0]
    t = 0 if np.linalg.norm(vec)<1.e-8 else np.clip(np.dot(point-segment[0],vec)/np.dot(vec,vec),0,1)
    return float(np.linalg.norm(point-segment[0]-t*vec))


def opening_rims(shape):
    """Recover each aperture from a single trimmed inner wire of a planar face."""
    rims=[]
    for face_number, raw in enumerate(indexed(shape,TopAbs_FACE),1):
        face=TopoDS.Face_s(raw)
        surf=BRepAdaptor_Surface(face,True)
        if surf.GetType()!=GeomAbs_Plane:
            continue
        direction=np.array(surf.Plane().Axis().Direction().Coord())
        outside=BRepTools.OuterWire_s(face)
        for raw_wire in indexed(face,TopAbs_WIRE):
            wire=TopoDS.Wire_s(raw_wire)
            if not outside.IsNull() and wire.IsSame(outside):
                continue
            arc=defaultdict(float)
            straight=[]
            for raw_edge in indexed(wire,TopAbs_EDGE):
                curve=BRepAdaptor_Curve(TopoDS.Edge_s(raw_edge))
                start, end=curve.FirstParameter(),curve.LastParameter()
                if curve.GetType()==GeomAbs_Circle:
                    circle=curve.Circle()
                    r=circle.Radius()
                    if not 2.2<=r<=16 or abs(np.dot(direction,circle.Axis().Direction().Coord()))<PARALLEL:
                        break
                    center=tuple(np.round(circle.Location().Coord(),6))
                    arc[(center,round(r,6))]+=abs(end-start)
                elif curve.GetType()==GeomAbs_Line:
                    delta=np.array(curve.Value(end).Coord())-curve.Value(start).Coord()
                    straight.append(delta)
                else:
                    break
            else:
                if len(arc)==1 and not straight and next(iter(arc.values()))>=6.27:
                    (center, radius),=arc
                    ends=np.array([center,center])
                    kind="圆孔"
                elif len(arc)==2 and len(straight)==2 and min(arc.values())>=3.11:
                    (a,radius),(b,rb)=sorted(arc)
                    ends=np.array([a,b])
                    span=ends[1]-ends[0]
                    if abs(radius-rb)>1.e-3 or np.linalg.norm(span)<.5 or any(
                            min(np.linalg.norm(d-span),np.linalg.norm(d+span))>.03 for d in straight):
                        continue
                    kind="长圆孔"
                else:
                    continue
                rims.append((face_number,direction,ends,radius,kind))
    return rims


def through_openings(rims):
    openings=[]
    for left,right in combinations(rims,2):
        ia,normal,a,r,k=left
        ib,other,b,rb,kb=right
        if ia==ib or k!=kb or abs(r-rb)>.001 or abs(normal@other)<PARALLEL:
            continue
        signed=float((b[0]-a[0])@normal)
        if abs(signed)<1.2 or np.max(np.abs(b-a-signed*normal))>.02:
            continue
        openings.append((ia,ib,normal,a,b,r,k))
    return openings


def hit(opening, origin, axis, shaft_length, radius):
    fa,fb,n,a,b,r,kind=opening
    proj=float(n@axis)
    if abs(proj)<PARALLEL or r<radius-FIT_MM or r-radius>max(1.5,.40*radius):
        return None
    t=[float((centers[0]-origin)@n)/proj for centers in (a,b)]
    offset = max(segment_distance(origin+dist*axis,centers) for dist,centers in zip(t,(a,b)))
    if offset>max(0,r-radius)+FIT_MM:
        return None
    lo,hi=sorted(t)
    overlap=min(shaft_length,hi)-max(0,lo)
    if overlap<max(1.0,.85*(hi-lo)):
        return None
    faces=(fa,fb) if t[0]<t[1] else (fb,fa)
    return {"mouth_faces": list(faces), "hole":kind,"radius_mm":r,"offset_mm":offset,
            "axial_overlap_mm":overlap,"interval_mm":[lo,hi]}


class Scanner:
    def __init__(self, graph):
        self.g=graph
        self.model=graph.model
        self.instance={self.model.ids[r["path"]]:r for r in self.model.rows}
        self.siblings=defaultdict(list)
        self.counts=Counter()
        self.plates=[]
        for i in self.g.ids:
            p=self.g.parts[i]
            self.siblings[p["parent"]].append(i)
            a,b,c=p["shape_mm"]
            if 1.2<=a<=.43*b and b>=24 and c>=53 and p["concentric_candidate"]:
                self.plates.append(self.g.look[i])
        self.plates=np.asarray(self.plates,dtype=int)

    @lru_cache(maxsize=1700)
    def planes_by_number(self,i):
        result={}
        for number,face,area,center,normal in self.g.planes(i):
            outward=-normal if face.Orientation()==TopAbs_REVERSED else normal
            result[number]=(face,area,center,outward)
        return result

    @lru_cache(maxsize=1000)
    def apertures(self,definition):
        return through_openings(opening_rims(self.model.shape(definition)))

    def placed_apertures(self,i):
        definition=self.instance[i]["definition"]
        orientation,center=self.model.pose(i)
        trans=self.model.axes@orientation
        shift=self.model.axes@(center-self.model.origin)
        return [(fa,fb,trans@n,a@trans.T+shift,b@trans.T+shift,r,k)
                for fa,fb,n,a,b,r,k in self.apertures(definition)]

    def project_faces(self,a,fa,b,fb,origin=None,axis=None):
        left,la,lc,ln=self.planes_by_number(a)[fa]
        right,ra,rc,rn=self.planes_by_number(b)[fb]
        cosine=float(-ln@rn)
        if cosine<PARALLEL:
            return None
        right_pivot=rc
        if axis is None:
            gap=float((rc-lc)@ln)
        else:
            ta=float(((lc-origin)@ln)/(axis@ln))
            tb=float(((rc-origin)@rn)/(axis@rn))
            gap=tb-ta
            right_pivot=origin+tb*axis
        x,y,z=np.cross(rn,-ln)
        skew=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
        theta=float(np.degrees(np.arctan2(np.linalg.norm([x,y,z]),cosine)))
        # Align nearly parallel faces at the measured shaft intercept for
        # a virtual projected overlap. The measured gap and tilt stay separate.
        rotation=np.eye(3)+skew+skew@skew/(1+cosine)
        shift=right_pivot-rotation@right_pivot-((right_pivot-lc)@ln)*ln
        projected=placed(right,rotation,shift)
        common=BRepAlgoAPI_Common(left,projected)
        if not common.IsDone() or common.Shape().IsNull():
            return None
        area,center=properties(common.Shape())
        if not 0.05<area<=min(la,ra)+.08:
            return None
        return {"parts":[a,b],"faces":[fa,fb],"gap_mm":round(gap,6),
                "normal_angle_deg":round(theta,6),"projected_area_mm2":round(area,6),
                "point_mm":np.round(center,6).tolist(),"method":"opposed_plane_projection"}

    @lru_cache(maxsize=4500)
    def patch(self,a,b,tolerance=TOUCH_MM):
        candidates=[]
        for j,(_,aa,ac,an) in self.planes_by_number(a).items():
            for k,(_,ba,bc,bn) in self.planes_by_number(b).items():
                if min(aa,ba)<4 or an@bn>-PARALLEL or abs((bc-ac)@an)>tolerance:
                    continue
                result=self.project_faces(a,j,b,k)
                if result and result["normal_angle_deg"]<=.001:
                    candidates.append(result)
        return max(candidates,key=lambda r:r["projected_area_mm2"],default=None)

    @lru_cache(maxsize=2500)
    def backings(self,plate):
        p=self.g.parts[plate]
        i=self.g.look[plate]
        positions=tuple(self.siblings[p["parent"]])
        candidates=[]
        for other in positions:
            if other==plate:
                continue
            a,b,c=self.g.parts[other]["shape_mm"]
            if a<20 or c<max(360,3.2*b):
                continue
            j=self.g.look[other]
            if np.max(intervals(self.g.lo[j],self.g.hi[j],self.g.lo[i],self.g.hi[i]))>TOUCH_MM:
                continue
            contact=self.patch(plate,other)
            if contact and contact["projected_area_mm2"]>=12:
                candidates.append(contact)
                continue
            nearby=self.g.nearest(plate,other)
            if not nearby or not nearby["inner"]:
                continue
            common=BRepAlgoAPI_Common(self.model.local_shape(plate),self.model.local_shape(other))
            if not common.IsDone() or common.Shape().IsNull():
                continue
            vol,center=properties(common.Shape(),volume=True)
            if 1.0<vol<=.30*p["volume"]:
                candidates.append({"parts":[plate,other],"faces":[],"gap_mm":0,
                                   "method":"local_solid_overlap","common_volume_mm3":round(vol,6),
                                   "point_mm":np.round(center,6).tolist()})
        return sorted(candidates,key=lambda r:(r["method"]!="opposed_plane_projection",-r.get("projected_area_mm2",0),
                                                -r.get("common_volume_mm3",0)))

    def ring(self,i,origin,axis,length,rad):
        p=self.g.parts[i]
        if p["shape_mm"][-1]>min(67,9*rad):
            return None
        accepted=[]
        for a,b,r,n in self.g.joint_cylinders(i,True):
            if not -FIT_MM<=r-rad<=max(.95,.18*rad):
                continue
            rel=np.array([a,b])-origin
            dist=rel@axis
            off=float(np.linalg.norm(rel-dist[:,None]*axis,axis=1).max())
            lap=min(length,float(dist.max()))-max(0,float(dist.min()))
            if off<=FIT_MM and lap>.25:
                accepted.append((off,lap,r,n))
        if not accepted:
            return None
        faces=self.planes_by_number(i)
        extent={}
        for sign in (-1,1):
            fits=[(float((center-origin)@axis),-area,j) for j,(_,area,center,n) in faces.items()
                  if sign*n@axis>PARALLEL and area>=4]
            if not fits:
                return None
            extent[sign]=min(fits,key=lambda t:(-sign*t[0],t[1]))
        lo,_,lf=extent[-1]
        hi,_,hf=extent[1]
        if not .45<hi-lo<=3.2*rad:
            return None
        off,lap,r,n=min(accepted,key=lambda t:(t[0],-t[1],t[3]))
        return {"instance":i,"rim_faces":[lf,hf],"interval_mm":[lo,hi],"bore_face":n,
                "radius_mm":r,"offset_mm":off,"axial_overlap_mm":lap}

    def plate_hits(self,bolt,origin,axis,length,rad):
        end=origin+length*axis
        window=np.array([np.minimum(origin,end)-rad-2,np.maximum(origin,end)+rad+2])
        nearby=self.plates[np.max(intervals(self.g.lo[self.plates],self.g.hi[self.plates],
                                             window[0],window[1]),axis=1)<FIT_MM]
        result=[]
        bolt_parent=self.g.parts[bolt]["parent"]
        for ix in nearby:
            i=self.g.ids[ix]
            if i==bolt or self.g.parts[i]["parent"]==bolt_parent:
                continue
            matches=[m for f in self.placed_apertures(i) if (m:=hit(f,origin,axis,length,rad))]
            if not matches:
                continue
            backing=self.backings(i)
            if backing:
                selected=min(matches,key=lambda v:(v["offset_mm"],-v["axial_overlap_mm"]))
                lo,hi=selected["interval_mm"]
                start=origin+(lo+.04)*axis
                probe=BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(*start),gp_Dir(*axis)),rad-.03,hi-lo-.08).Shape()
                common=BRepAlgoAPI_Common(self.model.local_shape(i),probe)
                if not common.IsDone():
                    continue
                obstruction=(abs(properties(common.Shape(),volume=True)[0]) if not common.Shape().IsNull() else 0)
                if obstruction>.05:
                    self.counts["shaft_bore_probe_blocked"]+=1
                    continue
                selected["shaft_probe_common_volume_mm3"]=round(obstruction,6)
                result.append(dict(selected,instance=i,backing=backing[0]))
        return sorted(result,key=lambda s:sum(s["interval_mm"]))

    def clamp(self,bolt,origin,axis,length,rad,shaft_face):
        head_planes=[(float((c-origin)@axis),number) for number,(_,area,c,n) in self.planes_by_number(bolt).items()
                     if n@axis>PARALLEL and area>.32*np.pi*rad**2
                     and -max(2.,.60*rad)<=float((c-origin)@axis)<=.2]
        if not head_planes:
            return None
        t_head,head_face=max(head_planes)
        rings=[a for i in self.siblings[self.g.parts[bolt]["parent"]]
               if i!=bolt and (a:=self.ring(i,origin,axis,length,rad))]
        if len(rings)<2:
            return None
        tail=sorted((a for a in rings if a["interval_mm"][0]>.48*length and
                     np.ptp(a["interval_mm"])>.72*rad),key=lambda a:-a["interval_mm"][1])
        if not tail:
            return None
        nut=tail[0]
        hosts=self.plate_hits(bolt,origin,axis,length,rad)
        if not hosts:
            return None
        best=[]
        for a,b in combinations(hosts,2):
            if a["backing"]["parts"][1]==b["backing"]["parts"][1]:
                continue
            if a["interval_mm"][0]<t_head-FIT_MM or b["interval_mm"][1]>nut["interval_mm"][0]+FIT_MM:
                continue
            if b["interval_mm"][0]<a["interval_mm"][1]-min(4.,.52*rad):
                continue
            layers=[dict(x,rim_faces=x["mouth_faces"]) for x in (a,b)]
            nuts_at_end=[r for r in rings if r!=nut and
                         t_head-FIT_MM<=r["interval_mm"][0]<=nut["interval_mm"][0]+FIT_MM]
            spacers=sorted(nuts_at_end+layers,key=lambda x:x["interval_mm"][0])
            chain=[{"instance":bolt,"rim_faces":[head_face,head_face],"interval_mm":[t_head,t_head]},*spacers,nut]
            contacts=[self.project_faces(p["instance"],p["rim_faces"][1],q["instance"],q["rim_faces"][0],origin,axis)
                      for p,q in zip(chain,chain[1:])]
            if any(x is None or x["projected_area_mm2"]<max(3,.12*rad**2) for x in contacts):
                continue
            confidence=sum(abs(x["gap_mm"]) for x in contacts)
            best.append((confidence,a,b,contacts,chain))
        if not best:
            self.counts["one_or_unclosed_clamp_candidates"]+=1
            return {"a":bolt,"b":hosts[0]["instance"],"kind":"shaft_stack","hosts":hosts[:2],"qualified":False,
                    "fasteners":[bolt,*sorted(r["instance"] for r in rings)],
                    "reason":"端板或承压面配对不完整", "source":"measured B-Rep features"}
        _,a,b,contacts,chain=min(best,key=lambda t:t[0])
        backing=[x["backing"] for x in (a,b)]
        physical=all(abs(x["gap_mm"])<=TOUCH_MM and x["normal_angle_deg"]<=.001 for x in contacts)
        self.counts["two_host_stacks"]+=1
        return {"a":bolt,"b":a["instance"],"kind":"shaft_stack","hosts":[a,b],"qualified":True,
                "fasteners":[bolt,*[r["instance"] for r in chain[1:] if r["instance"] not in (a["instance"],b["instance"])]],
                "shank":{"face":shaft_face,"radius_mm":rad,"axis_points_mm":
                         np.round([origin,origin+length*axis],6).tolist(),"length_mm":round(length,6)},
                "links":contacts,"same_plane":physical,"other_host_pairs":len(best)-1,
                "backing_faces":backing,"source":"B-Rep inner wires, shank cylinder and trimmed bearing faces"}

    def run(self, base_edges):
        new=[]
        for i in self.g.ids:
            p=self.g.parts[i]
            a,b,c=p["shape_mm"]
            if not p["concentric_candidate"] or not (12<=c<=108 and 3<=a<=b<=42):
                continue
            shafts=self.g.joint_cylinders(i,False)
            if not shafts:
                continue
            first,second,r,face=max(shafts,key=lambda s:np.linalg.norm(s[1]-s[0]))
            length=float(np.linalg.norm(first-second))
            if length<max(2.7*r,.47*c) or a<1.8*r or r>15:
                continue
            self.counts["shaft_candidates"]+=1
            for u,v in ((first,second),(second,first)):
                record=self.clamp(i,u,(v-u)/length,length,r,face)
                if record:
                    new.append(record)
                    break
        geometry_rows=self.g.types
        for e in base_edges:
            a,b=e["a"],e["b"]
            if e["kind"]!="frame_end" or geometry_rows.get(a)!="SKIN" or e["gap_mm"]>TOUCH_MM or e["inner"]:
                continue
            contact=self.patch(a,b)
            new.append({"a":a,"b":b,"kind":"panel_interface","side":e["side"],"nearest_gap_mm":e["gap_mm"],
                        "witness_point_mm":e["points"][0],"patch":contact,"basis_edge":e["edge_id"],
                        "source":"B-Rep opposing planar trim projection"})
        for n,e in enumerate(new,1):
            e["edge_id"]=f"V{n:05d}"
        return new
