"""Collect geometry-only relation witnesses for G+C. No engineering decisions here."""

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations

import numpy as np
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.GeomAbs import GeomAbs_Plane
from OCP.TopAbs import TopAbs_FACE
from OCP.TopoDS import TopoDS

from intake import indexed, placed, properties
from measure_geometry import gap
from recognize import PARAMETERS as P, casting_candidate, coords, geometry, upright_candidate


def intervals(lo, hi, a, b):
    return np.maximum(np.maximum(lo-b,a-hi),0)


class Graph:
    def __init__(self, model, packet):
        self.model = model
        self.packet = packet
        self.parts = packet["parts"]
        self.frame = packet["frame_box"]
        self.ids = sorted(i for i,p in self.parts.items() if "box" in p and p["valid"] and p["solids"]==1)
        self.look = {id:i for i,id in enumerate(self.ids)}
        boxes=np.array([self.parts[i]["box"] for i in self.ids])
        self.lo,self.hi = boxes[:,0],boxes[:,1]
        self.size = self.hi-self.lo
        self.center = (self.lo+self.hi)/2
        self.types={i:geometry(self.parts[i],self.frame)["rule"] for i in self.ids}
        self.span = np.ptp(self.frame,axis=0)
        self.counters=Counter()
        self.edges = []

    @lru_cache(None)
    def ancestry(self, ident):
        chain=[]
        p=self.parts[ident]["parent"]
        while p:
            chain.append(p)
            p=self.packet["tree"][p]
        return chain[::-1]

    def common_depth(self,a,b):
        return len(set(self.ancestry(a)) & set(self.ancestry(b)))

    @lru_cache(maxsize=2200)
    def planes(self, ident):
        s=self.model.local_shape(ident)
        planes=[]
        for i,raw in enumerate(indexed(s,TopAbs_FACE),1):
            face=TopoDS.Face_s(raw)
            adaptor=BRepAdaptor_Surface(face,True)
            if adaptor.GetType()!=GeomAbs_Plane:
                continue
            area,center=properties(face)
            if area>.05:
                norm=adaptor.Plane().Axis().Direction().Coord()
                planes.append((i,face,area,np.array(center),np.array(norm)))
        return planes

    def projection(self,a,b,side):
        ia,ib=self.look[a],self.look[b]
        end=(self.hi if side else self.lo)[ia,2]
        sign = 1 if side else -1
        parts=[p for p in self.planes(a) if abs(p[4][2])>.999999 and abs(p[3][2]-end)<.025]
        mates=[p for p in self.planes(b) if abs(p[4][2])>.999999 and
               -.04<=sign*(p[3][2]-end)<=P["port_search_mm"]]
        evidence=[]
        for ai,face,area,center,normal in parts:
            for bi,other,barea,bcenter,bnormal in mates:
                dz=bcenter[2]-center[2]
                target=placed(face,np.eye(3),np.array([0.,0.,dz]))
                measurement=BRepAlgoAPI_Common(target,other)
                if not measurement.IsDone() or measurement.Shape().IsNull():
                    continue
                face_area,xyz=properties(measurement.Shape())
                if face_area <= .02 or face_area > min(area,barea)+.05:
                    continue
                evidence.append({"face_on_member":ai,"face_on_casting":bi,"axial_gap_mm":round(abs(dz),6),
                                 "overlap_mm2":round(face_area,6),"witness_xyz_mm":np.round(xyz,6).tolist(),
                                 "member_face_center":np.round(center,6).tolist(),
                                 "cast_face_center":np.round(bcenter,6).tolist(),
                                 "face_normal":np.round(normal,6).tolist()})
        good=round(sum(x["overlap_mm2"] for x in evidence),5)
        distance=round(max((e["axial_gap_mm"] for e in evidence),default=0),6)
        return {"a":a,"b":b,"kind":"port","side":side,"gap_mm":distance,"area_mm2":good,
                "target_center":np.round(self.center[ib],5).tolist(),"qualified":good>=P["min_end_projection_mm2"],
                "faces":evidence,"source":"B-Rep terminal planar intersection"}

    def gather_ports(self):
        fittings=np.array([self.look[i] for i in self.ids if casting_candidate(self.parts[i],self.frame)])
        for a in self.ids:
            if not len(fittings) or not upright_candidate(self.parts[a],self.frame):
                continue
            ia=self.look[a]
            hor=intervals(self.lo[fittings,:2],self.hi[fittings,:2],self.lo[ia,:2],self.hi[ia,:2])
            for side in (0,1):
                limit=(self.hi if side else self.lo)[ia,2]
                dz=(self.lo if side else self.hi)[fittings,2]-limit
                allowed = (np.max(hor,axis=1)<5) & (np.abs(dz)<=P["port_search_mm"])
                allowed &= ((self.center[fittings,2]>self.center[ia,2])==bool(side))
                options=fittings[allowed]
                order=sorted(options,key=lambda i: abs(limit-self.center[i,2]))
                for ib in order[:8]:
                    try:
                        self.edges.append(self.projection(a,self.ids[ib],side))
                    except (ValueError,RuntimeError) as exc:
                        self.counters[f"projection:{type(exc).__name__}"]+=1
        self.counters["casting_candidates"]=len(fittings)

    @lru_cache(maxsize=600)
    def cylindrical(self, ident, inner, outer=False):
        row=self.model.by_path[self.model.paths[ident]]
        data=self.model.features[row["definition"]]
        ori,off=self.model.pose(ident)
        orientation=self.model.axes@ori
        offset=self.model.axes@(off-self.model.origin)
        picked=[]
        for i,f in enumerate(data["cylinders"],1):
            if f["inner"]!=inner or f["angle"]<2.5 or not 6<=f["r"]<=47:
                continue
            a,b = (orientation@np.asarray(f[k])+offset for k in ("a","b"))
            axis=(b-a)
            length=float(np.linalg.norm(axis))
            if length <= (P["long_rod_length_mm"]*.7 if outer else 5):
                continue
            if abs(axis[2])/length < .93:
                continue
            picked.append((a,b,f["r"],i))
        return picked

    def panel_candidates(self):
        panels=defaultdict(set)
        for i in self.ids:
            row=self.parts[i]
            a,b,c=row["shape_mm"]
            local=self.size[self.look[i]]
            # The folded blank can have flanges. Organic surface/volume thinness
            # prevents a nearby solid block from standing in for a door leaf.
            if local[2]<.38*self.span[2] or not 170<=b<=.93*self.span[1] or a>110:
                continue
            if 2*row["volume"]/max(1,row["area"])>11:
                continue
            for parent in self.ancestry(i)[2:]:
                panels[parent].add(i)
        return panels

    def find_panel(self,rod,fixture,panels):
        ancestors=set(self.ancestry(fixture)[2:])
        i=self.look[rod]
        for parent in reversed(self.ancestry(rod)[2:]):
            if parent not in ancestors:
                continue
            for p in sorted(panels[parent]):
                j=self.look[p]
                if p==rod or p==fixture:
                    continue
                d=intervals(self.lo[i],self.hi[i],self.lo[j],self.hi[j])
                axial=min(self.hi[i,2],self.hi[j,2])-max(self.lo[i,2],self.lo[j,2])
                if np.linalg.norm(d)<=65 and axial>=.52*self.size[i,2]:
                    return p
        return None

    def guide_matches(self):
        panels=self.panel_candidates()
        sizes=np.array([self.parts[i]["shape_mm"][2] for i in self.ids])
        fixtures=np.flatnonzero((sizes>=22)&(sizes<=170)&(self.size[:,2]>=12))
        for rod in self.ids:
            p=self.parts[rod]
            if p["shaft"]<7 or p["shape_mm"][1]>75 or p["shape_mm"][2]<P["long_rod_length_mm"]:
                continue
            a=self.look[rod]
            rods=self.cylindrical(rod,False,True)
            if not rods:
                continue
            boxes=intervals(self.lo[fixtures],self.hi[fixtures],self.lo[a],self.hi[a])
            for idx in fixtures[np.max(boxes,axis=1)<=7]:
                fixture=self.ids[idx]
                if self.common_depth(rod,fixture)<3:
                    continue
                panel=self.find_panel(rod,fixture,panels)
                if panel is None:
                    continue
                found=[]
                for u,v,rad,fi in rods:
                    shaft=v-u
                    length=np.linalg.norm(shaft)
                    direction=shaft/length
                    for q,s,br,gi in self.cylindrical(fixture,True):
                        radial=br-rad
                        if radial<-.09 or radial>5.0:
                            continue
                        offset=[w-u - ((w-u)@direction)*direction for w in (q,s)]
                        off=max(map(np.linalg.norm,offset))
                        along=sorted(float((w-u)@direction) for w in (q,s))
                        overlap=min(length,along[1])-max(0,along[0])
                        if off<=.35 and overlap>=max(5,.30*np.linalg.norm(s-q)):
                            xyz=u+direction*(max(0,along[0])+min(length,along[1]))/2
                            found.append((overlap,off,radial,xyz,rad,br,fi,gi,q,s))
                if found:
                    lap,off,radial,xyz,r,s,fi,gi,q,v=max(found,key=lambda t:(t[0],-abs(t[2]),-t[1]))
                    self.edges.append({"a":rod,"b":fixture,"kind":"axis_guide","axis_height":round(xyz[2],6),
                                       "source":"B-Rep cylinder axes and neutral assembly ancestry",
                                       "panel":panel,"radial_mm":round(radial,5),"axis_offset_mm":round(off,5),
                                       "axis_lap_mm":round(lap,5),"rod_radius_mm":round(r,5),"guide_radius_mm":round(s,5),
                                       "points":np.round([xyz,q,v],5).tolist(),"cylinder_indices":[fi,gi]})

    @lru_cache(maxsize=6000)
    def nearest(self,a,b):
        if self.parts[a]["faces"]*self.parts[b]["faces"]>70000:
            self.counters["distance_budget"]+=1
            return None
        self.counters["distance_checks"]+=1
        try:
            measured=gap(self.model.local_shape(a),self.model.local_shape(b))
        except (ValueError,RuntimeError,AssertionError) as exc:
            self.counters[f"distance:{type(exc).__name__}"]+=1
            return None
        return {"gap_mm":round(measured["distance_mm"],6),
                "points":[np.round(measured[f"point_{k}_mm"],5).tolist() for k in (1,2)],
                "inner":bool(measured["inner_solution"])}

    def frame_ends(self):
        members = np.array([self.look[i] for i in self.ids if self.types[i] in
                            {"TOP_RAIL","BOTTOM_RAIL","CROSS","BOTTOM_LONG","POST","CAST"} and
                            not (self.types[i]=="POST" and self.parts[i]["shaft"]>=7)])
        if len(members)==0:
            return
        for ident in self.ids:
            p=self.parts[ident]
            ia=self.look[ident]
            _,box,_,(_,_,z)=coords(p,self.frame)
            a,b,c = p["shape_mm"]
            thin=(self.types[ident]=="SKIN" and z>=.76 and box[1]>=.43*self.span[1] and box[0]<.23*self.span[1])
            segment=(self.types[ident]=="UNCLEAR" and 19<=a and c>=3.35*max(b,1) and 175<=c and
                     (z<.27 or z>.69))
            if not thin and not segment:
                continue
            kind="TOP_RAIL" if thin else None
            axis=1 if thin else int(np.argmax(self.size[ia]))
            options=members[(np.max(intervals(self.lo[members],self.hi[members],self.lo[ia],self.hi[ia]),axis=1)
                             <=P["frame_neighbor_search_mm"]) & (members!=ia)]
            for side in (0,1):
                end=(self.hi if side else self.lo)[ia,axis]
                take=options[(self.lo[options,axis]<=end+4)&(self.hi[options,axis]>=end-4)]
                take=[ib for ib in take if not kind or self.types[self.ids[ib]]==kind]
                take=sorted(take,key=lambda ib: np.linalg.norm(self.center[ib]-self.center[ia]))
                for ib in take[:5]:
                    b=self.ids[ib]
                    witness=self.nearest(ident,b)
                    if witness is None or witness["gap_mm"]>P["frame_neighbor_search_mm"]:
                        continue
                    if abs(witness["points"][0][axis]-end)>max(8,.07*self.size[ia,axis]):
                        continue
                    self.edges.append(dict(witness,a=ident,b=b,side=side,kind="frame_end",
                                           source="B-Rep shortest distance at measured end"))

    @lru_cache(maxsize=1400)
    def joint_cylinders(self, ident, inner):
        row=self.model.by_path[self.model.paths[ident]]
        props=self.model.features[row["definition"]]
        ori,off=self.model.pose(ident)
        orientation=self.model.axes@ori
        offset=self.model.axes@(off-self.model.origin)
        return [(orientation@np.asarray(c["a"])+offset,orientation@np.asarray(c["b"])+offset,c["r"],j)
                for j,c in enumerate(props["cylinders"],1)
                if c["inner"]==inner and 2.2<=c["r"]<=15.5 and c["angle"]>=2.65
                and np.linalg.norm(np.asarray(c["b"])-c["a"])>1.3]

    def beam_clamps(self):
        candidates=[]
        bolts=[]
        for ident in self.ids:
            p=self.parts[ident]
            a,b,c=p["shape_mm"]
            _, _, _, (_,_,z)=coords(p,self.frame)
            if not p["concentric_candidate"] or z<.73:
                continue
            if 8<=c<=82 and a>=3.7 and b<=38:
                bolts.append(ident)
            elif c>=68 and b>=21:
                candidates.append(self.look[ident])
        indices=np.array(candidates,dtype=int)
        if not len(indices):
            return
        for ident in bolts:
            i=self.look[ident]
            shaft=self.joint_cylinders(ident,False)
            if not shaft:
                continue
            options=indices[np.max(intervals(self.lo[indices],self.hi[indices],self.lo[i],self.hi[i]),axis=1)<=1.2]
            found=[]
            for j in options:
                neighbor=self.ids[j]
                p=self.parts[neighbor]
                a,b,c=p["shape_mm"]
                horizontal=abs(p["direction"][2])<=.25 and c>=190 and c>=3.1*max(1,b)
                if self.parts[neighbor]["faces"]>560:
                    continue
                for u,v,r,k in shaft:
                    length=np.linalg.norm(v-u)
                    if length<5 or length<.4*self.parts[ident]["shape_mm"][2]:
                        continue
                    direction=(v-u)/length
                    for q,s,rad,f in self.joint_cylinders(neighbor,True):
                        radial=rad-r
                        if not -.075<=radial<=2.3:
                            continue
                        axial=sorted(float((w-u)@direction) for w in (q,s))
                        residual=[w-u-((w-u)@direction)*direction for w in (q,s)]
                        off=max(map(np.linalg.norm,residual))
                        overlap=min(length,axial[1])-max(0,axial[0])
                        if off>.16 or overlap<1.3:
                            continue
                        mid=u+direction*(min(length,axial[1])+max(0,axial[0]))*.5
                        axis=int(np.argmax(np.abs(p["direction"])))
                        end=min(abs(mid[axis]-self.lo[j,axis]),abs(mid[axis]-self.hi[j,axis]))
                        long_end=bool(horizontal and end<=min(165,.23*c))
                        found.append({"a":ident,"b":neighbor,"kind":"beam_clamp","long_end":long_end,
                                      "source":"B-Rep coaxial terminal hole witness","points":np.round([mid,q,s],5).tolist(),
                                      "shaft_cylinder":k,"host_cylinder":f,"center_offset":round(off,5),
                                      "radial_clearance":round(radial,5),"axial_overlap":round(overlap,5)})
            if any(f["long_end"] for f in found):
                best={}
                for entry in sorted(found,key=lambda e:e["axial_overlap"],reverse=True):
                    best.setdefault((entry["b"],entry["shaft_cylinder"]),entry)
                self.edges.extend(best.values())

    def run(self):
        self.gather_ports()
        self.guide_matches()
        self.frame_ends()
        self.beam_clamps()
        for n,e in enumerate(self.edges,1):
            e["edge_id"]=f"E{n:05d}"
        counts=Counter(e["kind"] for e in self.edges)
        print("relations",counts,"measurement",self.counters,flush=True)
        return {"edges":self.edges,"counts":counts,"measurements":self.counters}
