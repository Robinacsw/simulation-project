# coding: utf-8
"""Copy retained CAD occurrences and reopen Parasolid in a separate job."""
import io
import json
import os
import traceback
from datetime import datetime

from SpaceClaim.Api.V252 import (Component, DesignBody, Document, ExportOptions,
                                IDesignBody, ImportOptions, ParasolidVersion, Part,
                                PartExportFormat, Task, WriteBlock)
from System import Enum

request_path = os.environ["FRAME_SIMPLIFICATION_REQUEST"]
work = os.path.dirname(request_path).replace("\\","/")
with io.open(request_path, encoding="utf-8") as stream:
    request = json.load(stream)
mode = request["mode"]
target = request.get("target")

def log(text):
    with io.open(work + "/spaceclaim.log", "a", encoding="utf-8") as f:
        f.write(unicode(datetime.now().isoformat()) + u" " + unicode(text) + u"\n")


def save(name, data):
    with io.open(work + "/" + name, "w", encoding="utf-8") as f:
        f.write(unicode(json.dumps(data, ensure_ascii=True)))


def name_token(value):
    return "".join(c if 32 <= ord(c) < 127 else "\\u%04x" % ord(c)
                   for c in unicode(value) if ord(c) != 0)


def caption(value):
    return u"".join(c for c in unicode(value)
                    if ord(c)>=32 and not 0xD800<=ord(c)<=0xDFFF).strip() or u"Component"


def xyz(point):
    return [float(getattr(point, axis))*1000 for axis in ("X","Y","Z")]


def metrics(obj, centers, areas):
    master = obj.Master
    shape = master.Shape
    transform = obj.TransformToMaster.Inverse
    extent = shape.GetBoundingBox(transform)
    tight = shape.GetBoundingBox(transform,True)
    volume = float(shape.Volume) * 1.0e9
    ident = master.Moniker.ToString()
    if ident not in areas:
        areas[ident] = sum(float(face.Area) for face in shape.Faces) * 1.0e6
    center = None
    if volume>0:
        if ident not in centers:
            centers[ident] = master.MassProperties.PrincipleAxes.Origin
        center = xyz(transform * centers[ident])
    parent= obj.Parent.Parent
    trace=name_token(parent.Name) if parent is not None and name_token(parent.Name).startswith("Q") else ""
    return {"name":name_token(obj.Name),"trace":trace,"box":xyz(extent.MinCorner)+xyz(extent.MaxCorner),
            "tight_box":xyz(tight.MinCorner)+xyz(tight.MaxCorner),
            "volume":volume,"faces":int(shape.Faces.Count),"center":center,"area":areas[ident]}


def inventory(doc):
    objs = list(doc.MainPart.GetDescendants[IDesignBody]())
    centers = {}
    areas = {}
    rows = []
    for n,obj in enumerate(objs,1):
        rows.append(metrics(obj, centers, areas))
        if n%2500==0:
            log("Measured %d body occurrences" % n)
    return objs,{"rows":rows,"count":len(rows),"solids":sum(r["volume"]>0 for r in rows),
                 "top_components":int(doc.MainPart.Components.Count)}


def assert_sample(actual, expected, index):
    assert actual["name"]==expected["name"], ("body name",index,actual["name"],expected["name"])
    assert actual["faces"]==expected["faces"], ("topology",index)
    assert all(abs(a-b)<=.005 for a,b in zip(actual["box"], expected["box"])), ("position",index)
    assert abs(actual["volume"]-expected["volume"])<=max(.003,abs(expected["volume"])*.00001), ("volume",index)
    assert (actual["center"] is None) == (expected["center"] is None), ("center type",index)
    if actual["center"] is not None:
        assert all(abs(a-b)<=.005 for a,b in zip(actual["center"], expected["center"])), ("center",index)


def source_lineage(obj):
    lineage = []
    part = obj.Parent
    while part.Parent is not None:
        occurrence = part.Parent
        title = occurrence.Name or getattr(part.Master,"Name",u"Component")
        lineage.append((occurrence.Moniker.ToString(), caption(title)))
        part = occurrence.Parent
        if len(lineage)>128:
            raise ValueError("Unexpected assembly depth")
    return list(reversed(lineage))


def assemble(doc, retained):
    groups = {():doc.MainPart}
    prototypes = {}
    placed = []
    for index,source in retained:
        path = ()
        for key,title in source_lineage(source):
            parent = groups[path]
            path = path + (key,)
            if path not in groups:
                child = Part.Create(doc,title)
                Component.Create(parent,child)
                groups[path] = child
        tag=u"Q%05d %s" % (index,caption(source.Name))
        ident=source.Master.Moniker.ToString()
        prototype=prototypes.get(ident)
        if prototype is None:
            template_tag=u"D%04d %s" % (len(prototypes)+1,caption(source.Name))
            prototype=Part.Create(doc,template_tag)
            DesignBody.Create(prototype,caption(source.Name),source.Master.Shape.Copy())
            prototypes[ident]=prototype
        result=Component.Create(groups[path],prototype)
        result.Placement=source.TransformToMaster.Inverse
        result.Name=tag
        if result is None:
            raise ValueError("Could not copy body %d" % index)
        placed.append(result)
        if len(placed)%2000==0:
            log("Copied %d retained bodies" % len(placed))
    return len(groups)-1, len(prototypes)


def clean():
    plan=request["plan"]
    log("Opening frozen CAD input")
    old=Document.Open(request["source"],ImportOptions.Create())
    objs,before=inventory(old)
    assert before["count"]==plan["source_bodies"] and before["solids"]==plan["source_solids"]
    remove=set(plan["removed_indices"])
    protected=set(plan["protected_body_indices"])
    assert not remove & protected and all(before["rows"][i-1]["volume"] > 0 for i in remove)
    for key,val in plan["before_samples"].items():
        assert_sample(before["rows"][int(key)-1],val,int(key))
    save("inventory.json",before)
    retained=[(i,obj) for i,obj in enumerate(objs,1) if i not in remove]
    result=Document.Create()
    group_count = [None]
    def copy_all():
        group_count[0] = assemble(result,retained)
    if WriteBlock.IsActive:
        copy_all()
    else:
        WriteBlock.ExecuteTask("Copy approved frame geometry", Task(copy_all))
    count=len(list(result.MainPart.GetDescendants[IDesignBody]()))
    assert count==len(retained),("Unexpected retained count",count,len(retained))
    opts=ExportOptions.Create()
    opts.ExportNames=True
    version_names=list(Enum.GetNames(ParasolidVersion))
    version="V33" if "V33" in version_names else "V27"
    opts.Parasolid.Version=getattr(ParasolidVersion,version)
    log("Exporting %d retained bodies in Parasolid %s" % (count,version))
    log("Shared body prototypes: %d" % group_count[0][1])
    assert not os.path.exists(target)
    result.MainPart.Export(PartExportFormat.ParasolidText,target,True,opts)
    assert os.path.getsize(target)>0
    save("export_receipt.json",{"label_count":count,"bytes":os.path.getsize(target),
         "version":version,"copied_component_occurrences":group_count[0][0],
         "shared_body_prototypes":group_count[0][1],
         "top_components":int(result.MainPart.Components.Count)})


def check():
    assert os.path.isfile(request["source"])
    log("Opening exported candidate in an independent SpaceClaim process")
    model=Document.Open(request["source"],ImportOptions.Create())
    _,after=inventory(model)
    save("inventory.json",after)
    log("Opened CAD input with %d bodies and %d solids" % (after["count"],after["solids"]))


try:
    if mode=="write":
        clean()
    else:
        assert mode=="read"
        check()
    log("COMPLETE")
except:
    log(traceback.format_exc())
    raise
