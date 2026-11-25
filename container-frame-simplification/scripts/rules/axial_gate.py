"""Limit through-clamp promotion to beam-end attachments with aligned members."""

from math import sqrt

stack_predict=predict
FINAL_REVISION="frame-EU-1.2"
MIN_ALIGNMENT=.8660254037844386


def unit(v):
    size=sqrt(sum(x*x for x in v))
    return [x/size for x in v] if size>1.e-8 else [0,0,0]


def predict(packet,prompt_sha):
    assert packet["rules"]==FINAL_REVISION
    parts=packet["base"]["parts"]
    passed=[]
    audit={}
    for edge in packet["relations"]:
        if edge["kind"]!="shaft_stack":
            passed.append(edge)
            continue
        if "shank" in edge:
            a,b=edge["shank"]["axis_points_mm"]
            shaft=unit([y-x for x,y in zip(a,b)])
        else:
            # A partial match may lack the exported detailed witness.
            shaft=unit(parts[edge["a"]]["direction"])
        aligned=[]
        for host in edge["hosts"]:
            ident=host["backing"]["parts"][1]
            beam=unit(parts[ident]["direction"])
            aligned.append({"member":ident,"abs_dot":round(abs(sum(s*t for s,t in zip(shaft,beam))),8)})
        count=sum(x["abs_dot"]>=MIN_ALIGNMENT for x in aligned)
        if count:
            revised=dict(edge,qualified=edge["qualified"] and count==2)
            passed.append(revised)
        audit[edge["edge_id"]]={"sides":aligned,"used":count>0,
                               "sufficient_for_keep":edge["qualified"] and count==2}
    result=stack_predict(dict(packet,rules=REVISION,relations=passed),prompt_sha)
    result["rules"]=FINAL_REVISION
    result["interfaces"]=packet["relations"]
    result["beam_end_alignment"]=audit
    return result
