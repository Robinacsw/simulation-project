"""Command line entry for anonymous frame recognition, review and Parasolid release."""

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command",required=True)
    scan = commands.add_parser("scan",help="measure and run anonymous rules on a STEP assembly")
    scan.add_argument("--step",type=Path,required=True)
    scan.add_argument("--configuration",required=True)
    scan.add_argument("--profile",type=Path,required=True)
    scan.add_argument("--out",type=Path,required=True)
    scan.add_argument("--mode",choices=("G","G+C"),default="G+C")
    scan.add_argument("--expected-leaves",type=int)
    scan.add_argument("--datum",type=Path)

    replay = commands.add_parser("replay",help="rerun frozen measures with updated project scope or shape-only mode")
    replay.add_argument("--analysis",type=Path,required=True)
    replay.add_argument("--profile",type=Path,required=True)
    replay.add_argument("--out",type=Path,required=True)
    replay.add_argument("--mode",choices=("G","G+C"))

    review = commands.add_parser("review",help="render a grouped CAD review of selected pending roles")
    review.add_argument("--analysis",type=Path,required=True)
    review.add_argument("--rule",action="append",required=True)
    review.add_argument("--out",type=Path,required=True)

    release = commands.add_parser("deliver",help="apply approved uniquely mapped deletions and reopen X_T")
    release.add_argument("--analysis",type=Path,required=True)
    release.add_argument("--geometry",type=Path,required=True)
    release.add_argument("--approval",type=Path,required=True)
    release.add_argument("--protection",type=Path)
    release.add_argument("--spaceclaim",type=Path,required=True)
    release.add_argument("--name",required=True,help="new X_T basename inside --out")
    release.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    if args.command == "scan":
        from analysis_bundle import scan_assembly
        scan_assembly(args.step,args.configuration,args.profile,args.out,args.mode,args.expected_leaves,args.datum)
    elif args.command == "replay":
        from analysis_bundle import replay as rerun
        rerun(args.analysis,args.profile,args.out,args.mode)
    elif args.command == "review":
        from review_image import render_review
        render_review(args.analysis,args.rule,args.out)
    else:
        from delivery import deliver
        deliver(args.analysis,args.geometry,args.approval,args.protection,args.spaceclaim,args.out,args.name)


if __name__=="__main__":
    main()
