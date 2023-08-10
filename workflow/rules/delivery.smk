import shutil

def symlink_dir(source, dest, *exclude):
    dest.mkdir(parents=True, exist_ok=True)
    for fn in source.glob("**/*"):
        fn_dest = dest / Path(*fn.parts[3:])
        if any([re.search(v, str(fn)) for v in exclude]):
            continue
        if fn.is_dir():
            fn_dest.mkdir(parents=True, exist_ok=True)
        else:
            if not fn_dest.exists():
                fn_dest.symlink_to(fn.absolute())

def deliver_analysis_targets():
    # Only deliver to root if nothing else is there
    use_root = True
    need = sorted(analyses.keys())
    for analysis in analyses:
        p = Path(deliverdir) / analysis / "analysis.xlsx"
        if p.exists():
            need = [v for v in need if v != analysis]
            use_root = False
        
    # Exlcude analyses linked in root of delivery dir
    root = Path(deliverdir) / "analysis.xlsx"
    if root.exists():
        have = root.resolve().parts[-2]
        need = [v for v in need if v != have]
        use_root = False

    # Define input and output files
    ret = {"input": expand([
        deseqdir + "/{v}/analysis.xlsx",
        gseadir + "/{v}/GSEA-{s}.xlsx"],
        v=need, s=gsea_sets)
    }
    if len(need) == 1 and use_root:
        ret["output"] = [deliverdir + "/Analysis/analysis.xlsx"]
        ret["output"] += [deliverdir + "/Analysis/GSEA/GSEA-{s}.xlsx" for s in gsea_sets]
    else:
        ret["output"] = expand([
            deliverdir + "/Analysis/{v}/analysis.xlsx",
            deliverdir + "/Analysis/{v}/GSEA/GSEA-{s}.xlsx"],
            v=need, s=gsea_sets)
    return(ret)


rule deliver_methods:
    input:
        methods = "resources/Methods.docx"
    output:
        methods = deliverdir + "/Methods.docx",
    run:
        shutil.copy(input.methods, output.methods)

rule deliver_qc:
    input:
        qc = rules.multiqc.output
    output:
        qc = deliverdir + "/QC.html"
    run:
        Path(output.qc).symlink_to(Path(input.qc[0]).resolve())

rule deliver_bamfiles:
    input:
        rules.deliver_qc.output,
        bams = rules.run_star.input
    output:
        bams = [
            deliverdir + "/BAMFiles/" + Path(v).name
            for v in rules.run_star.input]
    run:
        for ibam, obam in zip(input.bams, output.bams):
            Path(obam).symlink_to(Path(ibam).resolve())

rule deliver_counts:
    input:
        rules.deliver_qc.output,
        rules.deliver_bamfiles.output,
        rules.deliver_methods.output,
        counts = rules.run_counts.input
    output:
        counts = deliverdir + "/counts.csv"
    run:
        Path(output.counts).symlink_to(Path(input.counts[0]).resolve())

# Include DESeq and GSEA results
deliver_analysis = deliver_analysis_targets()
rule deliver_analysis:
    input:
        rules.deliver_qc.output,
        rules.deliver_methods.output,
        rules.deliver_bamfiles.output,
        rules.deliver_counts.output,
        deseq = ancient(deliver_analysis['input'])
    output:
        deseq = deliver_analysis['output'],
    run:
        for source, dest in zip(input.deseq, output.deseq):
            sdir = Path(source).parent
            ddir = Path(dest).parent
            symlink_dir(sdir, ddir, "/counts.csv")

