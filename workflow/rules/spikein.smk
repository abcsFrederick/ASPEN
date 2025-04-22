rule align2spikein:
    input:
        R1=rules.remove_BL.output.R1,
        R2=rules.remove_BL.output.R2,
    output:
        counts=join(RESULTSDIR,"spikein","{replicate}","{replicate}.counts"),
    params:
        replicate="{replicate}",
        workdir=RESULTSDIR,
        spikein=str(SPIKEIN).lower(),
        spikein_indexdir=SPIKEINDEXDIR,
        spikein_genome=SPIKEINGENOME,
        mem=getmemG("align2spikein")
    container: config["bwadocker"]
    threads: getthreads("align2spikein")
    shell:"""
set -exo pipefail
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then cd /lscratch/${{SLURM_JOB_ID}};else cd /dev/shm;fi
unset PYTHONPATH

echo -ne "{params.replicate}\\t" > {output.counts}
if [[ "{params.spikein}" == "true" ]];then
bwa-mem2 mem -t {threads} {params.spikein_indexdir}/{params.spikein_genome} {input.R1} {input.R2} | \
    samtools view -f 2 -q 30 | wc -l | awk '{{printf "%d\\n", $1/2}}' >> {output.counts}
else
echo -ne "0\\n" >> {output.counts}
fi
"""

localrules: compute_scaling_factors
rule compute_scaling_factors:
    input:
        counts=expand(join(RESULTSDIR,"spikein","{replicate}","{replicate}.counts"),replicate=REPLICATES),
    output:
        scaling_factors=join(RESULTSDIR,"spikein","scaling_factors.tsv"),
    params:
        replicates=REPLICATES,
        workdir=RESULTSDIR,
        scriptsdir=SCRIPTSDIR,
        spikein=str(SPIKEIN).lower(),
        script="_compute_downsampling_scaling_factors.py",
    container: config["masterdocker"]
    shell:"""
set -exo pipefail
TMPDIR="/lscratch/$SLURM_JOB_ID"
if [ ! -d $TMPDIR ];then
    TMPDIR="/dev/shm"
fi
unset PYTHONPATH
if [[ "{params.spikein}" == "true" ]];then
cat {input.counts} | python {params.scriptsdir}/{params.script} | sort > {output.scaling_factors}
else
cat {input.counts} | sort > ${{TMPDIR}}/allcounts.txt
while read replicateName count;do
    echo -ne "$replicateName\\t$count\\t1.0\\n"
done < ${{TMPDIR}}/allcounts.txt > {output.scaling_factors}
fi
"""

localrules: create_normalized_macs2_bigwigs
rule create_normalized_macs2_bigwigs:
    input:
        scaling_factors = join(RESULTSDIR,"spikein","scaling_factors.tsv"),
        bigwigFileLists = expand(join(RESULTSDIR,"peaks","macs2","{sample}.macs2.bigwigfiles"),sample=SAMPLES),
    output:
        expand(join(RESULTSDIR,"spikein","bigwig","{replicate}.macs2.norm.bw"),replicate=REPLICATES),
        # join(RESULTSDIR,"spikein","tmp")
    params:
        scriptsdir=SCRIPTSDIR,
        genomefile=GENOMEFILE,
        outdir=join(RESULTSDIR,"spikein","bigwig"),
        script1="_print_replicate_scaling_factor.py",
    container: config["ucscdocker"]
    shell:"""
set -exo pipefail
TMPDIR="/lscratch/$SLURM_JOB_ID"
if [ ! -d $TMPDIR ];then
    TMPDIR="/dev/shm"
fi
unset PYTHONPATH
cat {input.bigwigFileLists} > ${{TMPDIR}}/allbigwigs.txt
while read replicateName sampleName bigwig;do
    sf=$(cat {input.scaling_factors} | python {params.scriptsdir}/{params.script1} $replicateName)
    echo -ne "$replicateName\\t$sampleName\\t$bigwig\\t$sf\\n"
    bigWigToBedGraph $bigwig ${{TMPDIR}}/${{replicateName}}.bedGraph
    awk -v sf=$sf '{{print $1, $2, $3, $4 * sf}}' ${{TMPDIR}}/${{replicateName}}.bedGraph > ${{TMPDIR}}/${{replicateName}}.norm.bedGraph
    bedGraphToBigWig ${{TMPDIR}}/${{replicateName}}.norm.bedGraph {params.genomefile} {params.outdir}/${{replicateName}}.macs2.norm.bw
done < ${{TMPDIR}}/allbigwigs.txt > {output}
"""
