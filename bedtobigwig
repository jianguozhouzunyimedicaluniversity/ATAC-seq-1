rule bedtobigwig:
    input: 
        i1="bam/{sample}_sorted.bam",
        i2="macs2/{sample}_summits.bed"
    output:
        o1="macs2/{sample}_intersected.bam",
        o2="macs2/{sample}_intersected.bedgraph",
        o3="macs2/{sample}_sorted.bedgraph",
        o4="macs2/{sample}.bedgraph",
        o5="macs2/{sample}.bw"
    params: 
        err="macs2/{sample}_bedtobigwig.stderr", out="macs2/{sample}_bedtobigwig.stdout", name="bedtobigwig"
    message: 
        """--- bed to bigwig"""
    shell:
        "bedtools intersect -a {input.i1} -b {input.i2} > {output.o1}; bedtools genomecov -trackline -bg -ibam {output.o1} > {output.o2}; LC_COLLATE=C sort -k1,1 -k2,2n {output.o2} > {output.o3}; sed '$d' < {output.o3} > {output.o4}; bedGraphToBigWig {output.o4} /athena/khuranalab/scratch/kac2053/ref/chrom_sizes.txt {output.o5}"

