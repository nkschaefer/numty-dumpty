#! /usr/bin/env nextflow

process faidx_genome{
    time '30m'
    
    input:
    file(genome)
    
    output:
    tuple file(genome), file("${genome}.fai")
    
    script:
    """
    samtools faidx ${genome}
    """
}

def genome_basename(genome){
    def match = (genome =~ /(.*)\.(fa|fasta|FA|FASTA)(\.gz)?/)[0]
    return match[1]
}

process remove_MT_genome{
    time '30m'
    
    input:
    tuple file(genome), file(faidx), val(mito_name)
    
    output:
    file("*_noMT.fa.gz")

    script:
    excludeChrs = "^" + mito_name + "\$"
    if (params.exclude.length() > 0){
        excludeChrs = "^" + mito_name + "\$|" + params.exclude
    }
    genome_out = genome_basename(genome) + "_noMT.fa"
    """
    cat ${faidx} | cut -f1 | egrep -v "${excludeChrs}" | while read seqname; do
        samtools faidx ${genome} \${seqname} >> ${genome_out}
    done
    bgzip ${genome_out}
    """
}

process get_MT_genome{
    time '30m'
    
    input:
    tuple file(genome), file(faidx), val(mito_name)
    
    output:
    file("*_MT.fa")
    
    script:
    genome_out = genome_basename(genome) + "_MT.fa"
    """
    samtools faidx ${genome} ${mito_name} > ${genome_out}
    """
}

process bowtie2_index{
    cpus params.threads
    time '4h'
    
    input:
    file(genome_noMT)
    
    output:
    tuple file(genome_noMT), file("*.bt2")

    script:
    genome_out = genome_basename(genome_noMT)
    """
    bowtie2-build --threads ${params.threads} ${genome_noMT} ${genome_out}
    """
}

process simulate_reads{
    time '2h'
    
    input:
    file(genome_MT)
    
    output:
    file("*_mitosim.fq.gz")

    script:
    outbase = genome_basename(genome_MT)
    """
    art_illumina -na -ss NS50 -c ${params.nreads} -i ${genome_MT} -l 20 -o ${outbase}_mitosim
    gzip ${outbase}_mitosim.fq
    """
}

process map_reads{
    cpus params.threads
    time '4h'
    
    input: 
    tuple file(genome_noMT), file(genome_noMT_idx), file(MT_reads)
    
    output:
    tuple file("*_MTaln.bam"), file("*_MTaln.bam.bai")
     
    script:
    idxbase = ( genome_noMT_idx[0] =~ /(.*)\.(rev\.)?(\d+)\.bt2/ )[0][1]
    """
    # Map to reference and sort 
    bowtie2 -p ${params.threads} -x ${idxbase} -U ${MT_reads} | \
        samtools view -bh - | samtools sort -o ${idxbase}_MTaln.bam
    # Index BAM file
    samtools index ${idxbase}_MTaln.bam
    """
    
}

process call_peaks{
    time '2h'
    
    input:
    tuple file(genome), file(genome_idx)
    
    publishDir "${params.out}", mode: 'copy'

    output:
    file("*_numt.bed")
     
    script:
    genome_base = ( genome =~ /(.*)_noMT_MTaln\.bam/ )[0][1]
    """
    macs2 callpeak -t ${genome} --nomodel --nolambda \
        --keep-dup all -g ${params.macs2_effsize} -n ${genome_base} 
    mv ${genome_base}_peaks.narrowPeak ${genome_base}_numt.bed
    """
}

process mask_genome{
    time '30m'
    
    input:
    tuple file(genome), file(bed)
    
    output:
    file("*_numtmask.fa")

    script:
    genome_out = genome_basename(genome) + "_numtmask.fa"
    """
    if gzip -t ${genome} 2> /dev/null; then 
        gunzip -c ${genome} > unzipped.fa
        bedtools maskfasta -fi unzipped.fa -bed ${bed} -fo ${genome_out}
    else 
        bedtools maskfasta -fi ${genome} -bed ${bed} -fo ${genome_out}
    fi
    """
}

workflow mask_numts_fasta{
    take:
    genome
     
    main: 
    genome_indexed = faidx_genome(genome)
    genome_indexed_plusMT = genome_indexed.map{ genome, fai ->
        [genome, fai, params.chrM ]
    } 
    noMT_genome_indexed = remove_MT_genome(genome_indexed_plusMT) | bowtie2_index
    
    MT_genome = get_MT_genome(genome_indexed_plusMT)
    MT_reads = simulate_reads(MT_genome)
    
    peaksfile = noMT_genome_indexed.combine(MT_reads) | map_reads | call_peaks
    
    fa_masked = genome.combine(peaksfile) | mask_genome    
    
    emit:
    fa_masked
}

