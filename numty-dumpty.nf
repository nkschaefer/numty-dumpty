#! /usr/bin/env nextflow
nextflow.enable.dsl=2

include { rm_numts_annotation } from './workflows/numt_annotation.nf'
include { mask_numts_fasta } from './workflows/numt_fasta.nf'
include { genmap_blacklist_fasta } from './workflows/mappability_fasta.nf'

if (!params.out){
    error("Output directory is required.")
}
if (!params.genome){
    error("FASTA-format genome is required.")
}
if (!params.chrM){
    error("Mitochondrial sequence name is required.")
} 

process mask_both_beds{
    time '1h' 
    
    input:
    tuple file(numtmasked_fa), file(bl_bed)
    
    publishDir "${params.out}", mode: 'copy'

    output:
    file("*_numtmask_mapmask.fa.gz")
    file("*_numtmask.fa.gz")
    
    script:
    fnbase = ( numtmasked_fa =~ /(.*)_numtmask.fa/ )[0][1]
    genome_out = fnbase + "_numtmask_mapmask.fa"
    """
    bedtools maskfasta -fi ${numtmasked_fa} -bed ${bl_bed} -fo ${genome_out}
    bgzip ${numtmasked_fa}
    bgzip ${genome_out}
    """
}

process faidx{
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

process rm_short_scaffolds_annotation{
    time '2h'
    
    input:
    tuple file(genome), file(fai), file(annotation)
        
    output:    
    tuple file("*filt.fa"), file("*filt.{gtf,gff3}")
    
    script:
    genome_base = ( genome =~ /(.*)\.(fa|fna|fasta|FA|FASTA)(\.gz)?/ )[0][1]
    anno_match = ( annotation =~ /(.*)\.(gtf|gff3)(\.gz)?/ )[0]
    annotation_base = anno_match[1]
    annotation_ext = anno_match[2]
    """
    if [ -e ${genome_base}.filt.fa ]; then
        rm ${genome_base}.filt.fa
    fi
    cat ${fai} | cut -f1-2 | \
        awk '{if (\$2 >= ${params.min_len}){ print \$0; }}' | cut -f1 | 
        egrep -v "${params.exclude}" > keep_chrs.txt
    cat keep_chrs.txt | while read seq; do
        samtools faidx ${genome} \${seq} >> ${genome_base}_filt.fa 
    done
    ${baseDir}/scripts/filter_annotation_seqs.py -g ${annotation} \
        -s keep_chrs.txt > ${annotation_base}_filt.${annotation_ext}
    """
}

process rm_short_scaffolds_fa{
    time '2h'
    
    input:
    tuple file(genome), file(fai)
    
    output:
    file("*filt.fa")
    
    script:
    genome_base = ( genome =~ /(.*)\.(fa|fasta|FA|FASTA)(\.gz)?/ )[0][1]
    """
    if [ -e ${genome_base}.filt.fa ]; then
        rm ${genome_base}.filt.fa
    fi
    cat ${fai} | cut -f1-2 | \
        awk '{if (\$2 >= ${params.min_len}){ print \$0; }}' | cut -f1 | 
        egrep -v "${params.exclude}" | while read seq; do
        samtools faidx ${genome} \${seq} >> ${genome_base}_filt.fa
    done
    """
}

process compress_fa{
    time '30m'

    input:
    file(genome)
    
    publishDir "${params.out}", mode: 'copy'
    
    output:
    file("${genome}.gz")
    
    script:
    """
    bgzip ${genome}
    """
}

process compress_annotation{
    time '30m'

    input:
    file(annotation)
    
    publishDir "${params.out}", mode: 'copy'
    
    output:
    file("${annotation}.gz")
    
    script:
    """
    bgzip ${annotation}
    """
}

workflow{
    
    // index input FASTA 
    fa_indexed = faidx(Channel.fromPath(params.genome))
 
    // Toss out short sequences and ones marked excludable (i.e. 
    // random, chrUn, alt scaffolds)
    if (params.annotation){
        
        fa_ann_filt = fa_indexed.combine(Channel.fromPath(params.annotation)) \
            | rm_short_scaffolds_annotation
        rm_numts_annotation(fa_ann_filt)
        fa_filt = fa_ann_filt.map{ x ->
            x[0]
        }
        fa_ann_filt.map{ x -> 
            x[1]
        } | compress_annotation
    }
    else{
        fa_filt = rm_short_scaffolds_fa(fa_indexed)
    }
    
    fa_filt | compress_fa

    numtmasked_fa = mask_numts_fasta(fa_filt)
    bl_bed = genmap_blacklist_fasta(fa_filt)
    numtmasked_fa.combine(bl_bed) | mask_both_beds
    
}

