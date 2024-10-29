#! /usr/bin/env nextflow

process make_bed{
    cpus params.threads
    time '8h'
    
    input:
    file(genome)
    
    publishDir "${params.out}", mode: 'copy'
    
    output:
    file("*_lowMap_*.bed")
    
    script:
    genome_base = ( genome =~ /(.*)\.(fa|fasta|FA|FASTA)(\.gz)?/ )[0][1] 
    """
    echo $PATH
    echo ""
    if gzip -t ${genome} 2> /dev/null; then 
        gunzip -c ${genome} > unzipped.fa
        genmap index -F unzipped.fa -I genmap.idx
    else 
        genmap index -F ${genome} -I genmap.idx
    fi
    genmap map -T ${params.threads} -I genmap.idx -K ${params.bl_kmer} -bg -O genmap_out
    cat genmap_out.bedgraph | awk '{if (\$4 < 1){ print \$0; }}' | \
        bedtools merge -i stdin | \
        awk '{ if( \$3 - \$2 > ${params.bl_run} ){ print \$0; }}' \
        > ${genome_base}_lowMap_k${params.bl_kmer}_${params.bl_run}.bed
    """
}

workflow genmap_blacklist_fasta{
    take:
    genome
 
    main:
    blfile = make_bed(genome)
    
    emit:
    blfile
}
