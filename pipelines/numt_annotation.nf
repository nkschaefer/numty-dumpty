#! /usr/bin/env nextflow

process get_mitochondrial_genes{
    input:
    file(mito_gtf)
    
    output:
    file("mito_genes.txt")
    
    script:
    """
    ${baseDir}/scripts/get_mito_genes.py --gtf ${mito_gtf} > mito_genes.txt       
    """
}

process remove_mitochondrial_genes_gtf{
    input:
    tuple file(mito_gtf), file(mito_genes)

    output:
    file("*_nonumt.gtf") 
    
    script:
    filebase=( mito_gtf =~ /(.*)\.gtf(\.gz)?/ )[0][1]
    """
    ${baseDir}/scripts/rm_mito_genes.py -M ${mito_genes} -m ${params.chrM} -i ${params.gene_id} \
        -n ${params.gene_name} --gtf ${mito_gtf} > ${filebase}_nonumt.gtf
    """
}

process remove_mitochondrial_genes_gff3{
    input:
    tuple file(mito_gff3), file(mito_genes)

    output:
    file("*_nonumt.gff3")

    script:
    filebase=( mito_gff3 =~ /(.*)\.gff3(\.gz)?/ )[0][1]
    """
    ${baseDir}/scripts/rm_mito_genes.py -M ${mito_genes} -m ${params.chrM} -i ${params.gene_id} \
        -n ${params.gene_name} --gff3 ${mito_gff3} > ${filebase}_nonumt.gff3
    """

}

process lift_mt_annotations{
    input:
    tuple file(ref_annotation), file(query_annotation), file(ref_genome), file(query_genome)
    
    publishDir "${params.out}", mode: 'copy'
    
    output:
    file("*_nonumt_MTfix.g*.gz")

    script:
    match=( query_annotation =~ /(.*)\.(gff3|gtf)(\.gz)?/ )[0]
    filebase = match[1]
    file_ext = match[2]
    """
    samtools faidx ${query_genome} ${params.chrM} > query_chrM.fa
    liftoff -g ${ref_annotation} query_chrM.fa ${ref_genome} > liftoff.txt
    cat ${query_annotation} liftoff.txt > ${filebase}_MTfix.${file_ext}
    gzip ${filebase}_MTfix.${file_ext}
    """
}

def get_annotation_type(annotation){
    return (annotation.toString() =~ /(.*)\.(gtf|gff3)(\.gz)?/ )[0][2]
}

workflow rm_numts_annotation{
    take:
    anno_and_genome
 
    main: 
    if (params.cat){
        params.gene_id = "source_gene"
        params.gene_name = "source_gene_common_name" 
    }
    
    mito_genes = Channel.fromPath(baseDir + "/data/gencode.v46.chrM.gtf.gz")\
        | get_mitochondrial_genes
        
    def anno_type = get_annotation_type(params.annotation)    
    lift_data = null
    if (anno_type == "gtf"){
        anno_adj = anno_and_genome.map{ x ->
            x[1]
        }.combine(mito_genes) | 
            remove_mitochondrial_genes_gtf    
        lift_data = Channel.fromPath(baseDir + "/data/gencode.v46.chrM.gtf.gz")\
            .combine(anno_adj)
    }
    else if (anno_type == "gff3"){
        anno_adj = anno_and_genome.map{ x ->
            x[1]
        }.combine(mito_genes) | 
            remove_mitochondrial_genes_gff3
        lift_data = Channel.fromPath(baseDir + "/data/gencode.v46.chrM.gff3.gz")\
            .combine(anno_adj)
    }
    else{
        error("Unrecognized annotation type. Allowed types are gtf and gff3.")
    }
    
    anno_lift = lift_data.combine(Channel.fromPath(baseDir + "/data/hg38_chrM.fa.gz"))\
        .combine(anno_and_genome.map{ x -> x[0] }) | lift_mt_annotations
    
    emit:
    anno_lift
}
