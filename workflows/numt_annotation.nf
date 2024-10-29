#! /usr/bin/env nextflow

process get_mitochondrial_genes{
    time '15m'

    input:
    tuple file(mito_gtf), file(genome_gtf)
    
    output:
    file("mito_genes.txt")
    
    script:
    """
    # Get mitochondrial genes from GENCODE (human) annotation
    # Any genes with different (non-human Ensembl) IDs will be missed
    ${baseDir}/scripts/get_mito_genes.py --gtf ${mito_gtf} > mito_genes.txt 
    # Get mitochondrial genes from user-provided annotation
    # Genes in nuclear genome with matching names/IDs will be removed    
    ${baseDir}/scripts/get_mito_genes.py --gtf ${genome_gtf} -m ${params.chrM} \
        --id ${params.gene_id} --name ${params.gene_name} >> mito_genes.txt
    """
}

process remove_mitochondrial_genes_gtf{
    time '15m'

    input:
    tuple file(mito_gtf), file(mito_genes)

    output:
    file("*_nonumt.gtf") 
    
    script:
    filebase=( mito_gtf =~ /(.*)\.gtf(\.gz)?/ )[0][1]
    """
    ${baseDir}/scripts/rm_mito_genes.py -M ${mito_genes} -m ${params.chrM} -i ${params.gene_id} \
        -n ${params.gene_name} --gtf ${mito_gtf} --remove_all_MT > ${filebase}_nonumt.gtf
    """
}

process remove_mitochondrial_genes_nuclear_contigs_gtf{
    time '15m'

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
    time '15m'

    input:
    tuple file(mito_gff3), file(mito_genes)

    output:
    file("*_nonumt.gff3")

    script:
    filebase=( mito_gff3 =~ /(.*)\.gff3(\.gz)?/ )[0][1]
    """
    ${baseDir}/scripts/rm_mito_genes.py -M ${mito_genes} -m ${params.chrM} -i ${params.gene_id} \
        -n ${params.gene_name} --gff3 ${mito_gff3} --remove_all_MT > ${filebase}_nonumt.gff3
    """

}

process remove_mitochondrial_genes_nuclear_contigs_gff3{
    time '15m'

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
    time '30m'
    
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
    liftoff -g ${ref_annotation} query_chrM.fa ${ref_genome} |\
        grep -v "@" | grep -v "#" | grep -v aligning | grep -v lifting |\
        awk '{if ( NF >= 9){ print \$0; }}' > liftoff.txt
    cat ${query_annotation} liftoff.txt > ${filebase}_MTfix.${file_ext}
    gzip ${filebase}_MTfix.${file_ext}
    """
}

process copy_annotation{
    time '10m'
    
    input:
    file(input_annotation)
    
    publishDir "${params.out}", mode: 'copy'
    
    output:
    file("*_nonumt_MTfix.g*.gz")

    script:
    match=( input_annotation =~ /(.*)\.(gff3|gtf)(\.gz)?/ )[0]
    filebase = match[1]
    file_ext = match[2]
    """
    cp ${input_annotation} ${filebase}_MTfix.${file_ext}
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
        .combine(Channel.fromPath(params.annotation)) | get_mitochondrial_genes
        
    def anno_type = get_annotation_type(params.annotation)    
    lift_data = null
    anno_lift = null

    if (params.lift_human_MT){
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
    
    }
    else{
        if (anno_type == "gtf"){
            anno_adj = anno_and_genome.map{ x ->
                x[1]
            }.combine(mito_genes) | 
                remove_mitochondrial_genes_nuclear_contigs_gtf    
            lift_data = Channel.fromPath(baseDir + "/data/gencode.v46.chrM.gtf.gz")\
                .combine(anno_adj)
        }
        else if (anno_type == "gff3"){
            anno_adj = anno_and_genome.map{ x ->
                x[1]
            }.combine(mito_genes) | 
                remove_mitochondrial_genes_nuclear_contigs_gff3
            lift_data = Channel.fromPath(baseDir + "/data/gencode.v46.chrM.gff3.gz")\
                .combine(anno_adj)
        }
        else{
            error("Unrecognized annotation type. Allowed types are gtf and gff3.")
        }
        lift_data.view()
        anno_lift = lift_data.map{ x -> x[1] } | copy_annotation
    }
    
   
    emit:
    anno_lift
}
