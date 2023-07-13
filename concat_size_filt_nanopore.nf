nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "concat_ont-${date}"
params.fofn_fastq = "samples.tsv"
params.min_seq_size = 4000


fastqs = Channel.fromPath(params.fofn_fastq, checkIfExists: true)
    .splitCsv ( header: ['sample', 'fastq'], sep:'\t' )
    .map { row -> tuple(row.sample, (row.fastq =~ /.+_(\d+).fastq.gz/)[0][1], row.fastq) }


process size_filter {
    tag "${strain}_${fileNum}"
    //label 'btk'
    module 'minimap2/2.24--h7132678_1'
    input:
        tuple val(strain), val(fileNum), val(reads)

    output:
        tuple val(strain), val(fileNum), path("${strain}_${fileNum}.fasta.gz")

    script:
        """
        
        #iget $reads - | seqkit seq --min-len ${params.min_seq_size} \
        #    | gzip -c > ${strain}_${fileNum}.fasta.gz
	iget $reads ${strain}_${fileNum}.fasta.gz
        """
}

process concat_fasta {
    tag "$strain"
    publishDir "${params.outdir}/concatFasta",
        mode: 'move'
    
    input:
        tuple val(strain), path(reads)
    
    output:
        tuple val(strain), path("*.merged.fasta.gz")

    script:
        readList = reads.collect{it.toString()}
        """
        cat ${readList.sort().join(' ')} > ${strain}.merged.fasta.gz
        """
}


workflow {
    size_filter(
        fastqs
    )
    .groupTuple(by: [0])
    .map { it ->  [ it[0], it[2].flatten() ] }
    .set { ch_cat_fasta }
    concat_fasta(
        ch_cat_fasta
    )
}
