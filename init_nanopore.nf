nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "miniBtk-${date}"
params.fofn_fastq = "samples.tsv"
params.ecoli = "/lustre/scratch123/tol/teams/blaxter/projects/tol-nemotodes/Otipu/data/refGenomes/E_coli.fa"
params.fofn_summaries = "samples_summaries.tsv"

ecoli = Channel.fromPath(params.ecoli, checkIfExists: true).collect()

fastqs = Channel.fromPath(params.fofn_fastq, checkIfExists: true)
    .splitCsv ( header: ['sample', 'fastq'], sep:'\t' )
    .map { row -> tuple(row.sample, (row.fastq =~ /.+_(\d+).fastq.gz/)[0][1], row.fastq) }

summaries = Channel.fromPath(params.fofn_summaries, checkIfExists: true)
    .splitCsv ( header: ['sample', 'summary'], sep:'\t' )
    .map { row -> tuple(row.sample,  row.summary) }



process remove_ecoli {
    tag "${strain}_${fileNum}"
    label 'btk'

    input:
        tuple val(strain), val(fileNum), val(reads)
        path(assembly)

    output:
        tuple val(strain), val(fileNum), path("${strain}_${fileNum}.fasta.gz")

    script:
        """
        iget $reads - | minimap2 -ax map-ont -t ${task.cpus} $assembly - \
            | samtools fasta -f 4 - \
            | gzip -c > ${strain}_${fileNum}.fasta.gz
        """
}

process concat_fasta {
    tag "$strain"
    publishDir "${params.outdir}/concatFasta",
        mode: 'copy'
    
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

process pycoQC {
    tag "$strain"
    publishDir "${params.outdir}/pycoQC",
        mode: 'copy'
    
    input:
    tuple val(strain), val(summary_file)
    
    output:
    tuple val(strain), path("${strain}.html")

    script:
        """
        iget $summary_file
        pycoQC --summary_file *gz --html_outfile ${strain}.html
        """
}

workflow {
    // QC
    pycoQC(
        summaries
    )

    // Filter reads
    remove_ecoli(
        fastqs,
        ecoli
    )
    .groupTuple(by: [0])
    .map { it ->  [ it[0], it[2].flatten() ] }
    .set { ch_cat_fasta }
    concat_fasta(
        ch_cat_fasta
    )
}