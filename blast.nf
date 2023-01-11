nextflow.enable.dsl=2

params.query = "$HOME/sample.fa"
params.db = "$HOME/tools/blast-db/NR"
params.out = "out.blastn"
params.chunkSize = 200
params.format = "6"
params.blast = "blastp"
params.evalue = "0.0001"
// Validate inputs
db_name = file(params.db).name
db_path = file(params.db).parent

Channel
    .fromPath(params.query)
    .splitFasta(by: params.chunkSize)
    .set { fasta }

/*
 * Executes a exonerate job for each chunk emitted by the 'fasta' channel
 * and creates as output a channel named 'top_hits' emitting the resulting
 * BLAST matches
 */
process blast {
    input:
    file 'query.fa' 
    file db_path

    output:
    file 'blast.txt' 

    """
    ${params.blast} -query query.fa -db $db_path/$db_name -evalue ${params.evalue} -out blast.txt -outfmt "${params.format}" -num_threads $task.cpus
    """
}

workflow {
	blast (fasta, db_path).collectFile(name: "$params.out")
}
