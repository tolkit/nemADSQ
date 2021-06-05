nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "scaf_hic-${date}"
params.reads = "/home/ubuntu/oscheius/0-inputs/test.ccsf.fasta.gz"
params.assemblies = "/home/ubuntu/oscheius/0-inputs/test.fasta"
params.restriction_sites = "GATC"




reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(_test)?(_R1)?(_R2)?(\.subsamp)?(\.fastq)?(\.fq)?(\.gz)?$/, file.Name - ~/(_test)?(\.subsamp)?(\.fastq)?(\.fq)?(\.gz)?$/, file) }
assemblies = Channel.fromPath(params.assemblies, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(_hifi)?(\.flyemetav2.6_longest20G)?(\.hifiasm)?(\.canu)?(\.flye)?(\.wtdbg2)?(\.purged)?(\.tol)?(\.fa)?(\.fasta)?(\.gz)?$/, file.Name - ~/(\.fasta)?(\.gz)?$/, file) }


process bwa_index {
    tag "${assemName}"

    input:
      tuple val(strain), val(assemName), path(assembly)

    output:
      tuple val(strain), val(assemName), path("bwa")

    script:
      """
      mkdir bwa
      bwa index -p bwa/${assemName} $assembly
      """
}

process get_restriction_sites {
    tag "${assemName}"

    input:
      tuple val(strain), val(assemName), path(assembly)

    output:
      tuple val(strain), val(assemName), path("${assemName}_rest_site_positions.bed")

    script:
      """
      findRestSite --fasta $assembly \
        --searchPattern $params.restriction_sites \
        -o ${assemName}_rest_site_positions.bed
      """
}

process bwa_mem {
    tag "${assemName}_${read_ID}"

    input:
      tuple val(strain), val(assemName), path(indexBase), val(strain2), val(read_ID), path(readFile)
      

    output:
      tuple val(strain), val(assemName), path("${assemName}_${read_ID}.bam")

    script:
      """
      INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
      
      bwa mem -A1 -B4 -E50 -L0 \
        -t ${task.cpus} \
        \$INDEX \
        $readFile \
        | samtools view -@ 2 -bhS -o ${assemName}_${read_ID}.bam -
      """
}

process hic_matrix {
    tag "${assemName}"
    publishDir "$params.outdir/hic_QC", mode: 'copy'

    input:
      tuple val(strain), val(assemName), path(bam1), path(bam2), path(bed_sites) 

    output:
      tuple val(strain), val(assemName), path("${assemName}.h5"), emit: h5_matrix
      path("${assemName}_QC"), emit: folderQC

    script:
      """
      hicBuildMatrix --samFiles $bam1 $bam2 \
                 --binSize 10000 \
                 --restrictionSequence $params.restriction_sites \
                 --threads ${task.cpus} \
                 --inputBufferSize 100000 \
                 -o ${assemName}.h5 \
                 --QCfolder ${assemName}_QC
      """
}

process hic_correct {
    tag "${assemName}"

    input:
      tuple val(strain), val(assemName), path(h5_matrix)

    output:
      tuple val(strain), val(assemName), path("${assemName}_corrected.h5"), emit: h5_corr_matrix
      path("hic_diagnostic_${assemName}.png"), emit: hic_diagnostic

    script:
      """
      hicCorrectMatrix diagnostic_plot -m $h5_matrix -o hic_diagnostic_${assemName}.png
      hicCorrectMatrix correct -m $h5_matrix --filterThreshold -4 3 -o ${assemName}_corrected.h5
      """
}



process scaffold {
    tag "${assemName}"
    publishDir "$params.outdir/scaffold", mode: 'copy'

    input:
      tuple val(strain), val(assemName), path(h5_matrix), path(assembly)

    output:
      tuple val(strain), val(assemName), path("${assemName}_scaf")

    script:
      """
      if [ -f *.gz ]; then
            gunzip -c $assembly > assembly.fasta
        else
            ln -s $assembly assembly.fasta
      fi
      assemble -m $h5_matrix -o ${assemName}_scaf \
        --min_scaffold_length 100000 --bin_size 10000 \
        --misassembly_zscore_threshold -5 \
        --num_iterations 3 -f assembly.fasta
      rm ${assemName}_scaf/*h5 ${assemName}_scaf/*graphml
      bgzip ${assemName}_scaf/super_scaffolds.fa
      """
}


workflow {
    get_restriction_sites(assemblies)
    bwa_index(assemblies)
    bwa_mem(bwa_index.out.cross(reads)
        .map{it -> tuple(it[0][0], it[0][1], it[0][2], it[1][0], it[1][1], it[1][2]) } )
    hic_matrix(bwa_mem.out
      .groupTuple(by: [0,1])
      .join(get_restriction_sites.out, by: [0,1])
      .map{it -> tuple(it[0], it[1], it[2][0], it[2][1], it[3]) })
    hic_correct(hic_matrix.out.h5_matrix)
    scaffold(hic_correct.out.h5_corr_matrix.join(assemblies, by: [0,1]))
}
