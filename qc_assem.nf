nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "assemblyQC-${date}"
params.reads = "mini_DF5120.ccs.fasta.gz"
params.assemblies = "mini_DF5120.hifiasm.fasta.gz"
params.odb = 'nematoda_odb10'
params.telomere = 'TTAGGC'
params.busco2nigons = "gene2Nigon_busco20200927.tsv.gz"
params.min_occurr = 15
params.teloRepeatWindowSize = 10000
params.minimumGenesPerSequence = 15
params.minimumNigonFrac = 0.9
params.minimumFracAlignedTeloReads = 0.1
params.windowSizeQC = 5e5

reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(_filtered)?(\.telo)?(\.ccs)?(\.npr)?(\.fa)?(\.fasta)?(\.gz)?$/, file) }

fastFiles = Channel.fromPath(params.assemblies, checkIfExists: true)
assemblies = fastFiles.map { file -> tuple(file.Name - ~/(_filtered)?(\.hifiasm)?(\.flye)?(\.wtdbg2)?(\.canu_plus_flye)?(\.canu)?(\.purged)?(\.hic_scaff)?(\.fa)?(\.fasta)?(\.gz)?$/, file.Name - ~/(\.fa)?(\.fasta)?(\.gz)?$/, file) }

busco2nigons = Channel.fromPath(params.busco2nigons, checkIfExists: true).collect()

busco_dbs = Channel.of(params.odb.split(','))
geno_busco = assemblies.combine(busco_dbs)


process decompress_fasta {
    tag "${assembler}"

    input:
      tuple val(strain), val(assembler), path(assembly)

    output:
      tuple val(strain), val(assembler), path('assembly.fasta')

    script:
      """
      if [ -f *.gz ]; then
            gunzip -c $assembly > assembly.fasta
        else
            ln -s $assembly assembly.fasta
      fi
      """
}

process busco {
    tag "${assembler}_${busco_db}"
    publishDir "$params.outdir/busco", mode: 'copy'

    input:
      tuple val(strain), val(assembler), path(genome), val(busco_db)

    output:
      path "*single_copy_busco_sequences.{faa,fna}"
      path "${assembler}_${busco_db}_short_summary.txt"
      tuple val(assembler), path( "${assembler}_${busco_db}_full_table.tsv"), emit: busco_full

    script:
      """
      busco -c ${task.cpus} -l $busco_db -i $genome --out run_busco --mode geno
      awk 'BEGIN{FS="\\t";OFS=FS}(\$3 !~ /:/){print}' run_busco/run_*/full_table.tsv > ${assembler}_${busco_db}_full_table.tsv
      mv run_busco/short_summary* ${assembler}_${busco_db}_short_summary.txt
      #mv run_busco/run_*/full_table.tsv ${assembler}_${busco_db}_full_table.tsv
      for ext in .faa; do
        seqFile=${assembler}_${busco_db}_single_copy_busco_sequences\$ext
        for file in run_busco/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/*\$ext; do
          echo \">\$(basename \${file%\$ext})\" >> \$seqFile; tail -n +2 \$file >> \$seqFile;
        done
      done
      rm -rf run_busco/ busco_downloads/
      """
}

process red {
    tag "${assembler}"
    // publishDir "$params.outdir/red", mode: 'copy'

    input:
      tuple val(strain), val(assembler), path(assembly)

    output:
      tuple val(strain), val(assembler), path("Red_out/assembly.msk")

    script:
      """
      mkdir Red_dir Red_out
      mv $assembly Red_dir/assembly.fa

      Red -gnm Red_dir -msk Red_out/
      """
}

process red2bed {
    tag "${assembler}"
    label 'nemaQC'
    publishDir "$params.outdir/red", mode: 'copy'

    input:
      tuple val(strain), val(assembler), path(assembly)

    output:
      tuple val(assembler), path("${assembler}.red.bed.gz")

    script:
      """
      samtools faidx ${assembly}
      cut -f 1,2 ${assembly}.fai > ${assembler}.seqlen.tsv

      cat $assembly | soft_mask2bed | \
        perl -ne 's/ [^\\t]+\\t/\\t/; print' > ${assembler}.bed
      
      bedtools makewindows -g ${assembler}.seqlen.tsv -w ${params.teloRepeatWindowSize} | \
        bedtools coverage -a stdin -b ${assembler}.bed | \
        awk -F '\\t' 'BEGIN{OFS=FS}{print \$1, \$2, \$3, \$7, "repeats"}' | \
        gzip -c > ${assembler}.red.bed.gz
      rm ${assembler}.seqlen.tsv ${assembler}.bed
      """
}

process gc_by_windows {
    tag "${assembler}"
    label 'nemaQC'
    publishDir "$params.outdir/gc", mode: 'copy'

    input:
      tuple val(strain), val(assembler), path(assembly)

    output:
      tuple val(assembler), path("${assembler}.gc.bed.gz")

    script:
      """
      samtools faidx ${assembly}
      cut -f 1,2 ${assembly}.fai > ${assembler}.seqlen.tsv
      bedtools makewindows -g ${assembler}.seqlen.tsv \
        -w ${params.teloRepeatWindowSize} > windows.bed

      bedtools nuc -fi assembly.fasta -bed windows.bed | \
        cut -f 1-4 | tail -n+2 | \
        gzip -c > ${assembler}.gc.bed.gz
      
      rm assembly.fasta* ${assembler}.seqlen.tsv windows.bed
      """
}

process get_telomeric_reads {
    tag "${strain}"
    label 'nemaQC'
    publishDir "$params.outdir/teloReads", mode: 'copy'

    input:
      tuple val(strain), path(reads)
    output:
      tuple val(strain), path("${strain}.telo.fasta.gz")

    script:
      """
      zcat $reads | filter_telomeric_reads.py -r - -s ${params.telomere} \
        --times ${params.min_occurr} --out ${strain}.telo.fasta.gz
      """
}


process map_telomeric_reads {
    tag "${assemblies[1]}"
    publishDir "$params.outdir/teloMaps", mode: 'copy'
    label 'nemaQC'

    input:
      tuple val(reads), val(assemblies)

    output:
      tuple val("${assemblies[1]}"), path("${assemblies[1]}.teloMapped.paf.gz")

    script:
      """
      minimap2 ${assemblies[2]} ${reads[1]} | \
        gzip -c > ${assemblies[1]}.teloMapped.paf.gz
      """
}


process count_telomeric_repeat {
    tag "${assembler}"
    publishDir "$params.outdir/teloRepeatCounts", mode: 'copy'
    label 'nemaQC'

    input:
      tuple val(strain), val(assembler), path(assembly)
    
    output:
      tuple val(assembler), path( "${assembler}_teloRepeatCounts.tsv.gz")

    script:
      """
      samtools faidx ${assembly}
      cut -f 1,2 ${assembly}.fai > ${assembler}.seqlen.tsv
      seqkit locate -i --bed -M -G -p ${params.telomere} $assembly | \
        cut -f 1,2,3 | \
        sort -k1,1 -k2,2n > ${assembler}.teloRepeats.tsv
      bedtools makewindows -g ${assembler}.seqlen.tsv -w ${params.teloRepeatWindowSize} | \
        bedtools intersect -a stdin -b ${assembler}.teloRepeats.tsv -wa -wb | \
        bedtools groupby -i stdin -g 1,2,3 -c 1 -o count | \
        awk -F '\\t' 'BEGIN{OFS=FS}{print \$0, "telomeres"}' | \
        gzip -c > ${assembler}_teloRepeatCounts.tsv.gz
      rm ${assembler}.seqlen.tsv ${assembler}.teloRepeats.tsv
      """
}


process nematode_chromosome_QC {
    tag "${assembler}"
    publishDir "$params.outdir/nemaChromQC", mode: 'copy'
    label 'nemaQC'

    input:
      tuple val(assembler), path(buscoTable), path(teloMappedReads),
      path(teloRepeats), path(allRepeats), path(gcWindows)
      path(busco2nigons)
    
    output:
      path "${assembler}.{pdf,buscoString.txt,teloMappedBlocks.tsv}"
      path "${assembler}.chromQC.tsv" , emit: busco_full

    script:
      """
      nemaChromQC.R --assemblyName $assembler \
        --nigon $busco2nigons --busco $buscoTable \
        --teloMappedPaf $teloMappedReads \
        --teloRepeats $teloRepeats \
        --allRepeats $allRepeats \
        --gcFrac $gcWindows \
        --minimumGenesPerSequence $params.minimumGenesPerSequence \
        --minimumNigonFrac $params.minimumNigonFrac \
        --minimumFracAlignedTeloReads $params.minimumFracAlignedTeloReads \
        --windowSize $params.windowSizeQC
      """
}

process get_contiguity_stats {
    tag "all"
    publishDir "$params.outdir/", mode: 'copy'
    label 'nemaQC'

    input:
      path(assemblies)
    
    output:
      path "assemblies_contiguity_stats.tsv"

    script:
      """
      seqkit stats -a -T $assemblies > assemblies_contiguity_stats.tsv
      """
}


workflow {
    decompress_fasta(assemblies)
    busco(decompress_fasta.out.combine(busco_dbs))
    red(decompress_fasta.out) | red2bed
    gc_by_windows(decompress_fasta.out)
    count_telomeric_repeat(decompress_fasta.out)
    get_telomeric_reads(reads)
    map_telomeric_reads(get_telomeric_reads.out.cross(decompress_fasta.out))
    get_contiguity_stats(fastFiles.collect())
    nematode_chromosome_QC(busco.out.busco_full.join(map_telomeric_reads.out.join(count_telomeric_repeat.out.join(red2bed.out.join(gc_by_windows.out)))),
     busco2nigons)
}
