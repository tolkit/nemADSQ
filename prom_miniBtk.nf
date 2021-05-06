nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "miniBtk-${date}"
params.reads = "/home/ubuntu/oscheius/0-inputs/test.ccsf.fasta.gz"
params.blobtoolsPath = "${params.btkPath}/blobtools2/blobtools"
params.dmnd_db = "custom.dmnd"
params.max_target_seqs = 1000
params.evalue = 0.000001
params.taxid = 2613844
params.taxdump = "taxdump/"
params.taxrule = "bestsumorder"
params.btkFilterString = "bestsumorder_superkingdom--Keys=Bacteria"
params.kmer = '31'



// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/904/067/145/GCA_904067145.1_BOKI2/GCA_904067145.1_BOKI2_protein.faa.gz # Bursaphelenchus okinawaensis
// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_protein.faa.gz # Caenorhabditis elegans
// wget .../GCF_000005845.2_ASM584v2_protein.faa.gz # Escherichia coli str. K-12 substr. MG1655 
// wget .../GCF_002803535.1_ASM280353v1_protein.faa.gz # Ochrobactrum pituitosum
// wget .../GCF_004208635.1_ASM420863v1_protein.faa.gz # Leucobacter triazinivorans
// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/445/995/GCF_900445995.1_48290_B02/GCF_900445995.1_48290_B02_protein.faa.gz # Brevundimonas diminuta
// zcat dbs/custom/*.faa.gz | diamond makedb -p 8 -d custom --taxonmap prot.accession2taxid.gz --taxonnodes dbs/nodes.dmp

dmnd_db = Channel.fromPath(params.dmnd_db, checkIfExists: true).collect()
taxdump_db = Channel.fromPath(params.taxdump, checkIfExists: true).collect()
reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(_hifi)?(\.ccs)?(\.fa)?(\.fasta)?(\.gz)?$/, file) }

process decompress {
    tag "${strain}"

    input:
      tuple val(strain), path(reads)

    output:
      tuple val(strain), path("${strain}.decompReads.fasta")

    script:
      """
      if [ -f *.gz ]; then
            gunzip -c $reads > ${strain}.decompReads.fasta
        else
            ln -s $reads ${strain}.decompReads.fasta
      fi
      """
}

process shasta {
    tag "${strain}"

    input:
      tuple val(strain), path(reads)

    output:
      tuple val(strain), path("${strain}.shasta.fasta")

    script:
      """
      shasta-Linux-0.7.0 --input $reads \
                       --threads ${task.cpus}
      mv ShastaRun/Assembly.fasta ${strain}.shasta.fasta
      rm -r ShastaRun/
      """
}

process flye {
    tag "${strain}"
    publishDir "$params.outdir", mode: 'copy'

    input:
      tuple val(strain), path(reads)

    output:
      tuple val(strain), path("${strain}.flye.fasta")

    script:
      """
      flye --threads ${task.cpus} --nano-raw $reads --meta -o flyemeta
      mv flyemeta/assembly.fasta ${strain}.flye.fasta
      """
}

process mask_assembly {
    tag "${strain}"
    label 'btk'

    input:
      tuple val(strain), path(assembly)

    output:
      tuple val(strain), path("${strain}.masked.fasta")

    script:
      """
      windowmasker -in $assembly \
                      -infmt fasta \
                      -mk_counts \
                      -sformat obinary \
                      -out tmp.counts 2> log \
        && windowmasker -in $assembly \
                        -infmt fasta \
                        -ustat tmp.counts \
                        -dust T \
                        -outfmt fasta \
                        -out ${strain}.masked.fasta
      """
}


process chunk_assembly {
    tag "${strain}"
  label 'btk'
    input:
      tuple val(strain), path(assembly)

    output:
      tuple val(strain), path("${strain}.chunks.fasta")

    script:
      """
      chunk_fasta.py --in ${assembly} \
        --chunk 100000 --overlap 0 --max-chunks 20 \
        --out ${strain}.chunks.fasta
      """
}

process diamond_search {
    tag "${strain}"
    label 'btk'

    input:
      tuple val(strain), path(assembly)
      path(dmnd_db)

    output:
      tuple val(strain), path("${strain}.chunk.diamond.tsv")

    script:
      """
      diamond blastx \
            --query ${assembly} \
            --db $dmnd_db \
            --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --max-target-seqs $params.max_target_seqs \
            --max-hsps 1 \
            --evalue $params.evalue \
            --threads ${task.cpus} \
            > ${strain}.chunk.diamond.tsv
      """
}

process unchunk_hits {
    tag "${strain}"
    publishDir "$params.outdir/", mode: 'copy'
  label 'btk'
    input:
      tuple val(strain), path(diamondHits)

    output:
      tuple val(strain), path("${strain}.diamond.tsv")

    script:
      """
      unchunk_blast.py --in $diamondHits \
        --out ${strain}.diamond.tsv
      """
}

process subsample_reads {
    tag "${strain}"

    input:
      tuple val(strain), path(reads)

    output:
      tuple val(strain), path("${strain}.subsamp.fa.gz")

    script:
      """
      seqkit sample -p $params.sampRate $reads | gzip -c > ${strain}.subsamp.fa.gz
      """
}

process map_reads {
    tag "${strain}"
    label 'btk'

    input:
      tuple val(strain), path(reads), path(assembly)

    output:
      tuple val(strain), path("${strain}.bam")

    script:
      """
      minimap2 -a -k 19 -w 10 -I 10G -g 5000 -r 2000 -N 100 \
        --lj-min-ratio 0.5 -A 2 -B 5 -O 5,56 -E 4,1 -z 400,50 \
        --sam-hit-only -t ${task.cpus} ${assembly} \
        $reads | \
        samtools sort -@ ${task.cpus} -o ${strain}.bam
      """
}

process add_hits_and_coverage {
    tag "${strain}"
    publishDir "$params.outdir/btkDatasets", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(assembly), path(diamondHits), path(bam)
      path(taxdump_db)

    output:
      tuple val(strain), path("${strain}"), emit: blobDir
      path ("*.{png,svg}"), emit: blobplots

    script:
      """
      blobtools create \
            --fasta ${assembly} \
            --taxdump $taxdump_db \
            --taxid $params.taxid \
            $strain

      blobtools add \
            --hits ${diamondHits} \
            --taxdump $taxdump_db \
            --taxrule $params.taxrule \
            --cov ${bam}=reads \
            $strain
      blobtools view \
            --view blob \
            --param plotShape=circle \
            --param bestsumorder_phylum--Order=no-hit%2CNematoda%2CProteobacteria%2CActinobacteria \
            --format png \
            $strain
      """
}

process btk_static_images {
    tag "${strain}"
    publishDir "$params.outdir/blobplots", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(btkdir)

    output:
      path ("*.{png,svg}")

    script:
      """
      blobtools view \
            --view blob \
            --param plotShape=circle \
            --param bestsumorder_phylum--Order=no-hit%2CNematoda%2CProteobacteria%2CActinobacteria \
            --format png --format svg \
            $btkdir
      blobtools view \
            --view cumulative \
            --format png --format svg \
            $btkdir
      """
}


process filter_fasta {
    tag "${strain}"
    publishDir "$params.outdir/filteredData", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(btkdir), path(bam), path(assembly), path(reads)

    output:
      tuple val(strainName), path("$filtered_assemFile"), emit: filtered_assem
      tuple val(strainName), path("$filtered_readsFile"), emit: filtered_reads

    script:
    strainName = strain + "_filtered"
    btk_fltrd_assemFile = assembly.baseName - ~/(\.fasta)?(\.fa)?$/ + ".filtered.fasta"
    btk_fltrd_readsFile = reads.baseName - ~/(\.gz)?$/ + ".filtered.gz"
    filtered_assemFile = strainName + ".hifiasm.fasta"
    filtered_readsFile = strainName + ".ccs.fasta.gz"
      """
      blobtools filter \
        --query-string "${params.btkFilterString}" \
        --fasta $assembly \
        --fastq $reads \
        --cov $bam \
        $btkdir
      mv $btk_fltrd_readsFile $filtered_readsFile
      mv $btk_fltrd_assemFile $filtered_assemFile
      """
}

workflow raw_asses {
    take: reads
    main:
        decompress(reads)
        shasta(decompress.out) | mask_assembly | chunk_assembly
        diamond_search(chunk_assembly.out, dmnd_db) | unchunk_hits
        map_reads(decompress.out.join(shasta.out))
        add_hits_and_coverage(shasta.out.join(unchunk_hits.out.join(map_reads.out)), taxdump_db)
        filter_fasta(add_hits_and_coverage.out.blobDir.join(map_reads.out.join(shasta.out.join(reads))))
    emit:
        filter_fasta.out.filtered_reads
}

workflow fltd_asses {
    take: reads
    main:
        flye(reads) | mask_assembly | chunk_assembly
        diamond_search(chunk_assembly.out, dmnd_db) | unchunk_hits
        map_reads(reads.join(flye.out))
        add_hits_and_coverage(flye.out.join(unchunk_hits.out.join(map_reads.out)), taxdump_db)
    emit:
        add_hits_and_coverage.out.blobDir
}

workflow {
    raw_asses(reads)
    fltd_asses(raw_asses.out)
}

