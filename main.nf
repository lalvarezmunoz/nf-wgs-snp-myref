#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {bwa_index; bwa_mem} from "./modules/bwa_mem.nf"
include {flagstat} from "./modules/flagstat.nf"
include {sortbam} from "./modules/sortbam.nf"
include {markduplicates} from "./modules/markduplicates.nf"
include {indexbam} from "./modules/indexbam.nf"
include {freebayes} from "./modules/freebayes.nf"
include {filtervcf} from "./modules/filtervcf.nf"


workflow {

    // manifest reader, csv reader
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map {row -> tuple(row.id, file(row.R1), file(row.R2))}
        .set {raw_reads}


    // Create index based on reference FASTA
    bwa_index(params.genome_id, params.genome_fna)

    // bwa_mem
    // combine the reads with the bwa indexed channel
    raw_reads
        .combine(bwa_index.out.index)
        .set{bwa_align_input}
    // Align
    bwa_mem(bwa_align_input)
    
    // flagstats
    flagstat(bwa_mem.out.bam)

    // sort bam files
    sortbam(bwa_mem.out.bam)

    // mark duplicates with picard
    markduplicates(sortbam.out.bam)

    // index alignment
    indexbam(markduplicates.out.bam)

    // variant calling with freebayes
    freebayes(indexbam.out.indexedbam, params.genome_fna)

    // filter by quality
    filtervcf(freebayes.out.vcf)
    filtervcf.out.vcf.view()
}