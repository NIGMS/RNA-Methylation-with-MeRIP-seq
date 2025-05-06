//
// Convert BAM to normalised bigWig via bedGraph using BEDTools and UCSC
//

include { BEDTOOLS_GENOMECOV    } from '../../../modules/local/bedtools_genomecov'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { CUSTOM_GETCHROMSIZES  } from '../../../modules/nf-core/custom/getchromsizes/main'

workflow BAM_BEDGRAPH_BIGWIG_BEDTOOLS_UCSC {
    take:
    ch_bam_flagstat // channel: [ val(meta), [bam], [flagstat] ]
    ch_fasta  // channel: [ fasta ]

    main:

    ch_versions    = Channel.empty()
    ch_chrom_sizes = Channel.empty()

    //
    // Get chromosome sizes
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Create bedGraph coverage track
    //
    BEDTOOLS_GENOMECOV (
        ch_bam_flagstat
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    //
    // Create bigWig coverage tracks
    //
    UCSC_BEDGRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.bedgraph,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

    emit:
    bedgraph     = BEDTOOLS_GENOMECOV.out.bedgraph      // channel: [ val(meta), [ bedgraph ] ]
    scale_factor = BEDTOOLS_GENOMECOV.out.scale_factor  // channel: [ val(meta), [ txt ] ]

    bigwig       = UCSC_BEDGRAPHTOBIGWIG.out.bigwig     // channel: [ val(meta), [ bigwig ] ]

    versions     = ch_versions                          // channel: [ versions.yml ]
}