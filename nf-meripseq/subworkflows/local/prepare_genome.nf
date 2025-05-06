//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA       } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF         } from '../../modules/nf-core/gunzip'
include { UNTAR as UNTAR_RSEM_INDEX    } from '../../modules/nf-core/untar/main'
include { RSEM_PREPAREREFERENCE        } from '../../modules/nf-core/rsem/preparereference'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../modules/nf-core/rsem/preparereference'

workflow PREPARE_GENOME {
    take:
    fasta                //    path: path to genome fasta file
    gtf                  //    file: /path/to/genome.gtf
    rsem_index           //    file: /path/to/resm/index/

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], file(fasta, checkIfExists: true) ] ).gunzip.map{ it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta, checkIfExists: true))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //

    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], file(gtf, checkIfExists: true) ] ).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value(file(gtf, checkIfExists: true))
        }
    } else if (gff) {
        if (gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [:], file(gff, checkIfExists: true) ] ).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.value(file(gff, checkIfExists: true))
        }
        ch_gtf      = GFFREAD ( [ [:], ch_gff ] ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_gtf ).transcript_fasta
    ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)

    //
    // Uncompress RSEM index or generate from scratch if required
    //
    ch_rsem_index = Channel.empty()
        if (rsem_index) {
            if (rsem_index.endsWith('.tar.gz')) {
                ch_rsem_index = UNTAR_RSEM_INDEX ( [ [:], rsem_index ] ).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_RSEM_INDEX.out.versions)
            } else {
                ch_rsem_index = Channel.value(file(rsem_index))
            }
        } else {
            ch_rsem_index = RSEM_PREPAREREFERENCE ( ch_fasta, ch_gtf ).index
            ch_versions   = ch_versions.mix(RSEM_PREPAREREFERENCE.out.versions)
        }

    emit:
    fasta            = ch_fasta                  //    path: genome.fasta
    gtf              = ch_gtf                    //    path: genome.gtf
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    rsem_index       = ch_rsem_index             //    path: rsem/index/
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}