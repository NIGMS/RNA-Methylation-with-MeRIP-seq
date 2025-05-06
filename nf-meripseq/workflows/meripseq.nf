/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGN_STAR                                } from '../subworkflows/local/align_star'
include { PEAKS_CALL_MACS3_HOMER                    } from '../subworkflows/local/peaks_call_macs3_homer'
include { BED_CONSENSUS_MACS3                       } from '../subworkflows/local/bed_consensus_macs3'
include { QUANTIFY_RSEM                             } from '../subworkflows/local/quantify_rsem'
include { BAM_BEDGRAPH_BIGWIG_BEDTOOLS_UCSC         } from '../subworkflows/local/bam_bedgraph_bigwig_bedtools_ucsc'
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS      } from '../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'

include { MULTIQC                                   } from '../modules/nf-core/multiqc/main'
include { KHMER_UNIQUEKMERS                         } from '../modules/nf-core/khmer/uniquekmers/main'
include { HOMER_ANNOTATEPEAKS as ANNOT_PEAK         } from '../modules/nf-core/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as ANNOT_DIFF_PEAK    } from '../modules/nf-core/homer/annotatepeaks/main'
include { EXOMEPEAK2_PEAK as EXOMEPEAK2_SINGLE      } from '../modules/local/exomepeak2/peak'
include { EXOMEPEAK2_PEAK as EXOMEPEAK2_CONSENSUS   } from '../modules/local/exomepeak2/peak'
include { EXOMEPEAK2_DIFF                           } from '../modules/local/exomepeak2/diff/main'

include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText; validate_contrast;
          groupBamFilesByGroup; createContrastPairs } from '../subworkflows/local/utils_nfcore_meripseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Header files for MultiQC
ch_peak_count_header         = file("$projectDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header         = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_annotation_header    = file("$projectDir/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
ch_deseq2_pca_header        = Channel.value(file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true))
ch_deseq2_clustering_header = Channel.value(file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true))

workflow MERIPSEQ {

    take:
    ch_samplesheet         // channel: samplesheet read in from --input
    ch_versions            // channel: [ path(versions.yml) ]
    ch_fasta               // channel: path(genome.fa)
    ch_gtf                 // channel: path(genome.gtf)
    ch_star_rsem_index     // channel: path(rsem/index/)
    ch_transcript_fasta    // channel: path(transcript.fa)

    main:
    ch_multiqc_files = Channel.empty()

    // SUBWORKFLOW:  Run RNA-seq FASTQ preprocessing subworkflow
    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_samplesheet,                       // channel: [ val(meta), [ reads ] ]
        ch_fasta,                             // channel: /path/to/genome.fasta
        ch_transcript_fasta,                  // channel: /path/to/transcript.fasta
        ch_gtf,                               // channel: /path/to/genome.gtf
        null,                                 // ch_salmon_index,
        null,                                 // ch_sortmerna_index 
        null,                                 // ch_bbsplit_index   
        null,                                 // optional
        true,                                 // params.skip_bbsplit,
        params.skip_fastqc || params.skip_qc, // skip_fastqc 
        params.skip_trimming,                 // skip_trimming
        true,                                 // params.skip_umi_extract,
        true,                                 // !salmon_index_available,
        false,                                // make_sortmerna_index, boolean
        'trimgalore',                         // params.trimmer,
        params.min_trimmed_reads,             // integer: > 0
        params.save_trimmed,                  // 
        false,                                // params.remove_ribo_rna
        false,                                // params.with_umi,
        0,                                    // params.umi_discard_read,
        0.8,                                  // params.stranded_threshold,
        0.1,                                  // params.unstranded_threshold,
        true,                                 // params.skip_linting,
        false                                 // fastp_merge
    )
    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
    ch_versions                       = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)
    ch_strand_inferred_filtered_fastq = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
    ch_trim_read_count                = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_read_count

    //
    // MODULE: Run STAR alignment on trimmed reads
    //
    ch_genome_bam          = Channel.empty()
    ch_genome_bam_index    = Channel.empty()
    ch_star_log            = Channel.empty()
    ch_unaligned_sequences = Channel.empty()
    ch_transcriptome_bam   = Channel.empty()

    def is_aws_igenome = false
    if (params.fasta && params.gtf) {
        if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
            is_aws_igenome = true
        }
    }

    ALIGN_STAR (
        ch_strand_inferred_filtered_fastq,                  // Trimmed FASTQ files from Trim Galore!
        ch_star_rsem_index.map { [ [:], it ] },  // Transform to [meta2, path]
        ch_gtf.map { [ [:], it ] },         // Transform to [meta3, path]
        params.star_ignore_sjdbgtf ?: false,
        '',
        params.seq_center ?: '',
        is_aws_igenome,
        ch_fasta.map { [ [:], it ] }
    )

    ch_genome_bam          = ALIGN_STAR.out.bam
    ch_genome_bam_index    = ALIGN_STAR.out.bai
    ch_transcriptome_bam   = ALIGN_STAR.out.bam_transcript
    ch_star_log            = ALIGN_STAR.out.log_final
    ch_unaligned_sequences = ALIGN_STAR.out.fastq
    ch_multiqc_files = ch_multiqc_files.mix(ch_star_log.collect{it[1]})
    ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
    
    ch_multiqc_files = ch_multiqc_files
                .mix(ALIGN_STAR.out.stats.collect{it[1]})
                .mix(ALIGN_STAR.out.flagstat.collect{it[1]})
                .mix(ALIGN_STAR.out.idxstats.collect{it[1]})
    

    ch_genome_bam.join(ch_genome_bam_index, by: [0])
        .set {ch_genome_bam_bai }
    ch_genome_bam_bai
        .map {
            meta, bam, bai ->
                meta.control ? null : [ meta.id, [ bam ] , [ bai ] ]
        }
        .set { ch_control_bam_bai }
    ch_genome_bam_bai
        .map {
            meta, bam, bai ->
                meta.control ? [ meta.control, meta, [ bam ], [ bai ] ] : null
        }
        .combine(ch_control_bam_bai, by: 0)
        .map { it -> [ it[1] , it[2] + it[4], it[3] + it[5] ] }
        .set { ch_ip_control_bam_bai }


    //
    // SUBWORKFLOW: Normalised bigWig coverage tracks
    //
    BAM_BEDGRAPH_BIGWIG_BEDTOOLS_UCSC (
        ALIGN_STAR.out.bam.join(ALIGN_STAR.out.flagstat, by: [0]),
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_BEDGRAPH_BIGWIG_BEDTOOLS_UCSC.out.versions)


    //
    // SUBWORKFLOW: RNAseq RSEM gene counts for INPUT samples
    //
    ch_transcriptome_bam
        .map {
            meta, bam ->
                if (meta.containsKey('control') && meta.control) {
                    return null           
                } else {
                    return [ meta, bam ]  // Keep only control input samples
                }
        }
        .filter { it != null }            // Remove any null entries
        .set { ch_input_bam }


     QUANTIFY_RSEM (
        ch_strand_inferred_filtered_fastq,
        ch_star_rsem_index,
        ch_fasta.map { [ [:], it ] }
     )
   
    // You can now use the RSEM count outputs
    ch_versions = ch_versions.mix(QUANTIFY_RSEM.out.versions)
    

    //
    // MODULE: Calculute genome size with khmer
    //
    ch_macs_gsize     = Channel.empty()
    if (!params.macs_gsize) {
        KHMER_UNIQUEKMERS (
            ch_fasta,
            params.read_length
        )
        ch_macs_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it.text.trim() }
    }else{
        ch_macs_gsize     = params.macs_gsize
    }

    // Create channels: [ meta, ip_bam, control_bam ]
    ch_ip_control_bam_bai
        .map {
            meta, bams, bais ->
                [ meta , bams[0], bams[1] ]
        }
        .set { ch_ip_control_bam }

    //
    // SUBWORKFLOW: Call peaks with MACS3, annotate with HOMER and perform downstream QC
    //
    PEAKS_CALL_MACS3_HOMER (
        ch_ip_control_bam,
        ch_fasta,
        ch_gtf,
        ch_macs_gsize,
        "_peaks.annotatePeaks.txt",
        ch_peak_count_header,
        ch_frip_score_header,
        ch_peak_annotation_header,
        params.narrow_peak,
        params.skip_peak_annotation,
        params.skip_peak_qc
    )
    ch_versions = ch_versions.mix(PEAKS_CALL_MACS3_HOMER.out.versions)
    ch_multiqc_files = ch_multiqc_files
                .mix(PEAKS_CALL_MACS3_HOMER.out.frip_multiqc.collect{it[1]})
                .mix(PEAKS_CALL_MACS3_HOMER.out.peak_count_multiqc.collect{it[1]})
                .mix(PEAKS_CALL_MACS3_HOMER.out.plot_homer_annotatepeaks_tsv.collect())
    //
    //  Consensus peaks analysis
    //
    ch_subreadfeaturecounts_multiqc   = Channel.empty()
    ch_macs3_consensus_bed_lib        = Channel.empty()
    ch_macs3_consensus_txt_lib        = Channel.empty()
    ch_deseq2_pca_multiqc             = Channel.empty()
    ch_deseq2_clustering_multiqc      = Channel.empty()
    if (!params.skip_consensus_peaks) {
        // Create channels: [group, [ ip_bams ] ]
        ch_ip_control_bam
            .map {
                meta, ip_bam, control_bam ->
                   [ meta.group, ip_bam ]
            }
            .groupTuple()
            .set { ch_group_bams }

        BED_CONSENSUS_MACS3 (
            PEAKS_CALL_MACS3_HOMER.out.peaks,
            ch_group_bams,
            ch_fasta,
            ch_gtf,
            ch_deseq2_pca_header,
            ch_deseq2_clustering_header,
            params.narrow_peak,
            params.skip_peak_annotation,
            params.skip_deseq2_qc
        )
        ch_macs3_consensus_bed_lib       = BED_CONSENSUS_MACS3.out.consensus_bed
        ch_macs3_consensus_txt_lib       = BED_CONSENSUS_MACS3.out.consensus_txt
        ch_versions                      = ch_versions.mix(BED_CONSENSUS_MACS3.out.versions)
    }
    ch_multiqc_files = ch_multiqc_files
                .mix(BED_CONSENSUS_MACS3.out.featurecounts_summary.collect{it[1]})
                .mix(BED_CONSENSUS_MACS3.out.deseq2_qc_pca_multiqc.collect{it[1]})
                .mix(BED_CONSENSUS_MACS3.out.deseq2_qc_dists_multiqc.collect())
 
    //ExomePeak2 workflow
    if (!params.skip_exomepeak2) {
        def bsgenome = params.bsgenome ?: (params.genome && params.bsgenomes.containsKey(params.genome) ? params.bsgenomes[params.genome] : null)
        
        // Run exomePeak2 on individual samples
        if (!params.skip_exomepeak2_single) {
            ch_single_samples = ch_ip_control_bam.map { meta, ip, input ->
                [[id: meta.id, type: "single"], ip, input]
            }
            EXOMEPEAK2_SINGLE(
                ch_single_samples,
                ch_gtf,
                bsgenome
            )
        }

        // Run consensus peak calling for each group
        if (!params.skip_exomepeak2_consensus) {
            // Group BAM files by their group (reuse your existing function)
            ch_ip_control_bam_by_group = groupBamFilesByGroup(ch_ip_control_bam)
            ch_consensus_samples = ch_ip_control_bam_by_group.map { group_id, ip_bams, input_bams ->
                [[id: group_id, type: "consensus"], ip_bams, input_bams]
            }
            
            // Run consensus peak analysis
            EXOMEPEAK2_CONSENSUS(
                ch_consensus_samples,
                ch_gtf,
                bsgenome
            )
            ch_versions = ch_versions.mix(EXOMEPEAK2_CONSENSUS.out.versions)

            // Annotate peaks
            if (!params.skip_peak_annotation) {
                ANNOT_PEAK (
                    EXOMEPEAK2_CONSENSUS.out.bed,
                    ch_fasta,
                    ch_gtf
                )
                ch_versions = ch_versions.mix(ANNOT_PEAK.out.versions)
            }
        }
        
        // Run differential exomePeak2 analysis
        if (!params.skip_differential_peaks) {
            if (!params.contrast) {
                error("ERROR: Need --contrast to define groups in differential analysis.")
            } else {
                // Validate contrasts against sample sheet
                contrast_list = validate_contrast(params.contrast, ch_samplesheet)
                
                // Group BAM files by their group
                ch_ip_control_bam_grouped = groupBamFilesByGroup(ch_ip_control_bam)
                
                // Create contrast pairs for differential analysis
                ch_diff_contrasts = createContrastPairs(ch_ip_control_bam_grouped, contrast_list)
            }

            // Run differential analysis
            EXOMEPEAK2_DIFF(
                ch_diff_contrasts,
                ch_gtf,
                bsgenome
            )
            ch_versions = ch_versions.mix(EXOMEPEAK2_DIFF.out.versions)
            if (!params.skip_peak_annotation) {
                EXOMEPEAK2_DIFF.out.bed
                    .map { contrast, peak ->
                        // Create a meta map with id and single_end fields
                        def meta = [id: contrast, single_end: true]  // Adjust single_end based on your data
                        [ meta, peak ]
                    }
                    .set { ch_diff_peaks_for_annotation }

                ANNOT_DIFF_PEAK (
                    ch_diff_peaks_for_annotation,
                    ch_fasta,
                    ch_gtf
                )
                ch_versions = ch_versions.mix(ANNOT_DIFF_PEAK.out.versions)
            }
        }
    }


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'meripseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    MULTIQC (
        ch_multiqc_files.collect().reverse(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
