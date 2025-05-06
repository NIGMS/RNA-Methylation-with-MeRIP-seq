process EXOMEPEAK2_PEAK {
    tag "$meta.id"
    label 'process_high'

    // Conda environment with required packages
    conda "conda-forge::r-base conda-forge::r-optparse conda-forge::r-biocmanager bioconda::bioconductor-exomepeak2"

    // Container setup: Use local .sif or Docker image
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.singularity_image_path : params.docker_image_path }"

    input:
    tuple val(meta), path(ipbam), path(inputbam)  // IP and control BAM files
    path gtf                                        // GTF annotation file
    val bsgenome

    

    output:
    tuple val(meta), path("peaks.bed"),   emit: bed                   // BED format for peaks
    tuple val(meta), path("peaks.csv"),   emit: csv                   // Tabular results
    tuple val(meta), path("gc_fit.pdf"),     emit: gc, optional: true    // GC bias pdf, optional: true
    path "versions.yml",                            emit: versions              // Version info

    when:
    task.ext.when == null || task.ext.when  // Run unless explicitly disabled

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"  // Output prefix from meta.id
    def args   = task.ext.args   ?: ''      // Optional arguments from config   

    
    """
    # Run the exomePeak2 R script
    exomepeak2.r \\
        --bam_ip "${ipbam.join(',')}" \\
        --bam_input "${inputbam.join(',')}" \\
        --gtf "$gtf" \\
        --bsgenome "${bsgenome}" \\
        --outprefix "$prefix" \\
        --cores $task.cpus \\
        $args
    # Move outputs to working directory
    mv $prefix/peaks.bed .
    mv $prefix/peaks.csv .
    mv $prefix/gc_fit.pdf .

    # Generate versions.yml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/^.*R version //; s/ .*\$//')
        exomepeak2: \$(Rscript -e "library(exomePeak2); cat(as.character(packageVersion('exomePeak2')))")
    END_VERSIONS
    """
}