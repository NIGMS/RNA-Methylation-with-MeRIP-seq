process EXOMEPEAK2_DIFF {
    tag "${contrast}"
    label 'process_high'

    // Conda environment with required packages
    conda "conda-forge::r-base conda-forge::r-optparse conda-forge::r-biocmanager bioconda::bioconductor-exomepeak2"

    // Container setup: Use local .sif or Docker image
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.singularity_image_path : params.docker_image_path }"

    input:
        tuple val(contrast), val(group1), path(ip_bams1), path(input_bams1), val(group2), path(ip_bams2), path(input_bams2)
        path gtf
        val bsgenome


    output:
    tuple val(contrast), path("diffPeaks.bed"),   emit: bed                   // BED format for peaks
    tuple val(contrast), path("diffPeaks.csv"),   emit: csv                   // Tabular results
    tuple val(contrast), path("gc_fit.pdf"),  emit: gc, optional: true    // GC bias pdf, optional: true
    path "versions.yml",                      emit: versions              // Version info

    when:
    task.ext.when == null || task.ext.when  // Run unless explicitly disabled

    script:
    def prefix = task.ext.prefix ?: "${contrast}"  // Output prefix from contrast
    def args   = task.ext.args   ?: ''      // Optional arguments from config   
 

    def ip_bams1_paths = ip_bams1.join(',')
    def input_bams1_paths = input_bams1.join(',')
    def ip_bams2_paths = ip_bams2.join(',')
    def input_bams2_paths = input_bams2.join(',')

    """
    

    # Run the exomePeak2 R script
    exomepeak2.r \\
        --bam_ip "$ip_bams1_paths" \\
        --bam_input "$input_bams1_paths" \\
        --control_bam_ip "$ip_bams2_paths" \\
        --control_bam_input "$input_bams2_paths" \\
        --gtf "$gtf" \\
        --bsgenome "${bsgenome}" \\
        --outprefix "$prefix" \\
        --cores $task.cpus \\
        $args
       # Move outputs to working directory
    mv $prefix/diffPeaks.bed .
    mv $prefix/diffPeaks.csv .
    mv $prefix/gc_fit.pdf .

    # Generate versions.yml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/^.*R version //; s/ .*\$//')
        exomepeak2: \$(Rscript -e "library(exomePeak2); cat(as.character(packageVersion('exomePeak2')))")
    END_VERSIONS
    """
}