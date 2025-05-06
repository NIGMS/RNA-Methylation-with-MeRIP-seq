#!/usr/bin/env Rscript
user_lib <- Sys.getenv("R_LIBS_USER", unset = file.path(getwd(),"r_libs_user"))
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

# Load required libraries
suppressPackageStartupMessages({
    library(optparse)
    library(exomePeak2)
    library(GenomicFeatures)
    library(BiocParallel)
})

# Define command-line options
option_list <- list(
    make_option(c("--bam_ip"), type = "character", default = NULL, 
                help = "Path to IP BAM file", metavar = "FILE"),
    make_option(c("--bam_input"), type = "character", default = NULL, 
                help = "Path to control/input BAM file", metavar = "FILE"),
    make_option(c("--control_bam_ip"), type = "character", default = NULL, 
                help = "Path to IP BAM file", metavar = "FILE"),
    make_option(c("--control_bam_input"), type = "character", default = NULL, 
                help = "Path to control/input BAM file", metavar = "FILE"),
    make_option(c("--gtf"), type = "character", default = NULL, 
                help = "Path to GTF annotation file", metavar = "FILE"),
    make_option(c("--bsgenome"), type = "character", default = NULL, 
                help = "BSgenome package name (e.g., BSgenome.Hsapiens.UCSC.hg38)"),
    make_option(c("--outprefix"), type = "character", default = "exomepeak2", 
                help = "Prefix for output files [default: %default]", metavar = "PREFIX"),
    make_option(c("--cores"), type = "integer", default = 1, 
                help = "Number of CPU cores to use [default: %default]", metavar = "INT")
)


# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Set environment variables to help with parallel execution
Sys.setenv(R_PARALLEL_PORT = "random")
Sys.setenv(R_PARALLEL_SOCKET_TIMEOUT = "3600")

# Override the problematic featuresCounts function in exomePeak2
if (requireNamespace("exomepeak2", quietly = TRUE)) {
  # Create a wrapper for the original function
  original_featuresCounts <- exomepeak2:::featuresCounts
  
  # Define a fixed version that doesn't use SnowParam
  fixed_featuresCounts <- function(features,
                                   bam_dirs,
                                   strandness = c("unstrand", "1st_strand", "2nd_strand"),
                                   parallel = 1,
                                   yield_size = 5000000) {
    exist_indx <- file.exists(bam_dirs)
    
    if(any(!exist_indx)) {
      stop(paste0("Cannot find BAM file(s) of \n",
                  paste(bam_dirs[!exist_indx], collapse = ", "),
                  "\nplease re-check the input directories"))
    }
    
    strandness <- match.arg(strandness)
    
    # Only register MulticoreParam, avoid SnowParam
    BiocParallel::register(BiocParallel::SerialParam())
    if (parallel > 1) {
      suppressWarnings(BiocParallel::register(BiocParallel::MulticoreParam(workers = parallel)))
    }
    
    # Setup bam file list
    bam_lst = GenomicAlignments::BamFileList(file = bam_dirs, asMates = TRUE)
    GenomicAlignments::yieldSize(bam_lst) = yield_size
    
    # Count using summarizeOverlaps
    if(strandness == "1st_strand") {
      preprocess_func <- GenomicAlignments::readsReverse
    } else {
      preprocess_func <- NULL
    }
    
    se <- GenomicAlignments::summarizeOverlaps(
      features = features,
      reads = bam_lst,
      mode = "Union",
      inter.feature = FALSE,
      singleEnd = FALSE,
      preprocess.reads = preprocess_func,
      ignore.strand = strandness == "unstrand",
      fragments = TRUE
    )
    
    return(se)
  }
  
  # Replace the function in the exomepeak2 namespace
  tryCatch({
    unlockBinding("featuresCounts", asNamespace("exomepeak2"))
    assignInNamespace("featuresCounts", fixed_featuresCounts, ns="exomepeak2")
    lockBinding("featuresCounts", asNamespace("exomepeak2"))
    cat("Successfully patched exomepeak2::featuresCounts to avoid socket connections\n")
  }, error = function(e) {
    cat("Warning: Unable to patch exomepeak2::featuresCounts - ", e$message, "\n")
    cat("Will try to use environment settings instead\n")
  })
}


# Function to split space-separated paths into a character vector
parse_bam_paths <- function(path_string) {
    if (is.null(path_string)) return(NULL)
    # Split by space and remove any empty strings
    paths <- unlist(strsplit(path_string, ","))
    paths <- paths[paths != ""]
    
    # Verify all files exist
    for (path in paths) {
        if (!file.exists(path)) {
            stop(paste("BAM file not found:", path), call. = FALSE)
        }
    }
    
    return(paths)
}

# Parse BAM file paths
bam_ip <- parse_bam_paths(opt$bam_ip)
bam_input <- parse_bam_paths(opt$bam_input)
control_bam_ip <- parse_bam_paths(opt$control_bam_ip)
control_bam_input <- parse_bam_paths(opt$control_bam_input)

# Validate required inputs
if (is.null(opt$bam_ip) || is.null(opt$bam_input) || is.null(opt$gtf)) {
    print_help(opt)
    stop("Error: --bam_ip, --bam_input, and --gtf are required.", call. = FALSE)
}

# Print file information for debugging
cat("IP BAM files:", paste(bam_ip, collapse=", "), "\n")
cat("Input BAM files:", paste(bam_input, collapse=", "), "\n")
if (!is.null(control_bam_ip)) cat("Control IP BAM files:", paste(control_bam_ip, collapse=", "), "\n")
if (!is.null(control_bam_input)) cat("Control Input BAM files:", paste(control_bam_input, collapse=", "), "\n")

# Load BiocManager after ensuring it's installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", lib = user_lib, repos = "http://cran.rstudio.com/", quiet = TRUE)
}
suppressPackageStartupMessages(library(BiocManager))

# Install BSgenome package only if not NULL and not already installed
bsgenome <- NULL
if (!is.null(opt$bsgenome) && opt$bsgenome != "null") {
    if (!requireNamespace(opt$bsgenome, quietly = TRUE)) {
        cat("Installing ", opt$bsgenome, "...\n")
        BiocManager::install(opt$bsgenome, update = FALSE, ask = FALSE,lib = user_lib)
    }
    suppressPackageStartupMessages(library(opt$bsgenome, character.only = TRUE))
    bsgenome <- get(opt$bsgenome)
} else {
    message("Warning: --bsgenome is NULL or 'null'; skipping GC content correction.")
}

# set parallel workers
if (opt$cores > 1) {
    # Try to use MulticoreParam which doesn't rely on socket connections
    register(MulticoreParam(workers = opt$cores))
    cat("Using MulticoreParam with", opt$cores, "cores\n")
} else {
    register(SerialParam())
    cat("Using SerialParam (single core)\n")
}

# Run exomePeak2 peak calling
if(is.null(control_bam_ip) || is.null(control_bam_input)){  
    result <- exomePeak2(
        bam_ip          = bam_ip,      # IP sample BAM
        bam_input       = bam_input,   # Control/input BAM
        gff             = opt$gtf,         # Transcriptome annotation
        genome          = opt$bsgenome,    # BSgenome for GC bias correction
        parallel        = opt$cores,       # Use parallel processing 
        experiment_name = opt$outprefix    # Prefix for output files
    )
}else{
    result <- exomePeak2(
        bam_ip              = control_bam_ip,      # IP sample BAM
        bam_input           = control_bam_input,   # Control/input BAM
        bam_ip_treated      = bam_ip,              # IP sample BAM
        bam_input_treated   = bam_input,           # Control/input BAM
        gff                 = opt$gtf,                 # Transcriptome annotation
        genome              = opt$bsgenome,            # BSgenome for GC bias correction
        parallel            = opt$cores,               # Use parallel processing 
        experiment_name     = opt$outprefix            # Prefix for output files
    )
}

# Export results (BED and XLS are written automatically by exomePeak2)
# No additional export needed unless you want custom formats

cat("Peak calling completed.\n")
