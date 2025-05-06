process CHECKCONTRAST {
    input:
    tuple val(meta), val(reads)
    val contrast 

    output:
    tuple val(treatment_group), val(control_group)

    script:
    '''
    def groups = meta.collect { it['group'] }.toSet()
    def contrast_pairs = []
    def errors = []

    if (file(contrast).exists()) {
        file(contrast).text.readLines().each { line ->
            def pair = line.trim().split('_vs_')
            if (pair.size() == 2 && pair[0] in groups && pair[1] in groups) {
                contrast_pairs << pair
            } else {
                errors << "Invalid contrast: ${line}"
            }
        }
    } else {
        def pair = contrast.trim().split('_vs_')
        if (pair.size() == 2 && pair[0] in groups && pair[1] in groups) {
            contrast_pairs << pair
        } else {
            errors << "Invalid contrast format: ${contrast}"
        }
    }

    if (errors) {
        println "Contrast validation failed:"
        errors.each { println " - ${it}" }
        exit 1
    } else {
        contrast_pairs.each { pair ->
            println "Valid contrast: ${pair[0]} vs ${pair[1]}"
            emit(pair[0], pair[1])
        }
    }
    '''
}