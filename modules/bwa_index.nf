process BWA_INDEX {
    label 'process_low'
    input:
        path(reference)

    output:
        tuple path(reference), path("*.amb"), path('*.ann'), path("*.bwt"), path("*.pac"), path("*.sa")

    script:
        """
        bwa index ${reference} 
        """
}