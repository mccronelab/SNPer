process BWA_MEM {
    container "staphb/bwa:latest"

    input:
        path(reference)

    output:
        tuple path(reference), path("*.amb"), path('*.ann'), path("*.bwt"), path("*.pac"), path("*.sa")

    script:
        """
        bwa index ${reference} 
        """
}