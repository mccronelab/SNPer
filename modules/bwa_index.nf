process BWA_MEM {

    input:
        path(reference)

    output:
        tuple path(reference), path("*.amb"), path('*.ann'), path("*.bwt"), path("*.pac"), path("*.sa")

    script:
        """
        bwa index ${reference} 
        """
}