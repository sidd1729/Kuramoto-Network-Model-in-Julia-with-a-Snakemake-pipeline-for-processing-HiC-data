rule bam2p:
    input:
    
    shell:
        """
        bam2pairs 
        """