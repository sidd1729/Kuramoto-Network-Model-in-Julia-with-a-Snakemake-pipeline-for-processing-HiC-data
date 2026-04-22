FASTA_FILE_PATH = "temporaries/genome/danRer11.fa"

rule PC_saddle_calc:
    input:
        script = "workflow/scripts/PCA_saddle.py",
        cool_file = "temporaries/cools/{sample}.mcool"
    output:
        PC1 = "results/PCs/{resolution}/{sample}_{resolution}_PC1.feather",
        PC2 = "results/PCs/{resolution}/{sample}_{resolution}_PC2.feather",
        PC3 = "results/PCs/{resolution}/{sample}_{resolution}_PC3.feather",
        Saddleplot = "results/saddleplots/{resolution}/{sample}_{resolution}_saddleplot.png"
     
    params:
        resolution = "{resolution}",
        fasta_file_path = FASTA_FILE_PATH,
        out_base_dir_PCs = "results/PCs/{resolution}",
        out_base_dir_saddleplots = "results/saddleplots/{resolution}"
    conda:
        "cooltools"
    shell:
        """
        python {input.script} {input.cool_file} {params.resolution} {params.fasta_file_path} {params.out_base_dir_PCs} {params.out_base_dir_saddleplots}
        """ 

#cooler balance {input.cool_file}::/resolutions/{params.resolution}
    

