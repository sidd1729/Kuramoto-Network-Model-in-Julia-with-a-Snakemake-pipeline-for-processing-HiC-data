
rule cool_feather:
    input: 
        script = "workflow/scripts/cool_feather.py",
        cool_file = "temporaries/cools/{sample}.mcool"
    output:
        "temporaries/raw_feather/{resolution}/{sample}_{resolution}_{chromosome_no}.feather"
    params:
        resolution = "{resolution}",
        out_base_dir = "temporaries/raw_feather/{resolution}"
    conda:
        "cooltools"
    shell:
        """
        python {input.script} {input.cool_file} {params.resolution} {params.out_base_dir}
        """ 
        
