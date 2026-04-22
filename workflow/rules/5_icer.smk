#Add another rule parsing the input filename and storing in a temporary text file for later use, i.e storing final filename or naming the intermediate file

rule icer:
    input:
        script = "workflow/scripts/ICE_normalization.py",
        input_mat = "temporaries/raw_feather/{resolution}/{sample}_{resolution}_{chromosome_no}.feather"
    output:
        "temporaries/ice_feather/{resolution}/{sample}_{resolution}_{chromosome_no}_ICEd.feather"
    conda:
        "Kuramoto_scripts_env"
    shell:
        """
        python {input.script} {input.input_mat} {output}
        """