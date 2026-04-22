#WINDOW_SIZE = 5
rule insul_expected_calc:
    input: 
        script = "workflow/scripts/Insul_O_E_calc.py",
        cool_file = "temporaries/cools/{sample}.mcool"
    output:
        insulation = "results/insulation/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_insul.feather",
        expected = "temporaries/expected_feather/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_exp.feather",
        observed_minus_expected = "results/o_min_e/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_o_min_e.feather",
        observed_by_expected = "results/o_by_e/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_o_by_e.feather"
    params:
        resolution = "{resolution}",
        window_size = "{window_size}",
        out_base_dir_insul = "results/insulation/{resolution}",
        out_base_dir_expect = "temporaries/expected_feather/{resolution}",
        out_base_dir_o_min_e = "results/o_min_e/{resolution}",
        out_base_dir_o_by_e = "results/o_by_e/{resolution}"
    conda:
        "cooltools"
    shell:
        """
        cooler balance --force {input.cool_file}::resolutions/{params.resolution}
        python {input.script} {input.cool_file} {params.resolution} {params.window_size} {params.out_base_dir_insul} {params.out_base_dir_expect} {params.out_base_dir_o_min_e} {params.out_base_dir_o_by_e}
        """ 