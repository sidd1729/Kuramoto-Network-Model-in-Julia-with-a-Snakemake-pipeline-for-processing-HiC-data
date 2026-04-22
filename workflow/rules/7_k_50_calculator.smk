PARALLEL_WORKERS = 1
rule k_50_calculator:
    input:
        script = "workflow/scripts/K_50_calculator.R",
        kura_path  = "results/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_r_kura.feather",
        uni_path   = "results/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_r_uni.feather",
        k_path     = "results/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_k.feather"
    output:
        k_50_kura_path = "results/{resolution}/k_50/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_kura_k_50.feather",
        k_50_uni_path = "results/{resolution}/k_50/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_uni_k_50.feather",
        perpet_kura_path = "results/{resolution}/perpet_1/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_kura_perpet_1.feather",
        perpet_uni_path =  "results/{resolution}/perpet_1/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_uni_perpet_1.feather"
    params:
        parallel_workers = PARALLEL_WORKERS
    shell:
        """
        Rscript {input.script} {input.k_path} {input.kura_path} {input.uni_path} {params.parallel_workers} {output.k_50_kura_path} {output.k_50_uni_path} {output.perpet_kura_path} {output.perpet_uni_path}
        """
    