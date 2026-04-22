
PARAMETER_START = 0.00
PARAMETER_STOP = 0.060
PARAMETER_STEP = 0.0003
AVERAGING_NO = 100
SIMULATION_TIME = 500.0

rule kuramoto_calculator:
    input:
        script = "workflow/scripts/ICE_processor.jl",
        temp_iced_mat = "temporaries/ice_feather/{resolution}/{sample}_{resolution}_{chromosome_no}_ICEd.feather"
        #temp_iced_mat = lambda wildcards: f"temporaries/ice_feather/{wildcards.resolution}/{wildcards.sample}_{wildcards.resolution}_{wildcards.chromosome_no}_ICEd.feather"
    output:
        kura  = "results/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_r_kura.feather",
        uni   = "results/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_r_uni.feather",
        k     = "results/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}_k.feather"
    params:
        window_size = "{window_size}",
        parameter_start = PARAMETER_START,
        parameter_stop = PARAMETER_STOP,
        parameter_step = PARAMETER_STEP,
        averaging_no = AVERAGING_NO,
        simulation_time = SIMULATION_TIME,
        base_output_path = "results/{resolution}/{window_size}/{sample}_{resolution}_{chromosome_no}_{window_size}"
    shell:
        """
        julia {input.script} {input.temp_iced_mat} {params.window_size} {params.parameter_start} {params.parameter_stop} {params.parameter_step} {params.averaging_no} {params.simulation_time} {params.base_output_path}
        """