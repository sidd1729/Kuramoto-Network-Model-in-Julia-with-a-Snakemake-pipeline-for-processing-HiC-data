using CSV
using Arrow
using DataFrames
read = CSV.read("/home/ksslab/Programs/dcHiC/DifferentialResult/Regen_100kb_timepoint/fdr_result/differential.intra_sample_group.pcQnm.bedGraph",DataFrame;delim = '\t')

filtered_df = select(read,["chr","start","end","control_100kb","hpi_12_100kb","dpi_2_100kb","dpi_4_100kb","dpi_7_100kb"])
println(first(filtered_df,100))
Arrow.write("/home/ksslab/Programs/dcHiC/DifferentialResult/Regen_100kb_timepoint/fdr_result/differential.intra_sample_group.pcQnm.feather",filtered_df)