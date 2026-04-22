using CSV
using DataFrames

df = DataFrame(mat = ["control_50000.matrix","hpi_12_50000.matrix","dpi_2_50000.matrix","dpi_4_50000.matrix","dpi_7_50000.matrix"],
bed = ["control_50000_abs.bed","hpi_12_50000_abs.bed","dpi_2_50000_abs.bed","dpi_4_50000_abs.bed","dpi_7_50000_abs.bed"],
replicate_prefix = ["control_50kb_R","hpi_12_50kb_R","dpi_2_50kb_R","dpi_4_50kb_R","dpi_7_50kb_R"],
experiment_prefix = ["control_50kb","hpi_12_50kb","dpi_2_50kb","dpi_4_50kb","dpi_7_50kb"])

println(first(df,100))

CSV.write("/home/ksslab/Programs/dcHiC/Regeneration_50kb_input.txt",df;delim='\t',header=false)