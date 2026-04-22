using DataFrames
using Arrow
using BedgraphFiles
using Bedgraph
using FileIO
using CSV

df = DataFrame(Arrow.Table("results/Overall_table/Overall_PC1_500000.feather"))
println(first(df,5))
println(first(df[!,:chrom],5))
df[!,:chrom] .= map(x->"chr$(x)",df[!,:chrom])
println(first(df[!,:chrom],5))
println(first(df,5))
df[!,:start] .= replace!(df[!,:start],missing => 0)
df[!,:end] .= replace!(df[!,:end],missing => 0)
df[!,:control] .= replace(df[!,:control],missing => 0)
df[!,:"12hpi"] .= replace(df[!,:"12hpi"],missing => 0)
df[!,:"2dpi"] .= replace(df[!,:"2dpi"],missing => 0)
df[!,:"4dpi"] .= replace(df[!,:"4dpi"],missing => 0)
df[!,:"7dpi"] .= replace(df[!,:"7dpi"],missing => 0)

#println(size(df))
records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:control][i]) for i in 1:size(df)[1]]
bed_records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], i) for i in 1:size(df)[1]]
records_df = DataFrame(records)
bed_records_df = DataFrame(bed_records)
CSV.write("results/PC1_500kb_bedgraphs/tsv_control_PC1_500kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PC1_500kb_bedgraphs/t_control_PC1_500kb.bed",bed_records_df;delim='\t',header=false)
println(first(records_df,100))
sort!(records)
#println(records)
#sort!(records,by = x -> x.first)
println(propertynames(records))
save("results/PC1_500kb_bedgraphs/control_PC1_500kb.bedgraph",records)

records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:"12hpi"][i]) for i in 1:size(df)[1]]
records_df = DataFrame(records)
CSV.write("results/PC1_500kb_bedgraphs/tsv_12_hpi_PC1_500kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PC1_500kb_bedgraphs/t_12hpi_PC1_500kb.bed",bed_records_df;delim='\t',header=false)
sort!(records)
save("results/PC1_500kb_bedgraphs/12_hpi_PC1_500kb.bedgraph",records)

records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:"2dpi"][i]) for i in 1:size(df)[1]]
records_df = DataFrame(records)
CSV.write("results/PC1_500kb_bedgraphs/tsv_2_dpi_PC1_500kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PC1_500kb_bedgraphs/t_2dpi_PC1_500kb.bed",bed_records_df;delim='\t',header=false)
sort!(records)
save("results/PC1_500kb_bedgraphs/2_dpi_PC1_500kb.bedgraph",records)

records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:"4dpi"][i]) for i in 1:size(df)[1]]
records_df = DataFrame(records)
CSV.write("results/PC1_500kb_bedgraphs/tsv_4_dpi_PC1_500kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PC1_500kb_bedgraphs/t_4dpi_PC1_500kb.bed",bed_records_df;delim='\t',header=false)
sort!(records)
save("results/PC1_500kb_bedgraphs/4_dpi_PC1_500kb.bedgraph",records)

records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:"7dpi"][i]) for i in 1:size(df)[1]]
records_df = DataFrame(records)
CSV.write("results/PC1_500kb_bedgraphs/tsv_7_dpi_PC1_500kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PC1_500kb_bedgraphs/t_7dpi_PC1_500kb.bed",bed_records_df;delim='\t',header=false)

sort!(records)
save("results/PC1_500kb_bedgraphs/7_dpi_PC1_500kb.bedgraph",records)

df = DataFrame(bed_path = ["/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PC1_500kb_bedgraphs/tsv_control_PC1_500kb.bedgraph", "/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PC1_500kb_bedgraphs/tsv_12_hpi_PC1_500kb.bedgraph",
"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PC1_500kb_bedgraphs/tsv_2_dpi_PC1_500kb.bedgraph","/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PC1_500kb_bedgraphs/tsv_4_dpi_PC1_500kb.bedgraph",
"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PC1_500kb_bedgraphs/tsv_7_dpi_PC1_500kb.bedgraph"], pc_type = ["intra", "intra","intra","intra","intra"], rep_name = ["t_control_PC1_500kb", "t_12hpi_PC1_500kb","t_2dpi_PC1_500kb","t_4dpi_PC1_500kb","t_7dpi_PC1_500kb"], sample_name = ["t_control","t_12hpi","t_2dpi","t_4dpi","t_7dpi"])
CSV.write("results/PC1_500kb_bedgraphs/dcHiC_paths.tsv",df;delim='\t',header=false)