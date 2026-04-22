using DataFrames
using Arrow
using BedgraphFiles
using Bedgraph
using FileIO
using CSV

df = DataFrame(Arrow.Table("results/Overall_table/Overall_PC1_10000.feather"))
df_2 = DataFrame(Arrow.Table("results/Overall_table/Overall_PC2_10000.feather"))
df_3 = DataFrame(Arrow.Table("results/Overall_table/Overall_PC3_10000.feather")) 

println(first(df,5))
println(first(df[!,:chrom],5))
df[!,:chrom] .= map(x->"chr$(x)",df[!,:chrom])
println(first(df[!,:chrom],5))
df[!,:start] .= replace!(df[!,:start],missing => 0)
df[!,:end] .= replace!(df[!,:end],missing => 0)
df[!,:control] .= replace(df[!,:control],missing => 0)
df[!,:"12hpi"] .= replace(df[!,:"12hpi"],missing => 0)
df[!,:"2dpi"] .= replace(df[!,:"2dpi"],missing => 0)
df[!,:"4dpi"] .= replace(df[!,:"4dpi"],missing => 0)
df[!,:"7dpi"] .= replace(df[!,:"7dpi"],missing => 0)
println(first(df,5))

println(first(df_2,5))
println(first(df_2[!,:chrom],5))
df_2[!,:chrom] .= map(x->"chr$(x)",df_2[!,:chrom])
println(first(df_2[!,:chrom],5))
df_2[!,:start] .= replace!(df_2[!,:start],missing => 0)
df_2[!,:end] .= replace!(df_2[!,:end],missing => 0)
df_2[!,:control] .= replace(df_2[!,:control],missing => 0)
df_2[!,:"12hpi"] .= replace(df_2[!,:"12hpi"],missing => 0)
df_2[!,:"2dpi"] .= replace(df_2[!,:"2dpi"],missing => 0)
df_2[!,:"4dpi"] .= replace(df_2[!,:"4dpi"],missing => 0)
df_2[!,:"7dpi"] .= replace(df_2[!,:"7dpi"],missing => 0)
println(first(df_2,5))

println(first(df_3,5))
println(first(df_3[!,:chrom],5))
df_3[!,:chrom] .= map(x->"chr$(x)",df_3[!,:chrom])
println(first(df_3[!,:chrom],5))
df_3[!,:start] .= replace!(df_3[!,:start],missing => 0)
df_3[!,:end] .= replace!(df_3[!,:end],missing => 0)
df_3[!,:control] .= replace(df_3[!,:control],missing => 0)
df_3[!,:"12hpi"] .= replace(df_3[!,:"12hpi"],missing => 0)
df_3[!,:"2dpi"] .= replace(df_3[!,:"2dpi"],missing => 0)
df_3[!,:"4dpi"] .= replace(df_3[!,:"4dpi"],missing => 0)
df_3[!,:"7dpi"] .= replace(df_3[!,:"7dpi"],missing => 0)
println(first(df_3,5))

#println(size(df))'

#Control ---
records = []

#Chromosome 1 (PC2)
chrom_1_df = df_2[df_2[!,:chrom] .== "chr1",:]
for i in 1:size(chrom_1_df)[1]
push!(records,Bedgraph.Record(chrom_1_df[!,:chrom][i], chrom_1_df[!,:start][i], chrom_1_df[!,:end][i], chrom_1_df[!,:control][i]))
end

#Chromosome 2 (PC3)
chrom_2_df = df_3[df_3[!,:chrom] .== "chr2",:]
for i in 1:size(chrom_2_df)[1]
push!(records,Bedgraph.Record(chrom_2_df[!,:chrom][i], chrom_2_df[!,:start][i], chrom_2_df[!,:end][i], chrom_2_df[!,:control][i]))
end

#Chromosome 3 (PC1)
chrom_3_df = df[df[!,:chrom] .== "chr3",:]
for i in 1:size(chrom_3_df)[1]
push!(records,Bedgraph.Record(chrom_3_df[!,:chrom][i], chrom_3_df[!,:start][i], chrom_3_df[!,:end][i], chrom_3_df[!,:control][i]))
end

#Chromosome 4 (PC1)
chrom_4_df = df[df[!,:chrom] .== "chr4",:]
for i in 1:size(chrom_4_df)[1]
push!(records,Bedgraph.Record(chrom_4_df[!,:chrom][i], chrom_4_df[!,:start][i], chrom_4_df[!,:end][i], chrom_4_df[!,:control][i]))
end

#Chromosome 5 (PC3)
chrom_5_df = df_3[df_3[!,:chrom] .== "chr5",:]
for i in 1:size(chrom_5_df)[1]
push!(records,Bedgraph.Record(chrom_5_df[!,:chrom][i], chrom_5_df[!,:start][i], chrom_5_df[!,:end][i], chrom_5_df[!,:control][i]))
end
println(records)

#Chromosome 6 (PC3)
chrom_6_df = df_3[df_3[!,:chrom] .== "chr6",:]
for i in 1:size(chrom_6_df)[1]
push!(records,Bedgraph.Record(chrom_6_df[!,:chrom][i], chrom_6_df[!,:start][i], chrom_6_df[!,:end][i], chrom_6_df[!,:control][i]))
end

#Chromosome 7 (PC2)
chrom_7_df = df_2[df_2[!,:chrom] .== "chr7",:]
for i in 1:size(chrom_7_df)[1]
push!(records,Bedgraph.Record(chrom_7_df[!,:chrom][i], chrom_7_df[!,:start][i], chrom_7_df[!,:end][i], chrom_7_df[!,:control][i]))
end

#Chromosome 8 (PC3)
chrom_8_df = df_3[df_3[!,:chrom] .== "chr8",:]
for i in 1:size(chrom_8_df)[1]
push!(records,Bedgraph.Record(chrom_8_df[!,:chrom][i], chrom_8_df[!,:start][i], chrom_8_df[!,:end][i], chrom_8_df[!,:control][i]))
end

#Chromosome 9 (PC3)
chrom_9_df = df_3[df_3[!,:chrom] .== "chr9",:]
for i in 1:size(chrom_9_df)[1]
push!(records,Bedgraph.Record(chrom_9_df[!,:chrom][i], chrom_9_df[!,:start][i], chrom_9_df[!,:end][i], chrom_9_df[!,:control][i]))
end

#Chromosome 10 (PC2)
chrom_10_df = df_2[df_2[!,:chrom] .== "chr10",:]
for i in 1:size(chrom_10_df)[1]
push!(records,Bedgraph.Record(chrom_10_df[!,:chrom][i], chrom_10_df[!,:start][i], chrom_10_df[!,:end][i], chrom_10_df[!,:control][i]))
end

#Chromosome 11 (PC2)
chrom_11_df = df_2[df_2[!,:chrom] .== "chr11",:]
for i in 1:size(chrom_11_df)[1]
push!(records,Bedgraph.Record(chrom_11_df[!,:chrom][i], chrom_11_df[!,:start][i], chrom_11_df[!,:end][i], chrom_11_df[!,:control][i]))
end

#Chromosome 12 (PC3)
chrom_12_df = df_3[df_3[!,:chrom] .== "chr12",:]
for i in 1:size(chrom_12_df)[1]
push!(records,Bedgraph.Record(chrom_12_df[!,:chrom][i], chrom_12_df[!,:start][i], chrom_12_df[!,:end][i], chrom_12_df[!,:control][i]))
end

#Chromosome 13 (PC2)
chrom_13_df = df_2[df_2[!,:chrom] .== "chr13",:]
for i in 1:size(chrom_13_df)[1]
push!(records,Bedgraph.Record(chrom_13_df[!,:chrom][i], chrom_13_df[!,:start][i], chrom_13_df[!,:end][i], chrom_13_df[!,:control][i]))
end

#Chromosome 14 (PC3)
chrom_14_df = df_3[df_3[!,:chrom] .== "chr14",:]
for i in 1:size(chrom_14_df)[1]
push!(records,Bedgraph.Record(chrom_14_df[!,:chrom][i], chrom_14_df[!,:start][i], chrom_14_df[!,:end][i], chrom_14_df[!,:control][i]))
end

#Chromosome 15 (PC2)
chrom_15_df = df_2[df_2[!,:chrom] .== "chr15",:]
for i in 1:size(chrom_15_df)[1]
push!(records,Bedgraph.Record(chrom_15_df[!,:chrom][i], chrom_15_df[!,:start][i], chrom_15_df[!,:end][i], chrom_15_df[!,:control][i]))
end

#Chromosome 16 (PC3)
chrom_16_df = df_3[df_3[!,:chrom] .== "chr16",:]
for i in 1:size(chrom_16_df)[1]
push!(records,Bedgraph.Record(chrom_16_df[!,:chrom][i], chrom_16_df[!,:start][i], chrom_16_df[!,:end][i], chrom_16_df[!,:control][i]))
end

#Chromosome 17 (PC2)
chrom_17_df = df_2[df_2[!,:chrom] .== "chr17",:]
for i in 1:size(chrom_17_df)[1]
push!(records,Bedgraph.Record(chrom_17_df[!,:chrom][i], chrom_17_df[!,:start][i], chrom_17_df[!,:end][i], chrom_17_df[!,:control][i]))
end

#Chromosome 18 (PC2)
chrom_18_df = df_2[df_2[!,:chrom] .== "chr18",:]
for i in 1:size(chrom_18_df)[1]
push!(records,Bedgraph.Record(chrom_18_df[!,:chrom][i], chrom_18_df[!,:start][i], chrom_18_df[!,:end][i], chrom_18_df[!,:control][i]))
end

#Chromosome 19 (PC2)
chrom_19_df = df_2[df_2[!,:chrom] .== "chr19",:]
for i in 1:size(chrom_19_df)[1]
push!(records,Bedgraph.Record(chrom_19_df[!,:chrom][i], chrom_19_df[!,:start][i], chrom_19_df[!,:end][i], chrom_19_df[!,:control][i]))
end

#Chromosome 20 (PC2)
chrom_20_df = df_2[df_2[!,:chrom] .== "chr20",:]
for i in 1:size(chrom_20_df)[1]
push!(records,Bedgraph.Record(chrom_20_df[!,:chrom][i], chrom_20_df[!,:start][i], chrom_20_df[!,:end][i], chrom_20_df[!,:control][i]))
end

#Chromosome 21 (PC2)
chrom_21_df = df_2[df_2[!,:chrom] .== "chr21",:]
for i in 1:size(chrom_21_df)[1]
push!(records,Bedgraph.Record(chrom_21_df[!,:chrom][i], chrom_21_df[!,:start][i], chrom_21_df[!,:end][i], chrom_21_df[!,:control][i]))
end

#Chromosome 22 (PC2)
chrom_22_df = df_2[df_2[!,:chrom] .== "chr22",:]
for i in 1:size(chrom_22_df)[1]
push!(records,Bedgraph.Record(chrom_22_df[!,:chrom][i], chrom_22_df[!,:start][i], chrom_22_df[!,:end][i], chrom_22_df[!,:control][i]))
end

#Chromosome 23 (PC2)
chrom_23_df = df_2[df_2[!,:chrom] .== "chr23",:]
for i in 1:size(chrom_23_df)[1]
push!(records,Bedgraph.Record(chrom_23_df[!,:chrom][i], chrom_23_df[!,:start][i], chrom_23_df[!,:end][i], chrom_23_df[!,:control][i]))
end

#Chromosome 24 (PC3)
chrom_24_df = df_3[df_3[!,:chrom] .== "chr24",:]
for i in 1:size(chrom_24_df)[1]
push!(records,Bedgraph.Record(chrom_24_df[!,:chrom][i], chrom_24_df[!,:start][i], chrom_24_df[!,:end][i], chrom_24_df[!,:control][i]))
end

#Chromosome 25 (PC3)
chrom_25_df = df_3[df_3[!,:chrom] .== "chr25",:]
for i in 1:size(chrom_25_df)[1]
push!(records,Bedgraph.Record(chrom_25_df[!,:chrom][i], chrom_25_df[!,:start][i], chrom_25_df[!,:end][i], chrom_25_df[!,:control][i]))
end

println(records)
#=
records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:control][i]) for i in 1:size(df)[1]]
=#
bed_records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], i) for i in 1:size(df)[1]]
records_df = DataFrame(records)
bed_records_df = DataFrame(bed_records)
CSV.write("results/PCs_10kb_bedgraphs/tsv_control_PCs_10kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PCs_10kb_bedgraphs/t_control_PCs_10kb.bed",bed_records_df;delim='\t',header=false)
println(first(records_df,100))
sort!(records)
#println(records)
#sort!(records,by = x -> x.first)
println(propertynames(records))
save("results/PCs_10kb_bedgraphs/control_PCs_10kb.bedgraph",records)

#12 hpi ---
records = []

#Chromosome 1 (PC2)
chrom_1_df = df_2[df_2[!,:chrom] .== "chr1",:]
for i in 1:size(chrom_1_df)[1]
push!(records,Bedgraph.Record(chrom_1_df[!,:chrom][i], chrom_1_df[!,:start][i], chrom_1_df[!,:end][i], chrom_1_df[!,:"12hpi"][i]))
end

#Chromosome 2 (PC3)
chrom_2_df = df_3[df_3[!,:chrom] .== "chr2",:]
for i in 1:size(chrom_2_df)[1]
push!(records,Bedgraph.Record(chrom_2_df[!,:chrom][i], chrom_2_df[!,:start][i], chrom_2_df[!,:end][i], chrom_2_df[!,:"12hpi"][i]))
end

#Chromosome 3 (PC1)
chrom_3_df = df[df[!,:chrom] .== "chr3",:]
for i in 1:size(chrom_3_df)[1]
push!(records,Bedgraph.Record(chrom_3_df[!,:chrom][i], chrom_3_df[!,:start][i], chrom_3_df[!,:end][i], chrom_3_df[!,:"12hpi"][i]))
end

#Chromosome 4 (PC1)
chrom_4_df = df[df[!,:chrom] .== "chr4",:]
for i in 1:size(chrom_4_df)[1]
push!(records,Bedgraph.Record(chrom_4_df[!,:chrom][i], chrom_4_df[!,:start][i], chrom_4_df[!,:end][i], chrom_4_df[!,:"12hpi"][i]))
end

#Chromosome 5 (PC3)
chrom_5_df = df_3[df_3[!,:chrom] .== "chr5",:]
for i in 1:size(chrom_5_df)[1]
push!(records,Bedgraph.Record(chrom_5_df[!,:chrom][i], chrom_5_df[!,:start][i], chrom_5_df[!,:end][i], chrom_5_df[!,:"12hpi"][i]))
end
println(records)

#Chromosome 6 (PC3)
chrom_6_df = df_3[df_3[!,:chrom] .== "chr6",:]
for i in 1:size(chrom_6_df)[1]
push!(records,Bedgraph.Record(chrom_6_df[!,:chrom][i], chrom_6_df[!,:start][i], chrom_6_df[!,:end][i], chrom_6_df[!,:"12hpi"][i]))
end

#Chromosome 7 (PC2)
chrom_7_df = df_2[df_2[!,:chrom] .== "chr7",:]
for i in 1:size(chrom_7_df)[1]
push!(records,Bedgraph.Record(chrom_7_df[!,:chrom][i], chrom_7_df[!,:start][i], chrom_7_df[!,:end][i], chrom_7_df[!,:"12hpi"][i]))
end

#Chromosome 8 (PC3)
chrom_8_df = df_3[df_3[!,:chrom] .== "chr8",:]
for i in 1:size(chrom_8_df)[1]
push!(records,Bedgraph.Record(chrom_8_df[!,:chrom][i], chrom_8_df[!,:start][i], chrom_8_df[!,:end][i], chrom_8_df[!,:"12hpi"][i]))
end

#Chromosome 9 (PC3)
chrom_9_df = df_3[df_3[!,:chrom] .== "chr9",:]
for i in 1:size(chrom_9_df)[1]
push!(records,Bedgraph.Record(chrom_9_df[!,:chrom][i], chrom_9_df[!,:start][i], chrom_9_df[!,:end][i], chrom_9_df[!,:"12hpi"][i]))
end

#Chromosome 10 (PC2)
chrom_10_df = df_2[df_2[!,:chrom] .== "chr10",:]
for i in 1:size(chrom_10_df)[1]
push!(records,Bedgraph.Record(chrom_10_df[!,:chrom][i], chrom_10_df[!,:start][i], chrom_10_df[!,:end][i], chrom_10_df[!,:"12hpi"][i]))
end

#Chromosome 11 (PC2)
chrom_11_df = df_2[df_2[!,:chrom] .== "chr11",:]
for i in 1:size(chrom_11_df)[1]
push!(records,Bedgraph.Record(chrom_11_df[!,:chrom][i], chrom_11_df[!,:start][i], chrom_11_df[!,:end][i], chrom_11_df[!,:"12hpi"][i]))
end

#Chromosome 12 (PC3)
chrom_12_df = df_3[df_3[!,:chrom] .== "chr12",:]
for i in 1:size(chrom_12_df)[1]
push!(records,Bedgraph.Record(chrom_12_df[!,:chrom][i], chrom_12_df[!,:start][i], chrom_12_df[!,:end][i], chrom_12_df[!,:"12hpi"][i]))
end

#Chromosome 13 (PC2)
chrom_13_df = df_2[df_2[!,:chrom] .== "chr13",:]
for i in 1:size(chrom_13_df)[1]
push!(records,Bedgraph.Record(chrom_13_df[!,:chrom][i], chrom_13_df[!,:start][i], chrom_13_df[!,:end][i], chrom_13_df[!,:"12hpi"][i]))
end

#Chromosome 14 (PC3)
chrom_14_df = df_3[df_3[!,:chrom] .== "chr14",:]
for i in 1:size(chrom_14_df)[1]
push!(records,Bedgraph.Record(chrom_14_df[!,:chrom][i], chrom_14_df[!,:start][i], chrom_14_df[!,:end][i], chrom_14_df[!,:"12hpi"][i]))
end

#Chromosome 15 (PC2)
chrom_15_df = df_2[df_2[!,:chrom] .== "chr15",:]
for i in 1:size(chrom_15_df)[1]
push!(records,Bedgraph.Record(chrom_15_df[!,:chrom][i], chrom_15_df[!,:start][i], chrom_15_df[!,:end][i], chrom_15_df[!,:"12hpi"][i]))
end

#Chromosome 16 (PC3)
chrom_16_df = df_3[df_3[!,:chrom] .== "chr16",:]
for i in 1:size(chrom_16_df)[1]
push!(records,Bedgraph.Record(chrom_16_df[!,:chrom][i], chrom_16_df[!,:start][i], chrom_16_df[!,:end][i], chrom_16_df[!,:"12hpi"][i]))
end

#Chromosome 17 (PC2)
chrom_17_df = df_2[df_2[!,:chrom] .== "chr17",:]
for i in 1:size(chrom_17_df)[1]
push!(records,Bedgraph.Record(chrom_17_df[!,:chrom][i], chrom_17_df[!,:start][i], chrom_17_df[!,:end][i], chrom_17_df[!,:"12hpi"][i]))
end

#Chromosome 18 (PC2)
chrom_18_df = df_2[df_2[!,:chrom] .== "chr18",:]
for i in 1:size(chrom_18_df)[1]
push!(records,Bedgraph.Record(chrom_18_df[!,:chrom][i], chrom_18_df[!,:start][i], chrom_18_df[!,:end][i], chrom_18_df[!,:"12hpi"][i]))
end

#Chromosome 19 (PC2)
chrom_19_df = df_2[df_2[!,:chrom] .== "chr19",:]
for i in 1:size(chrom_19_df)[1]
push!(records,Bedgraph.Record(chrom_19_df[!,:chrom][i], chrom_19_df[!,:start][i], chrom_19_df[!,:end][i], chrom_19_df[!,:"12hpi"][i]))
end

#Chromosome 20 (PC2)
chrom_20_df = df_2[df_2[!,:chrom] .== "chr20",:]
for i in 1:size(chrom_20_df)[1]
push!(records,Bedgraph.Record(chrom_20_df[!,:chrom][i], chrom_20_df[!,:start][i], chrom_20_df[!,:end][i], chrom_20_df[!,:"12hpi"][i]))
end

#Chromosome 21 (PC2)
chrom_21_df = df_2[df_2[!,:chrom] .== "chr21",:]
for i in 1:size(chrom_21_df)[1]
push!(records,Bedgraph.Record(chrom_21_df[!,:chrom][i], chrom_21_df[!,:start][i], chrom_21_df[!,:end][i], chrom_21_df[!,:"12hpi"][i]))
end

#Chromosome 22 (PC2)
chrom_22_df = df_2[df_2[!,:chrom] .== "chr22",:]
for i in 1:size(chrom_22_df)[1]
push!(records,Bedgraph.Record(chrom_22_df[!,:chrom][i], chrom_22_df[!,:start][i], chrom_22_df[!,:end][i], chrom_22_df[!,:"12hpi"][i]))
end

#Chromosome 23 (PC2)
chrom_23_df = df_2[df_2[!,:chrom] .== "chr23",:]
for i in 1:size(chrom_23_df)[1]
push!(records,Bedgraph.Record(chrom_23_df[!,:chrom][i], chrom_23_df[!,:start][i], chrom_23_df[!,:end][i], chrom_23_df[!,:"12hpi"][i]))
end

#Chromosome 24 (PC3)
chrom_24_df = df_3[df_3[!,:chrom] .== "chr24",:]
for i in 1:size(chrom_24_df)[1]
push!(records,Bedgraph.Record(chrom_24_df[!,:chrom][i], chrom_24_df[!,:start][i], chrom_24_df[!,:end][i], chrom_24_df[!,:"12hpi"][i]))
end

#Chromosome 25 (PC3)
chrom_25_df = df_3[df_3[!,:chrom] .== "chr25",:]
for i in 1:size(chrom_25_df)[1]
push!(records,Bedgraph.Record(chrom_25_df[!,:chrom][i], chrom_25_df[!,:start][i], chrom_25_df[!,:end][i], chrom_25_df[!,:"12hpi"][i]))
end

println(records)
#=
records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:control][i]) for i in 1:size(df)[1]]
=#
bed_records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], i) for i in 1:size(df)[1]]
records_df = DataFrame(records)
bed_records_df = DataFrame(bed_records)
CSV.write("results/PCs_10kb_bedgraphs/tsv_12_hpi_PCs_10kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PCs_10kb_bedgraphs/t_12hpi_PCs_10kb.bed",bed_records_df;delim='\t',header=false)
println(first(records_df,100))
sort!(records)
#println(records)
#sort!(records,by = x -> x.first)
println(propertynames(records))
save("results/PCs_10kb_bedgraphs/12_hpi_PCs_10kb.bedgraph",records)

#2 dpi ---
records = []

#Chromosome 1 (PC2)
chrom_1_df = df_2[df_2[!,:chrom] .== "chr1",:]
for i in 1:size(chrom_1_df)[1]
push!(records,Bedgraph.Record(chrom_1_df[!,:chrom][i], chrom_1_df[!,:start][i], chrom_1_df[!,:end][i], chrom_1_df[!,:"2dpi"][i]))
end

#Chromosome 2 (PC3)
chrom_2_df = df_3[df_3[!,:chrom] .== "chr2",:]
for i in 1:size(chrom_2_df)[1]
push!(records,Bedgraph.Record(chrom_2_df[!,:chrom][i], chrom_2_df[!,:start][i], chrom_2_df[!,:end][i], chrom_2_df[!,:"2dpi"][i]))
end

#Chromosome 3 (PC1)
chrom_3_df = df[df[!,:chrom] .== "chr3",:]
for i in 1:size(chrom_3_df)[1]
push!(records,Bedgraph.Record(chrom_3_df[!,:chrom][i], chrom_3_df[!,:start][i], chrom_3_df[!,:end][i], chrom_3_df[!,:"2dpi"][i]))
end

#Chromosome 4 (PC1)
chrom_4_df = df[df[!,:chrom] .== "chr4",:]
for i in 1:size(chrom_4_df)[1]
push!(records,Bedgraph.Record(chrom_4_df[!,:chrom][i], chrom_4_df[!,:start][i], chrom_4_df[!,:end][i], chrom_4_df[!,:"2dpi"][i]))
end

#Chromosome 5 (PC3)
chrom_5_df = df_3[df_3[!,:chrom] .== "chr5",:]
for i in 1:size(chrom_5_df)[1]
push!(records,Bedgraph.Record(chrom_5_df[!,:chrom][i], chrom_5_df[!,:start][i], chrom_5_df[!,:end][i], chrom_5_df[!,:"2dpi"][i]))
end
println(records)

#Chromosome 6 (PC3)
chrom_6_df = df_3[df_3[!,:chrom] .== "chr6",:]
for i in 1:size(chrom_6_df)[1]
push!(records,Bedgraph.Record(chrom_6_df[!,:chrom][i], chrom_6_df[!,:start][i], chrom_6_df[!,:end][i], chrom_6_df[!,:"2dpi"][i]))
end

#Chromosome 7 (PC2)
chrom_7_df = df_2[df_2[!,:chrom] .== "chr7",:]
for i in 1:size(chrom_7_df)[1]
push!(records,Bedgraph.Record(chrom_7_df[!,:chrom][i], chrom_7_df[!,:start][i], chrom_7_df[!,:end][i], chrom_7_df[!,:"2dpi"][i]))
end

#Chromosome 8 (PC3)
chrom_8_df = df_3[df_3[!,:chrom] .== "chr8",:]
for i in 1:size(chrom_8_df)[1]
push!(records,Bedgraph.Record(chrom_8_df[!,:chrom][i], chrom_8_df[!,:start][i], chrom_8_df[!,:end][i], chrom_8_df[!,:"2dpi"][i]))
end

#Chromosome 9 (PC3)
chrom_9_df = df_3[df_3[!,:chrom] .== "chr9",:]
for i in 1:size(chrom_9_df)[1]
push!(records,Bedgraph.Record(chrom_9_df[!,:chrom][i], chrom_9_df[!,:start][i], chrom_9_df[!,:end][i], chrom_9_df[!,:"2dpi"][i]))
end

#Chromosome 10 (PC2)
chrom_10_df = df_2[df_2[!,:chrom] .== "chr10",:]
for i in 1:size(chrom_10_df)[1]
push!(records,Bedgraph.Record(chrom_10_df[!,:chrom][i], chrom_10_df[!,:start][i], chrom_10_df[!,:end][i], chrom_10_df[!,:"2dpi"][i]))
end

#Chromosome 11 (PC2)
chrom_11_df = df_2[df_2[!,:chrom] .== "chr11",:]
for i in 1:size(chrom_11_df)[1]
push!(records,Bedgraph.Record(chrom_11_df[!,:chrom][i], chrom_11_df[!,:start][i], chrom_11_df[!,:end][i], chrom_11_df[!,:"2dpi"][i]))
end

#Chromosome 12 (PC3)
chrom_12_df = df_3[df_3[!,:chrom] .== "chr12",:]
for i in 1:size(chrom_12_df)[1]
push!(records,Bedgraph.Record(chrom_12_df[!,:chrom][i], chrom_12_df[!,:start][i], chrom_12_df[!,:end][i], chrom_12_df[!,:"2dpi"][i]))
end

#Chromosome 13 (PC2)
chrom_13_df = df_2[df_2[!,:chrom] .== "chr13",:]
for i in 1:size(chrom_13_df)[1]
push!(records,Bedgraph.Record(chrom_13_df[!,:chrom][i], chrom_13_df[!,:start][i], chrom_13_df[!,:end][i], chrom_13_df[!,:"2dpi"][i]))
end

#Chromosome 14 (PC3)
chrom_14_df = df_3[df_3[!,:chrom] .== "chr14",:]
for i in 1:size(chrom_14_df)[1]
push!(records,Bedgraph.Record(chrom_14_df[!,:chrom][i], chrom_14_df[!,:start][i], chrom_14_df[!,:end][i], chrom_14_df[!,:"2dpi"][i]))
end

#Chromosome 15 (PC2)
chrom_15_df = df_2[df_2[!,:chrom] .== "chr15",:]
for i in 1:size(chrom_15_df)[1]
push!(records,Bedgraph.Record(chrom_15_df[!,:chrom][i], chrom_15_df[!,:start][i], chrom_15_df[!,:end][i], chrom_15_df[!,:"2dpi"][i]))
end

#Chromosome 16 (PC3)
chrom_16_df = df_3[df_3[!,:chrom] .== "chr16",:]
for i in 1:size(chrom_16_df)[1]
push!(records,Bedgraph.Record(chrom_16_df[!,:chrom][i], chrom_16_df[!,:start][i], chrom_16_df[!,:end][i], chrom_16_df[!,:"2dpi"][i]))
end

#Chromosome 17 (PC2)
chrom_17_df = df_2[df_2[!,:chrom] .== "chr17",:]
for i in 1:size(chrom_17_df)[1]
push!(records,Bedgraph.Record(chrom_17_df[!,:chrom][i], chrom_17_df[!,:start][i], chrom_17_df[!,:end][i], chrom_17_df[!,:"2dpi"][i]))
end

#Chromosome 18 (PC2)
chrom_18_df = df_2[df_2[!,:chrom] .== "chr18",:]
for i in 1:size(chrom_18_df)[1]
push!(records,Bedgraph.Record(chrom_18_df[!,:chrom][i], chrom_18_df[!,:start][i], chrom_18_df[!,:end][i], chrom_18_df[!,:"2dpi"][i]))
end

#Chromosome 19 (PC2)
chrom_19_df = df_2[df_2[!,:chrom] .== "chr19",:]
for i in 1:size(chrom_19_df)[1]
push!(records,Bedgraph.Record(chrom_19_df[!,:chrom][i], chrom_19_df[!,:start][i], chrom_19_df[!,:end][i], chrom_19_df[!,:"2dpi"][i]))
end

#Chromosome 20 (PC2)
chrom_20_df = df_2[df_2[!,:chrom] .== "chr20",:]
for i in 1:size(chrom_20_df)[1]
push!(records,Bedgraph.Record(chrom_20_df[!,:chrom][i], chrom_20_df[!,:start][i], chrom_20_df[!,:end][i], chrom_20_df[!,:"2dpi"][i]))
end

#Chromosome 21 (PC2)
chrom_21_df = df_2[df_2[!,:chrom] .== "chr21",:]
for i in 1:size(chrom_21_df)[1]
push!(records,Bedgraph.Record(chrom_21_df[!,:chrom][i], chrom_21_df[!,:start][i], chrom_21_df[!,:end][i], chrom_21_df[!,:"2dpi"][i]))
end

#Chromosome 22 (PC2)
chrom_22_df = df_2[df_2[!,:chrom] .== "chr22",:]
for i in 1:size(chrom_22_df)[1]
push!(records,Bedgraph.Record(chrom_22_df[!,:chrom][i], chrom_22_df[!,:start][i], chrom_22_df[!,:end][i], chrom_22_df[!,:"2dpi"][i]))
end

#Chromosome 23 (PC2)
chrom_23_df = df_2[df_2[!,:chrom] .== "chr23",:]
for i in 1:size(chrom_23_df)[1]
push!(records,Bedgraph.Record(chrom_23_df[!,:chrom][i], chrom_23_df[!,:start][i], chrom_23_df[!,:end][i], chrom_23_df[!,:"2dpi"][i]))
end

#Chromosome 24 (PC3)
chrom_24_df = df_3[df_3[!,:chrom] .== "chr24",:]
for i in 1:size(chrom_24_df)[1]
push!(records,Bedgraph.Record(chrom_24_df[!,:chrom][i], chrom_24_df[!,:start][i], chrom_24_df[!,:end][i], chrom_24_df[!,:"2dpi"][i]))
end

#Chromosome 25 (PC3)
chrom_25_df = df_3[df_3[!,:chrom] .== "chr25",:]
for i in 1:size(chrom_25_df)[1]
push!(records,Bedgraph.Record(chrom_25_df[!,:chrom][i], chrom_25_df[!,:start][i], chrom_25_df[!,:end][i], chrom_25_df[!,:"2dpi"][i]))
end

println(records)
#=
records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:control][i]) for i in 1:size(df)[1]]
=#
bed_records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], i) for i in 1:size(df)[1]]
records_df = DataFrame(records)
bed_records_df = DataFrame(bed_records)
CSV.write("results/PCs_10kb_bedgraphs/tsv_2_dpi_PCs_10kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PCs_10kb_bedgraphs/t_2dpi_PCs_10kb.bed",bed_records_df;delim='\t',header=false)
println(first(records_df,100))
sort!(records)
#println(records)
#sort!(records,by = x -> x.first)
println(propertynames(records))
save("results/PCs_10kb_bedgraphs/2_dpi_PCs_10kb.bedgraph",records)

#2 dpi ---
records = []

#Chromosome 1 (PC2)
chrom_1_df = df_2[df_2[!,:chrom] .== "chr1",:]
for i in 1:size(chrom_1_df)[1]
push!(records,Bedgraph.Record(chrom_1_df[!,:chrom][i], chrom_1_df[!,:start][i], chrom_1_df[!,:end][i], chrom_1_df[!,:"2dpi"][i]))
end

#Chromosome 2 (PC3)
chrom_2_df = df_3[df_3[!,:chrom] .== "chr2",:]
for i in 1:size(chrom_2_df)[1]
push!(records,Bedgraph.Record(chrom_2_df[!,:chrom][i], chrom_2_df[!,:start][i], chrom_2_df[!,:end][i], chrom_2_df[!,:"2dpi"][i]))
end

#Chromosome 3 (PC1)
chrom_3_df = df[df[!,:chrom] .== "chr3",:]
for i in 1:size(chrom_3_df)[1]
push!(records,Bedgraph.Record(chrom_3_df[!,:chrom][i], chrom_3_df[!,:start][i], chrom_3_df[!,:end][i], chrom_3_df[!,:"2dpi"][i]))
end

#Chromosome 4 (PC1)
chrom_4_df = df[df[!,:chrom] .== "chr4",:]
for i in 1:size(chrom_4_df)[1]
push!(records,Bedgraph.Record(chrom_4_df[!,:chrom][i], chrom_4_df[!,:start][i], chrom_4_df[!,:end][i], chrom_4_df[!,:"2dpi"][i]))
end

#Chromosome 5 (PC3)
chrom_5_df = df_3[df_3[!,:chrom] .== "chr5",:]
for i in 1:size(chrom_5_df)[1]
push!(records,Bedgraph.Record(chrom_5_df[!,:chrom][i], chrom_5_df[!,:start][i], chrom_5_df[!,:end][i], chrom_5_df[!,:"2dpi"][i]))
end
println(records)

#Chromosome 6 (PC3)
chrom_6_df = df_3[df_3[!,:chrom] .== "chr6",:]
for i in 1:size(chrom_6_df)[1]
push!(records,Bedgraph.Record(chrom_6_df[!,:chrom][i], chrom_6_df[!,:start][i], chrom_6_df[!,:end][i], chrom_6_df[!,:"2dpi"][i]))
end

#Chromosome 7 (PC2)
chrom_7_df = df_2[df_2[!,:chrom] .== "chr7",:]
for i in 1:size(chrom_7_df)[1]
push!(records,Bedgraph.Record(chrom_7_df[!,:chrom][i], chrom_7_df[!,:start][i], chrom_7_df[!,:end][i], chrom_7_df[!,:"2dpi"][i]))
end

#Chromosome 8 (PC3)
chrom_8_df = df_3[df_3[!,:chrom] .== "chr8",:]
for i in 1:size(chrom_8_df)[1]
push!(records,Bedgraph.Record(chrom_8_df[!,:chrom][i], chrom_8_df[!,:start][i], chrom_8_df[!,:end][i], chrom_8_df[!,:"2dpi"][i]))
end

#Chromosome 9 (PC3)
chrom_9_df = df_3[df_3[!,:chrom] .== "chr9",:]
for i in 1:size(chrom_9_df)[1]
push!(records,Bedgraph.Record(chrom_9_df[!,:chrom][i], chrom_9_df[!,:start][i], chrom_9_df[!,:end][i], chrom_9_df[!,:"2dpi"][i]))
end

#Chromosome 10 (PC2)
chrom_10_df = df_2[df_2[!,:chrom] .== "chr10",:]
for i in 1:size(chrom_10_df)[1]
push!(records,Bedgraph.Record(chrom_10_df[!,:chrom][i], chrom_10_df[!,:start][i], chrom_10_df[!,:end][i], chrom_10_df[!,:"2dpi"][i]))
end

#Chromosome 11 (PC2)
chrom_11_df = df_2[df_2[!,:chrom] .== "chr11",:]
for i in 1:size(chrom_11_df)[1]
push!(records,Bedgraph.Record(chrom_11_df[!,:chrom][i], chrom_11_df[!,:start][i], chrom_11_df[!,:end][i], chrom_11_df[!,:"2dpi"][i]))
end

#Chromosome 12 (PC3)
chrom_12_df = df_3[df_3[!,:chrom] .== "chr12",:]
for i in 1:size(chrom_12_df)[1]
push!(records,Bedgraph.Record(chrom_12_df[!,:chrom][i], chrom_12_df[!,:start][i], chrom_12_df[!,:end][i], chrom_12_df[!,:"2dpi"][i]))
end

#Chromosome 13 (PC2)
chrom_13_df = df_2[df_2[!,:chrom] .== "chr13",:]
for i in 1:size(chrom_13_df)[1]
push!(records,Bedgraph.Record(chrom_13_df[!,:chrom][i], chrom_13_df[!,:start][i], chrom_13_df[!,:end][i], chrom_13_df[!,:"2dpi"][i]))
end

#Chromosome 14 (PC3)
chrom_14_df = df_3[df_3[!,:chrom] .== "chr14",:]
for i in 1:size(chrom_14_df)[1]
push!(records,Bedgraph.Record(chrom_14_df[!,:chrom][i], chrom_14_df[!,:start][i], chrom_14_df[!,:end][i], chrom_14_df[!,:"2dpi"][i]))
end

#Chromosome 15 (PC2)
chrom_15_df = df_2[df_2[!,:chrom] .== "chr15",:]
for i in 1:size(chrom_15_df)[1]
push!(records,Bedgraph.Record(chrom_15_df[!,:chrom][i], chrom_15_df[!,:start][i], chrom_15_df[!,:end][i], chrom_15_df[!,:"2dpi"][i]))
end

#Chromosome 16 (PC3)
chrom_16_df = df_3[df_3[!,:chrom] .== "chr16",:]
for i in 1:size(chrom_16_df)[1]
push!(records,Bedgraph.Record(chrom_16_df[!,:chrom][i], chrom_16_df[!,:start][i], chrom_16_df[!,:end][i], chrom_16_df[!,:"2dpi"][i]))
end

#Chromosome 17 (PC2)
chrom_17_df = df_2[df_2[!,:chrom] .== "chr17",:]
for i in 1:size(chrom_17_df)[1]
push!(records,Bedgraph.Record(chrom_17_df[!,:chrom][i], chrom_17_df[!,:start][i], chrom_17_df[!,:end][i], chrom_17_df[!,:"2dpi"][i]))
end

#Chromosome 18 (PC2)
chrom_18_df = df_2[df_2[!,:chrom] .== "chr18",:]
for i in 1:size(chrom_18_df)[1]
push!(records,Bedgraph.Record(chrom_18_df[!,:chrom][i], chrom_18_df[!,:start][i], chrom_18_df[!,:end][i], chrom_18_df[!,:"2dpi"][i]))
end

#Chromosome 19 (PC2)
chrom_19_df = df_2[df_2[!,:chrom] .== "chr19",:]
for i in 1:size(chrom_19_df)[1]
push!(records,Bedgraph.Record(chrom_19_df[!,:chrom][i], chrom_19_df[!,:start][i], chrom_19_df[!,:end][i], chrom_19_df[!,:"2dpi"][i]))
end

#Chromosome 20 (PC2)
chrom_20_df = df_2[df_2[!,:chrom] .== "chr20",:]
for i in 1:size(chrom_20_df)[1]
push!(records,Bedgraph.Record(chrom_20_df[!,:chrom][i], chrom_20_df[!,:start][i], chrom_20_df[!,:end][i], chrom_20_df[!,:"2dpi"][i]))
end

#Chromosome 21 (PC2)
chrom_21_df = df_2[df_2[!,:chrom] .== "chr21",:]
for i in 1:size(chrom_21_df)[1]
push!(records,Bedgraph.Record(chrom_21_df[!,:chrom][i], chrom_21_df[!,:start][i], chrom_21_df[!,:end][i], chrom_21_df[!,:"2dpi"][i]))
end

#Chromosome 22 (PC2)
chrom_22_df = df_2[df_2[!,:chrom] .== "chr22",:]
for i in 1:size(chrom_22_df)[1]
push!(records,Bedgraph.Record(chrom_22_df[!,:chrom][i], chrom_22_df[!,:start][i], chrom_22_df[!,:end][i], chrom_22_df[!,:"2dpi"][i]))
end

#Chromosome 23 (PC2)
chrom_23_df = df_2[df_2[!,:chrom] .== "chr23",:]
for i in 1:size(chrom_23_df)[1]
push!(records,Bedgraph.Record(chrom_23_df[!,:chrom][i], chrom_23_df[!,:start][i], chrom_23_df[!,:end][i], chrom_23_df[!,:"2dpi"][i]))
end

#Chromosome 24 (PC3)
chrom_24_df = df_3[df_3[!,:chrom] .== "chr24",:]
for i in 1:size(chrom_24_df)[1]
push!(records,Bedgraph.Record(chrom_24_df[!,:chrom][i], chrom_24_df[!,:start][i], chrom_24_df[!,:end][i], chrom_24_df[!,:"2dpi"][i]))
end

#Chromosome 25 (PC3)
chrom_25_df = df_3[df_3[!,:chrom] .== "chr25",:]
for i in 1:size(chrom_25_df)[1]
push!(records,Bedgraph.Record(chrom_25_df[!,:chrom][i], chrom_25_df[!,:start][i], chrom_25_df[!,:end][i], chrom_25_df[!,:"2dpi"][i]))
end

println(records)
#=
records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:control][i]) for i in 1:size(df)[1]]
=#
bed_records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], i) for i in 1:size(df)[1]]
records_df = DataFrame(records)
bed_records_df = DataFrame(bed_records)
CSV.write("results/PCs_10kb_bedgraphs/tsv_2_dpi_PCs_10kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PCs_10kb_bedgraphs/t_2dpi_PCs_10kb.bed",bed_records_df;delim='\t',header=false)
println(first(records_df,100))
sort!(records)
#println(records)
#sort!(records,by = x -> x.first)
println(propertynames(records))
save("results/PCs_10kb_bedgraphs/2_dpi_PCs_10kb.bedgraph",records)

#4 dpi ---
records = []

#Chromosome 1 (PC2)
chrom_1_df = df_2[df_2[!,:chrom] .== "chr1",:]
for i in 1:size(chrom_1_df)[1]
push!(records,Bedgraph.Record(chrom_1_df[!,:chrom][i], chrom_1_df[!,:start][i], chrom_1_df[!,:end][i], chrom_1_df[!,:"4dpi"][i]))
end

#Chromosome 2 (PC3)
chrom_2_df = df_3[df_3[!,:chrom] .== "chr2",:]
for i in 1:size(chrom_2_df)[1]
push!(records,Bedgraph.Record(chrom_2_df[!,:chrom][i], chrom_2_df[!,:start][i], chrom_2_df[!,:end][i], chrom_2_df[!,:"4dpi"][i]))
end

#Chromosome 3 (PC1)
chrom_3_df = df[df[!,:chrom] .== "chr3",:]
for i in 1:size(chrom_3_df)[1]
push!(records,Bedgraph.Record(chrom_3_df[!,:chrom][i], chrom_3_df[!,:start][i], chrom_3_df[!,:end][i], chrom_3_df[!,:"4dpi"][i]))
end

#Chromosome 4 (PC1)
chrom_4_df = df[df[!,:chrom] .== "chr4",:]
for i in 1:size(chrom_4_df)[1]
push!(records,Bedgraph.Record(chrom_4_df[!,:chrom][i], chrom_4_df[!,:start][i], chrom_4_df[!,:end][i], chrom_4_df[!,:"4dpi"][i]))
end

#Chromosome 5 (PC3)
chrom_5_df = df_3[df_3[!,:chrom] .== "chr5",:]
for i in 1:size(chrom_5_df)[1]
push!(records,Bedgraph.Record(chrom_5_df[!,:chrom][i], chrom_5_df[!,:start][i], chrom_5_df[!,:end][i], chrom_5_df[!,:"4dpi"][i]))
end
println(records)

#Chromosome 6 (PC3)
chrom_6_df = df_3[df_3[!,:chrom] .== "chr6",:]
for i in 1:size(chrom_6_df)[1]
push!(records,Bedgraph.Record(chrom_6_df[!,:chrom][i], chrom_6_df[!,:start][i], chrom_6_df[!,:end][i], chrom_6_df[!,:"4dpi"][i]))
end

#Chromosome 7 (PC2)
chrom_7_df = df_2[df_2[!,:chrom] .== "chr7",:]
for i in 1:size(chrom_7_df)[1]
push!(records,Bedgraph.Record(chrom_7_df[!,:chrom][i], chrom_7_df[!,:start][i], chrom_7_df[!,:end][i], chrom_7_df[!,:"4dpi"][i]))
end

#Chromosome 8 (PC3)
chrom_8_df = df_3[df_3[!,:chrom] .== "chr8",:]
for i in 1:size(chrom_8_df)[1]
push!(records,Bedgraph.Record(chrom_8_df[!,:chrom][i], chrom_8_df[!,:start][i], chrom_8_df[!,:end][i], chrom_8_df[!,:"4dpi"][i]))
end

#Chromosome 9 (PC3)
chrom_9_df = df_3[df_3[!,:chrom] .== "chr9",:]
for i in 1:size(chrom_9_df)[1]
push!(records,Bedgraph.Record(chrom_9_df[!,:chrom][i], chrom_9_df[!,:start][i], chrom_9_df[!,:end][i], chrom_9_df[!,:"4dpi"][i]))
end

#Chromosome 10 (PC2)
chrom_10_df = df_2[df_2[!,:chrom] .== "chr10",:]
for i in 1:size(chrom_10_df)[1]
push!(records,Bedgraph.Record(chrom_10_df[!,:chrom][i], chrom_10_df[!,:start][i], chrom_10_df[!,:end][i], chrom_10_df[!,:"4dpi"][i]))
end

#Chromosome 11 (PC2)
chrom_11_df = df_2[df_2[!,:chrom] .== "chr11",:]
for i in 1:size(chrom_11_df)[1]
push!(records,Bedgraph.Record(chrom_11_df[!,:chrom][i], chrom_11_df[!,:start][i], chrom_11_df[!,:end][i], chrom_11_df[!,:"4dpi"][i]))
end

#Chromosome 12 (PC3)
chrom_12_df = df_3[df_3[!,:chrom] .== "chr12",:]
for i in 1:size(chrom_12_df)[1]
push!(records,Bedgraph.Record(chrom_12_df[!,:chrom][i], chrom_12_df[!,:start][i], chrom_12_df[!,:end][i], chrom_12_df[!,:"4dpi"][i]))
end

#Chromosome 13 (PC2)
chrom_13_df = df_2[df_2[!,:chrom] .== "chr13",:]
for i in 1:size(chrom_13_df)[1]
push!(records,Bedgraph.Record(chrom_13_df[!,:chrom][i], chrom_13_df[!,:start][i], chrom_13_df[!,:end][i], chrom_13_df[!,:"4dpi"][i]))
end

#Chromosome 14 (PC3)
chrom_14_df = df_3[df_3[!,:chrom] .== "chr14",:]
for i in 1:size(chrom_14_df)[1]
push!(records,Bedgraph.Record(chrom_14_df[!,:chrom][i], chrom_14_df[!,:start][i], chrom_14_df[!,:end][i], chrom_14_df[!,:"4dpi"][i]))
end

#Chromosome 15 (PC2)
chrom_15_df = df_2[df_2[!,:chrom] .== "chr15",:]
for i in 1:size(chrom_15_df)[1]
push!(records,Bedgraph.Record(chrom_15_df[!,:chrom][i], chrom_15_df[!,:start][i], chrom_15_df[!,:end][i], chrom_15_df[!,:"4dpi"][i]))
end

#Chromosome 16 (PC3)
chrom_16_df = df_3[df_3[!,:chrom] .== "chr16",:]
for i in 1:size(chrom_16_df)[1]
push!(records,Bedgraph.Record(chrom_16_df[!,:chrom][i], chrom_16_df[!,:start][i], chrom_16_df[!,:end][i], chrom_16_df[!,:"4dpi"][i]))
end

#Chromosome 17 (PC2)
chrom_17_df = df_2[df_2[!,:chrom] .== "chr17",:]
for i in 1:size(chrom_17_df)[1]
push!(records,Bedgraph.Record(chrom_17_df[!,:chrom][i], chrom_17_df[!,:start][i], chrom_17_df[!,:end][i], chrom_17_df[!,:"4dpi"][i]))
end

#Chromosome 18 (PC2)
chrom_18_df = df_2[df_2[!,:chrom] .== "chr18",:]
for i in 1:size(chrom_18_df)[1]
push!(records,Bedgraph.Record(chrom_18_df[!,:chrom][i], chrom_18_df[!,:start][i], chrom_18_df[!,:end][i], chrom_18_df[!,:"4dpi"][i]))
end

#Chromosome 19 (PC2)
chrom_19_df = df_2[df_2[!,:chrom] .== "chr19",:]
for i in 1:size(chrom_19_df)[1]
push!(records,Bedgraph.Record(chrom_19_df[!,:chrom][i], chrom_19_df[!,:start][i], chrom_19_df[!,:end][i], chrom_19_df[!,:"4dpi"][i]))
end

#Chromosome 20 (PC2)
chrom_20_df = df_2[df_2[!,:chrom] .== "chr20",:]
for i in 1:size(chrom_20_df)[1]
push!(records,Bedgraph.Record(chrom_20_df[!,:chrom][i], chrom_20_df[!,:start][i], chrom_20_df[!,:end][i], chrom_20_df[!,:"4dpi"][i]))
end

#Chromosome 21 (PC2)
chrom_21_df = df_2[df_2[!,:chrom] .== "chr21",:]
for i in 1:size(chrom_21_df)[1]
push!(records,Bedgraph.Record(chrom_21_df[!,:chrom][i], chrom_21_df[!,:start][i], chrom_21_df[!,:end][i], chrom_21_df[!,:"4dpi"][i]))
end

#Chromosome 22 (PC2)
chrom_22_df = df_2[df_2[!,:chrom] .== "chr22",:]
for i in 1:size(chrom_22_df)[1]
push!(records,Bedgraph.Record(chrom_22_df[!,:chrom][i], chrom_22_df[!,:start][i], chrom_22_df[!,:end][i], chrom_22_df[!,:"4dpi"][i]))
end

#Chromosome 23 (PC2)
chrom_23_df = df_2[df_2[!,:chrom] .== "chr23",:]
for i in 1:size(chrom_23_df)[1]
push!(records,Bedgraph.Record(chrom_23_df[!,:chrom][i], chrom_23_df[!,:start][i], chrom_23_df[!,:end][i], chrom_23_df[!,:"4dpi"][i]))
end

#Chromosome 24 (PC3)
chrom_24_df = df_3[df_3[!,:chrom] .== "chr24",:]
for i in 1:size(chrom_24_df)[1]
push!(records,Bedgraph.Record(chrom_24_df[!,:chrom][i], chrom_24_df[!,:start][i], chrom_24_df[!,:end][i], chrom_24_df[!,:"4dpi"][i]))
end

#Chromosome 25 (PC3)
chrom_25_df = df_3[df_3[!,:chrom] .== "chr25",:]
for i in 1:size(chrom_25_df)[1]
push!(records,Bedgraph.Record(chrom_25_df[!,:chrom][i], chrom_25_df[!,:start][i], chrom_25_df[!,:end][i], chrom_25_df[!,:"4dpi"][i]))
end

println(records)
#=
records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:control][i]) for i in 1:size(df)[1]]
=#
bed_records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], i) for i in 1:size(df)[1]]
records_df = DataFrame(records)
bed_records_df = DataFrame(bed_records)
CSV.write("results/PCs_10kb_bedgraphs/tsv_4_dpi_PCs_10kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PCs_10kb_bedgraphs/t_4dpi_PCs_10kb.bed",bed_records_df;delim='\t',header=false)
println(first(records_df,100))
sort!(records)
#println(records)
#sort!(records,by = x -> x.first)
println(propertynames(records))
save("results/PCs_10kb_bedgraphs/4_dpi_PCs_10kb.bedgraph",records)

#7 dpi ---
records = []

#Chromosome 1 (PC2)
chrom_1_df = df_2[df_2[!,:chrom] .== "chr1",:]
for i in 1:size(chrom_1_df)[1]
push!(records,Bedgraph.Record(chrom_1_df[!,:chrom][i], chrom_1_df[!,:start][i], chrom_1_df[!,:end][i], chrom_1_df[!,:"7dpi"][i]))
end

#Chromosome 2 (PC3)
chrom_2_df = df_3[df_3[!,:chrom] .== "chr2",:]
for i in 1:size(chrom_2_df)[1]
push!(records,Bedgraph.Record(chrom_2_df[!,:chrom][i], chrom_2_df[!,:start][i], chrom_2_df[!,:end][i], chrom_2_df[!,:"7dpi"][i]))
end

#Chromosome 3 (PC1)
chrom_3_df = df[df[!,:chrom] .== "chr3",:]
for i in 1:size(chrom_3_df)[1]
push!(records,Bedgraph.Record(chrom_3_df[!,:chrom][i], chrom_3_df[!,:start][i], chrom_3_df[!,:end][i], chrom_3_df[!,:"7dpi"][i]))
end

#Chromosome 4 (PC1)
chrom_4_df = df[df[!,:chrom] .== "chr4",:]
for i in 1:size(chrom_4_df)[1]
push!(records,Bedgraph.Record(chrom_4_df[!,:chrom][i], chrom_4_df[!,:start][i], chrom_4_df[!,:end][i], chrom_4_df[!,:"7dpi"][i]))
end

#Chromosome 5 (PC3)
chrom_5_df = df_3[df_3[!,:chrom] .== "chr5",:]
for i in 1:size(chrom_5_df)[1]
push!(records,Bedgraph.Record(chrom_5_df[!,:chrom][i], chrom_5_df[!,:start][i], chrom_5_df[!,:end][i], chrom_5_df[!,:"7dpi"][i]))
end
println(records)

#Chromosome 6 (PC3)
chrom_6_df = df_3[df_3[!,:chrom] .== "chr6",:]
for i in 1:size(chrom_6_df)[1]
push!(records,Bedgraph.Record(chrom_6_df[!,:chrom][i], chrom_6_df[!,:start][i], chrom_6_df[!,:end][i], chrom_6_df[!,:"7dpi"][i]))
end

#Chromosome 7 (PC2)
chrom_7_df = df_2[df_2[!,:chrom] .== "chr7",:]
for i in 1:size(chrom_7_df)[1]
push!(records,Bedgraph.Record(chrom_7_df[!,:chrom][i], chrom_7_df[!,:start][i], chrom_7_df[!,:end][i], chrom_7_df[!,:"7dpi"][i]))
end

#Chromosome 8 (PC3)
chrom_8_df = df_3[df_3[!,:chrom] .== "chr8",:]
for i in 1:size(chrom_8_df)[1]
push!(records,Bedgraph.Record(chrom_8_df[!,:chrom][i], chrom_8_df[!,:start][i], chrom_8_df[!,:end][i], chrom_8_df[!,:"7dpi"][i]))
end

#Chromosome 9 (PC3)
chrom_9_df = df_3[df_3[!,:chrom] .== "chr9",:]
for i in 1:size(chrom_9_df)[1]
push!(records,Bedgraph.Record(chrom_9_df[!,:chrom][i], chrom_9_df[!,:start][i], chrom_9_df[!,:end][i], chrom_9_df[!,:"7dpi"][i]))
end

#Chromosome 10 (PC2)
chrom_10_df = df_2[df_2[!,:chrom] .== "chr10",:]
for i in 1:size(chrom_10_df)[1]
push!(records,Bedgraph.Record(chrom_10_df[!,:chrom][i], chrom_10_df[!,:start][i], chrom_10_df[!,:end][i], chrom_10_df[!,:"7dpi"][i]))
end

#Chromosome 11 (PC2)
chrom_11_df = df_2[df_2[!,:chrom] .== "chr11",:]
for i in 1:size(chrom_11_df)[1]
push!(records,Bedgraph.Record(chrom_11_df[!,:chrom][i], chrom_11_df[!,:start][i], chrom_11_df[!,:end][i], chrom_11_df[!,:"7dpi"][i]))
end

#Chromosome 12 (PC3)
chrom_12_df = df_3[df_3[!,:chrom] .== "chr12",:]
for i in 1:size(chrom_12_df)[1]
push!(records,Bedgraph.Record(chrom_12_df[!,:chrom][i], chrom_12_df[!,:start][i], chrom_12_df[!,:end][i], chrom_12_df[!,:"7dpi"][i]))
end

#Chromosome 13 (PC2)
chrom_13_df = df_2[df_2[!,:chrom] .== "chr13",:]
for i in 1:size(chrom_13_df)[1]
push!(records,Bedgraph.Record(chrom_13_df[!,:chrom][i], chrom_13_df[!,:start][i], chrom_13_df[!,:end][i], chrom_13_df[!,:"7dpi"][i]))
end

#Chromosome 14 (PC3)
chrom_14_df = df_3[df_3[!,:chrom] .== "chr14",:]
for i in 1:size(chrom_14_df)[1]
push!(records,Bedgraph.Record(chrom_14_df[!,:chrom][i], chrom_14_df[!,:start][i], chrom_14_df[!,:end][i], chrom_14_df[!,:"7dpi"][i]))
end

#Chromosome 15 (PC2)
chrom_15_df = df_2[df_2[!,:chrom] .== "chr15",:]
for i in 1:size(chrom_15_df)[1]
push!(records,Bedgraph.Record(chrom_15_df[!,:chrom][i], chrom_15_df[!,:start][i], chrom_15_df[!,:end][i], chrom_15_df[!,:"7dpi"][i]))
end

#Chromosome 16 (PC3)
chrom_16_df = df_3[df_3[!,:chrom] .== "chr16",:]
for i in 1:size(chrom_16_df)[1]
push!(records,Bedgraph.Record(chrom_16_df[!,:chrom][i], chrom_16_df[!,:start][i], chrom_16_df[!,:end][i], chrom_16_df[!,:"7dpi"][i]))
end

#Chromosome 17 (PC2)
chrom_17_df = df_2[df_2[!,:chrom] .== "chr17",:]
for i in 1:size(chrom_17_df)[1]
push!(records,Bedgraph.Record(chrom_17_df[!,:chrom][i], chrom_17_df[!,:start][i], chrom_17_df[!,:end][i], chrom_17_df[!,:"7dpi"][i]))
end

#Chromosome 18 (PC2)
chrom_18_df = df_2[df_2[!,:chrom] .== "chr18",:]
for i in 1:size(chrom_18_df)[1]
push!(records,Bedgraph.Record(chrom_18_df[!,:chrom][i], chrom_18_df[!,:start][i], chrom_18_df[!,:end][i], chrom_18_df[!,:"7dpi"][i]))
end

#Chromosome 19 (PC2)
chrom_19_df = df_2[df_2[!,:chrom] .== "chr19",:]
for i in 1:size(chrom_19_df)[1]
push!(records,Bedgraph.Record(chrom_19_df[!,:chrom][i], chrom_19_df[!,:start][i], chrom_19_df[!,:end][i], chrom_19_df[!,:"7dpi"][i]))
end

#Chromosome 20 (PC2)
chrom_20_df = df_2[df_2[!,:chrom] .== "chr20",:]
for i in 1:size(chrom_20_df)[1]
push!(records,Bedgraph.Record(chrom_20_df[!,:chrom][i], chrom_20_df[!,:start][i], chrom_20_df[!,:end][i], chrom_20_df[!,:"7dpi"][i]))
end

#Chromosome 21 (PC2)
chrom_21_df = df_2[df_2[!,:chrom] .== "chr21",:]
for i in 1:size(chrom_21_df)[1]
push!(records,Bedgraph.Record(chrom_21_df[!,:chrom][i], chrom_21_df[!,:start][i], chrom_21_df[!,:end][i], chrom_21_df[!,:"7dpi"][i]))
end

#Chromosome 22 (PC2)
chrom_22_df = df_2[df_2[!,:chrom] .== "chr22",:]
for i in 1:size(chrom_22_df)[1]
push!(records,Bedgraph.Record(chrom_22_df[!,:chrom][i], chrom_22_df[!,:start][i], chrom_22_df[!,:end][i], chrom_22_df[!,:"7dpi"][i]))
end

#Chromosome 23 (PC2)
chrom_23_df = df_2[df_2[!,:chrom] .== "chr23",:]
for i in 1:size(chrom_23_df)[1]
push!(records,Bedgraph.Record(chrom_23_df[!,:chrom][i], chrom_23_df[!,:start][i], chrom_23_df[!,:end][i], chrom_23_df[!,:"7dpi"][i]))
end

#Chromosome 24 (PC3)
chrom_24_df = df_3[df_3[!,:chrom] .== "chr24",:]
for i in 1:size(chrom_24_df)[1]
push!(records,Bedgraph.Record(chrom_24_df[!,:chrom][i], chrom_24_df[!,:start][i], chrom_24_df[!,:end][i], chrom_24_df[!,:"7dpi"][i]))
end

#Chromosome 25 (PC3)
chrom_25_df = df_3[df_3[!,:chrom] .== "chr25",:]
for i in 1:size(chrom_25_df)[1]
push!(records,Bedgraph.Record(chrom_25_df[!,:chrom][i], chrom_25_df[!,:start][i], chrom_25_df[!,:end][i], chrom_25_df[!,:"7dpi"][i]))
end

println(records)
#=
records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], df[!,:control][i]) for i in 1:size(df)[1]]
=#
bed_records = [Bedgraph.Record(df[!,:chrom][i], df[!,:start][i], df[!,:end][i], i) for i in 1:size(df)[1]]
records_df = DataFrame(records)
bed_records_df = DataFrame(bed_records)
CSV.write("results/PCs_10kb_bedgraphs/tsv_7_dpi_PCs_10kb.bedgraph",records_df;delim='\t',header=false)
CSV.write("results/PCs_10kb_bedgraphs/t_7dpi_PCs_10kb.bed",bed_records_df;delim='\t',header=false)
println(first(records_df,100))
sort!(records)
#println(records)
#sort!(records,by = x -> x.first)
println(propertynames(records))
save("results/PCs_10kb_bedgraphs/7_dpi_PCs_10kb.bedgraph",records)

df = DataFrame(bed_path = ["/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs_10kb_bedgraphs/tsv_control_PCs_10kb.bedgraph", "/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs_10kb_bedgraphs/tsv_12_hpi_PCs_10kb.bedgraph",
"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs_10kb_bedgraphs/tsv_2_dpi_PCs_10kb.bedgraph","/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs_10kb_bedgraphs/tsv_4_dpi_PCs_10kb.bedgraph",
"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs_10kb_bedgraphs/tsv_7_dpi_PCs_10kb.bedgraph"], pc_type = ["intra", "intra","intra","intra","intra"], rep_name = ["t_control_PCs_10kb", "t_12hpi_PCs_10kb","t_2dpi_PCs_10kb","t_4dpi_PCs_10kb","t_7dpi_PCs_10kb"], sample_name = ["t_control","t_12hpi","t_2dpi","t_4dpi","t_7dpi"])
CSV.write("results/PCs_10kb_bedgraphs/dcHiC_paths.tsv",df;delim='\t',header=false)

#=
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
=#