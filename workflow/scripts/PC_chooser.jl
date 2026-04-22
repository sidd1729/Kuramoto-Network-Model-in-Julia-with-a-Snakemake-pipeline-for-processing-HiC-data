using Arrow
using DataFrames
using CSV
using Plots
using Measures
using Statistics

gcpt = CSV.read("/home/ksslab/Programs/dcHiC/danRer11_100000_goldenpathData/danRer11.GCpt.bedGraph",DataFrame,delim='\t',header=false)

cooltools_PC1 = DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Overall_table/Overall_PC1_100000.feather"))
cooltools_PC1[!,:start] .= replace!(cooltools_PC1[!,:start],missing => 0)
cooltools_PC1[!,:end] .= replace!(cooltools_PC1[!,:end],missing => 0)
cooltools_PC1[!,:control] .= replace(cooltools_PC1[!,:control],missing => 0)
cooltools_PC1[!,:"12hpi"] .= replace(cooltools_PC1[!,:"12hpi"],missing => 0)
cooltools_PC1[!,:"2dpi"] .= replace(cooltools_PC1[!,:"2dpi"],missing => 0)
cooltools_PC1[!,:"4dpi"] .= replace(cooltools_PC1[!,:"4dpi"],missing => 0)
cooltools_PC1[!,:"7dpi"] .= replace(cooltools_PC1[!,:"7dpi"],missing => 0)
cooltools_PC2 = DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Overall_table/Overall_PC2_100000.feather"))
cooltools_PC2[!,:start] .= replace!(cooltools_PC2[!,:start],missing => 0)
cooltools_PC2[!,:end] .= replace!(cooltools_PC2[!,:end],missing => 0)
cooltools_PC2[!,:control] .= replace(cooltools_PC2[!,:control],missing => 0)
cooltools_PC2[!,:"12hpi"] .= replace(cooltools_PC2[!,:"12hpi"],missing => 0)
cooltools_PC2[!,:"2dpi"] .= replace(cooltools_PC2[!,:"2dpi"],missing => 0)
cooltools_PC2[!,:"4dpi"] .= replace(cooltools_PC2[!,:"4dpi"],missing => 0)
cooltools_PC2[!,:"7dpi"] .= replace(cooltools_PC2[!,:"7dpi"],missing => 0)
cooltools_PC3 = DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Overall_table/Overall_PC3_100000.feather"))
cooltools_PC3[!,:start] .= replace!(cooltools_PC3[!,:start],missing => 0)
cooltools_PC3[!,:end] .= replace!(cooltools_PC3[!,:end],missing => 0)
cooltools_PC3[!,:control] .= replace(cooltools_PC3[!,:control],missing => 0)
cooltools_PC3[!,:"12hpi"] .= replace(cooltools_PC3[!,:"12hpi"],missing => 0)
cooltools_PC3[!,:"2dpi"] .= replace(cooltools_PC3[!,:"2dpi"],missing => 0)
cooltools_PC3[!,:"4dpi"] .= replace(cooltools_PC3[!,:"4dpi"],missing => 0)
cooltools_PC3[!,:"7dpi"] .= replace(cooltools_PC3[!,:"7dpi"],missing => 0)

data_name_vector = ["control","12hpi","2dpi","4dpi","7dpi"]
stacker_df = []

filtered_diff = CSV.read("/home/ksslab/Programs/dcHiC/DifferentialResult/Regen_100kb_timepoint/fdr_result/differential.intra_sample_group.Filtered.pcOri.bedGraph",DataFrame,delim='\t')
Arrow.write("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/diff_compartments/filtered_diff_compartments.feather",filtered_diff)

comb_bed = CSV.read("/home/ksslab/Programs/dcHiC/DifferentialResult/Regen_100kb_timepoint/fdr_result/differential.intra_sample_group.pcOri.bedGraph",DataFrame,delim='\t')
gene_density_df = DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Calculated_Gene_Density/danRer11_100000_gene_density.feather"))

for index in 1:5
    for chrom in 1:25
test_no = "$(chrom)"
test_chr = "chr$(test_no)"
index_to_look = index
cooltools_pc_1_df = cooltools_PC1[cooltools_PC1[!,:chrom] .== "$(test_no)",:]
cooltools_pc_2_df = cooltools_PC2[cooltools_PC2[!,:chrom] .== "$(test_no)",:]
cooltools_pc_3_df = cooltools_PC3[cooltools_PC3[!,:chrom] .== "$(test_no)",:]
pc_df = comb_bed[comb_bed[!,:chr] .== "$(test_chr)",:]#Dataframe[mask,:]
gc_df = gcpt[gcpt[!,:Column1] .== "$(test_chr)",:]
gden_df = gene_density_df[gene_density_df[!,:chrom] .== "$(test_chr)",:]

gd_track = gden_df[!,:gene_density]
gc_track = gc_df[!,:Column4]
gc_track_plot = plot(1:length(gc_track),gc_track,label = "GC content")
gd_track_plot = plot(1:length(gc_track),gd_track,label = "Gene Density",c =:brown)
pc_df_track = plot(1:length(pc_df[:,index_to_look+ 3]),pc_df[:,index_to_look+ 3],label = "dcHiC PC",c =:blue)
cooltools_pc1_track = plot(1:length(gc_track),cooltools_pc_1_df[:,index_to_look+3],label = "Cooltools PC1",c =:red)
cooltools_pc2_track = plot(1:length(gc_track),cooltools_pc_2_df[:,index_to_look+ 3],label = "Cooltools PC2",c =:green)
cooltools_pc3_track = plot(1:length(gc_track),cooltools_pc_3_df[:,index_to_look+ 3],label = "Cooltools PC3",c =:orange)
plot(gc_track_plot,gd_track_plot,pc_df_track,cooltools_pc1_track,cooltools_pc2_track,cooltools_pc3_track,link=:x, layout = Plots.grid(6, 1),size=(1000,1000),plot_title="PC choice of $(test_chr) of $(data_name_vector[index_to_look])",bottom_margin=[-7mm 0mm])#, widths = [0.5,0.5], heights=[0.5,0.5]))
savefig("/home/ksslab/Siddharth/Program_Directory/Jupyter/Julia/PC_choice/$(data_name_vector[index_to_look])_$(test_chr).png")
push!(stacker_df,["$(data_name_vector[index_to_look])","$(test_chr)",length(pc_df[:,index_to_look+3]) == length(gc_track) ? cor(pc_df[:,index_to_look+ 3],gc_track) : NaN ,cor(cooltools_pc_1_df[:,index_to_look+3],gc_track),cor(cooltools_pc_2_df[:,index_to_look+3],gc_track),cor(cooltools_pc_3_df[:,index_to_look+3],gc_track),length(pc_df[:,index_to_look+3]) == length(gd_track) ? cor(pc_df[:,index_to_look+ 3],gd_track) : NaN ,cor(cooltools_pc_1_df[:,index_to_look+3],gd_track),cor(cooltools_pc_2_df[:,index_to_look+3],gd_track),cor(cooltools_pc_3_df[:,index_to_look+3],gd_track)])
    end
end

# Source - https://stackoverflow.com/a/72957610
# Posted by Bogumił Kamiński, modified by community. See post 'Timeline' for change history
# Retrieved 2026-01-31, License - CC BY-SA 4.0

PC_info_df = DataFrame(mapreduce(permutedims, vcat, stacker_df), :auto)
rename!(PC_info_df,["timepoint","chr","dcHiC PC GC corr","cooltools PC1 GC corr","cooltools PC2 GC corr","cooltools PC3 GC corr","dcHiC PC GD corr","cooltools PC1 GD corr","cooltools PC2 GD corr","cooltools PC3 GD corr"])
CSV.write("/home/ksslab/Siddharth/Program_Directory/Jupyter/Julia/PC_choice/PC_choice_check/PC_choice_NaNs.csv",PC_info_df;sep = '\t')
Arrow.write("/home/ksslab/Siddharth/Program_Directory/Jupyter/Julia/PC_choice/PC_choice_check/PC_choice_NaNs.feather",PC_info_df)

PC_info_df[!,:"dcHiC PC GC corr"] .= replace!(PC_info_df[!,:"dcHiC PC GC corr"],NaN => 0)
PC_info_df[!,:"dcHiC PC GD corr"] .= replace!(PC_info_df[!,:"dcHiC PC GD corr"],NaN => 0)

PC_info_df. "dcHiC PC cscore" = 0.5.*(PC_info_df. "dcHiC PC GC corr" + PC_info_df. "dcHiC PC GD corr")
PC_info_df. "cooltools PC1 cscore" = 0.5.* (PC_info_df. "cooltools PC1 GC corr" + PC_info_df. "cooltools PC1 GD corr")
PC_info_df. "cooltools PC2 cscore" = 0.5.* (PC_info_df. "cooltools PC2 GC corr" + PC_info_df. "cooltools PC2 GD corr")
PC_info_df. "cooltools PC3 cscore" = 0.5.* (PC_info_df. "cooltools PC3 GC corr" + PC_info_df. "cooltools PC3 GD corr")

Chrom_score_col = []
Chrom_name_col = []
chromosomes = 1:25
for chr in chromosomes
   chrom = "chr$(chr)"
   push!(Chrom_name_col,chrom)
   chrom_filter_df = PC_info_df[PC_info_df[!,:"chr"] .== chrom,:] 
   push!(Chrom_score_col,(sum(chrom_filter_df. "dcHiC PC cscore" ),sum(chrom_filter_df. "cooltools PC1 cscore" ), sum(chrom_filter_df. "cooltools PC2 cscore" ), sum(chrom_filter_df. "cooltools PC3 cscore" ) ))
end
for i in 1 : size(PC_info_df)[1]-length(Chrom_name_col)
    push!(Chrom_name_col,"")
    push!(Chrom_score_col,"")
end
#println
PC_info_df. "Chromosome name" = Chrom_name_col
PC_info_df. "Chromosome score" = Chrom_score_col
CSV.write("/home/ksslab/Siddharth/Program_Directory/Jupyter/Julia/PC_choice/PC_choice_check/PC_choice_zeros.csv",PC_info_df;sep = '\t')
Arrow.write("/home/ksslab/Siddharth/Program_Directory/Jupyter/Julia/PC_choice/PC_choice_check/PC_choice_zeros.feather",PC_info_df)

#transform!(df, [:df_prop[3], :df_prop[4]] => ((x, y) -> x + y) => :Z)
#println(names(df))
#println(df. "dcHiC")