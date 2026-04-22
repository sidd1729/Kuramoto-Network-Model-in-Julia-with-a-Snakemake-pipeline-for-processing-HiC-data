#using Plots;
using Arrow
using LinearAlgebra
using DataFrames
using Combinatorics
using Measures
using Graphs
using GraphRecipes
using Printf
using Statistics
using OffsetArrays
using CairoMakie

using CSV
DLR = CSV.read("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/GM12878/1/results/ldr/chr1_ldr.bedGraph",DataFrame; header = 0)
Insulation = DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/GM12878/1/results/insul/4DNFITRVKRPA_50000_chr1_insul.feather"))
CTCF = CSV.read("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/GM12878/CTCF/GM12878_Chr1_CTCF_50kb.wig",DataFrame;delim='\t')
K4Me3 = CSV.read("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/GM12878/H3K4Me3/GM12878_Chr1_H3K4Me3_50kb.wig",DataFrame;delim='\t')
K27Ac = CSV.read("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/GM12878/H3K27Ac/GM12878_Chr1_H3K27Ac_50kb.wig",DataFrame;delim='\t')
K27Me3 = CSV.read("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/GM12878/H3K27Me3/GM12878_Chr1_H3K27Me3_50kb.wig",DataFrame;delim='\t')
PC1 = DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/4DNFITRVKRPA_50000_PC1.feather"))
PC2 = DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/4DNFITRVKRPA_50000_PC2.feather"))
PC3 = DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/4DNFITRVKRPA_50000_PC3.feather"))

PC3_sorted = PC3[1:length(Insulation[!,2]),:]
PC2_sorted = PC2[1:length(Insulation[!,2]),:]
PC1_sorted = PC1[1:length(Insulation[!,2]),:]

#push!(DLR, ["chr1",248950000,249000000,missing])

set_diff_start = setdiff(DLR[!,2],CTCF[!,2])
#set_diff_end = setdiff(DLR[!,3],CTCF[!,3])
for i in 1:length(set_diff_start)
    push!(CTCF,["chr1",set_diff_start[i],set_diff_start[i]+50000,0])
end

CTCF_sorted = sort(CTCF,:Column2)

deleteat!(CTCF_sorted, [4981])
CTCF_sorted

K4_setdiff = setdiff(CTCF_sorted[!,2],K4Me3[!,2])
K27_Ac_setdiff = setdiff(CTCF_sorted[!,2],K27Ac[!,2])
K27_Me3_setdiff = setdiff(CTCF_sorted[!,2],K27Me3[!,2])

for i in 1:length(K4_setdiff)
    push!(K4Me3,["chr1",K4_setdiff[i],K4_setdiff[i]+50000,0])
end
for i in 1:length(K27_Ac_setdiff)
    push!(K27Ac,["chr1",K27_Ac_setdiff[i],K27_Ac_setdiff[i]+50000,0])
end
for i in 1:length(K27_Me3_setdiff)
    push!(K27Me3,["chr1",K27_Me3_setdiff[i],K27_Me3_setdiff[i]+50000,0])
end

K4Me3_sorted = sort(K4Me3,:Column2)
deleteat!(K4Me3_sorted, [4981])
K27Ac_sorted = sort(K27Ac, :Column2)
deleteat!(K27Ac_sorted, [4981])
K27Me3_sorted = sort(K27Me3, :Column2)
deleteat!(K27Me3_sorted, [4981])

function z_normalization(arr)
    return map(x->(x-mean(arr))/std(arr),arr)
end

window_size = 5
df3 = Array(DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/GM12878/1/results/1_kura_k_50.feather")))
df4 = Array(DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/GM12878/1/results/1_kura_perpet_1.feather")))
trial = map(x-> x[1][1],df3)
trial_perpet = map(x-> x[1],df4)
trial_perpet = Int.(trial_perpet[trial_perpet .!= -1.0])
perpet_range = setdiff(1:length(trial),trial_perpet)
perpet_range_offset = (map(x -> x+floor((window_size+1)/2), perpet_range))

df5 = Matrix(DataFrame(Arrow.Table("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/GM12878/50000/ICEd/4DNFITRVKRPA_50000_ICEd_chr1.feather")))
df5[diagind(df5)] .= 0
log_df5 = log.(df5)
log_df5[diagind(log_df5)] .= 0 
bin_line = parse(Int,ARGS[1])
#gene_name = "cntf"#"gapdh"
plot_lims = (bin_line-15,bin_line+15)
#plot_1 = plot(perpet_range_offset,trial[perpet_range])
#vline!([bin_line] ,c =:black,label=false)
#plot_2 = heatmap(log.(df5),colorbar=false,c = cgrad(:Spectral,rev=true),ylims = plot_lims,xformatter=Returns(""))#,aspect_ratio=:equal)
#vline!([bin_line] ,c =:black,label=false)
#plot_3 = plot(1:length(Insulation[!,2]),K4Me3_sorted[!,2])
#plot_4 = plot(1:length(Insulation[!,2]),K27Ac_sorted[!,2])
#plot_5 = plot(1:length(Insulation[!,2]),K27Me3_sorted[!,2])
#plot(plot_1,plot_2,plot_3,plot_4,plot_5,link=:x,layout = Plots.grid(5, 1,heights=[0.10,0.6,0.10,0.10,0.10]),xlims = plot_lims,size=(800,800),plot_title="test",margin = 0mm,titlefont = font(3))#,widths = [0.5,0.5], heights=[0.3,0.7])
            #savefig("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/k_50_pairwise_plots/$(resolution)/k_50_pairwise_$(sample_1)_$(sample_2)_chr$(chromosome)_$(resolution).png")
            #savefig("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/k_50_pairwise_plots/$(resolution)/z_normed/k_50_pairwise_$(sample_1)_$(sample_2)_chr$(chromosome)_$(resolution)_z_normed.png")
#savefig("/home/ksslab/Siddharth/Program_Directory/Jupyter/Julia/GM12878/chr_1_fig_$(plot_lims[1])_$(plot_lims[2]).png")


f = Figure(size = (1000,1000)) #figure_padding = 5

ax = []
for i in 1:9, j in 1:1
    push!(ax,Axis(f[i, j]))
end

#Label(f[0, :], text = "First Supertitle", fontsize = 20)
#Label(f[-1, :], text = "Second Supertitle", fontsize = 30)
#Label(f[-2, :], text = "Third Supertitle", fontsize = 40)
lines!(f[1,1],perpet_range_offset,trial[perpet_range])
heatmap!(f[2,1],log_df5,colormap = Reverse(:Spectral))#,aspect_ratio = :equal)
lines!(f[3,1],1:length(Insulation[!,2]),K4Me3_sorted[!,4])
lines!(f[4,1],1:length(Insulation[!,2]),K27Ac_sorted[!,4])
lines!(f[5,1],1:length(Insulation[!,2]),K27Me3_sorted[!,4])
lines!(f[6,1],1:length(Insulation[!,2]),coalesce.(DLR[!,4],0))
lines!(f[7,1],1:length(Insulation[!,2]),coalesce.(PC1_sorted[!,1],0))
lines!(f[8,1],1:length(Insulation[!,2]),coalesce.(PC2_sorted[!,1],0))
lines!(f[9,1],1:length(Insulation[!,2]),coalesce.(PC3_sorted[!,1],0))
#linkxaxes!(f[1,1],f[2,1],f[3,1],f[4,1],f[5,1])
xlims!(ax[1],plot_lims[1],plot_lims[2])
vlines!(ax[1],bin_line,color = :black)
#ylims!(ax[1],plot_lims[1],plot_lims[2])
xlims!(ax[2],plot_lims[1],plot_lims[2])
ylims!(ax[2],plot_lims[1],plot_lims[2])
vlines!(ax[2],bin_line,color = :black)

xlims!(ax[3],plot_lims[1],plot_lims[2])
#ylims!(ax[3],plot_lims[1],plot_lims[2])
vlines!(ax[3],bin_line,color = :black)

xlims!(ax[4],plot_lims[1],plot_lims[2])
#ylims!(ax[4],plot_lims[1],plot_lims[2])
vlines!(ax[4],bin_line,color = :black)

xlims!(ax[5],plot_lims[1],plot_lims[2])
#ylims!(ax[5],plot_lims[1],plot_lims[2])
vlines!(ax[5],bin_line,color = :black)

xlims!(ax[6],plot_lims[1],plot_lims[2])
#ylims!(ax[5],plot_lims[1],plot_lims[2])
vlines!(ax[6],bin_line,color = :black)

xlims!(ax[7],plot_lims[1],plot_lims[2])
#ylims!(ax[5],plot_lims[1],plot_lims[2])
vlines!(ax[7],bin_line,color = :black)

xlims!(ax[8],plot_lims[1],plot_lims[2])
#ylims!(ax[5],plot_lims[1],plot_lims[2])
vlines!(ax[8],bin_line,color = :black)

xlims!(ax[9],plot_lims[1],plot_lims[2])
#ylims!(ax[5],plot_lims[1],plot_lims[2])
vlines!(ax[9],bin_line,color = :black)

rowgap!(f.layout,0)
rowsize!(f.layout,2,Auto(7))
#rowsize!(f.layout,Relative(1.0))
colgap!(f.layout,0)
#resize_to_layout!(f)
save("/home/ksslab/Siddharth/Program_Directory/Jupyter/Julia/GM12878/test_chr_1_fig_$(plot_lims[1])_$(plot_lims[2]).png",f)#,px_per_unit = 100)