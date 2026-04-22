using Arrow
using DataFrames
using Plots
resolution = 50000
#df = DataFrame(Arrow.Table("results/Overall_table/Overall_PC1_$(resolution).feather"))
df = DataFrame(Arrow.Table("results/Overall_table/Overall_PC1.feather"))

println(first(df,100))
for col in ["control","12hpi","2dpi","4dpi","7dpi"]
df[!,col][coalesce.(df[!,col] .< 0,true)] .= -1
df[!,col][coalesce.(df[!,col] .> 0,true)] .= +1
end
println(first(df,100))

#Arrow.write("results/Overall_table/Overall_compartment_$(resolution).feather",df)
#println(first(DataFrame(Arrow.Table("results/Overall_table/Overall_compartment_$(resolution).feather")),100))
Arrow.write("results/Overall_table/Overall_compartment.feather",df)
println(first(DataFrame(Arrow.Table("results/Overall_table/Overall_compartment.feather")),100))
heatmap(Matrix(df[!,["control","12hpi","2dpi","4dpi","7dpi"]]))
#savefig("/home/ksslab/Siddharth/Program_Directory/Jupyter/Julia/test.png")