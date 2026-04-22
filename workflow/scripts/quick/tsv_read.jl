using CSV
using DataFrames


df = CSV.read("/home/ksslab/Downloads/danRer11.refGene.gtf",DataFrame;delim = '\t')#,header = false)
println(first(df,100))