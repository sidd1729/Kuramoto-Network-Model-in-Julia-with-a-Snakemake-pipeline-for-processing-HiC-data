#Function for calculating Gene Density
using CSV
using DataFrames
using Arrow
content = ARGS[1]
chrom_sizes = ARGS[2]
tsv_reader = CSV.read(content,DataFrame;delim='\t')
chrom_sizes = CSV.read(chrom_sizes,DataFrame;delim='\t',header = false)
#tsv_reader[!,:txStart] .= parse.(Int, tsv_reader[!,:txStart])
println(first(chrom_sizes,25))
bin_size = parse(Int,ARGS[3])#100000
println(first(tsv_reader,100))
row_list = [] #["chrom","start","end","gene_density"]
for chrom in 1:25
    counter = 0
    bin_list  = []
    gene_density_list = []
    chrom_size = chrom_sizes[chrom_sizes[!,:Column1] .== "chr$(chrom)",:][!,:Column2][1]
    slice_to_consider = tsv_reader[tsv_reader[!,:chrom] .== "chr$(chrom)",:]
    println(chrom_size)
    for i in 0:bin_size:chrom_size
    push!(bin_list,(i,i+bin_size))
        #filtered = j.loc[((j['txStart'] > i) & (j['txEnd'] < i+bin_size)) | ((j['txStart'] > i) & (j['txStart'] < i+bin_size) & (j['txEnd'] < i+2*bin_size))]# Select entries lying in a bin, current ignoring those spanning between bins
    filtered = slice_to_consider[(slice_to_consider[!,:txStart] .> i) .& (slice_to_consider[!,:txStart] .< i+bin_size),:]
    unique_name_in_filtered =  unique(filtered[!,:name2])
        #print(filtered,i,i+bin_size)
    push!(gene_density_list,length(unique_name_in_filtered))
    push!(row_list,["chr$(chrom)",i,i+bin_size,length(unique_name_in_filtered)])
    end
    #println(gene_density_list)
end
gene_density_df = DataFrame(mapreduce(permutedims, vcat, row_list),:auto)
println(first(gene_density_df,5))
rename!(gene_density_df,["chrom","start","end","gene_density"])
println(names(gene_density_df))
println(first(gene_density_df,5))
CSV.write("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Calculated_Gene_Density/danRer11_$(bin_size)_gene_density.csv",gene_density_df)
Arrow.write("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Calculated_Gene_Density/danRer11_$(bin_size)_gene_density.feather",gene_density_df)
    #=
    for i in range(0,number_characters_in_fasta_file,bin_size)
        bin_list.append((i,i+bin_size))
        #filtered = j.loc[((j['txStart'] > i) & (j['txEnd'] < i+bin_size)) | ((j['txStart'] > i) & (j['txStart'] < i+bin_size) & (j['txEnd'] < i+2*bin_size))]# Select entries lying in a bin, current ignoring those spanning between bins
        filtered = j.loc[(j['txStart'] > i) & (j['txStart'] < i+bin_size)]
        unique_name_in_filtered =  set(filtered.loc[:,'name2'])
        #print(filtered,i,i+bin_size)
        gene_density_list.append(len(unique_name_in_filtered))
    end
    end
    #plot_list = [i for i in range(0,len(gene_density_list))]
    return (gene_density_list)
end
=#
#=
function RefSeq_Gene_Density_Calc(filename,bin_size ,number_characters_in_fasta_file)
    tsv_reader = CSV.read(content;delimiter='\t')
    counter = 0
    row_list = []
    for row in tsv_reader :
            #del row[len(row)-1]
            for i in range(0,len(row)):
                 row[i] = str(row[i])
            counter += 1
            row_list.append(row)
    j = pd.DataFrame(row_list)
    j.columns = j.iloc[0,:] 
    j.drop(0)
    j['txStart'] = pd.to_numeric(j['txStart'], downcast='integer', errors='coerce')
    j['txEnd'] = pd.to_numeric(j['txEnd'], downcast='integer', errors='coerce')
    bin_list = []
    gene_density_list = []
    for i in range(0,number_characters_in_fasta_file,bin_size):
        bin_list.append((i,i+bin_size))
        #filtered = j.loc[((j['txStart'] > i) & (j['txEnd'] < i+bin_size)) | ((j['txStart'] > i) & (j['txStart'] < i+bin_size) & (j['txEnd'] < i+2*bin_size))]# Select entries lying in a bin, current ignoring those spanning between bins
        filtered = j.loc[(j['txStart'] > i) & (j['txStart'] < i+bin_size)]
        unique_name_in_filtered =  set(filtered.loc[:,'name2'])
        #print(filtered,i,i+bin_size)
        gene_density_list.append(len(unique_name_in_filtered))
    end
    #plot_list = [i for i in range(0,len(gene_density_list))]
    return (gene_density_list)
end
=#