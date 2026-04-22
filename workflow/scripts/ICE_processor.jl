using Arrow
using DataFrames
using BenchmarkTools
using DifferentialEquations
using Distributions
using LinearAlgebra
using Statistics
using ProgressBars
#/home/ksslab/Siddharth/Temp_Dir/ICEd_mat

#function submatrix_reader(matrix::Matrix{Float64},submatrix_shape::Tuple{Int64, Int64},stride::Int8 = 1) #,stride_type = "diagonal"
function submatrix_reader(matrix,submatrix_shape,stride = 1,offset = (0,0),stride_type = "diagonal")
    submatrix_no = size(matrix)[1] - submatrix_shape[1] + 1
    submatrices=Matrix{Float64}[]
    for i in 1:submatrix_no
        submatrix = matrix[i+offset[1]:i+offset[1]+submatrix_shape[1]-1,i+offset[2]:i+offset[2]+submatrix_shape[2]-1]
        push!(submatrices,submatrix)
    end
    return submatrices
end

function kuramoto_preprocessor(path_var = ARGS[1],window_size = ARGS[2],parameter_start = ARGS[3],parameter_stop = ARGS[4],parameter_step = ARGS[5],averaging_no = ARGS[6],simulation_time = ARGS[7])
    mat_iced = Matrix(DataFrame(Arrow.Table(path_var)))
    window_size = parse(Int,window_size)
    submatrices = submatrix_reader(mat_iced,(window_size,window_size))
    par_start = parse(Float64,parameter_start)
    par_step = parse(Float64,parameter_step)
    par_stop = parse(Float64,parameter_stop)
    avg_no = parse(Int,averaging_no)
    sim_time = parse(Float64,simulation_time)
    return (submatrices,size(submatrices)[1],par_start:par_step:par_stop,avg_no,window_size,sim_time)
end

#Kuramoto Model Calculations

#Part 1 of model calculations
function vectorization_multiplication!(mult_matrix,inp_vector,temp_mat)
    @inbounds temp_mat .=  sum(mult_matrix .* (sin.((reshape(inp_vector,1,:)) .- (reshape(inp_vector,:,1)))),dims=2) #Using similar tricks to numpy broadcasting to have matrix multiplication isntead of loops and ideally improve efficiency
    return temp_mat
end

#Part 2 of model calculations
function broadcasting_Calculation!(du,omega,K,W,u,W_norm_fac,temp_mat)
    #if 0 in W_norm_fac
        #W_norm_fac .+= 0.0000001
    #end
    @inbounds temp_mat .=  (omega .+ (K .* (vectorization_multiplication!(W,u,temp_mat)./W_norm_fac)))
    return du.=temp_mat
end

#Intermediate calculations to get appropriate entries from weight matrix for kuramoto calculation
function W_norm_fac_calc(Matrix)
    W_norm_fac = sum(Matrix,dims = 2)
    #if 0 in W_norm_fac
        #W_norm_fac .+= 0.0000001
    #end
    return W_norm_fac
end

function W_overall_sum_calc(W_norm_fac)
    W_overall_sum = sum(W_norm_fac,dims = 1)
    return W_overall_sum[1]
end

#Combining part 1 and part 2 to get overall model caluclations
function kuramoto_network_model!(du,u,p,t)
  broadcasting_Calculation!(du,p.ω,p.K,p.W,u,p.W_norm_fac,p.temp_mat)#,p.temp_mat_broadcasting)
end

#Calculating Kuramoto Order Parameter
function kuramoto_order_parameter_calc_final(phase_solutions_mat)
    return mean(abs.(sum(exp.(im.*phase_solutions_mat),dims=1))/size(phase_solutions_mat,1))
end

#Calculating Universal Order Parameter ()

#Part 1 of universal order parameter calculations:
function mean_cos_calc(phase_solutions_mat)
    return Statistics.mean(cos.(reshape(phase_solutions_mat,size(phase_solutions_mat,1),1,size(phase_solutions_mat,2)) .- reshape(phase_solutions_mat,1,size(phase_solutions_mat,1),size(phase_solutions_mat,2))),dims=3)
end

#Part 2 of universal order parameter calculations:
function uni_param_fin_step(Adj_mat,mean_cos,W_overall_sum)
    return sum(Adj_mat .* @view mean_cos[:,:,1])/W_overall_sum
end

#Combined universal order parameter calculations
function universal_order_parameter_calculation_modified(phase_solutions_mat,Adj_mat,W_overall_sum)
        mean_cos = mean_cos_calc(phase_solutions_mat)
        uni_param  = uni_param_fin_step(Adj_mat,mean_cos,W_overall_sum)
        return uni_param
end

#Combining the Model calculations, Kuramoto order parameter and Universal order parameter calculations into 1 function:
function kuramoto_all_calc(u0,tspan,par)
    # Define the problem
    prob = ODEProblem(kuramoto_network_model!, u0, tspan, par)
    # Choose a solver, for stiff problems you can use AutoVern7(Rodas5()) #AutoTsit5(Rosenbrock23(autodiff = false))
    sol = solve(prob, AutoTsit5(Rosenbrock23(autodiff = false)),dt=0.01,saveat = tspan[2] - 40.0:10:tspan[2] ,save_everystep = false)#,abstol=1.49012e-8, reltol=1.49012e-8)
    #println(sol)
    phase_solutions_mat = convert(Array,sol)
    Kura_order = kuramoto_order_parameter_calc_final(phase_solutions_mat)
    Uni_order = universal_order_parameter_calculation_modified(phase_solutions_mat,par.W,par.W_overall_sum)
    return (Kura_order,Uni_order,sol)
end

function kuramoto_repeated_all_calc(index=8,range = 5:0.03:27,submatrix_list = submatrices, repetition_rate = 100,window_size = 5,end_time = 500.0,len_range = len_range)
len_range = len_range
Temp_storer_array = zeros(len_range)
Temp_storer_array_uni = zeros(len_range)
K_overall_list = range
n=window_size
tspan = (0.0, end_time)
for i in 1:repetition_rate
    for j in 1:len_range
    # Define initial conditions
    u0 = rand(Uniform(0,2*pi), n) 
    W = submatrix_list[index]
    #println(W)
    W[diagind(W)] .= 0
    #println(W)
    W_norm_fac = W_norm_fac_calc(W)
    W_overall_sum = W_overall_sum_calc(W_norm_fac)
    temp_mat = Array{Float64}(undef,n,1)
    par = (K = K_overall_list[j], n = n, ω = rand(Normal(1,0.1), n),W = W,W_norm_fac = W_norm_fac,W_overall_sum = W_overall_sum,temp_mat =temp_mat)
    Kura_order,Uni_order,sol = kuramoto_all_calc(u0,tspan,par)
    Temp_storer_array[j]+=Kura_order
    Temp_storer_array_uni[j]+=Uni_order
    end
end
    return (Temp_storer_array./repetition_rate,Temp_storer_array_uni./repetition_rate)
end

function Overall_calc_runner()
    submatrices,submatrices_size,par_range,avg_no,window_size,sim_time  = kuramoto_preprocessor()
    len_range = length(par_range)
    prealloc_arr_1 = Array{Float64,2}(undef,len_range,submatrices_size)
    prealloc_arr_2 = Array{Float64,2}(undef,len_range,submatrices_size)
    for k in tqdm(1:submatrices_size)
        prealloc_arr_1[1:len_range,k],prealloc_arr_2[1:len_range,k] = kuramoto_repeated_all_calc(k,par_range,submatrices,avg_no,window_size,sim_time,len_range)
    end
    return (prealloc_arr_1,prealloc_arr_2,par_range)
end

#println(kuramoto_preprocessor())
#bm = @benchmark Overall_calc_runner()
Out1 , Out2, par_range = Overall_calc_runner()

Arrow.write("$(ARGS[8])_r_kura.feather",DataFrame(Out1,:auto))
Arrow.write("$(ARGS[8])_r_uni.feather",DataFrame(Out2,:auto))
Arrow.write("$(ARGS[8])_k.feather",DataFrame(data = collect(par_range)))
#@show bm