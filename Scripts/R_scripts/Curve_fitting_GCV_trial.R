
library(DescTools) #Used for an O(n) search to verify correctness of binary search
library(psych) #Used for tr function for matrix trace
library(purrr) #Used for map function
#Implemented functions:

#GCV

#Implementation following tutorial in python from here:
#https://www.geeksforgeeks.org/machine-learning/mastering-generalized-cross-validation-gcv-theory-applications-and-best-practices/
GCV <- function(curve_fitting_function,output_extraction="",function_args,data_collection_list,input_values,observed_values,lower_lim=0,upper_lim=1,parameter_no=100){
  gcv_scores = c()
  working_interval<-seq(from=lower_lim,to=upper_lim,length.out = parameter_no)
  for (smooth in working_interval) {
    #print(smooth)
    #print(paste(deparse(substitute(curve_fitting_function)),function_args,sep=""))
    curve_output<-eval(parse(text =paste(deparse(substitute(curve_fitting_function)),function_args,sep="")))
    #print(curve_output)
    given_string <- paste("curve_output",output_extraction,sep="")
    parsed_string <- parse(text = given_string)
    output <- eval(parsed_string)
    if(typeof(observed_values) == "list"){
      observed_values =unlist(observed_values)
    }
    if (typeof(input_values) == "list"){
      input_values = unlist(input_values)
    }
    residual_sum_squares <- sum((output - observed_values)^2)
    mat_input = matrix(input_values) #Leaving as column matrix (dim : n x 1) so X^{T}X will be  a scalar as every column will depict a predictor variable
    X = mat_input #Values had to be in form of vector for correct operations rather than a list
    X_T = t(mat_input)
    #print(c(typeof(X),typeof(X_T)))
    X_T_X_sum_lambda_I_inv = solve( (X_T %*% X) + (smooth * diag(dim(X)[2])) ) #Diag creates the identity matrix
    hat_matrix = X %*% X_T_X_sum_lambda_I_inv %*% X_T
    trace_hat_matrix = tr(hat_matrix)
    gcv_score = residual_sum_squares/ (1- (trace_hat_matrix/dim(X)[1]))^2
    gcv_scores <-append(gcv_scores,gcv_score)
    
    optimal_smooth = working_interval[which.min(gcv_scores)]
  }
  return (optimal_smooth)
}

#My attempt at Binary Search:

#Below code for binary search based on implementation:
#From Stackoverflow answer by LC94
#https://stackoverflow.com/questions/20133344/find-closest-value-in-a-vector-with-binary-search
k_split_closest_search <- function(x,w,left_end=1,right_end = length(w)){
  midpoint <- floor((left_end + right_end) / 2)
  #print(midpoint)
  if (x < w[midpoint]){
    right_end = midpoint}
  else{
    left_end = midpoint
  }
  #if (length(w) == 2){
  if (right_end - left_end == 1){
    if (abs(w[left_end]-x) < abs(w[right_end]-x)){
      return (w[left_end])
    }
    else{
      return (w[right_end])
    }}
  else{
    return (k_split_closest_search(x,w,left_end,right_end))
  }
  
}

#Sorting Permutation hash generator:
hash_permutation_generator <- function(pre_sorted_item,post_sorted_item){
  pre_sorted_item <- map(pre_sorted_item,toString)
  #print(pre_sorted_item)
  post_sorted_item <- map(post_sorted_item,toString)
  #print(post_sorted_item)
  pre_sorted_hash <- c()
  for (val in pre_sorted_item) {
    pre_sorted_hash[val] <- match(val,pre_sorted_item)
  }
  post_sorted_hash <- c()
  for (val in post_sorted_item) {
    post_sorted_hash[val] <- match(val,post_sorted_item)
  }
  sorting_permutation_hash <- c()
  for (val in post_sorted_item){
    sorting_permutation_hash[post_sorted_hash[val]] <- pre_sorted_hash[val]
  }
  return (sorting_permutation_hash)
}

#Reading in Data (Replace appropriate data path)
K_trial <- as.list(read.delim("C:/Users/HP/Downloads/Siddharth/Jupyter/Python/temp_access_k", header = FALSE, sep = "\t",dec = "."))
R_trial <- as.list(read.delim("C:/Users/HP/Downloads/Siddharth/Jupyter/Python/temp_access_r", header = FALSE, sep = "\t",dec = "."))
K_trial[[length(K_trial)]] <- NULL
R_trial[[length(R_trial)]] <- NULL

#Performing Generalized Cross Validation:
smoothening_parameter <-GCV(lowess,"$y","(input_values,observed_values,f = smooth)",list(K_trial,R_trial),K_trial,R_trial,1/500,2/3,10000)
plot(K_trial,R_trial,pch=19, xlab='k', ylab='r')
lines(lowess(K_trial, R_trial),lwd = 3, col='red') 
#lines(lowess(K_trial, R_trial,smoothening_parameter),lwd = 3, col='red') #f = 1/1000 after that 1/50 used and then #1/10
#lowess_out <- lowess(K_trial, R_trial,f = smoothening_parameter)
lowess_out <- lowess(K_trial, R_trial,f = 0.2)

#Finding K_50
yan_attempt <- lowess_out$y
sorted_yan_attempt <- unlist(sort(yan_attempt))
sorted_actual_vals <- unlist(sort(unlist(R_trial)))

#Using a O(n) search comparing every element and finding minimum (without sorting)

#Actual indicates values in R_trial , lack of it indicates values in the lowess
Closest_vals <-Closest(yan_attempt,0.5)
Closest_vals_actual <- Closest(unlist(R_trial),0.5)

#Same things with which min explicitly shown:
index<- which.min(abs(yan_attempt - 0.5)) 
K50 <- K_trial[index]


Closest_indices <-match(Closest_vals,yan_attempt)
Closest_indices_actual <- match(Closest_vals_actual,R_trial)

K_50 <- lowess_out$x[Closest_indices]
K_50_actual <- K_trial[[Closest_indices_actual]]

#Binary search implementation
Bin_vals <- k_split_closest_search(0.5, sorted_yan_attempt)
Bin_Closest_vals_actual <-k_split_closest_search(0.5, sorted_actual_vals)

Bin_Closest_indices <- match(Bin_vals,sorted_yan_attempt)
Bin_Closest_indices_actual <- match(Bin_Closest_vals_actual,sorted_actual_vals)

#Mapping function to keep track of sorting, so that we can find K_50
loess_corress_hash <-hash_permutation_generator(yan_attempt,sorted_yan_attempt)
Actual_corress_hash <- hash_permutation_generator(R_trial,sorted_actual_vals)

K_50_Bin <- lowess_out$x[loess_corress_hash[Bin_Closest_indices]]
K_50_Bin_actual <- K_trial[[Actual_corress_hash[Bin_Closest_indices_actual]]]

output_list <- c(K_50_Bin,Bin_vals)
mapped <- map(output_list,print)
#print(output_list)

#Time comparisons
Fair_start_point <- unlist(R_trial)
time_taken_1_sum = 0
for (i in seq(from = 1, to = 100, length.out = 100)) {
  start.time <- Sys.time()
  sorted_actual_vals = sort(Fair_start_point)
  k_split_closest_search(0.5, sorted_actual_vals)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_taken_1 <- time.taken
  time_taken_1_sum <- time_taken_1_sum + time_taken_1
}
time_taken_1_avg <- time_taken_1_sum/100

#Closest is same as which abs(min)
time_taken_2_sum = 0
for (i in seq(from = 1, to = 100, length.out = 100)) {
  start.time <- Sys.time()
  Closest(0.5, Fair_start_point)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_taken_2 <- time.taken 
  time_taken_2_sum <- time_taken_2_sum + time_taken_2}
time_taken_2_avg <- time_taken_2_sum/100

#Absolute min suggested:
time_taken_4_sum <- 0
for (i in seq(from = 1, to = 100, length.out = 100)) {
  start.time <- Sys.time()
  desired_val <- Fair_start_point[which.min(abs(Fair_start_point-0.5))]
  #K_50_linear <- spl$x[desired_val]
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_taken_4 <- time.taken
  time_taken_4_sum <- time_taken_4_sum + time_taken_4
}
time_taken_4_avg <- time_taken_4_sum/100