

#https://scales.arabpsychology.com/stats/how-to-do-curve-fitting-in-r/#:~:text=The%20most%20common%20functions%20used%20for%20curve%20fitting,and%20the%20best%20fitting%20model%20can%20be%20chosen.
library(DescTools)
library(data.table)
library(MALDIquant)
library(origami)
library(psych) #Used for tr function for matrix trace
library(mgcv)
library(ggplot2)
library(origami)
library(assist)
#library(matlib)
#library(stringr)
#create data frame
df <- data.frame(x=1:15,
                 y=c(3, 14, 23, 25, 23, 15, 9, 5, 9, 13, 17, 24, 32, 36, 46))

#create a scatterplot of x vs. y
plot(df$x, df$y, pch=19, xlab='x', ylab='y')

#fit polynomial regression models up to degree 5
fit1 <- lm(y~x, data=df)
fit2 <- lm(y~poly(x,2,raw=TRUE), data=df)
fit3 <- lm(y~poly(x,3,raw=TRUE), data=df)
fit4 <- lm(y~poly(x,4,raw=TRUE), data=df)
fit5 <- lm(y~poly(x,5,raw=TRUE), data=df)

#create a scatterplot of x vs. y
plot(df$x, df$y, pch=19, xlab='x', ylab='y')

#define x-axis values
x_axis <- seq(1, 15, length=15)

#add curve of each model to plot
lines(x_axis, predict(fit1, data.frame(x=x_axis)), col='green')
lines(x_axis, predict(fit2, data.frame(x=x_axis)), col='red')
lines(x_axis, predict(fit3, data.frame(x=x_axis)), col='purple')
lines(x_axis, predict(fit4, data.frame(x=x_axis)), col='blue')
lines(x_axis, predict(fit5, data.frame(x=x_axis)), col='orange')

#calculated adjusted R-squared of each model
summary(fit1)$adj.r.squared
summary(fit2)$adj.r.squared
summary(fit3)$adj.r.squared
summary(fit4)$adj.r.squared
summary(fit5)$adj.r.squared

summary(fit4)

#https://www.geeksforgeeks.org/r-language/reading-files-in-r-programming/
  
#Reading in file:
#K_trial <- as.list(read.delim("C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/k_trial", header = FALSE, sep = "\t",dec = "."))
#R_trial <- as.list(read.delim("C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/r_trial", header = FALSE, sep = "\t",dec = "."))
K_trial <- as.list(read.delim("C:/Users/HP/Downloads/Siddharth/Jupyter/Python/temp_access_k", header = FALSE, sep = "\t",dec = "."))
#K_trial <- read.delim("C:/Users/HP/Downloads/Siddharth/Jupyter/Python/temp_access_k_7", header = FALSE, sep = "\t",dec = ".")
R_trial <- as.list(read.delim("C:/Users/HP/Downloads/Siddharth/Jupyter/Python/temp_access_r", header = FALSE, sep = "\t",dec = "."))
K_trial[[length(K_trial)]] <- NULL
R_trial[[length(R_trial)]] <- NULL
#K_trial
coupling_number = 1000
#assist::inc(unlist(R_trial),unlist(K_trial))
spl <- smooth.spline(unlist(K_trial),unlist(R_trial))
match_val <-spl$x[match(Closest(spl$y,0.5),spl$y)]
fit <- loess(unlist(R_trial) ~ unlist(K_trial))
plot(K_trial,R_trial,pch=19,main = "k-r for ICE-Kura order chr1 100kb index 5 of 10x10", xlab='k', ylab='r')
#lines(lowess(K_trial, R_trial,0.0179536),lwd = 3, col='red') #f = 1/1000 after that 1/50 used and then #1/10
lines(unlist(K_trial),fit$fitted,lwd = 3, col ='green')
abline(h = 0.5, col = "blue", lty = 2,lwd = 3)  # Horizontal line
abline(v = match_val, col = "red", lty = 3,lwd = 3)   # Vertical line
text("Predicted",x = match_val,y = max(unlist(R_trial)))
lines(spl,lwd = 3,col = 'orange')

time_taken_3_sum <- 0
for (i in seq(from = 1, to = 100, length.out = 100)) {
  start.time <- Sys.time()
  target_y <- 0.5
  root_fun <- function (x) predict(spl,x)$y - 0.5
  uniroot_val<-uniroot(root_fun,lower = 0.0, upper = 0.025)$root
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_taken_3 <- time.taken
  time_taken_3_sum <- time_taken_1_sum + time_taken_1
}
time_taken_3_avg <- time_taken_3_sum/100


predict(spl,uniroot_val)

#lines(gam_model,shade=TRUE,seWithMean=TRUE,scale=0)
lowess_out <- lowess(K_trial, R_trial,f = 0.0179536)
yan_attempt <- lowess_out$y

df = data.frame(unlist(R_trial),unlist(K_trial))
gam_model <- gam(unlist(R_trial) ~ s(unlist(K_trial)))
summary(gam_model)
x_seq <- seq(min(unlist(K_trial)),max(unlist(K_trial)),length.out = 500)
predictions <- predict(gam_model, data = data.frame(unlist(K_trial),unlist(R_trial)))
plot(gam_model,shade=TRUE,seWithMean=TRUE,scale=0)
lines(K_trial,R_trial,col = 'blue')

folds <- make_folds(data.frame(K_trial,R_trial))
#ggplot() +
  #geom_point(aes(x = K_trial, y = R_trial)) +
  #geom_line(data = data.frame(hp = new_data$hp, mpg = predictions$fit), 
            #aes(x = hp, y = mpg), color = "blue", size = 1) +
  #geom_ribbon(data = data.frame(hp = new_data$hp, fit = predictions$fit, 
                                #se = predictions$se.fit), aes(x = hp, 
                                                             # ymin = fit - 1.96 * se, 
                                                              #ymax = fit + 1.96 * se), alpha = 0.3) +
  #labs(title = "Generalized Additive Model (GAM) Fit for mpg vs. hp", 
       #x = "Horsepower", y = "Miles per Gallon") +
  #theme_minimal()
plot(unlist(K_trial),predictions)
sorted_yan_attempt <- unlist(sort(yan_attempt))
sorted_actual_vals <- unlist(sort(unlist(R_trial)))
sorted_yan_attempt[[545]]
match(0.5000777,yan_attempt)
which(abs(yan_attempt - 0.5) < 0.002)
lowess_out$x[[59]]

Closest_indices <-Closest(yan_attempt,0.5)
Closest_indices_actual <- Closest(unlist(R_trial),0.5)
match(Closest_indices,yan_attempt)
match(Closest_indices_actual,R_trial)
yan_attempt[[56]]
R_trial[[48]]
K_trial[[48]]
#function(){
match.closest(0.5,sorted_actual_vals)
#}
findInterval(0.5,sorted_yan_attempt)
findInterval(0.5,sorted_actual_vals)
sorted_actual_vals[[474]]

time_taken_5_sum = 0
for (i in seq(from = 1, to = 100, length.out = 100)) {
start.time <- Sys.time()
dt = data.table(sorted_actual_vals, val = sorted_actual_vals) # you'll see why val is needed in a sec
setattr(dt, "sorted", "sorted_actual_vals")
dt[J(0.5), roll = "nearest"]
end.time <- Sys.time()
time_taken_5 <- end.time - start.time
time_taken_5_sum <- time_taken_5 + time_taken_5_sum
#dt[J(0.5), .I, roll = "nearest", by = .EACHI]
}
time_taken_5_avg <- time_taken_5_sum/100
 # sorted_yan_attempt[is.nan(sorted_yan_attempt)]

x <- sorted_yan_attempt
left_end = 1
right_end <- length(x)
slice_1 <- x[left_end:(length(x)%/%2)]
slice_2 <- x[(length(x)%/%2):right_end]

#My attempt at Binary Search:
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
GCV(lowess,"$y","(input_values,observed_values,f = smooth)",list(K_trial,R_trial),K_trial,R_trial,1/500,2/3,10000)
start.time <- Sys.time()
k_split_closest_search(0.5, sorted_actual_vals)
end.time <- Sys.time()
time.taken <- end.time - start.time



#Above code for binary search based on below implementation:

#From Stackoverflow answer by LC94
#https://stackoverflow.com/questions/20133344/find-closest-value-in-a-vector-with-binary-search
NearestValueSearch = function(x, w){
  ## A simple binary search algo
  ## Assume the w vector is sorted so we can use binary search
  left = 1
  right = length(w)
  while(right - left > 1){
    middle = floor((left + right) / 2)
    if(x < w[middle]){
      right = middle
    }
    else{
      left = middle
    }
  }
  if(abs(x - w[right]) < abs(x - w[left])){
    return(right)
  }
  else{
    return(left)
  }
}


x = 4.5
w = c(1,2,4,6,7)
NearestValueSearch(x, w) # return 3
NearestValueSearch(0.5, sorted_actual_vals)
sorted_actual_vals[56]
unlist(R_trial)[254]

start.time <- Sys.time()
sorted_actual_vals = sorted_actual_vals <- unlist(sort(unlist(R_trial)))
k_split_closest_search(0.5, sorted_actual_vals)
end.time <- Sys.time()
time.taken <- end.time - start.time
time_taken_1 <- time.taken

start.time <- Sys.time()
Closest(0.5, sorted_actual_vals)
end.time <- Sys.time()
time.taken <- end.time - start.time
time_taken_2 <- time.taken

print(c(time_taken_1,time_taken_2))