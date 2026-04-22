
#Great Works #Stackoverflow answer of ah bon : https://stackoverflow.com/questions/19894365/running-r-script-from-python
#import rpy2.robjects as robjects
#robjects.r.source("C:/Users/HP/Downloads/Siddharth/Scripts/R_scripts/Curve_fitting_example.R", encoding="utf-8")



import subprocess
import sys
#subprocess.call(['RScript',"C:/Users/HP/Downloads/Siddharth/Scripts/R_scripts/Curve_fitting_example.R"])

#results = (subprocess.check_output(['RScript',"C:/Users/HP/Downloads/Siddharth/Scripts/R_scripts/Curve_fitting_GCV_trial.R"],text=True))
#results = (subprocess.check_output(['RScript',"C:/Users/HP/Downloads/Siddharth/Scripts/R_scripts/Curve_pipeline.R","C:/Users/HP/AppData/Local/Programs/Python/Python313/","C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/10k_coupling_k_vals_2_1000.pickle","C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/10k_coupling_trial_2_1000.pickle","C:/Users/HP/Downloads/plot.png","k-r for ICE-Kura order chr1 100kb index 5 of 10x10"],text=True))
results = (subprocess.check_output(['/usr/bin/Rscript',"/home/ksslab/Siddharth/Program_Directory/Scripts/R_scripts/Curve_pipeline.R","/home/ksslab/Programs/anaconda3/envs/Kuramoto_scripts_env/bin/python",sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]],text=True))
#print(type(results),type(str(results)),'\n')
results_list = results.split("\n")
del(results_list[len(results_list)-1])
#print(results,'\n',results[4:len(results)-2],'\n',results_list)
#print(list(map(lambda x : float(x[4:]), results_list)))
final_output = tuple(map(lambda x : float(x[4:]), results_list)) #Output in the form of K_50 value,R_value
print(str(final_output[0]) + ','+ str(final_output[1]))

#x = ( 1,2,3,4,5 )
#y =   ( 1,2,3,4,5 )
#result_test = subprocess.check_output(['RScript',"C:/Users/HP/Downloads/Siddharth/Scripts/R_scripts/Curve_pipeline.R","{}".format(x),"2"],text=True)
#print(result_test)
#subprocess.call("/usr/bin/Rscript --vanilla C:/Users/HP/Downloads/Siddharth/Scripts/R_scripts/Curve_fitting_example.R", shell=True)
'''
import pickle
def tuple_list_writer(filename,variable_name):
    file = open('{}'.format(filename), 'w')
    for t in variable_name:
      file.write(' '.join(str(s) for s in t) + '\n')
    file.close()

def tuple_list_reader(filename):
    file = open('{}'.format(filename), 'r')
    content = file.read().splitlines()
    tuple_list =[]
    file.close()
    for row in content:
        tuple_list.append(tuple(float(i) for i in row.split(' ')))
    return tuple_list
filename = "C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/10k_coupling_k_vals_2_1000.pickle"
tuple_list = tuple_list_reader("C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/10k_coupling_k_vals_2_1000")

tuple_list_obj = open(filename,'wb')
pickle.dump(tuple_list,tuple_list_obj)
tuple_list_obj.close()

tuple_list_obj = open(filename,'rb')
print(pickle.load(tuple_list_obj))
tuple_list_obj.close()
'''
'''
from rpy2.robjects import r

# Define R code as a string
r_code = """
K_trial <- as.list(read.delim("C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/k_trial", header = FALSE, sep = "\t",dec = "."))
R_trial <- as.list(read.delim("C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/r_trial", header = FALSE, sep = "\t",dec = "."))
K_trial[[length(K_trial)]] <- NULL
R_trial[[length(R_trial)]] <- NULL
#K_trial
plot(K_trial,R_trial,pch=19, xlab='x', ylab='y')
lines(lowess(K_trial, R_trial,f = 1/50), col='red')
"""

# Execute the R code and get the result
result = r(r_code)
print("R Output:", result)
'''
'''
read_delim = robjects.r['read.delim']
as_list = robjects.r['as.list']
length = robjects.r['length']
NULL = robjects.r['NULL']
FALSE = robjects.r['FALSE']
plot = robjects.r['plot']
lowess = robjects.r['lowes']
lines = robjects.r['lines']
K_trial = as_list(read_delim("C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/k_trial", header = FALSE, sep = "\t",dec = "."))
R_trial = as_list(read_delim("C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/r_trial", header = FALSE, sep = "\t",dec = "."))
K_trial[length(K_trial)] = NULL
R_trial[length(R_trial)] = NULL
#K_trial
plot(K_trial,R_trial,pch=19, xlab='x', ylab='y')
lines(lowess(K_trial, R_trial,f = 1/50), col='red')
'''