import re
from tqdm import tqdm
import pickle

def fasta_reader(filename):
    file = open("{}".format(filename),'r')
    content = file.read().splitlines()
    file.close()
    del content[0]
    total_len = 0
    for i in range(0,len(content)):
        content[i] = content[i].upper()
        total_len += len(content[i])
    return content,total_len

def file_char_counter(filename):
    file = open(filename,'r')
    content = file.read()
    file.close()
    len_N = len(re.findall(r'[N]',content))
    len_char = len(re.findall(r'[a-zA-Z]',content))
    return (len_N,len_char)

def general_pickler(filename,tbp_object=None,read = False,reader_fn = None):
    if read == True:
        tbp_object = reader_fn(filename)      
    with open(filename + ".pickle",'wb') as f:
        pickle.dump(tbp_object,f)

def bin_dependent_GC_content(content_list,total_len,bin_size,region_chosen = False,region_limits = (0,15000000),overlapping = True):
    bin_GC_storer = []
    total_undef = 0
    start_bin_no = 0 
    if overlapping == True:
        bin_step = 1
    else:
        bin_step = bin_size
    if region_chosen == False:
        bin_no = total_len-bin_size+1
        #print("{} is number of bins".format(bin_no))
    else:
        bin_no = (region_limits[1]-region_limits[0])-bin_size+1
        start_bin_no = region_limits[0]//bin_size
        #print("{} is number of bins".format(bin_no))
    stop_bin_no = start_bin_no + bin_no
    stop_region = region_limits[0] + bin_no
    #print("{},{} are start and end bin numbers".format(start_bin_no,stop_bin_no))
    fasta_line_len = len(content_list[0])
    line_no_per_bin = bin_size//fasta_line_len
    #print("{} is number of lines per bin".format(line_no_per_bin))
    if line_no_per_bin < 1:
        modulus_cycle = 0
        beginning_line_no = region_limits[0]//fasta_line_len 
        if start_bin_no %fasta_line_len == 0:
            modulus_cycle = -1
        current_line_no = beginning_line_no
        for i in tqdm(range(region_limits[0],stop_region,bin_step)):
            if i % fasta_line_len == 0:
                modulus_cycle +=1
                current_line_no = beginning_line_no + modulus_cycle
            bin_str=''
            #if i%fasta_line_len
            for j in (current_line_no,current_line_no+2):
            #for j in range((i*bin_size)//fasta_line_len,((i*bin_size)//fasta_line_len)+2):
            #for j in range(bin_size * i //fasta_line_len,(bin_size*i//fasta_line_len)+2):
            #for j in range(bin_size * (i - i%fasta_line_len),(bin_size * (i - i%fasta_line_len))+2):
                #print("{} is index being accessed at content_list".format(j))
                bin_str += content_list[j] #+'\t'
            #print("{} is string extracted to use in bin".format(bin_str))
            #print("{},{},{} are modulus of i with line length, quotient and i ".format(i%fasta_line_len,i//fasta_line_len,i))
            GC_content_str = bin_str[i%fasta_line_len:(i%fasta_line_len)+bin_size]
            #if i > 1203 and i < 1206:
            #print("{} is string for which we are calculating GC content".format(GC_content_str))
            len_GC = len(re.findall(r'[GC]',GC_content_str))
            len_AT = len(re.findall(r'[AT]',GC_content_str))
            len_others = bin_size - len_GC - len_AT
            total_undef += len_others
            if len_others == bin_size:
             GC_content = 0
            else:
             GC_content = len_GC/(len_GC+len_AT)*100
            bin_GC_storer.append(GC_content)
    else:
        beginning_line_no = region_limits[0]//fasta_line_len
        #print("{} is beginning line no".format(beginning_line_no))
        modulus_cycle = 0
        #print("{} is modulus cycle".format(modulus_cycle))
        if region_limits[0]%fasta_line_len == 0:
            modulus_cycle = -1
        current_line_no = beginning_line_no
        #print("{} is modulus cycle".format(modulus_cycle))
        #print("{} is current line no".format(current_line_no))
        slice_incrementer = 0
        for i in tqdm(range(region_limits[0],stop_region,bin_step)):
            if i % fasta_line_len == 0:
                modulus_cycle +=1
                slice_incrementer = 0
                current_line_no = beginning_line_no + modulus_cycle
            else:
                slice_incrementer+=1
            bin_str=''
            #print("{} is modulus cycle".format(modulus_cycle))
            #print("{} is current line no".format(current_line_no))
            for j in range(current_line_no, current_line_no + line_no_per_bin +1):
                bin_str += content_list[j]
                #print("{} is index being accessed at content_list".format(j))
            #print("{} is string extracted to use in bin".format(bin_str))
            #print("{},{},{} are modulus of i with line length, quotient and i ".format(i%fasta_line_len,i//fasta_line_len,i))
            #GC_content_str = bin_str[(i*bin_size)%fasta_line_len:((i*bin_size)%fasta_line_len) + bin_size]
            GC_content_str = bin_str[i % fasta_line_len: (i % fasta_line_len) + bin_size]
            #print("{} is string for which we are calculating GC content".format(GC_content_str))
            #if i > 1203 and i < 1206:
            #print(GC_content_str)
            len_GC = len(re.findall(r'[GC]',GC_content_str))
            len_AT = len(re.findall(r'[AT]',GC_content_str))
            len_others = bin_size - len_GC - len_AT
            total_undef += len_others
            if len_others == bin_size:
             GC_content = 0
            else:
             GC_content = len_GC/(len_GC+len_AT)*100
            bin_GC_storer.append(GC_content)
            #print(GC_content_str)
    return bin_GC_storer,total_undef
        

len_N,len_char = file_char_counter("/home/ksslab/Siddharth/Program_Directory/Data/hg38/chr1.fa")
content,total_len = fasta_reader("/home/ksslab/Siddharth/Program_Directory/Data/hg38/chr1.fa")
bin_GC_storer,total_undef = bin_dependent_GC_content(content,len_char,100000,region_chosen = True,region_limits = (0,248956425),overlapping = True)
general_pickler("/home/ksslab/Siddharth/Program_Directory/Data/hg38/100kbGC",bin_GC_storer)