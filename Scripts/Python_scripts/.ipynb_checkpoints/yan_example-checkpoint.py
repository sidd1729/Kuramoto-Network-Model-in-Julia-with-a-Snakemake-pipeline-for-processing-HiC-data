from multiprocessing import Pool, TimeoutError
import time
import os
from tqdm import tqdm

#def f(x):
    #return [y*y for y in x]
def f(x):
    return x*x
if __name__ == '__main__':
    start_time = time.perf_counter()   
    max_ = 2
    with Pool(processes=20) as p, tqdm(total=max_) as pbar:
            for result in p.imap(f, [1,2]):
                pbar.update()
                pbar.refresh()
    end_time = time.perf_counter()  
    print(start_time-end_time)