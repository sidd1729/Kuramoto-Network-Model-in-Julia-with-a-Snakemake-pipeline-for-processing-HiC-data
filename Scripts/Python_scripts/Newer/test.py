from cooler.create import ArrayLoader
import h5py
import cooler
import pandas as pd
h = h5py.File("cworld-test_hg19_C-40000-raw.hdf5", 'r')
heatmap = h['interactions']

# create some bins , using cooler-binnify or some other way
binsize = 40000
chromsizes = pd.read_csv(                                                            
    'hg19.reduced.chromsizes', 
    sep='\t', 
    names=['name', 'length']).set_index('name')['length']
bins = cooler.binnify(chromsizes, binsize)

# turn h5oy dataset (2D matrix) into a stream of sparse matrix chunks :
iterator = ArrayLoader(bins, heatmap, chunksize=int(1e6))

# load that into cooler:
cooler.create_cooler('output.40kb.cool', bins, iterator, dtypes={"count":"int"}, assembly="hg19")