import gzip
import bz2
import lzma

class RawFile(object):
    def __init__(self,filename):
        self.filename = filename
        if filename.endswith('.gz'):
            self.handle = gzip.open(filename,'rt')
        elif filename.endswith('bz2'):
            self.handle = bz2.open(filename,'rt')
        elif filename.endswith('xz'):
            self.handle = lzma.open(filenaem,'rt')
        else:
            self.handle = open(filename,'r')
    def __enter__(self):
        return self.handle
    def __exit__(self,dtype,value,traceback):
        self.handle.close()

def gen_ensem_id_map():
    '''
    Generates and Ensemble chromosome ID map. The chromosome
    names in ensemble are just numbers, this function prepends
    a 'chr' to the name and fixes the names for chrMt and chrX.
    '''
    ensem_map = {}
    def inc_range(start, end):
        '''
        inclusive range
        '''
        return range(start, end+1)
                            
    for i in inc_range(1,31):
        ensem_map[str(i)] = 'chr' + str(i)
        if 'MT' or 'X' not in ensem_map:
            ensem_map['MT'] = 'chrMt'
            ensem_map['X'] = 'chrX'
    return ensem_map   


