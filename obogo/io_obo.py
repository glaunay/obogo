from io import TextIOWrapper
import re
"""
[Term]
id: GO:0000001
name: mitochondrion inheritance
namespace: biological_process
def: "The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton." [GOC:mcc, PMID:10873824, PMID:11389764]
synonym: "mitochondrial inheritance" EXACT []
is_a: GO:0048308 ! organelle inheritance
is_a: GO:0048311 ! mitochondrion distribution

"""

class Buffer:
    """ Convenience Dict wrapper to transform obo text input into nx.add_nodes ***kwargs
        Basically, all dict values are lists
        Beware, the __getitem__ method will coherce any list dict value of len 1 into a scalar
    """
    def __init__(self):
        self.data = {}

    def __iter__(self):
        for k in self.data:
            yield k
    
    def __setitem__(self, k, v):
        if not k in self.data:
            self.data[k] = []
        self.data[k].append(v)
    
    def __getitem__(self, k):
        if not k in self.data:
            raise KeyError(f"key \"{k}\" not found in {self.data}")
        return self.data[k][0] if len(self.data[k]) == 1 else self.data[k]

    def __contains__(self, k):
        return k in self.data
    
    def clear(self):
        self.data = {}
    
    def __bool__(self):
        return bool(self.data)
    def __repr__(self):
        _ = { k : self.__getitem__(k) for k in self.data }
        return str(_)
    
    @property    
    def is_obsolete(self):
        return 'is_obsolete' in self.data
            
    def is_a_iter(self):
        if 'is_a' in self.data:
            return self.data['is_a']
        return []

    def consider_iter(self):
        if 'consider' in self.data:
            return self.data['consider']
        return []

    @property
    def nx_node_param(self):
        """ return dict view of the node with only k,val suitable for nx.node declaration
            effectively shadowing 'is_a' (used for edge) namespace fields
        """

        _ = { k if not k in {'id', 'is_a'} else '_id' if k == 'id' else '_is_a' : self.__getitem__(k)  for k in self if not k in { 'consider' } } 
        if 'consider' in self.data:
            _['consider'] = self.data['consider']
        
        return _
def obo_node_buffer_iter(fp:TextIOWrapper):
    buffer = Buffer()
    read_switch = False
    for line in fp:
        # update load bool
        line=line.rstrip()
        m = re.search(r'^\[(.+)\]$', line)
        if m:
            read_switch = m[1] == 'Term'
            continue
        # empty line is yield condition
        if re.match(r'^\s*$', line) and buffer:
          #  print(f"Poping {buffer}")
            yield buffer
            buffer.clear()
            continue
        if read_switch: # putting stuff into buffer         
            m = re.search(r'^(.+): "(.+)"[^"]+$', line)
            if not m:
                m = re.search(r'^(name): (.+)', line)
            if not m:
                m = re.search(r'^(.+): ([\S]+)', line)
                if not m:
                    raise ValueError(line)
            buffer[ m[1] ] = m[2]

    # pop trailer
    if buffer:
        yield buffer
    fp.close()

