from networkx.classes.digraph import DiGraph
import sys, time
from tree import reader
if __name__ == "__main__":
    
    start = time.time()
    kpo = False
    G:DiGraph  = reader(sys.argv[1], keep_obsolete=kpo)
    end = time.time()
    
    print(f"DiGraph of {len(G)} nodes and {len(G.edges())} edges loaded in {end - start} [obsol_keep = {kpo}]")
    #for u, v in G.edges():
    #    print (f"{u} => {v}")
