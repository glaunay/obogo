from .tree import GO_tree
from uniprot_redis.store.schemas import UniprotAC, UniprotDatum
import scipy.stats as stats 
from typing import Literal, Union, Iterator

SortCrit = Literal['pvalue', 'count', 'bkfq']
OraNormalizer = Literal["background", "measured"]

class ORA_error(Exception):
    pass
    
def score_ora_tree(tree:GO_tree, delta_prot:Union[ Iterator[UniprotDatum], Iterator[UniprotAC]],
                         norm:OraNormalizer="background"):
    """
    Compute pvalue of Fisher Test on every Go term
    Parameters : tree, a GO_tree object where both measured and background protein sets have been percolated
                 delta_prot, an iterator over the "abundant" elements
    The Backgound population is defined by "norm" paramter: 
                            - 'background' will use the whole proteome as background population
                            - 'measured' will use the experimentally measured proteome as background population
    Implicitly background protein sets is a superset of the measured which is a superset of the delta_prot
    Returns:

    """
    delta, N, pop_key = ora_validator(tree, delta_prot, norm)

    for n in tree.concrete_nodes():
        _ = _node_ora(n, delta, N, pop_key)
        if not _ is None:
            yield _

def ora_validator(tree, delta_prot, norm):
    """ coherce delta into strings, set the backgroud pop and convert it into strings
    """
    if not tree.ora_rdy:
        raise ORA_error("GO term tree not ready")
    inputs_as_string = False
    omega_bkg, omega_mea = tree.uniprot_omega
    N = omega_bkg if norm == "background" else omega_mea
    # consume iterator into persitent list of uniprotDatum and convert to uniprotACs
    delta = set([ _ if type(_) is str else _.id for _ in delta_prot ])
    N     = set([ _.id for _ in N ])
    
    return delta, N, f"perc_{norm}"

def compute_node_ora(tree:GO_tree, delta_prot:Union[ Iterator[UniprotDatum], Iterator[UniprotAC]], go_id:str,
                         norm:OraNormalizer="background"):
    delta, N , pop_key = ora_validator(tree, delta_prot, norm)
    n = tree.get_go_node(go_id)
    return _node_ora(n, delta, N, pop_key)

def _node_ora(n:dict, delta:set[UniprotAC], N:set[UniprotAC], pop_key):
    
    # proteins abundant and path member
    n_path = set( [ _.id for _ in n[pop_key] ] )
    n_not_path = N - n_path
    s11 = delta & n_path

    if not s11:
        return None
    
    delta_0  = N - delta
    # proteins abundant and  not path member
    s12 = delta   & n_not_path
    s21 = delta_0 & n_path
    s22 = delta_0 & n_not_path
    
    """     | Path   | not(Path)
    --------------------------s
    Delta_+ |   s11   |   s12
    -------------------------
    Delta_0 |   s21   |   s22
        
    """
    
    odd_ratio, p_value = stats.fisher_exact( [ ( len( s11 ) , len( s12 ) ),
                                                ( len( s21 ) , len( s22 ) )
                                            ])
    return  ( n['_id'], n['name'], len(s11), odd_ratio, p_value,\
                    [ ( len( s11 ) , len( s12 ) ),
                        ( len( s21 ) , len( s22 ) )
                ])