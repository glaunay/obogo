import networkx as nx
from typing import Iterator, Union, NewType, Literal, get_args
from networkx.classes.digraph import DiGraph
from .io_obo import obo_node_buffer_iter
from uniprot_redis.store.schemas import UniprotDatum, UniprotAC
from pyproteinsext.uniprot import Entry as Uniprot

from .type_checkers import literal_arg_checker
import sys

NodeID   = NewType("NodeID", str)
NodeName = NewType("NodeName", str)
ProteinsType  = Literal["background", "measured"]
PercolateType = Literal["background", "measured", "both"]

import re


class GO_tree(nx.DiGraph):
    def __init__(self):
        super().__init__()
        self.protein_load_status = (False, False)
        self.percolated_status   = (False, False)
        self.uniprot_omega       = (set(), set()) # background uniprot data, measured uniprot data
        self.names_index = {}
    def add_node(self, id, **params):
        """
        nx add_node wrapper to account for node aliases
        """

        super(GO_tree, self).add_node(id, **params)
        if 'name' in params:
            self.names_index[ params['name'] ] = id
        
        if 'alt_id' in params:
            alt_ids = params['alt_id'] if type(params['alt_id']) is list else [ params['alt_id'] ]        
            real_n = self.nodes[id]
            for alt_id in alt_ids:
                super(GO_tree, self).add_node(alt_id, alias_to = real_n, _id = alt_id)
    
    def get_go_node(self, node_id:Union[NodeID, NodeName], many_to_consider=False)->Union[dict, list[dict]]:
        """
        nx G.node wrapper for the automatic forward of alias go_id to concrete node 
        
        [GL comments 01/11/2023]
        Invalid entry policy and the get_go_node method
        2 types of invalid GO identifiers::
        1) The obsolete ones, eg:
        [Term]
        id: GO:0000005
        name: obsolete ribosomal chaperone activity
        ...
        is_obsolete: true
        consider: GO:0042254
        consider: GO:0044183
        consider: GO:0051082

        Bascially means that GO:0000005 was removed and splited/mapped onto SEVERAL valid GO terms
        It makes sense to silently forward this one-to-many application solely in the
        load_protein methods. This corresponds to protein anntoted with an outdated GO versions, for
        which we dont wanna loose this oboslete piece of information. In this situation 
        get_go_node('GO:0000005') will return the list of 'consider' nodes [get_go_node('GO:0042254'), ...]
        In all other cases access to an obsolete entry will raise a warning along w/ a None return value providing the possible repalcement
        terms in the warning message. This behaviour is controled by the "many_to_consider" boolean parameter.
        Hence, Obsolete not should not carry any protein list.

        2) the removed ones, eg:
        [Term]
        id: GO:0000003
        name: reproduction
        ...
        alt_id: GO:0019952
        alt_id: GO:0050876
        ...

        Here id: GO:0000003 is valid but GO:0019952 and GO:0050876 were merged into GO:0000003
        No problem anticipated here, ee only create ghost node with  GO:0019952 or GO:0050876 as id 
        and reference to the actual node GO:0000003.

        2a) 
        [Term]
        id: GO:0003783
        name: obsolete barbed-end actin capping activity
        ...      
        is_obsolete: true
        replaced_by: GO:0051016

        Seems like the above reciprocal, the term itself no longer used and replaced by a single 
        term. We silently forward from the deprecated to the replacement one. 
        """

        isGOterm = lambda value: True if re.match(r'^GO:[0-9]+', value) else False

        if not isGOterm(node_id):
            if not node_id in self.names_index:
                raise KeyError(f"This value \"{node_id}\" is not a GO identifier or a GO name")
            node_id = self.names_index[node_id]

        _ = self.nodes[node_id]
        
        if "alias_to" in _:
            return _["alias_to"]
        
        if not 'is_obsolete' in _: 
            return _
        
        if 'replaced_by' in _:
            return self.nodes[ _['replaced_by'] ]
        
        if not 'consider' in _:
            sys.stderr.write(f"{node_id} is an obsolete GO term with no valid term to replace it\n")
            return _
        if not many_to_consider:
            sys.stderr.write(f"{node_id} is an obsolete GO term, please consider \"{ ','.join(_['consider']) }\"\n")
            return None
        return [ self.get_go_node(valid_go_id) for valid_go_id in _['consider'] ]
    
 
    def view_go_node(self, go_id:str)->dict:
        """ returns string representation of concrete node where proteins values are simply represented by their uniprotID """
        n = self.get_go_node(go_id)
        return { k : [ u.id for u in v ] if k in ['background', 'measured', 'perc_measured', 'perc_background'] else v for k,v in n.items() }
        

    def concrete_nodes(self)->Iterator[dict]:
        """ node iterator over guaranteed concrete nodes only """
        for k in self.nodes:
            d = self.maybe_concrete(k)
            if d is None:
                continue
            yield d

    def maybe_concrete(self, node_id:NodeID)->Union[dict, None]:
        """ Boolean evaluation of queried node as a node that is not obsolete nor and alias
        """
        _ = self.nodes[node_id]
        if not 'is_obsolete' in _ and not 'alias_to' in _:
            return _
        return None
        
    def successors(self, node_id:NodeID,  many_to_consider = False):
        """ Wrapping for alias forwarding 
            Do we make it possible to go down from a deprecated GO term ?
            It is implemented but should not be used. 
            As, in pratice successors are meant to be accessed from a process that starts from valid go term leaves.
        """
        go_node = self.get_go_node(node_id, many_to_consider)
        if not type(go_node) is list :
            for s in super(GO_tree, self).successors(go_node['_id']):
                yield s
        # This will be slow ^^
        else:
            for g in go_node:
                s = super(GO_tree, self).successors(node_id)
                for _ in s:
                    yield _

    @literal_arg_checker
    def load_proteins(self, k:ProteinsType, protein_coll:Iterator[ Union[UniprotDatum, Uniprot] ]) -> None:
        """
        Iterate through a uniprot datum collection and attach UniprotDatum to 
        corresponding concerte node to background or measured attribute
        we eventually forward obsolete entry to their many valid nodes
        """
        self.clear_proteins(k)
        for unip_datum in protein_coll:
            for go_datum in unip_datum.go:
                try:
                    go_node = self.get_go_node(go_datum.id, many_to_consider=True)
                    go_node = [ go_node ] if not type(go_node) is list else go_node
                    for g in go_node:
                        if not k in g:
                            g[k] = set()
                        g[k].add(unip_datum)
                except KeyError as e:
                    print(f"{go_datum.id} not found")
                    print(e)

        self.protein_load_status = (True, self.protein_load_status[1]) if k == "background" \
        else (self.protein_load_status[0], True)

    @literal_arg_checker
    def clear_proteins(self, k: ProteinsType):
        """ clear all the protein lists: "perc_" and "classic" of the specified type: "background" or "measured" """
        for node_dic in self.nodes.values():
            _ = node_dic.pop(k, None)
            if f"perc_{k}" in node_dic:
                _ = node_dic.pop(f"perc_{k}", None)
        self.protein_load_status = (False, self.protein_load_status[1]) if k == "background" \
        else (self.protein_load_status[0], False)

    @property
    def root_ids(self)->set[NodeID]:
        r = set()
        for n in self.concrete_nodes():
            if not '_is_a' in n:
                r.add(n['_id'])
        return r

    @literal_arg_checker
    def get_proteins(self, node_id:NodeID, k:ProteinsType="background", deep=True)->set[UniprotDatum]:
        """ Get all UniprotDatum attached to subtree rooted at provided node id """

        def _get_proteins(node_id, deep)->set[UniprotDatum]:
            curr_node = self.get_go_node(node_id)
            results   = curr_node[k] if k in curr_node else set()
            if not deep:
                return results
            
            for next_node_id in self.successors(node_id):
                results |= _get_proteins(next_node_id, deep)
            return results
            
        return _get_proteins(node_id, deep)
    
    @property
    def leave_ids(self)->Iterator[NodeID]:
        for d in self.concrete_nodes():
            bot = True
            for _ in self.successors(d['_id']):
                bot = False
                break
            if bot:
                yield d['_id']
    
    @property
    def ora_rdy(self):
        return self.protein_load_status[0] and self.protein_load_status[1] and \
               self.percolated_status[0]   and self.percolated_status[1]

    @literal_arg_checker
    def percolate(self, percol_type:PercolateType="both"):
        """ makes the perc_background/measure attribute of one go_node the union of its descendants
            A parent node will pbbly be visited more than as all its descendant must bubble up their proteins
        """
        #assert(self.precolate_rdy)
        def uniprot_upercolator(go_id:str, buffer_bkg:set[UniprotDatum],  buffer_mea:set[UniprotDatum]):
            """ makes the proteins attribute of one go_node the union of its descendants
                A parent node will pbbly be visited more than as all its descendant must bubble up their proteins
            """
            #print(f"Visiting {go_id}")
            #print(f"Target tree curr members:{tree_tgt.nodes}")
            # Check if current node has already been visited
            n_dict = self.nodes[go_id]
            if percol_type in ["both", "background"]:
                if 'perc_background' in n_dict:
                    n_dict['perc_background'] = n_dict['perc_background'] | buffer_bkg
                else :
                    n_dict['perc_background'] = n_dict['background'] | buffer_bkg if 'background' in n_dict \
                                                else buffer_bkg
        
            if percol_type in ["both", "measured"]:
                if 'perc_measured' in n_dict:
                    n_dict['perc_measured'] = n_dict['perc_measured'] | buffer_mea
                else:
                    n_dict['perc_measured'] = n_dict['measured']      | buffer_mea if 'measured' in n_dict \
                                                else buffer_mea
             
            for go_id_pred in self.predecessors(go_id):                        
                uniprot_upercolator( go_id_pred, 
                    n_dict['perc_background'].copy() if percol_type in ["both", "background"] else set(),
                    n_dict['perc_measured'].copy()   if percol_type in ["both", "measured"]  else set()
                )
               
        for go_id in self.leave_ids:
            uniprot_upercolator(go_id, set(), set())
        
        for root_id in self.root_ids:
            if percol_type in ["both", "background"]:
                root_bkg = self.get_go_node(root_id)['perc_background']
                self.uniprot_omega = ( self.uniprot_omega[0] | root_bkg, self.uniprot_omega[1])  
            if percol_type in ["both", "measured"]:
                root_mea = self.get_go_node(root_id)['perc_measured']
                self.uniprot_omega = ( self.uniprot_omega[0],  self.uniprot_omega[1] | root_mea)  
        
        self.percolated_status= ( True if percol_type in ["both", "background"] else self.percolated_status[0], 
                                  True if percol_type in ["both", "measured"]   else self.percolated_status[1] )

# MAybe move to io
def reader(obo_file_path:str, keep_obsolete=True, ns=None)-> DiGraph :
    G = GO_tree() #nx.DiGraph()

    for node_buffer in obo_node_buffer_iter(open(obo_file_path, 'r')):
        if not keep_obsolete and node_buffer.is_obsolete:
            continue
        G.add_node(node_buffer['id'], **node_buffer.nx_node_param)
        for node_parent_id in node_buffer.is_a_iter():
            G.add_edge(node_parent_id, node_buffer['id'], type='is_a')

    return G

