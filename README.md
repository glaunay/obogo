# obogo: yet another GeneOntology reader/crawler !!!

## installation
`pip install obogo`
### prerequistes
You will need to grab the latest Gene Ontology description in `.obo` format. Download it from the [official GO site](https://geneontology.org/docs/download-ontology/).
 This is a reasonable format, where each GO term is encoded as a paragraph of key:value lines:
```txt
[Term]
id: GO:0000001
name: mitochondrion inheritance
namespace: biological_process
def: "The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton." [GOC:mcc, PMID:10873824, PMID:11389764]
synonym: "mitochondrial inheritance" EXACT []
is_a: GO:0048308 ! organelle inheritance
is_a: GO:0048311 ! mitochondrion distribution
```

Obvisously, `is_a` is the main relationship of interest to build the DAG of GO terms. But `consider`, `replaced_by` and `alias_to` properties are also handled by this package.

### Build the GO DAG structure
Simply read from a flat `.obo` file
```python
from obogo import create_tree_from_obo
obogo_tree = create_tree_from_obo('../data/go-basic.obo')
```

You can query a go term by a name or its GO identifier
```python
obogo_tree.view_go_node('GO:1903507')
obogo_tree.view_go_node('biological process')
```
### Build a protein collection
In order to set the background population of proteins for each GO term, you will need to build a collection of uniprot data containers. In obogo, these are called `UniprotDatum` and can be directly created from a [uniprot proteome xml reference file](https://www.uniprot.org/proteomes?query=*) (eg: [E.coli K12](https://www.uniprot.org/proteomes?facets=proteome_type:1&query=(organism_id:83333))).
Supposed we downloaded the above mentioned E.Coli K12 proteome xml file named `uniprotkb_proteome_UP000000625.xml`, the collection of uniprot containers can be buildt this way:

```python
from uniprot_redis.store.mockup import UniprotStoreDummy
from uniprot_redis.wrapper import Collector
#populate store
store = UniprotStoreDummy()
load_data = store.load_uniprot_xml(file="path/to/uniprotkb_proteome_UP000000625.xml")
store.save_collection('ecoli_K12', load_data)
#retrieve collection
my_collection = Collector(store, 'ecoli_K12')
```
This collection is iterable, slice-able or get-able via uniprot Accession numbers.
```python
print(my_collection['P19636'])
for uniprot_datum in my_collection:
    print(uniprot_datum.id, uniprot_datum.go)
```


### Assign whole proteome to the GO structure
Each protein of the collection now has to be attached to the GO terms that are described in its `UniprotDatum` `go` field (see above).
```python
obogo_tree.load_proteins('background', my_collection)
```
The `obogo_tree` represents each GO term as a straight newtworkx node. The `load_proteins` call will set the  for each node the value of their 'background' key to a list of UniprotDatum.

In most ORA analysis the population of proteins attached to a given node is the union of the proteins attached to its descendant ("GO annotation goes up": "any specific GO term implicitly carries the meaning of a less specific"). Hence, an additional operation is required to propagate protein populations up the tree.
```python
obogo_tree.percolate(percol_type="background")
```
This currently takes around 2mn for E.Coli proteome.

###### NB: At this stage, serializing the obogo_tree could be handy
```python
import pickle
pik_fpath = "obogo_ecoliK12.pik"
with open(pik_fpath, "wb") as fp:
    pickle.dump(obogo_tree, fp)
```
### Load the experimental protein set
For this tutorial, we will create a dummy collection of experimental proteins based on a slice of 1200 protein from the proteome and load it into obogo_tree. Note that this time, it is loaded using the `'measured'` argument. Then, we also propagate this additional protein population up the tree.
```python
obo_tree.load_proteins('measured', my_collection[1000:2200])
obo_tree.percolate(percol_type='measured')
```

### Define sample protein set
We now define a subset of `measured` proteins as "of interest" (aka: over-abundant).

```python
sample = my_collection[1100:1180]
```
#### Compute ORA of the GO terms within the sample
For a particular GO term
```python
from obogo.statistics import compute_node_ora
print( compute_node_ora(obo_tree, sample, 'GO:0006811') )
print( compute_node_ora(obo_tree, sample, 'metal ion transport') )
print( compute_node_ora(obo_tree, sample, 'GO:0006811', norm='measured') )
```
The sample parameter can also be a straight Uniprot AC iterator (eg: `['P02930', 'P03819', 'P0A910']`).
The `norm`` parameter controls the reference population for the Fisher statistic:
- 'background' : the whole proteome (default)
- 'measured'   : the proteins of the experiment

The returned value is a tuple of the form:

`(GO_identfier, GO_name, total_sample protein carrying the GO term, Fisher test log_odd, Fisher test pvalue, contingency table)`

A similar operation can be applied to the entiere tree, where a generator of `score` tuples will be returned
```python
from obogo.statistics import score_ora_tree
for go_score in score_ora_tree(obo_tree, sample):
    print(go_score)
```
