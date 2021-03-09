# GO_term_grapher
Graphs logfold across gene oncology terms. One script will graph genes across a term passed in as a string, the other will graph several terms from a list of terms passed in as an input file.
# USAGE: 
* To compare individual genes within a GO term
python ./graph_genes_within_goterm.py <GO_TERM(string)> <first deseq2> <second deseq2>
* To compare several GO terms
python ./graph_several_goterms.py <filename for a list of GO terms> <first deseq2> <second deseq2>
  
I've attached some small, made up data to run a few examples
python ./graph_genes_within_goterm.py GO:0042393 dummydata1 dummydata2
python ./graph_genes_within_goterm.py GO:00044012 dummydata3 dummydata4
python ./graph_several_goterms.py several_goterms.txt dummydata3 dummydata4

If you include a goterm and it doesn't have any genes associated with it, the graph will segfault. 
This was one of my first projects and is about a year old. Hoping to revise it this summer.
