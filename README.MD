This project for COP5725-Fall 2023 includes the reimplementation of the BD-Anchor in Java mentioned in the paper "Text Indexing for Long Paterns: Anchors are All you Need" by Lorraine A. K. Ayad, Grigorios Loukides, and Solon P. Pissis.

* We have retained the original implementation but erased the core algorithms, then utilized our own re-implmeentation.

In order to run, first run the pre-install.sh shell script from BDA-Index_I or BDA-Index_II, this is provided by the original author to set up the environments.
Then, run command "make -f Makefile.64-bit.gcc". This is also provided by the original author's github page
In order to run the program, first move into psascan and sparsePhi of BDA-Index_I or BDA-Index_II and run the makefile for each folder, then return back to BDA-Index_I or BDA-Index_II and use command in the following form: 
Usage: 
./index <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size> <index_filename>

<text_file> - name of input text file
<ell> - minimum size of pattern to consider searching for within text. 
<pattern_file> - name of input file containing patterns
<output_filename> - name of output file where pattern occurrences will be placed.
<ram_use> - ram usage for external SA and LCP
<block_size> - size of block size b to use.
<index_filename> - name of the index file to be used if they exist otherwise to be created.

Example:
 $ ./bda-index_I ./data/text 3 ./data/patterns out 150 10 index
 $ ./bda-index_II ./data/text 3 ./data/patterns out 150 10 index

 (Info taken from original author's github page: https://github.com/lorrainea/BDA-index/tree/main)