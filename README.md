# Genome-Assembler

A basic assembler to assemble a genome from a set of error-prone reads. This project is based on the [Genome Assembly Programming Challenge](https://www.coursera.org/learn/assembling-genomes) capstone project on Coursera.

The assembler work as follows:
* It breaks given reads into k-mers (with k = 20)
* Build De Bruijn graph from these k-mers
* Remove tips from De Bruijn graph
* Remove bubbles from De Bruijn graph to account for sequencing errors
* Generate contigs of the original genome by finding maximal non-branching paths and isolated cycles in the De Bruijn graph
