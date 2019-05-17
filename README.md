# Generalized network dismantling

This Python 2.7 code implements the Generalized Network Dismantling algorithm by Ren, Gleinig, Helbing, and Antulov-Fantulina, PNAS 116(14):6554-6559 (2019). It tries to follow the description in the paper rather than the published code by the authors.

The program can be run as: "python gnd.py network.adj [degree/unit]" where network.adj is an edgelist representation of the network and the second argument specifies if you want a degree-weighted or unit cost.

Output is as follows (per row): [accumulated deletion cost] [deleted node id] [size of largest connected component]

The randomness in the output comes from the code breaking ties by randomness. To make the output independent of the labeling of nodes.
