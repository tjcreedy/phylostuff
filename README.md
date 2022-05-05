# Phylostuff

This repository comprises various scripts for preparing, processing and manipulating phylogenetic 
trees and associated data. These are summarised here, see the documentation of each script using 
the help argument.

**phylabel.R** relabels the tips and/or nodes of a phylogenetic tree according to other supplied 
data, such as taxonomy and metadata

**treedentify.R** infers taxonomic classification for unknown terminals on a tree based on 
node taxonomisation performed by phylabel.R

**partitioner.py** expands basic partition tables as output by catfasta2phyml.pl into gene, codon 
and/or direction-based partitions in RAxML or nexus formats

**droptips.R** drops tips from a phylogeny

**phynames.R** reports the tip labels of a phylogenetic tree

**grafttree.R** grafts one phylogeny to a specified node in a host phylogeny with Open Tree of 
Life node IDs

**phylofuncs.R** provides a variety of R functions for working with phylogenetic trees in R

All of the above scripts apart from phylofuncs.R are designed for CLI usage.

## Installation

Download the specific script you need or git clone the repository. The R scripts require an up to 
date installation of R. 

The CLI R scripts are designed for Linux although will probably work on Mac. They are untested on 
Windows. They require the R libraries ape, getopt, geiger, taxize and plyr, you can install these
easily by running:
```
Rscript -e 'x<-c("ape", "getopt", "geiger", "taxize", "plyr");install.packages(x[!x %in% installed.packages()[,"Package"]], repos = "https://cloud.r-project.org")'
```

phylofuncs.R may require other libraries.

### Automated installation on Ubuntu Linux

The install.sh and uninstall.sh scripts provided can be used to automatically install the 
necessary R libraries and make the R scripts available on the PATH. This is very much a barebones 
method and may not work on many operating systems and configurations - it is only tested on Ubuntu 
Linux.

Clone the repository to a sensible location, e.g. ~/software/ for a local installation or /opt/ 
for installation for all users.
```
git clone https://github.com/tjcreedy/phylostuff.git
```

Run the installation script. Running as sudo will attempt to create symlinks in /usr/local/bin/,
otherwise the script will add a line to your .bashrc file to add the directory to your PATH.
```
cd phylostuff
bash install.sh
# or
sudo bash install.sh
```

To unistall, you **must** be in the phylostuff directory
```
bash uninstall.sh
# or
sudo bash uninstall.sh 
```
