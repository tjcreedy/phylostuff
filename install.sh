#! /usr/bin/env bash
Rscript -e 'x<-c("ape", "getopt", "geiger", "taxize", "plyr");install.packages(x[!x %in% installed.packages()[,"Package"]], repos = "https://cloud.r-project.org")'
if (( $EUID == 0 )); then
    find *.R -maxdepth 1 -executable -type f | while read l; do ln -s $PWD/$l /usr/local/bin/$l; echo "Created symlink for $l"; done
else
    echo "export PATH=$PATH:$PWD" >> ~/.bashrc
    source ~/.bashrc
    echo "PATHs installed, you may need to log out and in again for this to take effect"
fi

