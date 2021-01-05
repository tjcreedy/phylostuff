#! /usr/bin/env bash
if (( $EUID == 0 )); then
        find *.R -maxdepth 1 -executable -type f | while read l; do rm /usr/local/bin/$l; echo "Removed symlink for $l"; done
else
    echo "Run as root to remove symlinks"
fi

