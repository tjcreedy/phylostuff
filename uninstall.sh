#! /usr/bin/env bash
if (( $EUID == 0 )); then
   find *.R -maxdepth 1 -executable -type f | while read l; do rm /usr/local/bin/$l; echo "Removed symlink for $l" done
else
   sed -i "/^export .*${PWD}$/d" ~/.bashrc
   echo "PATHs removed, you may need to log out and in again for this to take effect"
fi

