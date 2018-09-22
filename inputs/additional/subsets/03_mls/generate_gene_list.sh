#! /bin/bash

cd generated_files/under_version_control/genes

grep -i macrolide * | grep flash_key | awk -F'flash_key:' '{print $2}' | awk -F'|' '{print $1}' | sort -f
