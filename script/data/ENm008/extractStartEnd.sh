#!/bin/bash - 

sed -e '1,1d' Nature_2010_ContactMap_ENm008_K562.txt  | awk '{print $1}' | awk -F: '{print $2}' | awk -F- '{print $1"\t"$2}' > ENm008_K562_StartEnd.txt
sed -e '1,1d' Nature_2010_ContactMap_ENm008_GM12878.txt  | awk '{print $1}' | awk -F: '{print $2}' | awk -F- '{print $1"\t"$2}' > ENm008_GM12878_StartEnd.txt

