#!/bin/sh

#make .sample file to use as reference panel in shapeit
echo "sample population group sex" > $1
awk '{ print $1, $1, "POP", "1" }' $2 >> $1
