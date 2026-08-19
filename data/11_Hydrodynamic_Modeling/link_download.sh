#!/bin/sh
for file in $(cat $1); do
    echo "Processing $file"
    wget $file
done
