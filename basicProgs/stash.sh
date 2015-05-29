#!/bin/bash 

usage() 
# Usage statement for when things go wrong 
{ 
    echo "stash - stash a copy of the program, expects a folder named .stash 
usage:
    stash <fileToStash.ext> 
Stashes the current file with a date and time stamp in a hidden folder .stash. 
Currently .stash must already exist, comments are not supported but may be 
coming soon ;) " 1>&2
}
# Spit usage when no arguments are given
if [ $# -lt 1 ]; then
    usage 
    exit 255
fi

# Check if there is a stash, make one if there is not
export DIRECTORY=".stash"
if [ ! -d "$DIRECTORY" ]; then
    mkdir "$DIRECTORY" 
fi

# Check that the file given is a valid file, spit usage otherwise
export fileToStash="$1" 
if [ ! -s "${fileToStash}" ]; then 
    echo "ERROR: can not read fileToStash.ext file: $fileToStash" 1>&2 
    usage
fi

cd .stash
echo ${fileToStash} > aRandomFileName1
date >> aRandomFileName1 
cat aRandomFileName1 | tr -d ' \t\n\r\f' > aRandomFileName2
echo  >> aRandomFileName2
file="aRandomFileName2"
while IFS= read -r line
do 
    cp "${fileToStash}" "$line"
done < "$file" 
rm aRandomFileName1
rm aRandomFileName2
echo "The file ${fileToStash} was succesfully stashed" 
