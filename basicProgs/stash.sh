#!/bin/bash 

usage() 
# Usage statement for when things go wrong 
{ 
    echo "stash - stash a copy of the program in a hidden folder named .stash.  
usage:
    stash <fileToStash.ext> "Comment inside quotes, white space is fine"  
Stashes the current file with a date and time stamp in a hidden folder .stash. 
Helpful for private/local source code management " 1>&2
}
# Spit usage when no arguments are given
if [ $# -lt 1 ]; then
    usage 
    exit 255
fi

# Check that the file given is a valid file, spit usage otherwise
export fileToStash="$1" 
if [ ! -s "${fileToStash}" ]; then 
    echo "ERROR: can not read fileToStash.ext file: $fileToStash" 1>&2 
    usage
    exit 255
fi

# Check if there is a stash, make one if there is not
export DIRECTORY=".stash"
if [ ! -d "$DIRECTORY" ]; then
    mkdir "$DIRECTORY" 
fi

# Go to the stash
cd .stash

# Prepare the file name, uses two dummy files in the process
echo ${fileToStash} > aRandomFileName1
date >> aRandomFileName1 
cat aRandomFileName1 | tr -d ' \t\n\r\f' | tr ':' '-' > aRandomFileName2
echo  >> aRandomFileName2
# Read the file name from the dummy file
file="aRandomFileName2"
while IFS= read -r line
do
    # Paste the comment if it exists
    if [ $# -gt 1 ]; then
        echo $2 > "$line"
    fi
    # Paste the rest of the file 
    cat "../${fileToStash}" >>  "$line"
done < "$file" 

# Remove the dummy files
rm aRandomFileName1
rm aRandomFileName2

# Tell the user that everything is ok
echo "The file ${fileToStash} was succesfully stashed" 
if [ $# -gt 1 ]; then
    echo "The comment provided for this stash is \"$2\"" 
fi
