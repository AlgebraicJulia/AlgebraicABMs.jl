#!/bin/bash

DEFAULT_REPO='AlgebraicABMs'
DEFAULT_UUID='5A5E3447-9604-46E6-8D91-CB86F5F51721' 

usage="This script is for initializing the template with the new repository name and UUID. Please provide the new repository name and UUID in that order. The repository name cannot be 'Test.'\n
Example:\n
./init.sh ${DEFAULT_REPO} ${DEFAULT_UUID}"

REPO=${1:-"${PWD##*/}"}
UUID=${2:-$(uuidgen)}

# set to lowercase
UUID=${UUID,,}

if [ ! $REPO ] || [ "$REPO" = 'Test' ] || [ ! $UUID ]; then
  echo ""
  printf "$usage" 
  exit 1
fi

read -p "By continuing, the following substitutions will be made:

REPO: $DEFAULT_REPO => $REPO
UUID: $DEFAULT_UUID => $UUID

Are you sure? [y/N]" -n 1 -r -s
echo

if [[ ! $REPLY =~ ^[Yy]$ ]]
then 
  exit 1
fi

echo "Doing the thing..."

# get version
unameOut="$(uname -s)"

case "${unameOut}" in
  Linux*)  
    git grep -l $DEFAULT_REPO | xargs sed -i "s/${DEFAULT_REPO}/${REPO}/g";
    git grep -l $DEFAULT_UUID | xargs sed -i "s/${DEFAULT_UUID}/${UUID}/g";;
  Darwin*) 
    git grep -l $DEFAULT_REPO | xargs sed -i '' -e "s/${DEFAULT_REPO}/${REPO}/g";
    git grep -l $DEFAULT_UUID | xargs sed -i '' -e "s/${DEFAULT_UUID}/${UUID}/g";;
  *)       
    echo UNKNOWN:${unameOut};; 
esac

# rename 
if [[ $REPO == *.jl ]]
then
  mv src/$DEFAULT_REPO.jl src/$REPO
else
  mv src/$DEFAULT_REPO.jl src/$REPO.jl
fi

read -p "Would you like this script to add, commit, and push the new changes? [y/N]" -n 1 -r -s
echo

if [[ ! $REPLY =~ ^[Yy]$ ]]
then
  exit 1
fi

git commit -am "Initialized $REPO"
git push
