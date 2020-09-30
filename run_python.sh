#!/bin/bash
set -e

#TODO don't use this script (Find better way to package dependencies)
cd $(dirname "$0")
source ./.venv/bin/activate
python $1 $2