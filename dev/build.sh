#!/bin/bash

echo "Rebuilding CHERAB extension modules (in place)..."
python setup.py build_ext --inplace $1 $2 $3 $4 $5

echo

echo "Rebuilding CHERAB F2PY extension modules (in place)..."
python setup_f2py.py build_ext --inplace $1 $2 $3 $4 $5

