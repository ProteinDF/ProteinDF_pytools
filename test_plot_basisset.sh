#!/bin/bash

export PDF_HOME=${HOME}/local/intel/ProteinDF.openmpi
export PATH=${PDF_HOME}/bin:${PATH}

export PYTHONPATH=${HOME}/local/intel/ProteinDF.openmpi/lib/python2.7/site-packages/


${PDF_HOME}/bin/pdf-plot-basisset.py "${1}"

