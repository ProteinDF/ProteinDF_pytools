#!/bin/bash

# ======================================================================
# doctest
# ======================================================================
DOCTESTS=" \
  basis2 \
  basisset \
  gauparam \
  matrix \
  orbinfo \
  pdfcommon \
  pdfgraph \
  pdfmath \
  pdfparam \
  pdfsim \
  process \
  qmsim \
  vector \
  "

for i in ${DOCTESTS}; do
    echo ">>>> doctest: ${i}"
    python -m proteindf_tools.${i}
    echo "done."
done


# ======================================================================
# unit test
# ======================================================================
UNITTESTS="\
  basis2 \
  matrix \
  vector \
  basisset \
  basisset_parser \
  pdfmath \
  functions \
  orbinfo \
  pdfparam_basic \
  pdfparam_model \
  pdfparam_hdf5 \
  pdfarchive_sqlite3 \
  process \
  pdfcommon \
  gauparam \
  poputils \
  pdfsim \
  package \
  pdfgraph_base \
  pdfgraph_extra \
  report \
  pdf_test_eri \
  plot_scripts \
  "

for i in ${UNITTESTS}; do
    echo "test ${i}"
    python -m unittest tests.test_${i} 2>&1 | tee out.test_${i}
    echo "done."
    echo
done
