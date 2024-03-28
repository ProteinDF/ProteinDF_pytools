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
  "

for i in ${UNITTESTS}; do
    echo "test ${i}"
    python -m unittest tests.test_${i} 2>&1 | tee out.test_${i}
    echo "done."
    echo
done
