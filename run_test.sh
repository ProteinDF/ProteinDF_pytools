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
  plot_elevel \
  plot_scripts \
  archive_scripts \
  input_scripts \
  info_scripts \
  pop_scripts \
  report_run_scripts \
  tools_small \
  roks_run_types \
  pdfparam_object \
  open_archive \
  h5_dash_d_basic \
  h5_dash_d_rest \
  "

for i in ${UNITTESTS}; do
    echo "test ${i}"
    python -m unittest tests.test_${i} 2>&1 | tee out.test_${i}
    echo "done."
    echo
done
