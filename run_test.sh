#!/bin/bash

TESTS="\
  basis2 \
  "

for i in ${TESTS}; do
    echo "test ${i}"
    python -m unittest tests.test_${i} 2>&1 | tee out.test_${i}
    echo "done."
    echo
done


