#!/bin/bash

DATA_DIR=$1

# Set the coverage threshold
min_coverage=3

# Run Docker with Bedtools, executing all commands in one session
docker run --rm -v "$(pwd):$(pwd)" -w "$(pwd)" -it quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6 /bin/bash -c "
    gzip -dc ${DATA_DIR}/data_coverage.per-base.bed.gz | \
    egrep -v 'HLA|decoy|random|alt|chrUn|chrEBV' | \
    awk -v OFS='\t' -v min_coverage=${min_coverage} '\$4 >= min_coverage { print }' | \
    bedtools merge -d 1 -c 4 -o mean -i - > ${DATA_DIR}/x3_coverage.bed;
"

echo "✅ Bedtools analysis completed successfully!"
