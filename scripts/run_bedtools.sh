#!/bin/bash

DATA_DIR=$1
min_coverage=$2

# Run Docker with Bedtools, executing all commands in one session
docker run --rm -v "$(pwd):$(pwd)" -v ${DATA_DIR}:${DATA_DIR} -w "$(pwd)" -it quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6 /bin/bash -c "
    gzip -dc ${DATA_DIR}/data_coverage.per-base.bed.gz | \
    egrep -v 'HLA|decoy|random|alt|chrUn|chrEBV' | \
    awk -v OFS='\t' -v min_coverage=${min_coverage} '\$4 >= min_coverage { print }' | \
    bedtools merge -d 1 -c 4 -o mean -i - > ${DATA_DIR}/above_min_coverage.bed;
"

echo "✅ Bedtools analysis completed successfully!"
