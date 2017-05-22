#!/bin/bash
version="1.4";

echo "============================================================================="
echo "check versions:"
echo ""
java -version
echo ""
tblastn -version
echo ""
echo "============================================================================="
echo ""

java -jar GeMoMa-${version}.jar CLI Extractor a=test_data/annotation.gff g=test_data/reference.fasta outdir=results/precomputed/

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""

java -jar GeMoMa-${version}.jar CLI GeMoMa t=test_data/truncated_tblastn.txt c=test_data/precomputed_cds-parts.fasta a=test_data/precomputed_assignment.txt tg=test_data/contig.fasta  outdir=results/precomputed/