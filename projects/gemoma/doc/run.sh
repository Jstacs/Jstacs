#!/bin/bash
jar=$(eval "ls GeMoMa-*.jar")

# This script allows to run GeMoMa with minimal work from command line. 
#
# A simple example without RNA-seq
# using tblastn is
# ./run.sh tblastn tests/gemoma/target-fragment.fasta tests/gemoma/ref-annotation.gff tests/gemoma/ref-fragment.fasta results/sw-tblastn
#
# or using mmseqs is
# ./run.sh mmseqs tests/gemoma/target-fragment.fasta tests/gemoma/ref-annotation.gff tests/gemoma/ref-fragment.fasta results/sw-mmseqs
#
# If you like to run GeMoMa with RNA-seq data you have to run the script with the following parameters
# ./run.sh <search> <reference annotation> <reference geneome> <out-dir> <target genome> <library type> <mapped reads>
# where library type is one of {FR_UNSTRANDED, FR_FIRST_STRAND, FR_SECOND_STRAND} and mapped reads is a SAM or BAM file.
#
# In both cases, the final prediction is located in ${out}/filtered_predictions.gff.

#parameters
threads=1;

search=$1
target=$2
annotation=$3
reference=$4
out=$5

if [ $# -ne 5 ]; then
	echo "GeMoMa using RNA-seq data: library type=" $6 "mapped reads=" $7
	lib=$6;
	reads=$7;
else 
	echo "GeMoMa without RNA-seq data"
fi



echo "============================================================================="
echo ""
echo "check versions:"
echo ""

java -version

echo ""
if [ "$search" == "tblastn" ]
then
	tblastn -version
	sort=false;
	score="Trust"
else
	echo "mmseqs"
	mmseqs version
	sort=true;
	score="ReAlign"
fi

echo ""
echo "============================================================================="
echo ""

if [ $# -ne 5 ]; then
	echo "ERE:"
	echo ""

	java -jar $jar CLI ERE c=true s=${lib} m=${reads} outdir=${out}

	echo ""
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	echo ""
	#no Denoise so far
fi
echo "Extractor:"
echo ""

java -jar $jar CLI Extractor a=${annotation} g=${reference} p=true Ambiguity=AMBIGUOUS outdir=${out}

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""

if [ "$search" == "tblastn" ]
then
	echo "makeblastdb:"
	echo ""
	
	makeblastdb -out ${out}/blastdb -hash_index -in ${target} -title "target" -dbtype nucl
	
	echo ""
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	echo ""
	echo "tblastn:"
	echo ""
	
	tblastn -query ${out}/cds-parts.fasta -db ${out}/blastdb -evalue 100.0 -out ${out}/search.txt -outfmt "6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles" -db_gencode 1 -matrix BLOSUM62 -seg no -word_size 3 -comp_based_stats F -gapopen 11 -gapextend 1
else
	echo "create mmseqs dbs:"
	echo ""
	
	mmseqs createdb ${target} ${out}/mmseqsdb --dont-split-seq-by-len -v 2
	mkdir ${out}/ref
	mmseqs createdb ${out}/cds-parts.fasta ${out}/ref/mmseqsdb -v 2
	
	echo ""
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	echo ""
	echo "mmseqs:"
	echo ""

	mmseqs search ${out}/ref/mmseqsdb ${out}/mmseqsdb ${out}/ref/mmseqsdb_align.out ${out}/ref/mmseqsdb_tmp -e 100.0 --threads ${threads} -s 8.5 -a --comp-bias-corr 0 --max-seqs 500 --mask 0 --orf-start-mode 1 -v 2

	echo ""
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	echo ""
	echo "convertalis:"
	echo ""

	mmseqs convertalis ${out}/ref/mmseqsdb ${out}/mmseqsdb ${out}/ref/mmseqsdb_align.out ${out}/search.txt --threads ${threads} --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen" -v 2
fi

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo "GeMoMa:"
echo ""

if [ $# -eq 5 ]; then
	java -jar $jar CLI GeMoMa s=${out}/search.txt c=${out}/cds-parts.fasta a=${out}/assignment.tabular q=${out}/proteins.fasta t=${target} sort=${sort} Score=${score} outdir=${out}
else
	if [ ${lib} == "FR_UNSTRANDED" ]; then
echo "UNSTRANDED"
		java -jar $jar CLI GeMoMa s=${out}/search.txt c=${out}/cds-parts.fasta a=${out}/assignment.tabular q=${out}/proteins.fasta t=${target} sort=${sort} Score=${score} outdir=${out} i=${out}/introns.gff coverage=UNSTRANDED coverage_unstranded=${out}/coverage.bedgraph
	else
echo "STRANDED"
		java -jar $jar CLI GeMoMa s=${out}/search.txt c=${out}/cds-parts.fasta a=${out}/assignment.tabular q=${out}/proteins.fasta t=${target} sort=${sort} Score=${score} outdir=${out} i=${out}/introns.gff coverage=STRANDED coverage_forward=${out}/coverage_forward.bedgraph coverage_reverse=${out}/coverage_reverse.bedgraph
	fi
fi

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo "GAF:"
echo ""

java -jar $jar CLI GAF g=${out}/predicted_annotation.gff outdir=${out}

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo "AnnotationFinalizer:"
echo ""

java -jar $jar CLI AnnotationFinalizer a=${out}/filtered_predictions.gff outdir=${out} rename=NO