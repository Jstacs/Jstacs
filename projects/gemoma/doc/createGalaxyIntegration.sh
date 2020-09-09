#!/bin/bash

# This script allows to create the XML files needed for an Galaxy integration as we typically use it.
# Especially the VM arguments might be adjusted for larger genomes, ...
#
# The script has a single parameter which is the version of GeMoMa.

jar=$(eval "ls GeMoMa-*.jar")

java -jar $jar --create
java -jar $jar --create GeMoMa -Xms5G -Xmx40G
java -jar $jar --create GeMoMaPipeline -Xms5G -Xmx40G
java -jar $jar --create GAF -Xms5G -Xmx40G


cp Extractor.xml Extractor.xml-default
cp ERE.xml ERE.xml-default
mv GeMoMaPipeline.xml GeMoMaPipeline.xml-default
grep -v restart GeMoMaPipeline.xml-default | grep -v "_PATH" > GeMoMaPipeline.xml

function change {
	sed -i "s/<\/command>/\n\t\t#if \$$3_ps_$2:\n\t\t\t\$$2\n\t\t\#end if\n<\/command>/g" $1
	sed -i "s/<\/outputs>/<data format=\"$4\" name=\"$2\" label=\"\#if str\(\$$3_jobname\) == \'\' then \$tool\.name + \' on \' + \$on_string else \$$3_jobname\#: $2\">\n\t<filter>$3_ps_$2<\/filter>\n<\/data>\n<\/outputs>/g" $1
}


array=("proteins" "cds" "genomic")
for a in "${array[@]}"
do
	change Extractor.xml $a Extractor fasta
done

array=("predicted_proteins" "predicted_CDSs" "predicted_genomic_regions")
for a in "${array[@]}"
do
	change GeMoMaPipeline.xml $a GeMoMa_pipeline fasta
done

change ERE.xml coverage Extract_RNA_seq_Evidence bedgraph
sed -i "s/coverage<\/filter>/coverage and Extract_RNA_seq_Evidence_ps_Stranded=='FR_UNSTRANDED'<\/filter>/g" ERE.xml
sed -i "s/\$coverage/\#if \$Extract_RNA_seq_Evidence_ps_Stranded=='FR_UNSTRANDED':\n\t\t\t\t\$coverage\n\t\t\t\#else\n\t\t\t\t\$coverage_forward \$coverage_reverse\n\t\t\t\#end if/g" ERE.xml
sed -i "s/<\/outputs>/<data format=\"bedgraph\" name=\"coverage_forward\" label=\"\#if str(\$Extract_RNA_seq_Evidence_jobname) == \'\' then \$tool.name + \' on \' + \$on_string else \$Extract_RNA_seq_Evidence_jobname\#: coverage (forward strand)\">\n\t<filter>Extract_RNA_seq_Evidence_ps_coverage and Extract_RNA_seq_Evidence_ps_Stranded!=\'FR_UNSTRANDED\'<\/filter>\n<\/data>\n<data format=\"bedgraph\" name=\"coverage_reverse\" label=\"\#if str(\$Extract_RNA_seq_Evidence_jobname) == \'\' then \$tool.name + \' on \' + \$on_string else \$Extract_RNA_seq_Evidence_jobname\#: coverage (reverse strand)\">\n\t<filter>Extract_RNA_seq_Evidence_ps_coverage and Extract_RNA_seq_Evidence_ps_Stranded!=\'FR_UNSTRANDED\'<\/filter>\n<\/data>\n<\/outputs>/g" ERE.xml


sed -i "s/<command>/<command>\n\t/g" *.xml