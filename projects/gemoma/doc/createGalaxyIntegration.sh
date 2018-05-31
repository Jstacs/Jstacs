#!/bin/bash

# This script allows to create the XML files need for an Galaxy integration as we typically use it.
# Especially the VM arguments might be adjusted for larger genomes, ...
#
# The script has a single parameter which is the version of GeMoMa.

java -jar GeMoMa-$1.jar --create
java -jar GeMoMa-$1.jar --create GeMoMa -Xms5G -Xmx40G

cp Extractor.xml Extractor.xml-default

sed -i 's/<\/command>/\n#if \$Extractor_ps_proteins:\n\t\$proteins\n\#end if\n#if \$Extractor_ps_cds:\n\t\$cds\n\#end if\n<\/command>/g' Extractor.xml
sed -i "s/<\/outputs>/<data format=\"fasta\" name=\"proteins\" label=\"\#if str\(\$Extractor_jobname\) == \'\' then \$tool\.name + \' on \' + \$on_string else \$Extractor_jobname\#: proteins\">\n\t<filter>Extractor_ps_proteins<\/filter>\n<\/data>\n\n<data format=\"fasta\" name=\"cds\" label=\"\#if str\(\$Extractor_jobname\) == \'\' then \$tool\.name + \' on \' + \$on_string else \$Extractor_jobname\#: cds\">\n<filter>Extractor_ps_cds<\/filter>\n<\/data>\n<\/outputs>/g" Extractor.xml

cp ERE.xml ERE.xml-default

sed -i 's/<\/command>/\n#if \$Extract_RNA_seq_Evidence_ps_coverage_output:\n\t\$coverage\n\#end if\n<\/command>/g' ERE.xml
sed -i "s/<\/outputs>/<data format=\"bedgraph\" name=\"coverage\" label=\"\#if str\(\$Extract_RNA_seq_Evidence_jobname\) == \'\' then \$tool\.name + \' on \' + \$on_string else \$Extract_RNA_seq_Evidence_jobname\#: coverage\">\n\t<filter>Extract_RNA_seq_Evidence_ps_coverage_output<\/filter>\n<\/data>\n<\/outputs>/g" ERE.xml
