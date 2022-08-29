#!/bin/bash
# Run this script, if you want to polish substitutions based on mapping of reads to the corrected assembly. You need to map reads to corrected assembly before you start this bash script.
#Please install igvtools first: https://software.broadinstitute.org/software/igv/igvtools_commandline

#set output folder
outputFolder="outputTest/"

#set path to corrected ONT assembly
assembly="outputTest/correctedTALEs.fa"

#set path to TALEcorrection.jar
pathToTALECorrection="TALEcorrection.jar"

#set path to sorted bam file with mapped reads on corrected assembly
bamFile="inputTest/nanoSampledReads_mappedToassembly.correctedTALEs.sorted.bam"

#set path to igvtools output folder
pathFiles=${outputFolder}igvtools_count_output/
mkdir -p $pathFiles

####igvtools count##### 
echo "... run igvtools count ..."
bash ${outputFolder}igvtools.count.Substitution.script.sh ${assembly} ${bamFile} ${pathFiles}

cat ${pathFiles}out.igvtools.count.*.wig > ${pathFiles}all.wig

###polish corrected####
echo "... start substitution polishing ..."
java -jar ${pathToTALECorrection} polish c=${assembly} i=${pathFiles}all.wig outdir=${outputFolder}
