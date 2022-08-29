#!/bin/bash
# Run this script, if you want to correct TALEs in ONT assemblies.
# You can run correction of the example genome or uncomment parameters outputFolder=$1 (line 7) and uncorrectedAssembly=$2 (line 11) to correct your own assembly. Please change the pathToHMMs="HMMs/Xoo/" to the matching HMM pathovar.

#set output folder 
outputFolder="outputTest/" 
#outputFolder=$1

#set path to uncorrected ONT assembly
uncorrectedAssembly="inputTest/uncorrectedAssembly.fasta" 
#uncorrectedAssembly=$2

#set path to TALEcorrection.jar
pathToTALECorrection="TALEcorrection.jar"

#set path to folder with TALE HMMs
pathToHMMs="HMMs/Xoo/" #or Xoc: #pathToHMMs="HMMs/Xoc/" #or custom: #pathToHMMs="HMMs/custom/"

#### run nHMMER ####
echo "...start nhmmer..."
parts=("N-terminus.10bpRepeat1" "repeat" "C-terminus")
for p in ${parts[@]}; do	
	nhmmer ${pathToHMMs}${p}.hmm ${uncorrectedAssembly} >${outputFolder}out_nhmmer.${p}.txt
done

#### correction of TALEs #####
echo "...start correction of TALEs..."
java -jar ${pathToTALECorrection} correct s=${uncorrectedAssembly} n=${outputFolder}out_nhmmer.${parts[0]}.txt r=${outputFolder}out_nhmmer.${parts[1]}.txt c=${outputFolder}out_nhmmer.${parts[2]}.txt outdir=${outputFolder}

#### if you want to polish substitutions #####
echo "...finished..."
echo "If you want to polish substitutions: map reads against corrected assembly and continue with runSubPolish.sh"
