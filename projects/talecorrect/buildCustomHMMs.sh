#!/bin/bash
# Run this script, if you want to build custom HMMs for TALEs of your specified pathovar 
# If you have an own training data set with TALEs of that pathovar, please uncomment line 26 pathToTALEsFasta=$1 #($1 should contain path to a file with TALE DNA sequences in fastA format)
# Alternatively you can use TALE sequences stored in AnnoTALE, so you need to change line 16 customSpeciesName="Xcc" to the matching pathovar. Maybe you should load current data set of All_TALE_DNA_sequences.fasta from AnnoTALE first.
# Save the changes and run this script and find your custom HMMs in the folder "HMMs/custom".

### Please install HMMER, see: http://hmmer.org/documentation.html and clustalo, see: http://www.clustal.org/omega/ ###

#set output folder  
pathToHMMs="HMMs/custom/"

#set path to input TALE DNA sequences (You can load current TALE sequences from AnnoTALE: http://www.jstacs.de/index.php/AnnoTALE)
pathToAnnoTALE_AllTALEsFasta=${pathToHMMs}All_TALE_DNA_sequences.fasta

#set customSpeciesName, only if you want to use TALE sequences from AnnoTALE
customSpeciesName="Xcc"

#Delete pseudo genes in AnnoTALE TALE sequences
sed '/Pseudo/{N;d;}' ${pathToAnnoTALE_AllTALEsFasta} >${pathToAnnoTALE_AllTALEsFasta%.fasta}_withoutPseudo.fasta

#filter for ${customSpeciesName} in AnnoTALE TALE sequences
grep -A1 ${customSpeciesName} ${pathToAnnoTALE_AllTALEsFasta%.fasta}_withoutPseudo.fasta | grep -v "\-\-" - >${pathToHMMs}${customSpeciesName}_TALEsAnnoTALE.fasta

#set path to TALE DNA sequences, that will be used as training data set to build TALE HMMs 
pathToTALEsFasta=${pathToHMMs}${customSpeciesName}_TALEsAnnoTALE.fasta # use this path, in case you want to use TALE sequences from AnnoTALE 
#pathToTALEsFasta=$1 #uncomment this line to provide a path to a custom TALE data set in fastA format, that should be used for hmmbuild

pathToAnnoTALE="AnnoTALEcli-1.5.jar" #Path to command line version of AnnoTALE. You can download current release here: http://www.jstacs.de/index.php/AnnoTALE

#### run AnnoTALE to split TALEs in parts N-terminus, C-terminus and repeats ####
echo "...run AnnoTALE..."
java -jar ${pathToAnnoTALE} analyze t=${pathToTALEsFasta} outdir=${pathToHMMs}output

pathToTALEsParts=${pathToHMMs}output/TALE_DNA_parts.fasta

#### prepare TALE dna parts files for hmmbuild ####
echo "...prepare TALE DNA parts files..."
java -jar  TALEcorrection.jar prepare t=${pathToTALEsParts} outdir=${pathToHMMs}output

#### run clustalo and hmmbuild ####
echo "...run clustalo and hmmbuild..."
parts=("N-terminus.10bpRepeat1" "repeat" "C-terminus")
for pattern in ${parts[@]}; do	
	#MSA with clustalo
	clustalo -i ${pathToTALEsParts%.fasta}.${pattern}.fasta -o ${pathToTALEsParts%.fasta}.${pattern}.clustalo.clw --outfmt=clu 
	
	#build HMMs
	hmmbuild ${pathToHMMs}${pattern}.hmm ${pathToTALEsParts%.fasta}.${pattern}.clustalo.clw
done
