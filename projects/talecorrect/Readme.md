# General information

This packages contains classes related to the correction of transcription activator-like effector (TALE) genes in ONT-based assemblies of *Xanthomonas* bacteria.

The implementation of our new approach "TALEcorrection" consists of three sub-tools:
- *PrepareTALEparts* implements prepartion of TALE DNA sequences to build TALE HMMs for TALE parts N-terminus, repeat and C-terminus,
- *CorrectTALESequences* implements correction of TALEs in ONT assemblies based on TALE HMMs,
- *PolishTALESubstitutions* implemtents substitution polishing based on mapping of reads to the corrected assembly.
All three tools are implemented in separate Java classes, and are combined in a main class TALEcorrection.

To correct TALE-containing ONT assemblies with "TALEcorrection", we provide bash scripts
- *buildCustomHMMs.sh* for building custom hidden Markov models (HMMs, from HMMer) from FastA files of TALE DNA sequences as, for instance, generated by the "Load and View TALE Classes" tools of [AnnoTALE](https://www.jstacs.de/index.php/AnnoTALE),
- *runTALEcorrection.sh* for running the correction of TALE genes on an ONT assembly using custom or pre-defined HMMs,
- *runSubPolish.sh* for running the substitution polishing of previously corrected TALE sequences based on mapped ONT reads (as BAM file).
The scripts can be configured using internal variables according to your purpose. 

For details see the contained "help.txt" file.

# Requirements

The TALE correction scripts require the external tools HMMER, clustalo and igvtools:
- HMMER, see: http://hmmer.org/documentation.html or install using conda: `conda install -c bioconda hmmer`
- clustalo, see: http://www.clustal.org/omega/ or install using conda: `conda install -c bioconda clustalo`
- igvtools, see https://software.broadinstitute.org/software/igv/igvtools_commandline or install using conda: `conda install -c bioconda igvtools`

In addition, the Java code needs a Java Runtime Environment.

The recommended way of installing all dependencies is within a separate conda environment:

`conda create -n talecorrect hmmer clustalo igvtools openjdk`
    

followed by

`conda activate talecorrect`


# Binary version packed with test data

We provide an [archive](https://www.jstacs.de/downloads/TALECorrection_scripts.zip) with pre-compiled Java code, all scripts, and test data.
You can test the correct installation of dependencies and provided scripts by calling those scripts in the succession


    bash buildCustomHMMs.sh
    bash runTALEcorrection.sh
    bash runSubPolish.sh


After all scripts terminated successfully, the directory *outputTest* should contain (among other files) a FastA file *correctedTALEs.fa* containing the complete input assembly (cf. inputTest/uncorrectedAssembly.fasta) but with TALE sequences corrected and polished.
Please note that corrected nucleotides are output in lower case, while all unaltered nucleotides are output in upper case.

# Tutorial

Here, we provide a brief tutorial of our TALE correction scripts a) with pre-defined HMMs and b) with custom HMMs. This assumes that you are working under Linux (or have a bash, grep and cat installed) have already installed all dependencies and downloaded the [TALE correction archive](https://www.jstacs.de/downloads/TALECorrection_scripts.zip) with binaries.

## Using pre-defined HMMs

Within the archive, we provide pre-defined HMMs for *Xanthomonas oryzae* pv. *oryzae* (*Xoo*) and *Xanthomonas oryzae* pv. *oryzicola* (*Xoc*), for correcting TALEs in genome assemblies of other *Xanthomonas* species, see the next section.

Open the script *runTALEcorrection.sh* in your favourite text editor and adapt the variables in the header to your needs. Specifically, you may want to adapt `uncorrectedAssembly` to the path to your assembly FastA and, if you work with *Xoc* instead of *Xoo*, `pathToHMMs` to `HMMs/Xoc`.

Run the script by calling

`bash runTALEcorrection.sh`

Afterwards, a file *correctedTALEs.fa* should be present in the directory *outputTest*, which contains the sequence of your assembly with the TALE sequences corrected. Please note that corrected nucleotides are output in lower case, while all unaltered nucleotides are output in upper case.
In addition, this script generated further files that are required for the next, optional, step of substitution polishing.

If you map ONT reads to the corrected assembly (*correctedTALEs.fa*), the nucleotides corrected in the previous step can further be polished based on experimental evidence of nucleotides.
To do so, open the script *runSubPolish.sh* in a text editor and adjust the variables in the header. Specifically, you may want to adapt `bamFile` to the path of the BAM file containing your mapped reads.

Run the script by calling

`bash runSubPolish.sh`

Afterwards, a file *polishedSubstitutions.fa* should be present in the directory *outputTest*, which contains the assembly with TALE sequences corrected and polished.

## Building custom HMMS for other strains

For building custom HMMs for your *Xanthomonas* species of interest, we provide an additional script *buildCustomHMMs.sh*.
This script may either use TALE sequences from [AnnoTALE](https://www.jstacs.de/index.php/AnnoTALE) or a custom FastA file containing the DNA sequences of your own collection of TALEs (DNA sequences of TALEs on a single line). Updated versions of the file *HMMs/custom/All_TALE_DNA_sequences.fasta* can be obtained from [AnnoTALE](https://www.jstacs.de/index.php/AnnoTALE) using the "Load and View TALE Classes" tool and replacing this file with the downloaded version.

Open the script *buildCustomHMMs.sh* in your favourite text editor and adapt the variables in the header to your needs. Specifically, you may want to change `pathToAnnoTALE_AllTALEsFasta` if you use a custom file with TALE DNA sequences and/or `customSpeciesName` if you want to filter sequences based on the FastA header.

### a) Custom species, provided TALE sequences

For instance, if you would like to build custom HMMs for *Xanthomonas translucens* based on the existing *All_TALE_DNA_sequences.fasta* file, you would need to adapt the definition of `customSpeciesName` to

`customSpeciesName="Xt"`

### b) Custom TALE sequences

Alternatively, if you would like to build custom HMMs from your own file (assumed to be stored in *path/to/myTALEs.fasta*), you would need to adapt the definition of `pathToTALEsFasta` to

`pathToTALEsFasta=path/to/myTALEs.fasta`.

You can comment out all lines before except the line where the output folder is set: `pathToHMMs="HMMs/custom/"`.

### Run the script after option a) or b)

Then, run the script by calling

`bash buildCustomHMMs.sh`

Afterwards, the directory *HMMs/custom* should contain three HMM files *C-terminus.hmm*, *N-terminus.10bpRepeat1.hmm* and *repeat.hmm* with your custom HMMs.

To use your custom HMMs for TALE correction, you now need to configure *runTALEcorrection.sh* accordingly.

Open the script *runTALEcorrection.sh* in your favourite text editor and modify the variable `pathToHMMs` to

`pathToHMMs="HMMs/custom/"`

Besides, you can follow the description in the previous section "Using pre-defined HMMs" for the remaining steps.


# Citation and contact

If you use TALEcorrection, please refer to our [preprint](https://doi.org/10.1101/2022.08.17.504259):

A.Erkes, R. Grove, M. Žarković, S. Krautwurst, R. Koebnik, R. D. Morgan, G. G. Wilson, M. Hölzer, M. Marz, J. Boch, J. Grau. Assembling highly repetitive Xanthomonas TALomes using Oxford Nanopore sequencing. *bioRxiv*, 2022, doi: [10.1101/2022.08.17.504259](https://doi.org/10.1101/2022.08.17.504259)

If you experience problems using TALEcorrection, please [contact](mailto:grau@informatik.uni-halle.de) us.
