This packages contains classes related to the correction of transcription activator-like effectors (TALEs) in ONT assemblies of *Xanthomonas* bacteria.

The implementation of our new approach "TALEcorrection" consists of three sub-tools: 
- *PrepareTALEparts* implements prepartion of TALE DNA sequences to build TALE HMMs for TALE parts N-terminus, repeat and C-terminus.
- *CorrectTALESequences* implements correction of TALEs in ONT assemblies based on TALE HMMs.
- *PolishTALESubstitutions* implemtents substitution polishing based on mapping of reads to the corrected assembly.

To correct TALE-containing ONT assemblies with "TALEcorrection", please use the bash scripts provided under: XXX and configure the scripts according to your purpose. For details read the contained "help.txt" file.

First, it is necessary to install the following tools: HMMER, clustalo and igvtools:
- HMMER, see: http://hmmer.org/documentation.html or with conda: conda install -c bioconda hmmer 
- clustalo, see: http://www.clustal.org/omega/ or with conda: conda install -c bioconda clustalo 
- igvtools, see https://software.broadinstitute.org/software/igv/igvtools_commandline or with conda: conda install -c bioconda igvtools 

If you use TALEcorrection, please refer to preprint:

A.Erkes, R. Grove, M. Žarković, S. Krautwurst, R. Koebnik, R. D. Morgan, G. G. Wilson, M. Hölzer, M. Marz, J. Boch, J. Grau. Assembling highly repetitive Xanthomonas TALomes using Oxford Nanopore sequencing. *bioRxiv*, 2022.

If you experience problems using TALEcorrection, please contact_ us.

.. _contact: mailto:grau@informatik.uni-halle.de
