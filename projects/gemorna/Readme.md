# GeMoRNA source code repository

GeMoRNA reconstructs genes and transcript models from mapped RNA-seq reads (in coordinate-sorted BAM format) and reports these in GFF format.
It is intended as a companion for the homology-based gene prediction program [GeMoMa](https://www.jstacs.de/index.php/GeMoMa).
In a typical workflow, predictions of transcript models may be obtained from GeMoRNA for a collection of BAM files individually and subsequently merged using the GeMoMa Annotation Filter (GAF). Optionally, homology-based gene prediction may be performed using GeMoMa and the resulting GFF files may be merged using the Merge tool of GeMoRNA.

# Binary version

A binary version (Java required) of GeMoRNA is available from the [GeMoRNA webpage](https://www.jstacs.de/index.php/GeMoRNA).
There, command line parameters of the GeMoRNA tools are described in detail.

# Source code structure

The main class of GeMoRNA is `TranscriptPrediction`, which handles the individual sub-tasks and provides parallelization. Relevant data structures are defined in classes `ReadGraph`and `SplicingGraph`. Auxiliary methods for predicting CDSs using the longest ORF heuristic are implemented in `PredictCDSFromGFF`. A tool for merging RNA-seq based GeMoRNA and homology-based GeMoMa prediction is implemented in `MergeGeMoMaGeMoRNA`.
