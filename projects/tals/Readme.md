This packages contains classes related to the prediction of transcription activator-like effectors (TALEs) from *Xanthomonas* bacteria.

The implementation of our recent approach "PrediTALE" is distributed over different sub-packages. Specifically

- sub-package *linear* contains the base classes implementing the model and objective function: *LF0Conditional* implements the specificities at position 0 of a TALE target box, while the remaining specificities are implemented in *LFSpecificity_parallel_cond9C*. *LFPosition_mixture* implements the position-dependent term. *LFModularConditional9C* implements the model wrapper, combining specificities, position-dependent term, and group-specific linear transformation. *MSDFunction* implements the objective function for learning model parameter, including the penalty term.
- sub-package *prediction* contains the applications for predicting TALE targets in sequence data: *QuickTBSPredictionTool* implements the scanning routine in method *run*, while *PrediTALECli* and *PrediTALEGalaxy* provide the renderers for different user interfaces. The XML file *preditale_quantitative_PBM.xml* contains an XML representation of the model parameters, which is loaded by *QuickTBSPredictionTool*.
- sub-package *rnaseq* contain the implementation of *DerTALE*, which allows for filtering genome-wide predictions of TALE targets using RNA-seq data, and an R script for visualizing DerTALE output.
- sub-package *training* contains the main class for learning PrediTALE models from data, where the main class for training is *PrediTALE_Training*.

Further information about PrediTALE and DerTALE from a user perspective may be found at http://jstacs.de/index.php/PrediTALE.
