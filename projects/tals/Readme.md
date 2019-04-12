This packages contains classes related to the prediction of transcription activator-like effectors (TALEs) from *Xanthomonas* bacteria.

The implementation of our recent approach "PrediTALE" is distributed over different sub-packages. Specifically

- sub-package *linear* contains the base classes implementing the model and objective function
- sub-package *prediction* contains the applications for predicting TALE targets in sequence data
- sub-package *rnaseq* contain the implementation of DerTALE, which allows for filtering genome-wide predictions of TALE targets using RNA-seq data
- sub-package *training* contains the main class for learning PrediTALE models from data

Further information about PrediTALE and DerTALE from a user perspective may be found at http://jstacs.de/index.php/PrediTALE.
