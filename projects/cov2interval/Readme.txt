This package allows to convert and visualize coverage files obtained from samtools depth (e.g. Aegilops_markgrafii_KP2012106_Julius_2D.cov).

First, the coverage file is parsed using small java program:
java projects.cov2interval.Cov2Interval projects/cov2interval/Aegilops_markgrafii_KP2012106_Julius_2D.cov projects/cov2interval/Julius-remapping.txt 1 > Ae_mark_Julius_2D.tsv

Second, this output can be visualized using the R script:
Rscript projects/cov2interval/plot_profile.R Ae_mark_Julius_2D.pdf profile - - Ae_mark_Julius_2D.tsv
The R script will create one page per chromosome contained in the input.
If you like to include the position of the centromer in the plots, you need to speciffy the second and third parameter accordingly:
Rscript projects/cov2interval/plot_profile.R Ae_mark_Julius_2D.pdf profile projects/cov2interval/centromer.txt Julius Ae_mark_Julius_2D_with_centromer.tsv
The R script also allows to plot several coverage files in the same plot. If you like to do so, you need to provide several coverage files.

Feel free to adapt the R script according to you needs.