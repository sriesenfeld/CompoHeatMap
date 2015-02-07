# SJ Riesenfeld

R code that uses ggplot2 for creating nice heatmaps composed with
other bar plots.

Features:

(1) Uses ggplot2, which means the plots look good and can be adapted
by adding additional ggplot2 functions.

(2) Rows and columns in the heatmap can be split visually by line
segments by specifying a fractional height at which to cut the
dendrograms clustering them. Clustering can be done in the same
command with creating the plot, or separately.

(3) Able to create complex composed plots with many types of marginal
information (e.g., sample origin, size, etc) represented by bar plots
that align correctly with the rows of the heatmap.

Code:

compoHeatMap.R: Contains all functionality. (Not much error checking.)

test.CHmap.R: Series of examples/demonstrations.

See tests.outdir directory for the pdfs output by test.CHmap.R.

Known bug: An empty Rplots.pdf is created by the test code when it is
run using Rscript (i.e., outside the interpreter). This is a known
issue associated with ggplot2. If I can resolve it, I will update
here.

