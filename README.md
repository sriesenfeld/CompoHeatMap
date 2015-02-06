# SJ Riesenfeld

R code that uses ggplot2 for creating nice heatmaps composed with
other barplots. Rows and columns in the heatmaps can be split visually
with line segments by specifying a fractional height at which to cut
dendrograms clustering them. 

compoHeatMap.R: Contains all functionality. (Not much error checking.)

test.CHmap.R: Series of examples/demonstrations.

Known bug: For some reason, an empty Rplots.pdf is created by the test
code. This is a known issue associated with ggplot2. Once I resolve
it, I will update here.

