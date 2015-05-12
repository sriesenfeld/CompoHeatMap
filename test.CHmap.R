#!/usr/bin/Rscript
### SJ Riesenfeld
## Tests/demos the code in compoHeatMap.R
source("compoHeatMap.R")

out.dir = "tests.outdir"
if (!file.exists(out.dir)) {
    dir.create(out.dir)
}

## Define the colors corresponding to the cell types (or other types)
ctypes=paste0("type_", 1:4)
ct.colors=c("#ffff33", "#ff9933", "#99ccff", "#33cc33")
names(ct.colors)=ctypes
## Sample colors
samples=c("date_1", "date_2")
sample.colors=c("grey", "black")
names(sample.colors)=samples
## Define data matrix for heatmap
test.m=rbind(
    c(1,1,10,15,2),
    c(1,1,20,15,2),
    c(1,1,1,1,1),
    c(2,1,1,2,1),
    c(20,1,2,2,3),
    c(20,2,1,1,1),
    c(10,10,1,1,5),
    c(1,1,20,15,10)
)
rownames(test.m)=paste0("clu_", 1:nrow(test.m))
colnames(test.m)=paste0("gene_", 1:ncol(test.m))

## Define data for row cluster size => will be used to create a vertical barplot
clu.sz=c(100, 10, 30, 50, 200, 150, 300, 75)
names(clu.sz)=rownames(test.m)
clu.sz.colors=color.palette(c("white", "purple"), 20)(20)

## Define data for for row conditions => will be used to create a colorstack
row.conditions=c(rep("type_1", 2), rep("type_2", 2), rep("type_3", 3), rep("type_4",1))
names(row.conditions)=rownames(test.m)

## Define data for row sample => will be used to create a colorstack
row.samples=c("date_1", "date_2", "date_1", rep("date_2", nrow(test.m)-3))
names(row.samples)=rownames(test.m)

## Test building each plot component and composing them all in one command.
## The supporting components of the plot can be accessed in the returned list.
p.complete.l=create.gg.hmap.w.barps(hmap.dat=test.m,
    barp.dat.v=clu.sz,
    colorstack.dat.l=list(row.conditions, row.samples),
    colorstack.fill.l=list(ct.colors, sample.colors),
    widths=c(3,1,1,15), plot.title="Testing complete plot construction in 1 command",
    leg.title.l=c("Condition", "Date"),
    print.plot=FALSE,
    ## heatmap-specific arguments follow
    dend.c.ord=ncol(test.m):1,
    cut.frac.c.h=0.9,
    hmap.col=get.hmap.col(range.val=range(test.m), mid.val=5),
    leg.title="Intensity",
    norm.rows=TRUE)
out.f=paste0(out.dir, "/out.test1.main.pdf")
pdf(out.f, w=7,h=8)
print(p.complete.l$p.compo)
trash=dev.off()
print(paste0("Printed 1-command complete plot to file ", out.f))

## Test building each plot component and composing them all in one
## command, mapping the continuous cluster size data to a colorstack
## instead of a barplot. The supporting components of the plot can be
## accessed in the returned list.
p.complete.nobarp.l=create.gg.hmap.w.barps(hmap.dat=test.m,
    barp.dat.v=NULL,
    colorstack.dat.l=list(clu.sz, row.conditions, row.samples),
    colorstack.fill.l=list(clu.sz.colors, ct.colors, sample.colors),
    discrete.v=c(FALSE, TRUE, TRUE),
    widths=c(1,1,1,15), heights=c(2.5,10),
    plot.title="Testing complete plot construction in 1 command\n(continuous colorstack instead of a bar plot)",
    leg.title.l=c("Cluster Size", "Condition", "Date"), leg.direction="vertical",
    print.plot=FALSE,
    ## heatmap-specific arguments follow
    dend.c.ord=ncol(test.m):1,
    cut.frac.c.h=0.9,
    hmap.col=get.hmap.col(range.val=range(test.m), mid.val=5),
    leg.title="Intensity",
    norm.rows=TRUE)
out.f=paste0(out.dir, "/out.test2.main.pdf")
pdf(out.f, , w=6.5,h=8.5)
print(p.complete.nobarp.l$p.compo)
trash=dev.off()
print(paste0("Printed 1-command complete plot (second example) to file ", out.f))

## NOW TEST BUILDING EACH PLOT PIECE SEPARATELY

## Build and plot dendrograms from test.m with gg.dendro()
out.f=paste0(out.dir, "/out.test3.gg.dendro.pdf")
pdf(out.f)
gg.dendro.l=gg.dendro(dat=test.m, p=5, no.labels=FALSE)
print(gg.dendro.l$gg.p+ggtitle("Row Dendrogram (p=5)"))
## test keeping reverse of ordering of rows as much as possible
gg.dendro.l=gg.dendro(dat=test.m, dend.ord=nrow(test.m):1, no.labels=FALSE)
print(gg.dendro.l$gg.p+ggtitle("Row Dendrogram, default reverse order for leaves"))
## test different ways of measuring distance (via get.hc ())
gg.dendro.l=gg.dendro(dat=test.m, norm.rows=T, asinh.transf=TRUE, no.labels=FALSE)
print(gg.dendro.l$gg.p+ggtitle("Row Dendrogram From Normalized, Asinh-Transformed Rows"))
gg.dendro.l=gg.dendro(dat=test.m, norm.rows=T, no.labels=FALSE)
print(gg.dendro.l$gg.p+ggtitle("Row Dendrogram from Normalized Rows"))
## test column dendrogram
gg.dendro.col.l=gg.dendro(dat=t(test.m), no.labels=FALSE)
print(gg.dendro.col.l$gg.p + ggtitle("Basic column dendrogram"))
trash=dev.off()
print(paste0("Printed dendrogram plots to file ", out.f))

## Test heatmaps with rows split by cutting row dendrogram at 50% of
## height and columns split by cutting column dendrogram at 90% of
## height; also test various coloring strategies.
spl.v=get.splits(gg.dendro.l$dend, 0.5)
cl.v=get.splits(gg.dendro.col.l$dend, 0.9)
out.f=paste0(out.dir,"/out.test4.hmap.splits.colors.pdf")
pdf(out.f, w=5, h=7)
p.hm=gg.hmap(dat=test.m[gg.dendro.l$ord, ], fix.ord=TRUE, no.x.labels=FALSE, no.y.labels=FALSE,
    r.splits=spl.v, c.splits=cl.v, leg.title="Intensity")
print(p.hm+ggtitle("Heatmap with row and column splits"))
p.hm=gg.hmap(dat=test.m[gg.dendro.l$ord, ], fix.ord=TRUE, no.x.labels=FALSE, no.y.labels=FALSE,
    r.splits=spl.v, c.splits=cl.v,
    hmap.col=get.hmap.col(range.val=range(test.m), mid.val=5),
    leg.title="Intensity")
print(p.hm+ggtitle("Heatmap with skewed colors"))
p.hm=gg.hmap(dat=test.m[gg.dendro.l$ord, ], fix.ord=TRUE, no.x.labels=FALSE, no.y.labels=FALSE,
    r.splits=spl.v, c.splits=cl.v, cap=c(1,10), leg.title="Intensity")
print(p.hm+ggtitle("Heatmap with values capped"))
## In one command, make a nice heatmap as above, but with default
## ordering of columns reversed (i.e., that ordering is used whenever
## the clustering does not completely determine the column order).
p.hm.l=gg.hmap.via.dendro(dat=test.m, dend.c.ord=ncol(test.m):1, cut.frac.c.h=0.9,
    norm.rows=TRUE, hmap.col=get.hmap.col(range.val=range(test.m), mid.val=5), leg.title="Intensity")
print(p.hm.l$gg.hm + ggtitle("Heatmap in 1 cmd, skewed colors,\ndefault reverse order for columns"))
## print the unlabeled dendrograms as well, just to see them
print(p.hm.l$gg.dend.r + ggtitle("Row dendrogram (no labels)"))
print(p.hm.l$gg.dend.c + ggtitle("Column dendrogram (no labels)"))
## Test doing the hierarchical clustering separately (could use this
## approach to do the clustering differently from default, but here we
## are just testing that we get the same overall output) and then
## using it in the heatmap.
### Get hclust for rows
hc.rows=get.hc(dat=test.m, norm.rows=TRUE)
## Get hclust for columns; reverse order columns by default
hc.cols=get.hc(dat=t(test.m), dend.ord=ncol(test.m):1)
## Now create plots with already computed hc object, include labels to check ordering
p.hm.l=gg.hmap.via.dendro(dat=test.m, hc.r=hc.rows, hc.c=hc.cols, cut.frac.c.h=0.9,
    r.d.labels=TRUE, c.d.labels=TRUE, leg.title="Intensity",
    hmap.col=get.hmap.col(range.val=range(test.m), mid.val=5))
print(p.hm.l$gg.hm + ggtitle("Same heatmap,\nrow/col clustering done separately"))
print(p.hm.l$gg.dend.r + ggtitle ("Row dendrogram"))
print(p.hm.l$gg.dend.c + ggtitle ("Column dendrogram"))
trash=dev.off()
print(paste0("Printed splits-demo plots to file ", out.f))


## Test cluster size bar plot or colorstack, with the splits added
p.clusz.barp=gg.barp(clu.sz[p.hm.l$ord.r], y.axis.lab="Cluster Size", x.labs.only=TRUE, splits=p.hm.l$r.splits)
p.clusz.colstk=gg.colorstack(clu.sz[p.hm.l$ord.r], cols=NULL, discrete=FALSE,
                    blank=T, leg.title="Cluster Size", splits=p.hm.l$r.splits)
p.clusz.colstk.bw=gg.colorstack(clu.sz[p.hm.l$ord.r], cols=clu.sz.colors, discrete=FALSE,
    blank=T, leg.title="Cluster Size", splits=p.hm.l$r.splits)
out.f=paste0(out.dir, "/out.test5.bars_and_continuous_colorstacks.pdf")
pdf(out.f, w=5, h=7)
print(p.clusz.barp + ggtitle("Testing cluster size bar plot"))
print(p.clusz.colstk + ggtitle("Testing cluster size colorstack (default colors)"))
print(p.clusz.colstk.bw + ggtitle("Testing cluster size colorstack (other colors)"))
trash=dev.off()
print(paste0("Printed continuous-data bar plot and colorstack to file ", out.f))

## Test plotting the discrete colorstacks for condition and batch, with the splits added
p.cond=gg.colorstack(row.conditions[p.hm.l$ord.r], cols=ct.colors,blank=T, leg.title="Condition", splits=p.hm.l$r.splits)
p.samp=gg.colorstack(row.samples[p.hm.l$ord.r], cols=sample.colors, blank=T, leg.title="Date", splits=p.hm.l$r.splits)
p.samp.def=gg.colorstack(row.samples[p.hm.l$ord.r], cols=NULL, blank=T, leg.title="Date", splits=p.hm.l$r.splits)
out.f=paste0(out.dir, "/out.test6.discrete_colorstacks.pdf")
pdf(out.f, w=5, h=7)
print(p.cond+ggtitle("Testing condition colorstack"))
print(p.samp+ggtitle("Testing batch colorstack"))
print(p.samp.def+ggtitle("Testing batch colorstack (default colors)"))
trash=dev.off()
print(paste0("Printed discrete-data colorstack plots to file ", out.f))

## Test composing premade heatmap, bar, and colorstack plots
p.compose.all=compose.gg.hmap.w.barps(p.hm.l$gg.hm, p.clusz.barp, gg.colorstack.l=list(p.cond, p.samp),
    widths=c(3,1,1,15), plot.title="Testing composing all plot components")
p.compose.colorstacks=compose.gg.hmap.w.barps(p.hm.l$gg.hm, gg.barp=NULL, gg.colorstack.l=list(p.cond, p.samp),
    widths=c(1,1,15), plot.title="Testing composing heatmap and colorstacks")
p.hm.barp=compose.gg.hmap.w.barps(p.hm.l$gg.hm, p.clusz.barp,
    plot.title="Testing composing heatmap and bar plot", leg.top=TRUE)
out.f=paste0(out.dir, "/out.test7.hmap.compos.w.bars.colorstacks.pdf")
pdf(out.f)
print(p.compose.all)
print(p.hm.barp)
print(p.compose.colorstacks)
trash=dev.off()
print(paste0("Printed composed plots to file ", out.f))

