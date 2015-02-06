#!/usr/bin/Rscript
### SJ Riesenfeld
## Tests/demos the code in vis.hmap.plus.20150129.R"
source("/ahg/regevdata/users/riesen/Rscripts/vis.hmap.plus.20150129.R")
# source("/ahg/regevdata/users/riesen/CyTOF/Analyses/R_code_2014/utils/hmap.cytof.20150130.R")

out.dir = "test.hmaps.outdir"
dir.create(out.dir)

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
    print.plot=TRUE,
    ## heatmap-specific arguments follow
    dend.c.ord=ncol(test.m):1,
    cut.frac.c.h=0.9,
    hmap.col=get.hmap.col(range.val=range(test.m), mid.val=5),
    leg.title.l=c("Condition", "Date"),
    leg.title="Intensity",
    norm.rows=TRUE)
pdf(paste0(out.dir, "/out.test.complete.build.pdf"), w=7,h=8)
print(p.complete.l$p.compo)
dev.off()


## NOW TEST EACH PLOT PIECE SEPARATELY

## Build and plot dendrograms from test.m with gg.dendro()
pdf(paste0(out.dir, "/out.test.gg.dendro.pdf"))
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
gg.dendro.col.l=gg.dendro(dat=t(test.m), no.labels=FALSE)
print(gg.dendro.col.l$gg.p + ggtitle("Column dendrogram"))
trash=dev.off()

## Test heatmaps with rows split by cutting row dendrogram at 50% of
## height and columns split by cutting column dendrogram at 90% of
## height; also test various coloring strategies.
spl.v=get.splits(gg.dendro.l$dend, 0.5)
cl.v=get.splits(gg.dendro.col.l$dend, 0.9)
pdf(paste0(out.dir,"/out.test.hmap.splits.colors.pdf"), w=5, h=7)
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
trash=dev.off()

## In one command, make a nice heatmap as above, but with default
## ordering of columns reversed (i.e., that ordering is used whenever
## the clustering does not completely determine the column order).
p.hm.l=gg.hmap.via.dendro(dat=test.m, dend.c.ord=ncol(test.m):1, cut.frac.c.h=0.9,
    norm.rows=TRUE, hmap.col=get.hmap.col(range.val=range(test.m), mid.val=5), leg.title="Intensity")
pdf(paste0(out.dir, "/out.test.hmap.via.dendr.pdf"), w=5, h=7)
print(p.hm.l$gg.hm + ggtitle("Heatmap, skewed colors,\ndefault reverse order for columns"))
## print the unlabeled dendrograms as well, just to see them
print(p.hm.l$gg.dend.r + ggtitle("Row dendrogram (no labels)"))
print(p.hm.l$gg.dend.c + ggtitle("Column dendrogram (no labels)"))
trash=dev.off()

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
pdf(paste0(out.dir, "/out.test.hmap.via.sep.hc.pdf"), w=5, h=7)
print(p.hm.l$gg.hm + ggtitle("Same heatmap,\nrow/col clustering done separately"))
print(p.hm.l$gg.dend.r + ggtitle ("Row dendrogram"))
print(p.hm.l$gg.dend.c + ggtitle ("Column dendrogram"))
trash=dev.off()

## Test cluster size bar plot alone
p.barp=gg.barp(clu.sz[p.hm.l$ord.r], y.axis.lab="Cluster Size", x.labs.only=TRUE)
pdf(paste0(out.dir, "/out.test.bar.clu.sz.pdf"), w=5, h=7)
print(p.barp + ggtitle("Testing bar plot"))
trash=dev.off()

## Test plotting cluster size bar plot with heatmap
p.hm.barp=compose.gg.hmap.w.barps(p.hm.l$gg.hm, p.barp, # out.f=out.f, w=7, h=7,
    plot.title="Heatmap plus bar plot", leg.top=TRUE)
pdf(paste0(out.dir, "/out.test.hmap.w.bar.clu.sz.pdf"))
print(p.hm.barp)
dev.off()

## Test plotting just the colorstacks for condition and batch
p.cond=gg.colorstack(row.conditions[p.hm.l$ord.r], cols=ct.colors,blank=T, leg.title="Condition")
p.samp=gg.colorstack(row.samples[p.hm.l$ord.r], cols=sample.colors, blank=T, leg.title="Date")
pdf(paste0(out.dir, "/out.test.colorstacks.pdf"), w=5, h=7)
print(p.cond+ggtitle("Testing condition colorstack"))
print(p.samp+ggtitle("Testing batch colorstack"))
trash=dev.off()

## Test composing premade heatmap, bar, and colorstack plots
p.compose.all=compose.gg.hmap.w.barps(p.hm.l$gg.hm, p.barp, gg.colorstack.l=list(p.cond, p.samp),
    widths=c(3,1,1,15), plot.title="Testing composing all plot components")
p.compose.colorstacks=compose.gg.hmap.w.barps(p.hm.l$gg.hm, gg.barp=NULL, gg.colorstack.l=list(p.cond, p.samp),
    widths=c(1,1,15), plot.title="Testing composing heatmap and colorstacks")
pdf(paste0(out.dir, "/out.test.hmap.w.bars.colorstacks.pdf"))
print(p.compose.all)
print(p.compose.colorstacks)
dev.off()
