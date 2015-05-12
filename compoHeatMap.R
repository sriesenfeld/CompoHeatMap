## SJ Riesenfeld
##
## Functions for implementing better heatmaps that are composed with
## row-aligned bar plots. In particular, see create.gg.hmap.w.barps().

## NOTE FOR FUTURE: PLAN TO ADD FUNCTIONALITY TO BE ABLE TO BIND
## COLORSTACKS AND HEATMAP HORIZONTALLY, OR WILL EXPLAIN WHY THIS
## DOESN'T WORK (I MAY HAVE TRIED AND FAILED IN THE PAST).

library("ggplot2")
library("reshape2") # for melt
library("gridExtra")
library("ggdendro")
library("scales") # for "squish", "muted"
## library("RColorBrewer") # Used previously for color options

mag <- function(v, p=2) {  # magnitude of vector, using p-norm (Euclidean dist, if p=2)
    dist(rbind(v, rep(0, length(v))), method="minkowski", p=p)
}

normalize <- function(v, p=2) {  # normalize vector by dividing by its p-norm
    v/mag(v, p=p)
}

## Wrapper function for colorRampPalette based on
## http://stackoverflow.com/questions/13327326/r-image-function-in-r
## It allows for the definition of the number of intermediate colors
## between the main colors.  Using this option, one can stretch out
## colors that should predominate the palette spectrum. Additional
## arguments of colorRampPalette can also be added regarding the type
## and bias of the subsequent interpolation.
### ARGS:
## steps: integer.
### n.steps.between: NULL or integer.
## ...: Unspecified arguments sent to colorRampPalette().
### RETURNS:
## a color palette function, as returned by colorRampPalette.
### Usage:
## Compare pal.1 <- colorRampPalette(c("blue", "cyan", "yellow",
### "red"), bias=1)
## with
### pal.2 <- color.palette(c("blue", "cyan", "yellow", "red"),
### n.steps.between=c(10,1,10))
color.palette <- function(steps, n.steps.between=NULL, ...){
    if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
    if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
    fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[,fill.steps] <- col2rgb(steps)
    for(i in which(n.steps.between>0)){
        col.start=RGB[,fill.steps[i]]
        col.end=RGB[,fill.steps[i+1]]
        for(j in seq(3)){
            vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]
            RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
        }
    }
    new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}

## Get good colors for use in the heatmap. A wrapper function for
## color.palette(). Allows steps to be rescaled so that middle color
## corresponds to a given value in the range.
### ARGS:
## steps: vector of colors.
### n.steps.between: integer vector of #steps between each color given
## range.val: numeric vector of length 2 giving lower and upper limits
## of range of values that will be plotted;
### mid.val: numeric value in range.val that should be represented by
### the middle index (rounding up) in the steps vector; ignored if
### range.val is NULL.
## n.steps.final: the number of colors desired in the output vector.
### ...: Unspecified arguments are sent to color.palette().
## RETURNS:
### a vector of length n.steps.final.
get.hmap.col <- function(steps=c("blue", "cyan", "yellow", "red"), n.steps.between=c(9,1,10),
                         range.val=NULL, mid.val=NULL, n.steps.final=30,...) {
    if (!is.null(range.val)) {
        if (is.null(mid.val)) {
            mid.val=(range.val[2]-range.val[1])/2.0
        }
        mid.index=ceiling(length(n.steps.between)/2)
        ## fractional steps in the low vs. high ranges
        frac.steps.low=n.steps.between[1:mid.index]/sum(n.steps.between[1:mid.index])
        frac.steps.high=n.steps.between[(mid.index+1):length(n.steps.between)]/sum(n.steps.between[(mid.index+1):length(n.steps.between)])
        ## fraction of actual values in the low vs. high ranges
        frac.low=(mid.val-range.val[1])/(range.val[2]-range.val[1])
        frac.high=(range.val[2]-mid.val)/(range.val[2]-range.val[1])
        ## Get the right resolution and scale:
        ## n.steps is the total number of steps that will be used
        n.steps=max(255, ceiling(10^abs(log10(min(frac.low*frac.steps.low)))), ceiling(10^abs(log10(min(frac.high*frac.steps.high)))))
        ## sum(frac.high*frac.steps.high)+sum(frac.low*frac.steps.low) == 1
        n.steps.low=round(frac.low*frac.steps.low*n.steps)
        n.steps.high=round(frac.high*frac.steps.high*n.steps)
        n.steps.between=c(n.steps.low, n.steps.high)
    }
    hmcol=color.palette(steps=steps, n.steps.between=n.steps.between, ...)(n.steps.final)
    return(hmcol)
}

## Heatmap functionality via ggplot2 (more adjustable than heatmap.2),
## adapted from
## http://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/.
## Segmentation lines are inspired by:
## http://stackoverflow.com/questions/13258454/marking-specific-tiles-in-geom-tile-geom-raster.
### ARGS:
## dat: numeric matrix to be heatmapped, with named rows and columns
## (no adjustments to the data are done here except capping).
### cap: numeric vector of two values for limits of plotted values.
## hmap.col: vector of colors to use for the heatmap; must set to NULL
## if a discrete scale should be selected by ggplot2 (in which case
## "discrete" should also be set to TRUE).
### discrete: Boolean T/F, whether or not heatmap data is discrete; if
### fill colors are to be selected by ggplot2, set hmap.col=NULL.
## leg.labels: NULL or character vector of labels for discrete color
## guide; ignored unless discrete==TRUE.
### leg.labels: NULL or vector of values giving order for discrete color
### guide; ignored unless discrete==TRUE.
## r.splits/c.splits: NULL or integer vector specifying where lines
## dividing rows/columns should be drawn (indexing rows/columns from
## bottom/left).
### fix.ord: Boolean, whether or not to fix the row and column order as
### provided in the matrix; otherwise, ggplot may reorder.
## leg.title: string giving name of legend (e.g., values being
## plotted).
### base.size: base font size to use in theme.
## axis.text.size: the base font size for the axis labels.
### split.l.size: width of lines used for splitting rows/columns.
## splits.col: the color for the row/col split lines.
### no.x.labels/no.y.labels: Boolean T/F, whether or not to show the
### x-axis/y-axis labels.
## leg.direction: "vertical" or "horizontal"; direction of legend
### RETURNS:
## a ggplot2 plot object of a heatmap.
gg.hmap <- function(dat, cap=NULL,
                    hmap.col=get.hmap.col(), discrete=FALSE,
                    leg.labels=NULL, leg.breaks=NULL,
                    r.splits=NULL, c.splits=NULL,
                    fix.ord=TRUE,
                    leg.title="", split.l.size=1, splits.col="black",
                    base.size=14, axis.text.size=14,
                    na.color="gray75",
                    no.x.labels=FALSE, no.y.labels=FALSE, leg.direction="horizontal") {
    barwidth=8
    barheight=NULL
    if (leg.direction=="vertical") {
        barwidth=NULL
        barheight=6
    }
    if (is.null(cap)) {
        cap=c(min(dat, na.rm=TRUE), max(dat, na.rm=TRUE))
    }
    dat.df = melt(dat)
    colnames(dat.df)=c("ID", "Gene", "Value") ## names of the melted dataframe (will be used in the legend)
    if (fix.ord) {
        dat.df$ID=factor(dat.df$ID, levels=rownames(dat))
        dat.df$Gene=factor(dat.df$Gene, levels=colnames(dat))
    }
    if (!discrete) {
        p = ggplot(dat.df, aes(x=Gene, y=ID)) +
            ## geom_tile(aes(fill = Value), colour = "lightgray") + ## getting rid of "colour" doesn't remove the lines
            geom_tile(aes(color=Value), fill=NA) + geom_tile(aes(fill=Value), color=NA) + ## to get ride of lines, need to plot both color and fill
                scale_fill_gradientn(colours=hmap.col, limits=cap, oob=squish, na.value=na.color,
                                     guide=guide_colorbar(title=leg.title, direction=leg.direction, title.position="top", barwidth=barwidth, barheight=barheight)) +
                                         scale_colour_gradientn(colours=hmap.col, limits=cap, oob=squish, na.value=na.color,
                                                                guide=FALSE)
    } else {
        if (is.null(leg.labels)) {
            leg.labels=waiver()
        }
        if (is.null(leg.breaks)) {
            leg.breaks=waiver()
        }
        if(!is.factor(dat.df$Value)) {
            dat.df$Value=factor(dat.df$Value)
        }
        p = ggplot(dat.df, aes(x=Gene, y=ID)) +
            geom_tile(aes(fill = Value)) #, colour = "lightgray")
        if (!is.null(hmap.col)) {
            p=p+scale_fill_manual(values=hmap.col, limits=cap, na.value=na.color, labels=leg.labels, breaks=leg.breaks)
        } else {
            p=p+scale_fill_discrete(na.value=na.color, labels=leg.labels, breaks=leg.breaks)
        }
        p = p + guides(fill=guide_legend(title=leg.title, direction=leg.direction, title.position="top"))
    }
    if (! is.null(r.splits)) {
        p =p + geom_segment(data=data.frame("splits"=r.splits), size=split.l.size, colour=splits.col,
            x=0.5, xend=ncol(dat)+0.5, aes(y=splits+0.5, yend=splits+0.5))
    }
    if (! is.null(c.splits)) {
        p = p + geom_segment(data=data.frame("splits"=c.splits), size=split.l.size, colour=splits.col,
            y=0.5, yend=nrow(dat)+0.5, aes(x=splits+0.5, xend=splits+0.5))
    }
    p = p + theme_gray(base_size=base.size)
    p = p + theme_dendro()
    ## NOTE: theme_dendro() sets most ggplot theme options to blank,
    ## i.e. blank theme elements for panel grid, panel background,
    ## axis title, axis text, axis line and axis ticks.
    p = p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
    if (!no.x.labels) {
        if (!no.y.labels) {
            p =  p + theme(axis.text=element_text(colour="#696969"),
                axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5, size=axis.text.size),
                axis.text.y = element_text(size=axis.text.size, hjust=1))
        } else {
            p =  p + theme(axis.text=element_text(colour="#696969"),
                axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5, size=axis.text.size))
        }
    } else if (! no.y.labels) {
        p =  p + theme(axis.text=element_text(colour="#696969"),
            axis.text.y = element_text(size=axis.text.size, hjust=1))
    }
    p=p + theme(legend.position="bottom")
    return(p)
}

## Compute a dendrogram using hierarchical clustering.
### ARGS:
## dat: a numeric matrix with m rows and n columns from which a row
## dendrogram should be constructed (by first computing a distance
## matrix).
### dist.mat: a numeric distance matrix, as would be computed by
### dist() from dat; ignored if "dat" is provided.
## p: an integer for the "minkowski" distance;
### asinh.transf: Boolean T/F, whether or not to transform dat by asinh
### before computing distances;
## norm.rows: Boolean T/F, whether or not to normalize rows before
## computing distances (done after asinh transformation, if
## specified);
### dend.ord: NULL or a numeric vector specifying an index ordering of
### the rows that should be used as much as possible, while respecting
### the dendrogram computed.
## RETURNS: a hierarchical clustering object or NULL if one cannot be
## computed due to presence of NAs in distances.
get.hc <- function(dat=NULL, dist.mat=NULL, p=2, asinh.transf=FALSE, norm.rows=FALSE, dend.ord=NULL) {
    if (!is.null(dat)) {
        if (asinh.transf) {
            dat=asinh(dat)
        }
        if (norm.rows) {
            dat=t(apply(dat, 1, normalize, p=p))
        }
        d=dist(dat, method="minkowski", p=p)
    } else {
        if (is.null(dist.mat)) {
            stop("get.hc: specify either dat or dist.mat")
        }
        d=as.dist(dist.mat)
    }
    if (TRUE %in% is.na(d)) {
        return(NULL)
    } else {
        hc=hclust(d)
        if (!is.null(dend.ord)) {
            hc=as.hclust(reorder(as.dendrogram(hc), dend.ord))
        }
        return(hc)
    }
}

## Make a ggplot of a dendrogram.
### ARGS:
## hc: NULL or object returned by hclust.
### no.labels: Boolean T/F, whether or not to show leaf labels.
## base.size: Font size for leaf labels.
### dat: NULL or a dissimilarity matrix with rows and columns labeled
### by leaves. If hc==NULL, dat must be specified.
## ...: Unspecified arguments are passed, along with dat, to get.hc().
### RETURNS:
## a list of three objects: "gg.p": a ggplot2 plot object of the
## dendrogram computed (vertical, leaves at bottom); "dend": the
## computed dendrogram object; and "ord": a permutation of the leaf
## indices (i.e., row indices of dat) to give the order in which they
## appear in the dendrogram (from left).
gg.dendro <- function(hc=NULL, no.labels=FALSE, base.size=14, dat=NULL, ...) {
    if (is.null(hc)) {
        hc=get.hc(dat=dat,...)
    }
    if (is.null(hc)) {
        stop("Cannot compute dendrogram due to NAs in distances computed")
    }
    ddata=dendro_data(hc, type="rectangle")
    p = ggplot(segment(ddata)) +
        geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
            theme_dendro()
    if (!no.labels) {
        p = p+ scale_x_discrete(labels=ddata$labels$label, expand=c(0,0.5)) +
            theme(axis.text=element_text(colour="#696969"),
                  axis.text.x = element_text(angle=-90, hjust = 0, size=base.size))
    }
    ret.l=list("gg.p"=p, "dend"=as.dendrogram(hc), "ord"=hc$ord)
    return(ret.l)
}

## Compute a numeric vector specifying the leaf splits (from left)
## created by cutting the dendrogram at a given fractional height.
### ARGS:
## dend: dendrogram.
### cut.h: cut height as fraction of height of tree (from the leaves,
### so smaller fractions tend to give more splits).
## RETURNS: a vector of positions for the line segments between leaf
## groups split by cut.  The position counts leaves from left to right
## on a vertical representation of the dendrogram.
get.splits <- function(dend, frac.cut.h) {
    cut.h=frac.cut.h*attributes(dend)$height
    leaf.groups=lapply(cut(dend, cut.h)$lower, labels)
    temp=unlist(lapply(leaf.groups, length))
    leaf.breaks=sapply(1:length(temp), function(i) {sum(temp[1:i])})
                                        # return(leaf.breaks[-length(leaf.breaks)])
    return(c(0,leaf.breaks)) ## always include empty and complete
}

## Change a given ggplot2 theme so that the y-axis elements are blank.
theme.only.x <- function(default.theme=theme_bw(), axis.text.size=14) {
    return(default.theme + theme(axis.text=element_text(colour="#696969", size=axis.text.size),
                                 axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(),
                                 panel.grid=element_blank(), panel.background=element_blank()))
}
## Change a given ggplot2 theme so that the x-axis elements are blank.
theme.only.y <- function(default.theme=theme_bw(), axis.text.size=14) {
    return(default.theme + theme(axis.text=element_text(colour="#696969", size=axis.text.size),
                                 axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),
                                 panel.grid=element_blank(), panel.background=element_blank()))

}


## Makes a ggplot heatmap after ordering the rows and columns by
## computing dendrograms.
### ARGS:
## dat: numeric matrix.
### hc.r/hc.c: NULL or hclust objects for the rows/columns; if NULL,
### one will be computed if dend.r/dend.c==TRUE.
## r.dist.mat/c.dist.mat: numeric distance matrix for rows/columns of
## heatmap, as would be computed by dist() from dat; ignored if "dat"
## is provided.
### dend.r/dend.c: Boolean T/F; if T, compute row/column dendrogram
### from hc.r/hc.c, if specified, or dat.
## dend.r.ord/dend.c.ord: NULL or numeric vector of row/column indices
## specifying order of rows/columns to use as much as possible, while
## respecting row/column dendrogram computed.
### cut.frac.r.h/cut.frac.c.h: NULL or real number in [0,1] specifying
### where the dendrogram should be cut (as a fraction of its height)
### to create row/column splits.
## r.d.labels/c.d.labels: Boolean T/F; if T, show leaf labels on
## row/column dendrograms.
### dend.base.size: font base size for dendrograms.
## p: exponent for "minkowski" distance for dendrogram
## computations.
### asinh.transf: Boolean T/F; if T, transform data by asinh only for
### computing the dendrograms.
## norm.rows/norm.cols: Boolean T/F; if T, normalize rows/columns
## only for computing row/column dendrogram.
### ...: Unspecified arguments are passed to gg.hmap().
## RETURNS:
### a list of 7 components: "gg.hm": ggplot2 plot object of a heatmap;
### "gg.dend.r"/"gg.dend.c": NULL or ggplot2 plot object of dendrogram
### computed on the rows/columns; "ord.r"=order of rows in heatmap;
### "ord.c"=order of columns in heatmap; "r.splits"/"c.splits": NULL
### or vector describing where rows/columns are split by lines.
gg.hmap.via.dendro <- function(dat, hc.r=NULL, hc.c=NULL,
                               r.dist.mat=NULL, c.dist.mat=NULL,
                               dend.r=TRUE, dend.c=TRUE,
                               dend.r.ord=NULL, dend.c.ord=NULL,
                               cut.frac.r.h=0.5, cut.frac.c.h=0.5,
                               r.d.labels=FALSE, c.d.labels=FALSE,
                               dend.base.size=14, p=2, asinh.transf=FALSE,
                               norm.rows=FALSE, norm.cols=FALSE,
                               ...) {
    ret.l=list("gg.dend.r"=NULL, "gg.dend.c"=NULL)
    if (dend.r) {
        if (is.null(hc.r) && is.null(r.dist.mat)) {
            hc.dat=dat
        } else {
            hc.dat=NULL
        }
        gg.dendro.r.l=gg.dendro(dat=hc.dat, dist.mat=r.dist.mat, hc=hc.r, no.labels=!r.d.labels, base.size=dend.base.size,
            dend.ord=dend.r.ord, p=p, asinh.transf=asinh.transf, norm.rows=norm.rows)
        dat=dat[gg.dendro.r.l$ord, ]
        ret.l$ord.r=gg.dendro.r.l$ord
        ret.l$gg.dend.r=gg.dendro.r.l$gg.p
    } else {
        ret.l$ord.r=1:nrow(dat)
    }
    if (dend.c) {
        if (is.null(hc.c) && is.null(c.dist.mat)) {
            hc.dat=t(dat)
        } else {
            hc.dat=NULL
        }
        gg.dendro.c.l=gg.dendro(dat=hc.dat, dist.mat=c.dist.mat, hc=hc.c, no.labels=!c.d.labels, base.size=dend.base.size,
            dend.ord=dend.c.ord, p=p, asinh.transf=asinh.transf, norm.rows=norm.cols)
        dat=dat[, gg.dendro.c.l$ord]
        ret.l$ord.c=gg.dendro.c.l$ord
        ret.l$gg.dend.c=gg.dendro.c.l$gg.p
    } else {
        ret.l$ord.c=1:ncol(dat)
    }
    r.splits.v=NULL
    if (!is.null(cut.frac.r.h)) {
        r.splits.v=get.splits(gg.dendro.r.l$dend, cut.frac.r.h)
    }
    c.splits.v=NULL
    if (!is.null(cut.frac.c.h)) {
        c.splits.v=get.splits(gg.dendro.c.l$dend, cut.frac.c.h)
    }
    if (!is.null(cut.frac.r.h)) {
        p.hm=gg.hmap(dat, fix.ord=TRUE, no.x.labels=FALSE, no.y.labels=FALSE,
            r.splits=r.splits.v, c.splits=c.splits.v,...)
    } else {
        p.hm=gg.hmap(dat, fix.ord=TRUE, no.x.labels=FALSE, no.y.labels=FALSE,...)
    }
    ret.l$r.splits=r.splits.v
    ret.l$c.splits=c.splits.v
    ret.l$gg.hm=p.hm
    return(ret.l)
}

## Create simple bar ggplot, typically for use with other ggplot
## objects.
### ARGS:
## dat.v: named numeric vector.
### fix.ord: Boolean T/F; if T, fix the order as in the vector v;
### otherwise ggplot may reorder the bars.
## flip: Boolean T/F; if T, make a vertical plot; otherwise,
## horizontal.
### rev: Boolean T/F; if T, reverse the x-axis.
## fill.c: NULL or vector of the same length as dat.v of values to be
## mapped to the fill attribute in the barplot. If fill.c is a factor
## or discrete, set discrete==TRUE.
### cols: NULL or a vector of colors, to be used in the mapping of
### fill.c. If (discrete==TRUE), vector should be named by the unique
### entries in fill.c. If (fill.c==TRUE and cols==NULL), ggplot2 will
### select colors.
## default.fill: color to use as fill colour for the bars (ignored
## unless fill.c==NULL).
### discrete: Boolean T/F; if, T use discrete color scale; otherwise,
### continuous; ignored if fill.c==NULL.
### leg.title: title of the legend.
### y.axis.lab: title of y axis.
## blank: Boolean T/F; if T, use theme_dendro() to make many plot
## elements blank (ignored if x.labs.only==TRUE or if
## y.labs.only==TRUE).
### x.labs.only/y.labs.only: Boolean T/F; if T, keep only the
### x-axis/y-axis tick labels; y.labs.only is ignored if
### x.labs.only==TRUE.
## base.size: base font size to use in theme.
### axis.text.size: the base font size for the axis labels.
## splits: NULL or integer vector specifying where lines dividing bars
## should be drawn (indexing rows/columns from bottom/left).
### split.l.size: width of lines used for splitting rows/columns.
## splits.col: the color for the split lines.
### leg.direction: "vertical" or "horizontal"; direction of legend
## RETURNS:
### a ggplot2 plot object.
gg.barp <- function(dat.v, fix.ord=TRUE, flip=TRUE, rev=TRUE,
                    fill.c=NULL, cols=NULL, default.fill="#888888", discrete=TRUE,
                    leg.title="", y.axis.lab="", ## e.g., "size" or "condition"
                    blank=FALSE, x.labs.only=TRUE, y.labs.only=FALSE,
                    base.size=14, axis.text.size=14,
                    splits=NULL, split.l.size=1, splits.col="black", leg.direction="horizontal") {
    barwidth=8
    barheight=NULL
    if (leg.direction=="vertical") {
        barwidth=NULL
        barheight=6
    }
    if (fix.ord) {
        id=factor(names(dat.v), levels=names(dat.v))
    } else {
        id = factor(names(dat.v))
    }
    df=data.frame("ID"=id, "Value"=dat.v)
    if (!is.null(fill.c)) {
        df=cbind(df, "Type"=fill.c)
        p = ggplot(df, aes(x=ID, y=Value)) +
            geom_bar(stat="identity", width=1, aes(fill=Type))
        if (discrete) {
            if (is.null(cols)) {
                p = p + scale_fill_discrete()
            } else {
                p = p + scale_fill_manual(values=cols)
            }
            p = p + guides(fill=guide_legend(title=leg.title, direction=leg.direction, title.position="top"))
        } else {
            if (!is.null(cols)) {
                p = p + scale_fill_gradientn(colours=cols)
            }

            p = p + guides(fill=guide_colorbar(title=leg.title, direction=leg.direction, title.position="top", barwidth=barwidth, barheight=barheight))
        }
    } else {
        p = ggplot(df, aes(x=ID, y=Value)) +
            geom_bar(stat="identity", width=1, fill=default.fill) #width=0.99
    }
    if (! is.null(splits)) {
        if (!rev) {
            ymax=max(dat.v)
        } else {
            ymax=-max(dat.v)
        }
         p = p + geom_segment(data=data.frame("splits"=splits), size=split.l.size, colour=splits.col,
             y=0, yend=ymax, aes(x=splits+0.5, xend=splits+0.5))
    }
    p = p +ylab(y.axis.lab)
    p = p + scale_x_discrete(expand=c(0,0))
    p = p +  theme_grey(base_size = base.size) + theme(axis.text=element_text(colour="#696969", size=axis.text.size))
    expon=min(2,floor(log10(max(dat.v))))
    max.break=(ceiling(max(dat.v) / 10^expon))*(10^expon)
    if (rev) {
        p = p + scale_y_continuous(trans="reverse", expand=c(0,0),
            breaks=c(0, round(max.break/2), max.break), limits=c(max.break,0))
    } else {
        p = p + scale_y_continuous(expand=c(0,0),
            breaks=c(0, round(max.break/2), max.break), limits=c(0,max.break))
    }
    if (x.labs.only) {
        p = p+theme.only.x(theme_grey(base_size=base.size), axis.text.size=axis.text.size) + geom_vline(xintercept=0.5, linetype="dashed")
    } else if (y.labs.only) {
        p = p+theme.only.y(theme_grey(base_size=base.size), axis.text.size=axis.text.size) #+ geom_hline(yintercept=0, linetype="dashed")
    } else if (blank) {
        p = p+theme_dendro()
    }
    if (flip) {
        p = p + coord_flip()
    }
    return(p)
}

## Wrapper function for gg.barp(). Use when (1) the vector of data is
## categorical or discrete, i.e., only "fill.c", rather than "dat.v",
## is meaningful, or (2) when a 1D heatmap-type plot is desired,
## rather than a bar plot, to show continuous numeric data. Creates a
## ggplot2 bar plot object where all bars are the same height (which
## we call a "colorstack"). See gg.barp() specs for details.
gg.colorstack <- function(fill.c,rev=FALSE, x.labs.only=FALSE, y.labs.only=FALSE, blank=TRUE,...) {
    dat.v=rep(1, length(fill.c))
    names(dat.v)=names(fill.c)
    p=gg.barp(dat.v, fill.c=fill.c, rev=rev, blank=blank, x.labs.only=x.labs.only, y.labs.only=y.labs.only,...)
    p=p+theme(legend.position="bottom")
    return(p)
}


## Extract legend from a ggplot2 plot object.
## https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(gg.p) {
    tmp = ggplot_gtable(ggplot_build(gg.p))
    leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend = tmp$grobs[[leg]]
    return(legend)
}

## Creates a ggplot with plot margins, panel margins, and panel
## borders set so that is easier to use as a component of a larger
## plot (e.g., via grid.arrange()).
### ARGS:
## gg.p: a ggplot2 plot object.
### mar: a numeric vector of length 4 specifying the plot margins, in
### lines, in order top, right, bottom, left.
## RETURNS:
### a ggplot2 plot object.
ggspace <- function(gg.p, mar=c(0.1, 0.1, 0.1, 0.1)) {
    p=gg.p+theme(plot.margin=unit(mar,"lines"), panel.margin = unit(0.1, "lines"),
        panel.border=element_blank()) + labs(x=NULL, y=NULL)
    return(p)
}


## Binds the grobs in the unspecified arguments into a table that makes sense.
### ARGS:
## ...: Unspecfied arguments are grobs to be put in the table.
### c.bind: Boolean T/F; if T, column bind bgrobs; else row bind.
## size: how to choose the size of the table; should really be
## size="max" but there is a bug in ggplot2, so this needs to be set
## to "last"; sometimes it is appropriate to set size="first".
### widths: a vector of widths for the panels in the output table
### (used only if c.bind==TRUE).
## heights: a vector of heights for the panels in the output table
## (used only if c.bind==FALSE).
### RETURNS:
## a grob table.
bind.grobs <- function(..., c.bind=TRUE, size="last", widths=NULL, heights=NULL) {
    big.grob=NULL
    for (x in list(...)) {
        if (is.null(big.grob)) {
            big.grob=x
            next
        }
        if (c.bind) {
            big.grob=gtable:::cbind_gtable(big.grob, x, size=size)
        } else {
            big.grob=gtable:::rbind_gtable(big.grob, x, size=size)
        }
    }
    if (c.bind && !is.null(widths)) {
        panels=grepl("panel", big.grob$layout$name) # locate the panels in the gtable layout
        r.edge.panels=big.grob$layout$r[panels]
        ## assign new (relative) widths to the panels
        big.grob$widths[r.edge.panels] <- sapply(widths, function(x) {return(list(unit(x, "null")))})
    }
    if (!c.bind && (!is.null(heights))) {
        panels=grepl("panel", big.grob$layout$name) # locate the panels in the gtable layout
        t.edge.panels=big.grob$layout$t[panels]
        ## assign new (relative) heights to the panels
        big.grob$heights[t.edge.panels] <- sapply(heights, function(x) {return(list(unit(x, "null")))})
    }
    return(big.grob)
}

## Creates (using grid.arrange()) a composition plot based on the
## input ggplot objects containing a heatmap, bar plot (which can be
## NULL but is not by default), and optional colorstacks.
### ARGS:
## gg.hm: ggplot plot object containing a heatmap, assumed to be
## created with gg.hmap().
### gg.barp: NULL or ggplot plot object containing a bar plot, assumed
### to be created with gg.barp(). The number of rows in the heatmap
### and the bar plot are assumed to be the same.
## gg.colorstack.l: NULL or a list of ggplot objects containing
## colorstacks, each assumed to created with gg.colorstack() and
## containing a number of rows equal to the number of rows in gg.hm.
### widths: numerical vector giving the relative widths for the plots,
### in this order: barplot, colorstacks, heatmap; the length of widths
### should be equal to sum(c(gg.hm, gg.barp,
### unlist(gg.colorstack.l))!=NULL).
## heights: numerical vector of length 2, giving the relative heights
## for the legends versus the plots, or of length 3, such that the
## first element is the relative height of the title.
### plot.title: title for composed plot.
## leg.top: Boolean T/F; if T, put legend on top; else on bottom.
### RETURNS:
## a ggplot object with composed plot
compose.gg.hmap.w.barps <- function(gg.hm, gg.barp, gg.colorstack.l=NULL,
                                    widths=c(2,10), heights=c(1.5,10),
                                    plot.title="", leg.top=FALSE, leg.direction="horizontal") {
    leg.hm=g_legend(gg.hm)
    grob.hm=ggplotGrob(ggspace(
        gg.hm + theme(legend.position="none", axis.text.y = element_blank()),
        mar=c(1,1,1,0.4)))
    if (!is.null(gg.barp)) {
        grob.barp=ggplotGrob(ggspace(gg.barp, mar=c(1,0.1,1,1)))
    }
    if (!is.null(gg.colorstack.l)) {
        leg.colorstack.l=lapply(gg.colorstack.l, g_legend)
        grob.colorstack.l=lapply(gg.colorstack.l, function(p) {
            ggplotGrob(ggspace( p + theme(legend.position="none"), mar=c(1,0.1,1,0.4)))})
        if (!is.null(gg.barp)) {
            grobs.plots=bind.grobs(grob.barp, do.call(bind.grobs, grob.colorstack.l), grob.hm, widths=widths)
        } else {
            grobs.plots=bind.grobs(do.call(bind.grobs, grob.colorstack.l), grob.hm, widths=widths)
        }
        grobs.legs=bind.grobs(do.call(bind.grobs, leg.colorstack.l), leg.hm, size="first")
    } else if (!is.null(gg.barp)) {
        grobs.plots=bind.grobs(grob.barp, grob.hm, widths=widths)
        grobs.legs=leg.hm
    } else {
        grobs.plots=grob.hm
        grobs.legs=leg.hm
    }
    if (length(heights)==2) {
        if (plot.title!="") {
            heights=c(0.1*sum(heights), heights)
        } else {
            heights=c(0,heights)
        }
    } ## now heights has length 3, with legend heights in position 2

    ## Creating empty ggplot with annotation is an alternative to
    ## calling textGrob below -- tried it to avoid querying graphical
    ## parameters and hence creating empty Rplots.pdf file, but one is
    ## still created. May be due to arrangeGrob().
    ### p.title=ggplot(data.frame()) + annotate("text", label=plot.title, x=0.5, y=0.5) + theme_dendro()
    ### titleGrob=ggplotGrob(ggspace(p.title))
    if (leg.top) {
        p=arrangeGrob(textGrob(plot.title, just="center", gp=gpar(cex=1.25)),
            grobs.legs, grobs.plots,
            nrow=3, heights=unit(heights, "null"))
    } else {
        p=arrangeGrob(textGrob(plot.title, just="center", gp=gpar(cex=1.25)),
            grobs.plots, grobs.legs,
            nrow=3, heights = unit(c(heights[1], rev(heights[2:length(heights)])), "null"))
    }
    return(p)
}


## Creates ggplot objects for heatmap, barplot, and colorstacks from
## given data and then feeds them to composed.gg.hmap.w.barps() to get
## a composed plot.
### ARGS:
## hmap.dat: data matrix for heatmaps that will get fed to
## gg.hmap.via.dendro().
### barp.dat.v: NULL or data vector for barplot that will get fed to
### gg.barp().
## colorstack.dat.l: NULL or list of vectors for colorstack plots to
## be fed individually to gg.colorstack().
### colorstack.fill.l: NULL or a list (of same length as
### colorstack.dat.l) of vectors of colors (some entries of the list
### may instead be null), each named by the unique entries in the
### corresponding vector in colorstack.dat.l.
## discrete.v: NULL or a Boolean T/F vector describing whether each
## vector in colorstack.dat.l should be plotted as discrete (T) or
## continuous (F) data. If NULL, assumes all are discrete.
### leg.title.l: NULL or a list (of same length as colorstack.dat.l)
### of titles for the legends for the colorstacks.
## widths: numerical vector giving the relative widths for the plots,
## in this order: barplot, colorstacks, heatmap; the length of widths
## should be equal to sum(list(hmap.dat, barp.dat.v,
## unlist(colorstack.dat.l))!=NULL).
### heights: numerical vector of length 2 giving the relative heights
### for the legends versus the plots.
## plot.title: title for composed plot.
### leg.top: Boolean T/F; if T, put legend on top; else on bottom.
## y.labs.only: Boolean T/F; if T, keep the y-axis tick labels for
## the colorstack plots.
### split.l.size: width of lines used for splitting rows/columns.
## splits.col: the color for the row/col split lines.
### base.size: base font size to use in theme.
## axis.text.size: the base font size for the axis labels.
### print.plot: Boolean T/F; if T, composed plot is printed as a side
### effect.
## leg.direction: "vertical" or "horizontal"; direction of legends.
### ...: Unspecified arguments are passed, along with hmap.dat, to
### gg.hmap.via.dendro().
## RETURNS:
### list with 4 named components: "p.compo", the ggplot object of the
### composed plot, "p.hm.l", containing the list returned by
### gg.hmap.via.dendro(); "p.barp", a ggplot object returned by
### gg.barp(); and "p.colorstack.l", the list of ggplot objects
### creacted by all calls to gg.colorstack().
create.gg.hmap.w.barps <- function(hmap.dat, barp.dat.v,
                                   colorstack.dat.l=NULL, colorstack.fill.l=NULL,
                                   discrete.v=NULL, leg.title.l=NULL,
                                   widths=c(2,10), heights=c(1,10),
                                   plot.title="", leg.top=FALSE,
                                   y.labs.only=FALSE, leg.direction="horizontal",
                                   splits.col="black",split.l.size=1,
                                   print.plot=FALSE, base.size=14, axis.text.size=14,...) {
    p.hm.l=gg.hmap.via.dendro(dat=hmap.dat, base.size=base.size, axis.text.size=axis.text.size, split.l.size=split.l.size, leg.direction=leg.direction,...)
    if (!is.null(barp.dat.v)) {
        p.barp=gg.barp(barp.dat.v[p.hm.l$ord.r], x.labs.only=TRUE, splits=p.hm.l$r.splits, splits.col=splits.col, split.l.size=split.l.size,
            base.size=base.size, axis.text.size=axis.text.size)
    } else {
        p.barp=NULL
    }
    p.colorstack.l=NULL
    if (!is.null(colorstack.dat.l)) {
        if (is.null(discrete.v)) {
            discrete.v=rep(TRUE, length(colorstack.dat.l))
        }
        p.colorstack.l=list()
        for (i in 1:length(colorstack.dat.l)) {
            p.colorstack.l[[i]]=gg.colorstack(colorstack.dat.l[[i]][p.hm.l$ord.r],
                              cols=colorstack.fill.l[[i]], discrete=discrete.v[i],
                              y.labs.only=y.labs.only, leg.title=leg.title.l[[i]], leg.direction=leg.direction,
                              splits=p.hm.l$r.splits, splits.col=splits.col, split.l.size=split.l.size, base.size=base.size, axis.text.size=axis.text.size)
        }
    }
    p.compo=compose.gg.hmap.w.barps(p.hm.l$gg.hm, p.barp, gg.colorstack.l=p.colorstack.l,
        widths=widths, heights=heights, plot.title=plot.title, leg.top=leg.top, leg.direction=leg.direction)
    if(print.plot) {
        print(p.compo)
    }
    return(list("p.compo"=p.compo, "p.hm.l"=p.hm.l, "p.barp"=p.barp, "p.colorstack.l"=p.colorstack.l))
}

## Annotate a plot, in particular, outside the plot area.
### ARGS:
## gg.p: ggplot plot object.
### txt: text for annotation
## x/y: x/y coordinates for annotation; (0,0) is bottom left.
### hjust/vjust: horizontal and vertical justification of text.
## fontface: style of font (e.g., "italic", "plain", "bold").
### fontsize: size of font.
## col: color of text.
### ...: Unspecified arguments are passed to textGrob().
## RETURNS:
### the annotated plot (for printing)
plot.w.footnote <- function(gg.p, txt,
                            x=0, y=0, hjust=-0.1, vjust=-0.5,
                            fontface="italic", fontsize=11, col="#696969",...) {
    return(arrangeGrob(gg.p, sub = textGrob(txt, x=x, y=y, hjust=hjust, vjust=vjust,
                                 gp = gpar(fontface = fontface, fontsize = fontsize, col=col)),...))
}
