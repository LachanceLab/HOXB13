#    Haplostrips
#
#    Copyright 2016 Davide Marnetto
#
#    This file is part of Haplostrips.
#
#    Haplostrips is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Haplostrips is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

cat(">>> USING CUSTOM plotter_general.r <<<\n")
write("USING_CUSTOM_PLOTTER", file="haplostrips_custom_marker.txt")

rasterize <<- TRUE 

paste0 = function(...){paste(...,sep='')}

cluster = function(mat) {
	matd=dist(t(mat),method='manhattan')
	if (any(is.na(matd))) {stop("some distances between haplotypes could not be computed.")}
	return(hclust(matd,method='single'))
	}

plot_one = function(mat,meta,refpop,pops,hapsorder,ref,color_palette,tree,matwidth,colderived) {
	if (is.null(ref)) {
		refmat=mat[,rownames(meta)[meta$pop==refpop]]
		#ref = do.call(pmax,refmat)
		ref = round(apply(refmat,1,function(x){mean(c(x,0.501))}))
	}
	else {
		ref = mat[,ref]
	}
	distsmat=mat
	afs=rowMeans(distsmat,na.rm=T)
	for(i in seq(length(afs))) {
		distsmat[i,is.na(distsmat[i,])]=afs[i]
	}

	dists = apply(distsmat,2,function(x){sum(abs(x-ref),na.rm=T)})
	pdists=dists+as.numeric(factor(meta[colnames(mat),'pop'], levels=pops, ordered=TRUE))/(10*length(pops))
	
	## --- BEGIN: explicit haplotype order override (by columns) ---
	## If env var HAPLOSTRIPS_HAP_ORDER points to a file with hap IDs (one per line),
	## reorder the *columns* of `mat` to match exactly that list (extras, if any,
	## are appended after in their original order). Prevent further re-sorting.

	hap_order_file <- Sys.getenv("HAPLOSTRIPS_HAP_ORDER", unset = "")
	if (nzchar(hap_order_file)) {
	  wanted <- readLines(hap_order_file, warn = FALSE)
	  # Drop blank lines / comments and trim whitespace
	  wanted <- trimws(wanted)
	  wanted <- wanted[nchar(wanted) > 0 & !grepl("^#", wanted)]

	  missing <- setdiff(wanted, colnames(mat))
	  if (length(missing)) {
	    stop("Haplostrips: IDs in hap_order not found among matrix columns: ",
	         paste(missing, collapse = ", "))
	  }

	  # Columns present but not listed; we’ll append them after the wanted order
	  extras <- setdiff(colnames(mat), wanted)

	  # Optional strict mode: fail if extras exist
	  strict <- tolower(Sys.getenv("HAPLOSTRIPS_HAP_ORDER_STRICT", "false")) %in% c("1","true","yes")
	  if (strict && length(extras)) {
	    stop("Haplostrips: matrix contains hap IDs not in order list (strict mode): ",
	         paste(extras, collapse = ", "))
	  }

	  newcols <- c(wanted, extras)
	  mat <- mat[, newcols, drop = FALSE]

	  # Disable any further sorting/clustering to preserve the explicit order
	  hapsorder <- 0
	  if (tree != "none") tree <- "none"
	}
	## --- END: explicit haplotype order override ---


	if (tree=='row' & hapsorder!=1) {
		cat("Warning: no clustering and sorting (1) but tree requested: omitting tree\n")
		tree='none'
	}
	if (hapsorder==1) {
		#clustering and then sorting for distance to ref
		mat = mat[,order(pdists)]
		dend=as.dendrogram(cluster(mat))
		sorteddend=reorder(dend,pdists[order(pdists)],agglo.FUN=min)
		sorting=rev(sorteddend)
	}
	else if (hapsorder==2) {
		#sorting for distance to ref
		mat = mat[,order(pdists)]
		sorting=FALSE
		
	}
	else if (hapsorder==3) {
        	#pops and then clustering
		permut=vector();
		for (i in pops) {
			popmat=mat[,rownames(meta)[meta$pop==i]]
			if (length(popmat)>3) {
				poptree=cluster(popmat)
				permut=append(permut,rownames(meta)[meta$pop==i][poptree$order])
			} else {
				permut=append(permut,rownames(meta)[meta$pop==i])
			}
		}
		mat=mat[,permut]
		sorting=FALSE
        }
	else if (hapsorder==0) {
		#no sorting
		sorting=FALSE
	}
	else {
		stop(hapsorder," is not a valid value for reorder option.") 
	}
	
	thispops = meta[colnames(mat),'pop']
	fact_pops=factor(thispops, levels=pops, ordered=TRUE)
	rsidecols = color_palette[as.integer(fact_pops)]
	if (colderived) {
		colorref=ref
		colorref[which(is.na(colorref))]=0
		mat=mat*(colorref+1)
	}
	cex.snpnames= min(4,8*matwidth/nrow(mat))
	legendcols=floor(matwidth/(0.24*(max(nchar(levels(fact_pops)))+4)))
	legendrows=ceiling(length(levels(fact_pops))/legendcols)
	h=heatmap.2(t(mat), dendrogram=tree, lhei=c(0.6+legendrows,40), lwid=c(0.5,matwidth/3),
			Rowv=sorting, Colv=NULL, scale='none', trace='none', key=F,
			col=c('white','black','grey50'), breaks=c(0,0.5,1.5,2),na.color='grey',
			labRow='', RowSideColors=rsidecols, cexCol= cex.snpnames,
			margins=c(nchar(rownames(mat)[nrow(mat)])*cex.snpnames*0.5+0.8,0.8))
	legend(0,1,legend=levels(fact_pops), ncol=legendcols, fill=color_palette,bty='n',cex=1.8,xjust=0,yjust=1)
	h$carpet[h$carpet==2]=1
	h$carpet[is.na(h$carpet)]='.'
	list(dists,h)
}



plot_haplotypes = function(	haps_in,
				happdf_fname,
				poptable,
				add_poptable=NULL,
				interval=NULL,
				spops=NULL,
				private_cutoff=.05,
				hapsorder=1,
				ref=NULL,
				colors=NULL,
				tables=FALSE,
				tree='none',
				missing=NULL,
				colderived=FALSE,
				debug=FALSE) {

	# private_cutoff is a filter as follows - remove any variant that
		# has a frequency < private_cutoff in *all* populations or 
		# has a frequency > 1-private_cutoff in *all* populations
	idcols = c("CHROM","POS","REF","ALT")
	cat("reading haps input\n")
	hapsinput = read.delim(haps_in,check.names=F,as.is=T,na.strings=".")
	names(hapsinput)[names(hapsinput)=="#CHROM"] <- "CHROM"

	# Now we can restrict to a smaller region of interest
	if (debug) {cat("subsetting snps\n")}
	if (!is.null(interval)) {
		interv_coords = strsplit(interval,'-',fixed=T)
		start = as.numeric(interv_coords[[1]][2])
		stop = as.numeric(interv_coords[[1]][3])
		hapsinput = subset(hapsinput,(CHROM == interv_coords[[1]][1]) & (POS >= start) & (POS < stop))
	}

	ainner_mat = hapsinput[,setdiff(names(hapsinput),idcols)]
	row.names(ainner_mat) = hapsinput[['POS']]
	popt = read.delim(poptable,as.is=T)

	# add an extra population table
	if (debug) {cat("adding an extra population table\n")}
	if (!is.null(add_poptable)) {
		toadd=read.table(add_poptable, header=FALSE, as.is=T)
		if (dim(toadd)[2]==2) {
			toadd$super_pop=NA}
		else if (dim(toadd)[2]<2 | dim(toadd)[2]>3) {
			stop(paste0(add_poptable,": unexpected file format (2/3 columns required)"))
		}
		colnames(toadd)=c("sample","pop","super_pop")
		cat("adding population metadata:",paste(unique(toadd$pop)),'\n')
		popt = rbind(popt[!popt$sample %in% toadd$sample,],toadd)
	}
	
	#bind meta and matrix
	if (debug) {cat("binding meta and matrix\n")}
	meta = data.frame(colid=colnames(ainner_mat),stringsAsFactors=F)
	meta$sample=substr(meta$colid,1,nchar(as.character(meta$colid))-2)
	meta = merge(meta,popt,by='sample',all.x=T)

	row.names(meta) = meta$colid
	nasamples = unique(meta$sample[which(is.na(meta$pop))])
	meta$pop[which(is.na(meta$pop))]='unknown'
	meta$super_pop[which(is.na(meta$super_pop))]='unknown'
	if (length(nasamples)>0) {
		hsamples=paste(head(nasamples,n=10),collapse=", ")
		if (length(nasamples)>10) {hsamples=paste(hsamples,", ...")}
		cat("missing population metadata for",length(nasamples),"samples, setting population to unknown (",hsamples,")")
	}
	#verify references and populations
	if (debug) {cat("verifying references and populations\n")}
	if (!is.null(spops)) {
		pops=vector()
		spops=strsplit(spops,",")[[1]]
		for (i in spops) {
			newpops=unique(meta$pop[meta$super_pop==i | meta$pop==i])
			if (length(newpops)==0) {
				stop(i," population is not present in input.")
			}
			pops=append(pops,newpops)
		}
		pops=unique(pops)
	}
	else {pops=unique(meta$pop)
	}
	cat("included populations:",pops,'\n')
	refpop=pops[1]
	if (!is.null(ref)) {
		if (!(ref %in% meta$colid)) {
			stop(ref," reference haplotype is not present in input.") }
		else {
			cat("reference haplotype: ",ref,'\n') }
	} else {
			cat("reference population: ",refpop,'\n') 
	}
	
	#keep just wanted populations
	if (debug) {cat("subsetting populations\n")}
	thismat = ainner_mat[,rownames(meta)[meta$pop %in% pops]]
	if (!is.null(ref)) {
		if (!ref %in% colnames(thismat)) {
			stop(ref," reference haplotype is not in the populations you selected.")
		}
	}


	if (missing!='True') {
		#remove sites with missing data
		cat(sprintf("sites excluded due to missing/low quality data: %s out of %s total\n",sum(!complete.cases(thismat)),nrow(thismat)))
		thismat=thismat[complete.cases(thismat),]
	}

	# Now filter rare private variants - private in *all* populations
	minfreqs = rep(2,nrow(thismat))   #set a high frequency that disappear after the first loop
	maxfreqs = rep(-1,nrow(thismat))  #set a low frequency that disappear after the first loop
	
	for (pop in pops)	{   # private means within pop
		variantfreqs = rowMeans(thismat[,meta$colid[meta$pop==pop]],na.rm=T)
		minfreqs = pmin(variantfreqs,minfreqs,na.rm=T)
		maxfreqs = pmax(variantfreqs,maxfreqs,na.rm=T)
	}
	private =  (minfreqs > 1-as.numeric(private_cutoff)) | (maxfreqs < as.numeric(private_cutoff))
	cat(sprintf("sites excluded due to low private frequency: %s out of %s total\n", sum(private),nrow(thismat)))
	thismat = thismat[!private,]

	## --- BEGIN: hardcoded haplotype order override ---
	## Your order file: one hap ID per line, e.g. "18-1-2-0319_2"
	hap_order_path <- normalizePath("~/Documents/HOXB13/hoxb13/ReorderforHaplostrips/hap_order.txt",
	                                mustWork = TRUE)

	wanted <- readLines(hap_order_path, warn = FALSE)
	wanted <- trimws(wanted)
	wanted <- wanted[nchar(wanted) > 0 & !grepl("^#", wanted)]

	# Check they exist *after* any population/missing/private filtering
	missing_cols <- setdiff(wanted, colnames(thismat))
	if (length(missing_cols)) {
	  stop("Haplostrips order file contains IDs not present in the matrix after filtering: ",
	       paste(missing_cols, collapse = ", "))
	}

	# Reorder and subset to exactly the listed columns (drop everything else)
	thismat <- thismat[, wanted, drop = FALSE]

	# Keep metadata in the same order (and subset to the same columns)
	meta <- meta[wanted, , drop = FALSE]

	# Force 'no sorting' and 'no tree' downstream so the order stays intact
	hapsorder <- 0
	tree <- "none"

	cat("Applied explicit haplotype order from: ", hap_order_path, "\n")
	cat("First columns now: ", paste(head(colnames(thismat), 6), collapse = ", "), "\n")
	## --- END: hardcoded haplotype order override ---

	
	#describe the matrix that will be plotted
	cat("sites to plot:",nrow(thismat),'\n')
	thiswidth=max(10,min(24,nrow(thismat)/20+9))
	cat("haps to plot:",ncol(thismat),'\n')
	uniq = unique(as.data.frame(t(thismat)))
	cat("unique haps:",nrow(uniq),'\n')
	if (nrow(thismat)<2 | ncol(thismat)<2) {
		stop("matrix too small to be plotted.")
	}
	# define palette
	if (is.null(colors)) {
		colors=c('#060633','#E6F0FF',"#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
	}
	else {
		colors=strsplit(colors,",")[[1]]
	}
	if (length(pops)>length(colors)) {
		 stop("more populations (",length(pops),") than colors available (",length(colors),"): use -p to choose your populations of interest or -C to define a longer color palette.")
	}
	# aaaaaand now plot!
	plotname=happdf_fname
	pdf(plotname,width=thiswidth,height=18)
	par(xpd=TRUE,mar=c(0,5,0,0))
	p = plot_one(thismat,meta,refpop,pops,hapsorder=hapsorder,ref=ref,color_palette=colors,tree,thiswidth,colderived)
	dev.off()
	cat("file created: ",plotname,'\n')
	plotname=gsub(".pdf", ".dist.pdf",happdf_fname)
	pdf(plotname,width=7,height=7)
	plot(sort(p[[1]]),1:length(p[[1]]),
		ylab='Haplotypes',xlab='Number of differences',
		cex.lab=1.25,cex.axis=0.75,cex=1,pch=15)

	dev.off()
	cat("file created: ",plotname,'\n')

	if (tables) {
		tabdists_df=data.frame(	haplotype=rev(colnames(p[[2]]$carpet)), 
					population=meta[rev(colnames(p[[2]]$carpet)),'pop'],
					distance=p[[1]][rev(colnames(p[[2]]$carpet))],
					distance_over_SNPs=p[[1]][rev(colnames(p[[2]]$carpet))]/dim(p[[2]]$carpet)[1])

		write.table(tabdists_df,file=gsub(".pdf", ".distances_tab",happdf_fname),
				quote=F,row.names=F,col.names=T,sep="\t")
		write.table(t(p[[2]]$carpet)[rev(1:dim(p[[2]]$carpet)[2]),],file=gsub(".pdf", ".mat",happdf_fname),
				quote=F,row.names=T,col.names=NA,sep="\t")
		cat("file created: ",gsub(".pdf", ".distances_tab",happdf_fname),'\n',"file created: ",gsub(".pdf", ".mat",happdf_fname),'\n')
		}

}


#   The following function 'heatmap.2' is taken and modified from gplots 2.11.3,
#   which is under license GPL-2, on date 07/12/2016.
#   Below is reported the copyright notice and warranty disclaimer of the original package:

#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

heatmap.2=function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (is.null(Rowv) || is.na(Rowv)) 
            Rowv <- FALSE
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (is.null(Colv) || is.na(Colv)) 
            Colv <- FALSE
        else if (Colv == "Rowv" && !isTRUE(Rowv)) 
            Colv <- FALSE
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE, useRaster=rasterize)
	if (dendrogram=='none') {
	axis(2, at=seq(0,1,length.out=5), labels=rev(seq(0,1,length.out=5)*nr),cex.axis=2.5,las=2,col="white",col.ticks=par("fg"))
	}
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE, useRaster=rasterize)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, useRaster=rasterize, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!is.null(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}
