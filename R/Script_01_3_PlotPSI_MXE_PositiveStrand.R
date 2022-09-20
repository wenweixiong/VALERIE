#' @title Plot percent spliced-in (PSI) values for mutually exclusive exons (MXE) located on positive strand
#'
#' @description
#' \code{PlotPSI.MXE.Pos} computes percent spliced-in (PSI) at each genomic coordinate for mutually exclusive exons (MXE) located on positive (forward) strand.
#'
#' @details
#' This function computes percent spliced-in (PSI) at each genomic coordinate for mutually exclusive exons (MXE) located on positive (forward) strand. Formula for computing PSI is number of reads with non-N CIGAR operation divided by the total number of reads. Total number of reads is the sum of reads with non-N CIGAR operation and reads with N-CIGAR operation
#'
#' @param tran_id Character string. Splicing event nomenclature.
#' @param Bam Character string. Path to folder where the BAM files and their corresponding index files are located.
#' @param BamPheno object of class data.frame. Mandatory columns are \code{bam.file.name} and \code{cell.type}. \code{bam.file.name} column indicates BAM file names as per that found in the \code{Bam} folder. \code{cell.type} column indicates the cell group names.
#' @param cell.types Character string. Cell types to plot. Should be the same number of cell groups or less than the \code{cell.type} column of the \code{BamPheno} argument.
#' @param min.coverage Numeric value. Coverage (Total reads) threshold below which the PSI value of the genomic coordinate is annotate as missing value, i.e. no coverage.
#' @param cons.exon.cutoff Numeric value. Limit the number of bases to plot for the constitutive exons. This allow users to focus the plots on the alternative exon.
#' @param method Character string. Statistical test to compare the PSI values across the different cell types. \code{"wilcox"}, \code{"t.test"}, \code{"ks"}, \code{"ad"}, and \code{"dts"} available for 2-group comparison. \code{"ANOVA"} and \code{"kw"} available for 3- or more group comparison. \code{"ks"}, \code{"ad"}, \code{"dts"}, and \code{"kw"}, represent Kolmogorov–Smirnov, Anderson-Darling, DTS, and Kruskal-Wallis test, respectively.
#' @param method.adj Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param sig.pval Numeric value. Adjust p-value, below which, the p-value is considered statistically significant.
#' @param cell.types.colors Character string. Legend colors for each cell type. Should be of same length as \code{cell.types} argument. To use ggplot2 default color scheme, please specify \code{"ggplot.default"}.
#' @param plot.title Character string. Main title for plot. Examples are gene ID, gene names, splicing ID etc..
#' @param plot.width Numeric value. Width of plot.
#' @param plot.height Numeric value. Height of plot.
#' @param plot.out Character string. Path to folder to output plot.
#' @param track Logical. If set to \code{TRUE} (default), a process of reading in the BAM files, which is the rate-limiting step, will be tracked on the console.
#' @param nboots Numeric value. When \code{method} set to \code{"dts"}, the number of bootstrap iterations for computing the p-value.
#' @param show.mean.ci Logical value. If set to \code{TRUE}, the 95percent confidence interval of the per-cell group mean PSI values will not be shown. Default is \code{FALSE}.
#' @export
#' @return A plot in PDF format located in the folder specified by \code{plot.out} argument.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import IRanges
#' @import Rsamtools
#' @import ggplot2
#' @import pheatmap
#' @import ggplotify
#' @import ggpubr
#' @import scales
#' @importFrom reshape2 dcast
#' @import grDevices
#' @import kSamples
#' @import twosamples
#' @examples
#' \donttest{
#' # DO NOT RUN
#' PlotPSI.MXE.Pos(
#'  tran_id="chr7:128748573:128748804:+@chr7:128754262:128754455
#'           :+@chr7:128754529:128754722:+@chr7:128758871:128759037",
#'  Bam="/Users/BAM/",
#'  BamPheno=BamPheno,
#'  cell.types=c("iPSC", "Endoderm"),
#'  min.coverage=10,
#'  cons.exon.cutoff=100,
#'  method="ks",
#'  method.adj="bonferroni",
#'  cell.types.colors="ggplot.default",
#'  plot.title="SNRPN",
#'  plot.width=5,
#'  plot.height=8,
#'  plot.out=paste(tempdir(), "Plot.pdf", sep="")
#'  )
#'  }


PlotPSI.MXE.Pos <- function(tran_id, Bam, BamPheno, cell.types, min.coverage, cons.exon.cutoff, method, method.adj, sig.pval=0.10, cell.types.colors, plot.title, plot.width, plot.height, plot.out, track=TRUE, nboots=2000, show.mean.ci=TRUE) {
        
    #tran_id <- "chr2:271866:271939:+@chr2:272037:272150:+@chr2:272192:272305:+@chr2:275140:275201"
    #Bam <- "/Users/seanwen/Documents/VALERIE/VALERIE/Dataset/BAM/"
    #BamPheno <- read.table("/Users/seanwen/Documents/VALERIE/VALERIE/Dataset/BAM_PhenoData.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    #cell.types <- c("iPSC", "Endoderm")
    #min.coverage <- 10
    #cons.exon.cutoff <- 100
    #method <- "ks"
    #method.adj <- "fdr"
    #cell.types.colors <- "ggplot.default"
    #plot.title <- "ACP1"
    #plot.width <- 5
    #plot.height <- 8
    #plot.out <- "/Users/seanwen/Documents/VALERIE/VALERIE/Dataset/Plots/ACP1.pdf"
    #track <- TRUE
    #sig.pval <- 0.10
    
    ##########################################################################

    # Determine cell group order
    cell.types <- factor(cell.types, levels=cell.types)

    ##########################################################################
    ############################# PLOT COLORS ################################
    ##########################################################################
    
    if(cell.types.colors[1]=="ggplot.default") {

        gg_color_hue <- function(n) {
          hues = seq(15, 375, length = n + 1)
          hcl(h = hues, l = 65, c = 100)[1:n]
        }
        n = length(cell.types)

        cell.types.colors <- gg_color_hue(n)

    } else {
        
        cell.types.colors <- cell.types.colors
    
    }
    
    ##########################################################################
    ############################# TRIM EXON ##################################
    ##########################################################################
    
    # 5' constitutive exon
    . <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]][1]
    chr <- strsplit(., split=":", fixed=TRUE)[[1]][1]
    start <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][2])
    end <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][3])
    exon.length <- (end - start) + 1
    
    if(exon.length > cons.exon.cutoff) {
        
        start.new <- end - (cons.exon.cutoff - 1)
        exon.1 <- paste(chr, start.new, end, sep=":")
        
    } else {
        
        exon.1 <- paste(chr, start, end, sep=":")
        
    }
    
    # 5' alt. exon (Do nothing)
    exon.2 <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]][2]
    
    # 3' alt. exon (Do nothing)
    exon.3 <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]][3]
    
    # 3' constitutive exon
    . <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]][4]
    chr <- strsplit(., split=":", fixed=TRUE)[[1]][1]
    start <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][2])
    end <- as.numeric(strsplit(., split=":", fixed=TRUE)[[1]][3])
    exon.length <- (end - start) + 1
    
    if(exon.length > cons.exon.cutoff) {
        
        end.new <- start + (cons.exon.cutoff - 1)
        exon.4 <- paste(chr, start, end.new, sep=":")
        
    } else {
        
        exon.4 <- paste(chr, start, end, sep=":")
        
    }
    
    # Merge
    tran_id <- paste(exon.1, exon.2, exon.3, exon.4, sep=":+@")

    ##########################################################################
    ############################# COMPUTE PSI ################################
    ##########################################################################

    # Retrieve BAM files to analyse
        # Retrieve all files within directory
        files <- list.files(Bam)

        # Retrieve non-index files
        files <- grep(".bam$", files, value=TRUE)
        
        # Subset BAM files present in sample metadata
        overlap <- intersect(files, BamPheno$bam.file.name)
        BamPheno <- BamPheno[which(BamPheno$bam.file.name %in% overlap), ]
        files <- files[which(files %in% overlap)]

    # Retrieve cell types to analyse
    BamPheno <- BamPheno[which(BamPheno$cell.type %in% cell.types), ]
    files <- files[which(files %in% BamPheno$bam.file.name)]

    BamPheno$cell.type <- factor(BamPheno$cell.type, levels=cell.types)
    
    # Check if header contains chr
        # Read example file
        bamfile.GA <- readGAlignments(paste(Bam, BamPheno$bam.file.name[1], sep="/"))
        
        # Retrieve header
        header <- names(coverage(bamfile.GA))
        
        # Check if header contains chr
        header <- grepl("^chr", header[1])
                
    # Specify event coordinates
        # chr
        exons <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]][1]
        chr <- strsplit(exons, split=":", fixed=TRUE)[[1]][1]
    
        if(length(header)==TRUE) {
            
            chr <- chr
            
        } else {
            
            chr <- gsub("chr", "", chrs)
            
        }
        
        # Retrieve start, exon length
        exons <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]]
        
        starts <- NULL
        widths <- NULL
        
        for(i in 1:length(exons)) {
                        
            start <- as.numeric(strsplit(exons[i], split=":", fixed=TRUE)[[1]][2])
            end <- as.numeric(strsplit(exons[i], split=":", fixed=TRUE)[[1]][3])
            width <- (end - start) + 1
            
            starts[i] <- start
            widths[i] <- width
            
        }
                        
        # Create GRanges object
        gr <- GRanges(seqnames=chr, ranges=IRanges(start=starts, width=widths))
                
    # Read-in BAM for specificied coordinates
    print("Reading in BAM files...")
    
    bamfileGA.list <- vector(mode="list", length=length(files))
    
    if(track==TRUE) {
        
        pb <- txtProgressBar(1, length(files), style=3)
        
    }
    
    for(i in 1:length(files)) {
        
        # Read file
        bamfileGA.list[[i]] <- readGAlignments(file=paste(Bam, files[i], sep="/"), index=paste(Bam, files[i], sep="/"), param=ScanBamParam(which=gr))
                
        # Track progress
        if(track==TRUE) {
            
            setTxtProgressBar(pb, i)
            
        }
                
    }
                            
    # Specify per-base coordinates
    exons <- strsplit(tran_id, split=":+@", fixed=TRUE)[[1]]
    
    range.list <- list()
    
    for(i in 1:length(exons)) {
                    
        start <- as.numeric(strsplit(exons[i], split=":", fixed=TRUE)[[1]][2])
        end <- as.numeric(strsplit(exons[i], split=":", fixed=TRUE)[[1]][3])
        range.list[[i]] <- seq(start, end)
        
    }
    
    coord <- unlist(range.list)
                    
    # Compute PSI
    print("Computing PSI...")
    
    psi.list <- list()
    
    for(i in 1:length(files)) {

        # Retrieve GAalignment object
        bamfile.GA <- bamfileGA.list[[i]]

        # Retrieve read counts
        read.counts <- as.vector(coverage(bamfile.GA)[[chr]])[coord]

        # Retrieve read + skipped counts
        all.counts <- as.vector(coverage(granges(bamfile.GA))[[chr]])[coord]
        
        # Set threshold for coverage
        all.counts[which(all.counts < min.coverage)] <- NA

        # Compute PSI
        psi <- read.counts/all.counts

        # Save as data frame
        psi <- data.frame("bam.file.name"=files[i], "chr.coord"=paste(chr, coord, sep=":"), "chr"=chr, "coord"=coord, "psi"=psi, stringsAsFactors=FALSE)

        # Save PSI in list
        psi.list[[i]] <- psi
        
        # Remove BAM file
        remove(bamfile.GA)
            
    }
    
    df <- do.call(rbind.data.frame, psi.list)
    
    ##########################################################################
    ######################### PREPARE TO PLOT ################################
    ##########################################################################
    
    print("Plotting...")
    
    # Annotate with sample metadata
    df <- join(df, BamPheno, by="bam.file.name", type="left")

    # Annotate constitutive, alt. exons
    df.1 <- data.frame("chr.coord"=paste(chr, range.list[[1]], sep=":"),
                       "exon.type"="5' Cons. exon",
                       stringsAsFactors=FALSE
                       )
    df.2 <- data.frame("chr.coord"=paste(chr, range.list[[2]], sep=":"),
                       "exon.type"="5' Alt. exon",
                       stringsAsFactors=FALSE
                       )
    df.3 <- data.frame("chr.coord"=paste(chr, range.list[[3]], sep=":"),
                       "exon.type"="3' Alt. exon",
                       stringsAsFactors=FALSE
                       )
    df.4 <- data.frame("chr.coord"=paste(chr, range.list[[4]], sep=":"),
                      "exon.type"="3' Cons. exon",
                      stringsAsFactors=FALSE
                      )
                      
    df.merged <- rbind.data.frame(df.1, df.2, df.3, df.4)

    df <- join(df, df.merged, by="chr.coord", type="left")
    
    # Set factor levels
    df$cell.type <- factor(df$cell.type, levels=cell.types)
    df$chr.coord <- factor(df$chr.coord, levels=unique(df$chr.coord))
    df$exon.type <- factor(df$exon.type,
                           levels=c("5' Cons. exon", "5' Alt. exon", "3' Alt. exon", "3' Cons. exon")
                           )

    # Compute mean, 95% CI PSI for each base
    .list <- list()

    for(i in 1:length(cell.types)) {

        # Subset relevant cell type
        df.small <- df[which(df$cell.type==cell.types[i]), ]
        
        # Mean
        ave <- tapply(df.small$psi, df.small$chr.coord, function(x) {
            
                y <- x[!is.na(x)]
                mean(y, na.rm=TRUE)
            
               })
        
        
        # Error
        error <- tapply(df.small$psi, df.small$chr.coord, function(x) {
                
                y <- x[!is.na(x)]
                qt(0.975, df=length(y)-1)*sd(y)/sqrt(length(y))
                
                })
        
        # CI
        ci.lower <- ave - error
        ci.higher <- ave + error
        
        # Save into data frame
        .list[[i]] <- data.frame("chr.coord"=levels(df$chr.coord),
                              "mean"=ave,
                              "ci.lower"=ci.lower,
                              "ci.higher"=ci.higher,
                              "cell.type"=cell.types[i]
                              )
    }

    df.mean <- do.call(rbind.data.frame, .list)
    
    df.mean$ci.higher[df.mean$ci.higher > 1] <- 1
    df.mean$ci.lower[df.mean$ci.lower < 0] <- 0
    
    df.mean$chr.coord <- factor(df.mean$chr.coord, levels=unique(df.mean$chr.coord))

    # Compute collapse mean by alt. exon PSI
    df.small <- df[which(df$exon.type=="5' Alt. exon"), ]

    ave <- tapply(df.small$psi, df.small$bam.file.name, function(x) {
        
            y <- x[!is.na(x)]
            mean(y, na.rm=TRUE)
        
           })

    df.ave.alt <- data.frame("bam.file.name"=names(ave),
                             "psi"=ave,
                             stringsAsFactors=FALSE)
    row.names(df.ave.alt) <- NULL

    df.ave.alt <- join(df.ave.alt, BamPheno, by="bam.file.name", type="left")

    df.ave.alt <- df.ave.alt[order(df.ave.alt$cell.type, -df.ave.alt$psi), ]
    
    # Compute p-values
    p.val <- NULL

    if(method=="wilcox") {
        
        coords <- levels(df$chr.coord)
        
        for(i in 1:length(coords)) {
            
            df.small <- df[which(df$chr.coord==coords[i]), ]
            
            non.na <- tapply(df.small$psi, df.small$cell.type, function(x) {sum(!is.na(x))})
            non.na <- non.na[which(non.na < 3)]
            non.na <- length(non.na)
            
            if(non.na != 0) {
                
                p.val[i] <- NA
                
            } else {
                
            p.val[i] <- wilcox.test(psi ~ cell.type, na.omit(df.small))$p.value
            
            }
            
        }
        
        ###############################
        
    } else if(method=="t.test") {
        
        coords <- levels(df$chr.coord)
        
        for(i in 1:length(coords)) {
            
            df.small <- df[which(df$chr.coord==coords[i]), ]
            
            non.na <- tapply(df.small$psi, df.small$cell.type, function(x) {sum(!is.na(x))})
            non.na <- non.na[which(non.na < 3)]
            non.na <- length(non.na)
            
            if(non.na != 0) {
                
                p.val[i] <- NA
                
            } else {
            
            p.val[i] <- t.test(psi ~ cell.type, na.omit(df.small))$p.value
            
            }
            
        }
        
        ###############################
        
    } else if(method=="ks") {

        coords <- levels(df$chr.coord)
        
        for(i in 1:length(coords)) {
            
            df.small <- df[which(df$chr.coord==coords[i]), ]
            
            non.na <- tapply(df.small$psi, df.small$cell.type, function(x) {sum(!is.na(x))})
            non.na <- non.na[which(non.na < 3)]
            non.na <- length(non.na)
            
            if(non.na != 0) {
                
                p.val[i] <- NA
                
            } else {
                
                x <- df.small[which(df.small$cell.type==cell.types[1]), "psi"]
                y <- df.small[which(df.small$cell.type==cell.types[2]), "psi"]
                
                x <- na.omit(x)
                y <- na.omit(y)
                
                p.val[i] <- ks.test(x, y)$p.value
                
                }
            
        }
        
        ###############################
        
    } else if(method=="ad") {

        coords <- levels(df$chr.coord)
        
        for(i in 1:length(coords)) {
            
            df.small <- df[which(df$chr.coord==coords[i]), ]
            
            non.na <- tapply(df.small$psi, df.small$cell.type, function(x) {sum(!is.na(x))})
            non.na <- non.na[which(non.na < 3)]
            non.na <- length(non.na)
            
            if(non.na != 0) {
                
                p.val[i] <- NA
                
            } else {
                
                x <- df.small[which(df.small$cell.type==cell.types[1]), "psi"]
                y <- df.small[which(df.small$cell.type==cell.types[2]), "psi"]
            
                x <- na.omit(x)
                y <- na.omit(y)
                
                error.check <- tryCatch(ad.test(x, y), error=function(err) "Error")
                
                if(error.check[1] == "Error") {
                
                    p.val[i] <- 1
                    
                    
                } else {
                    
                    p.val[i] <- ad.test(x, y, method="asymptotic")$ad[1,3]
                
                }
                
            }
            
        }
        
        ###############################
    
    } else if(method=="dts") {

        coords <- levels(df$chr.coord)
        
        for(i in 1:length(coords)) {
            
            df.small <- df[which(df$chr.coord==coords[i]), ]
            
            non.na <- tapply(df.small$psi, df.small$cell.type, function(x) {sum(!is.na(x))})
            non.na <- non.na[which(non.na < 3)]
            non.na <- length(non.na)
            
            if(non.na != 0) {
                
                p.val[i] <- NA
                
            } else {
                
                x <- df.small[which(df.small$cell.type==cell.types[1]), "psi"]
                y <- df.small[which(df.small$cell.type==cell.types[2]), "psi"]
                
                x <- na.omit(x)
                y <- na.omit(y)
            
                p.val[i] <- dts_test(x, y, nboots=nboots)[2]
                
                }
            
        }
        
        ###############################
        
    } else if(method=="ANOVA") {
        
        coords <- levels(df$chr.coord)
        
        for(i in 1:length(coords)) {
            
            df.small <- df[which(df$chr.coord==coords[i]), ]
            
            non.na <- tapply(df.small$psi, df.small$cell.type, function(x) {sum(!is.na(x))})
            non.na <- non.na[which(non.na < 3)]
            non.na <- length(non.na)
            
            if(non.na != 0) {
                
                p.val[i] <- NA
                
            } else {
                
            p.val[i] <- summary(aov(psi ~ cell.type, na.omit(df.small)))[[1]][["Pr(>F)"]][1]
            
            }
            
        }
        
        ###############################
    
    } else if(method=="kw") {
        
        coords <- levels(df$chr.coord)
        
        for(i in 1:length(coords)) {
            
            df.small <- df[which(df$chr.coord==coords[i]), ]
            
            non.na <- tapply(df.small$psi, df.small$cell.type, function(x) {sum(!is.na(x))})
            non.na <- non.na[which(non.na < 3)]
            non.na <- length(non.na)
            
            if(non.na != 0) {
                
                p.val[i] <- NA
                
            } else {
                
            p.val[i] <- kruskal.test(psi ~ cell.type, na.omit(df.small))$p.value
            
            }
            
        }
        
        ###############################
        
    }
    
    if(method != "dts") {
        
        p.val.adj <- p.adjust(p.val, method=method.adj, n=length(p.val))
        
    } else if(method == "dts"){
        
        p.val.adj <- p.adjust(p.val, method="none", n=length(p.val))
        
    }
    
    p.val.adj[is.na(p.val.adj)] <- NA
    set.seed(1)
    p.val.adj[p.val.adj == 0 & !is.na(p.val.adj)] <- runif(n=length(p.val.adj[p.val.adj == 0 & !is.na(p.val.adj)]),
                                       min=2.0e-16,
                                       max=2.2e-16
                                       )
    p.val.adj.transformed <- -log10(p.val.adj)

    df.pval <- data.frame("chr.coord"=levels(df$chr.coord),
                          "p.val.adj.transformed"=p.val.adj.transformed
                          )
    
    df.pval$chr.coord <- factor(df.pval$chr.coord, levels=unique(df.pval$chr.coord))

    ##########################################################################
    ######################### HEATMAP: PSI PLOT ##############################
    ##########################################################################

    # Create expression matrix
    df.exp <- dcast(df,  bam.file.name ~ chr.coord, value.var="psi")
    row.names(df.exp) <- df.exp$bam.file.name
    df.exp$bam.file.name <- NULL

    # Reorder by mean alt. exon expression, cell type
    df.exp <- df.exp[df.ave.alt$bam.file.name, ]
    table(row.names(df.exp)==df.ave.alt$bam.file.name)
    
    # Define column labels
        # Create data frame
        . <- data.frame("chr.coord"=names(df.exp))
        . <- join(., unique(df[,c("chr.coord", "exon.type")]), by="chr.coord", type="left")

        # Create row names
        row.names(.) <- .$chr.coord
        .$chr.coord <- NULL
        names(.) <- "exon.type"

        # Save as new object
        annotation.col.lab <- .

    # Define row labels
        # Create data frame
        . <- data.frame("bam.file.name"=row.names(df.exp))
        . <- join(., unique(df[,c("bam.file.name", "cell.type")]), by="bam.file.name", type="left")

        # Create row names
        row.names(.) <- .$bam.file.name
        .$bam.file.name <- NULL
        names(.) <- "cell.type"

        # Save as new object
        annotation.row.lab <- .

    # Define column and colors
        # Columnm color
        annotation.col.color <- c("black", "orange", "brown", "black")
        names(annotation.col.color) <- c("5' Cons. exon", "5' Alt. exon", "3' Alt. exon", "3' Cons. exon")
        
        # Row colors
        annotation.row.color <- cell.types.colors
        names(annotation.row.color) <- cell.types

        # Create color list
        annotation.row.col.color <- list("exon.type"=annotation.col.color, "cell.type"=annotation.row.color)

    # Color range
    color <- grDevices::colorRampPalette(c("yellow", "white", "blue"))((20))

    # Plot
    p1 <- pheatmap(as.matrix(df.exp), cluster_cols=FALSE, cluster_rows=FALSE, scale="column", show_rownames=FALSE, show_colnames=FALSE, color=color, annotation_col=annotation.col.lab, annotation_row=annotation.row.lab, annotation_colors=annotation.row.col.color, border_color=NA, fontsize=10, main=plot.title, annotation_names_row=FALSE, annotation_names_col=FALSE, silent=TRUE)

    # Convert pheatmap class to grob class
    p1 <- as.grob(p1)

    ##########################################################################
    ######################### LINE PLOT: MEAN PSI ############################
    ##########################################################################
    
    # Check maximum coloum z-score
    z.max <- max(apply(df.exp, 2, scale), na.rm=TRUE)
    z.min <- min(apply(df.exp, 2, scale), na.rm=TRUE)
    z.max.abs <- max(abs(c(z.max, z.min)))
    
    if(z.max.abs <= 10) {

        accuracy <- 0.01

        } else {

        accuracy <- 0.001

    }
    
    # Plot
    if(show.mean.ci==TRUE){
        
        p2 <- ggplot(df.mean, aes(x=chr.coord, y=mean, group = cell.type, col=cell.type, fill=cell.type)) +
                geom_line() +
                geom_ribbon(aes(ymin=ci.lower, ymax=ci.higher), alpha=0.2, colour = NA) +
                labs(y="Mean (PSI)", x=NULL) +
                scale_y_continuous(labels=scales::number_format(accuracy=accuracy), limits=c(0, 1), position="right") +
                scale_fill_manual(values=cell.types.colors) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border=element_blank(),
                      axis.line.y.right = element_line(color="black"),
                      axis.title.y.right=element_text(size=11, angle=0, vjust=0.5, margin = margin(t = 0, r = 0, b = 0, l = 20)),
                      axis.text=element_text(size=13),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      legend.position="none",
                      plot.margin = margin(0.5, 0.38, 0, 0.7, "cm")
                      )

    } else {
        
        p2 <- ggplot(df.mean, aes(x=chr.coord, y=mean, group = cell.type, col=cell.type, fill=cell.type)) +
                geom_line() +
                #geom_ribbon(aes(ymin=ci.lower, ymax=ci.higher), alpha=0.2, colour = NA) +
                labs(y="Mean (PSI)", x=NULL) +
                scale_y_continuous(labels=scales::number_format(accuracy=accuracy), limits=c(0, 1), position="right") +
                scale_fill_manual(values=cell.types.colors) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border=element_blank(),
                      axis.line.y.right = element_line(color="black"),
                      axis.title.y.right=element_text(size=11, angle=0, vjust=0.5, margin = margin(t = 0, r = 0, b = 0, l = 20)),
                      axis.text=element_text(size=13),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      legend.position="none",
                      plot.margin = margin(0.5, 0.38, 0, 0.7, "cm")
                      )
        
    }
    
    ##########################################################################
    ######################### LINE PLOT: P-VALUES ############################
    ##########################################################################
    
    if(max(df.pval$p.val.adj.transformed, na.rm=TRUE) > 10 & z.max.abs <= 10) {

        accuracy <- 0.1

    } else if(max(df.pval$p.val.adj.transformed, na.rm=TRUE) < 10 & z.max.abs <= 10) {

        accuracy <- 0.01

    } else if(max(df.pval$p.val.adj.transformed, na.rm=TRUE) > 10 & z.max.abs > 10) {
        
        accuracy <- 0.01
        
    } else if(max(df.pval$p.val.adj.transformed, na.rm=TRUE) < 10 & z.max.abs > 10) {
        
        accuracy <- 0.001
        
    }
        
    # Plot
    p3 <- ggplot(data=df.pval, aes(x=chr.coord, y=p.val.adj.transformed, group=1)) +
            geom_line() +
            geom_hline(yintercept=-log10(sig.pval), col="red", linetype="dashed") +
            labs(x=NULL, y="-log10(p-value)") +
            scale_y_continuous(labels=scales::number_format(accuracy=accuracy),
                               limits=c(0, max(df.pval$p.val.adj.transformed)), position="right") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border=element_blank(),
                  axis.line.y.right = element_line(color="black"),
                  axis.title.y.right=element_text(size=11, angle=0, vjust=0.5),
                  axis.text=element_text(size=13),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  plot.margin = margin(0.5, 0.16, 0, 0.7, "cm")
                  )


    ##########################################################################
    ############################# FINAL PLOT #################################
    ##########################################################################

    # Arrange plots
    plot.final <- ggarrange(p1, p2, p3, ncol=1, nrow=3, widths=0.25)

    # Save plot
    ggsave(plot.out, plot.final, device="pdf", width=plot.width, height=plot.height)


}

