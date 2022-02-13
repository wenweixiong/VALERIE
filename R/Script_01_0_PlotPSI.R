#' @title Plot percent spliced-in (PSI) values for exon-level alternative splicing events
#'
#' @description
#' \code{PlotPSI} computes percent spliced-in (PSI) at each genomic coordinate for exon-level alternative splicing events.
#'
#' @details
#' This function computes percent spliced-in (PSI) at each genomic coordinate for exon-level alternative splicing events, namely skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS). Formula for computing PSI is number of reads with non-N CIGAR operation divided by the total number of reads. Total number of reads is the sum of reads with non-N CIGAR operation and reads with N-CIGAR operation
#'
#' @param tran_id Character string. Splicing event nomenclature.
#' @param event.type Character string. Specify \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"} or \code{"A3SS"}.
#' @param strand Character string. Specify \code{"positive"} or \code{"negative"} to indicate forward or negative strand, respectively.
#' @param Bam Character string. Path to folder where the BAM files and their corresponding index files are located.
#' @param BamPheno object of class data.frame. Mandatory columns are \code{bam.file.name} and \code{cell.type}. \code{bam.file.name} column indicates BAM file names as per that found in the \code{Bam} folder. \code{cell.type} column indicates the cell group names.
#' @param cell.types Character string. Cell types to plot. Should be the same number of cell groups or less than the \code{cell.type} column of the \code{BamPheno} argument.
#' @param min.coverage Numeric value. Coverage (Total reads) threshold below which the PSI value of the genomic coordinate is annotate as missing value, i.e. no coverage.
#' @param cons.exon.cutoff Numeric value. Limit the number of bases to plot for the constitutive exons. This allow users to focus the plots on the alternative exon.
#' @param method Character string. Statistical test to compare the PSI values across the different cell types. \code{"wilcox"}, \code{"t.test"}, \code{"ks"}, and \code{"ad"} available for 2-group comparison. \code{"ANOVA"} and \code{"kw"} available for 3- or more group comparison. \code{"ks"}, \code{"ad"}, and \code{"kw"}, represent Kolmogorov–Smirnov, Anderson-Darling, and Kruskal-Wallis test, respectively.
#' @param method.adj Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param sig.pval Numeric value. Adjust p-value, below which, the p-value is considered statistically significant.
#' @param cell.types.colors Character string. Legend colors for each cell type. Should be of same length as \code{cell.types} argument. To use ggplot2 default color scheme, please specify \code{"ggplot.default"}.
#' @param plot.title Character string. Main title for plot. Examples are gene ID, gene names, splicing ID etc..
#' @param plot.width Numeric value. Width of plot.
#' @param plot.height Numeric value. Height of plot.
#' @param plot.out Character string. Path to folder to output plot.
#' @param track Logical. If set to \code{TRUE} (default), a process of reading in the BAM files, which is the rate-limiting step, will be tracked on the console.
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
#' # Read sample metadata
#' path_to_file <- system.file("extdata", "BAM_PhenoData_Small.txt", package="VALERIE")
#' BamPheno <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
#' head(BamPheno)
#'
#' # Plot
#' PlotPSI(
#'  tran_id="chr18:82554580:82554750:+@chr18:82561778:82561855:+@chr18:82572825:82572926",
#'  event.type="SE",
#'  strand="positive",
#'  Bam=system.file("extdata/BAM", package="VALERIE"),
#'  BamPheno=BamPheno,
#'  cell.types=c("Ctrl", "EAE"),
#'  min.coverage=10,
#'  cons.exon.cutoff=100,
#'  method="t.test",
#'  method.adj="bonferroni",
#'  cell.types.colors="ggplot.default",
#'  plot.title="Mbp",
#'  plot.width=5,
#'  plot.height=8,
#'  plot.out=paste(tempdir(), "Plot.pdf", sep="")
#'  )

PlotPSI <- function(tran_id, event.type, strand, Bam, BamPheno, cell.types, min.coverage, cons.exon.cutoff, method, method.adj, sig.pval=0.10, cell.types.colors, plot.title, plot.width, plot.height, plot.out, track=TRUE, nboots=2000) {
        
    if(event.type=="SE" & strand=="positive") {
        
        PlotPSI.SE.Pos(tran_id=tran_id, Bam=Bam, BamPheno=BamPheno,
                       cell.types=cell.types, min.coverage=min.coverage,
                       cons.exon.cutoff=cons.exon.cutoff, method=method,
                       method.adj=method.adj, sig.pval, cell.types.colors=cell.types.colors,
                       plot.title=plot.title, plot.width=plot.width, plot.height=plot.height,
                       plot.out=plot.out,
                       track=track, nboots=nboots)
        
    } else if(event.type=="SE" & strand=="negative") {
        
        PlotPSI.SE.Neg(tran_id=tran_id, Bam=Bam, BamPheno=BamPheno,
                       cell.types=cell.types, min.coverage=min.coverage,
                       cons.exon.cutoff=cons.exon.cutoff, method=method,
                       method.adj=method.adj, sig.pval, cell.types.colors=cell.types.colors,
                       plot.title=plot.title, plot.width=plot.width, plot.height=plot.height,
                       plot.out=plot.out,
                       track=track, nboots=nboots)
                       
    } else if(event.type=="MXE" & strand=="positive") {
        
        PlotPSI.MXE.Pos(tran_id=tran_id, Bam=Bam, BamPheno=BamPheno,
                       cell.types=cell.types, min.coverage=min.coverage,
                       cons.exon.cutoff=cons.exon.cutoff, method=method,
                       method.adj=method.adj, sig.pval, cell.types.colors=cell.types.colors,
                       plot.title=plot.title, plot.width=plot.width, plot.height=plot.height,
                       plot.out=plot.out,
                       track=track, nboots=nboots)
                       
    } else if(event.type=="MXE" & strand=="negative") {
        
        PlotPSI.MXE.Neg(tran_id=tran_id, Bam=Bam, BamPheno=BamPheno,
                       cell.types=cell.types, min.coverage=min.coverage,
                       cons.exon.cutoff=cons.exon.cutoff, method=method,
                       method.adj=method.adj, sig.pval, cell.types.colors=cell.types.colors,
                       plot.title=plot.title, plot.width=plot.width, plot.height=plot.height,
                       plot.out=plot.out,
                       track=track, nboots=nboots)

    } else if(event.type=="RI" & strand=="positive") {
        
        PlotPSI.RI.Pos(tran_id=tran_id, Bam=Bam, BamPheno=BamPheno,
                       cell.types=cell.types, min.coverage=min.coverage,
                       cons.exon.cutoff=cons.exon.cutoff, method=method,
                       method.adj=method.adj, sig.pval, cell.types.colors=cell.types.colors,
                       plot.title=plot.title, plot.width=plot.width, plot.height=plot.height,
                       plot.out=plot.out,
                       track=track, nboots=nboots)

    } else if(event.type=="RI" & strand=="negative") {
        
        PlotPSI.RI.Neg(tran_id=tran_id, Bam=Bam, BamPheno=BamPheno,
                       cell.types=cell.types, min.coverage=min.coverage,
                       cons.exon.cutoff=cons.exon.cutoff, method=method,
                       method.adj=method.adj, sig.pval, cell.types.colors=cell.types.colors,
                       plot.title=plot.title, plot.width=plot.width, plot.height=plot.height,
                       plot.out=plot.out,
                       track=track, nboots=nboots)

    } else if(event.type=="A5SS" & strand=="positive") {
        
        PlotPSI.A5SS.Pos(tran_id=tran_id, Bam=Bam, BamPheno=BamPheno,
                       cell.types=cell.types, min.coverage=min.coverage,
                       cons.exon.cutoff=cons.exon.cutoff, method=method,
                       method.adj=method.adj, sig.pval, cell.types.colors=cell.types.colors,
                       plot.title=plot.title, plot.width=plot.width, plot.height=plot.height,
                       plot.out=plot.out,
                       track=track, nboots=nboots)

    } else if(event.type=="A5SS" & strand=="negative") {
        
        PlotPSI.A5SS.Neg(tran_id=tran_id, Bam=Bam, BamPheno=BamPheno,
                       cell.types=cell.types, min.coverage=min.coverage,
                       cons.exon.cutoff=cons.exon.cutoff, method=method,
                       method.adj=method.adj, sig.pval, cell.types.colors=cell.types.colors,
                       plot.title=plot.title, plot.width=plot.width, plot.height=plot.height,
                       plot.out=plot.out,
                       track=track, nboots=nboots)
                       
    } else if(event.type=="A3SS" & strand=="positive") {
        
        PlotPSI.A3SS.Pos(tran_id=tran_id, Bam=Bam, BamPheno=BamPheno,
                       cell.types=cell.types, min.coverage=min.coverage,
                       cons.exon.cutoff=cons.exon.cutoff, method=method,
                       method.adj=method.adj, sig.pval, cell.types.colors=cell.types.colors,
                       plot.title=plot.title, plot.width=plot.width, plot.height=plot.height,
                       plot.out=plot.out,
                       track=track, nboots=nboots)

    } else if(event.type=="A3SS" & strand=="negative") {
        
        PlotPSI.A3SS.Neg(tran_id=tran_id, Bam=Bam, BamPheno=BamPheno,
                       cell.types=cell.types, min.coverage=min.coverage,
                       cons.exon.cutoff=cons.exon.cutoff, method=method,
                       method.adj=method.adj, sig.pval, cell.types.colors=cell.types.colors,
                       plot.title=plot.title, plot.width=plot.width, plot.height=plot.height,
                       plot.out=plot.out,
                       track=track, nboots=nboots)

    }

}
