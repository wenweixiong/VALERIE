#' @title Percent spliced-in (PSI) visualization for alternative splicing events
#'
#' @description
#' \code{PlotPSI} visualizes percent spliced-in (PSI) for each genomic coordinate for alternative splicing events across two groups of single cells.
#'
#' @details
#' This function visualizes the percent spliced-in (PSI) at each genomic coordinate encompassing the alternative exon and its flanking constitutive exons for each single cell in the form of a heatmap. The PSI mean for the respective groups are also display in the form of a line graph to summarize the PSI distributions of the respective groups. Pair-wise comparison of PSI at each genomic coordinate is performed using either the parametric (student t-test/ANOVA) or non-parameteric (Wilcoxon rank-sum/Kruskal-Wallis) test. The p-values can be adjusted for multiple testing using the \code{p.adjust} function.
#'
#' @param object Object of class rehab generated using \code{ComputePSI}.
#' @param SampleInfo Tab-delimited file describing the naming and grouping of the single cells. First column should contain the names of the binary alignment map (BAM) files. Second column indicates the grouping for each single cell, namely Group1, Group2, etc. Third column indicates the group names. Example file provided in extdata directory of the package.
#' @param ExonInfo Tab-delimited file describing the alternative splicing events. First column contains the alternative splicing nomenclature as per BRIE (Huang et al, Genome Biology, 2019) or MISO (Katz et al, Nature Methods, 2010). Second column indicates the type of alternative splicing event, namely SE, MXE, RI, A5SS, and A3SS. Third column contains the gene name or any personal notation. Example file provided in extdata directory of the package.
#' @param statistical.test Method for comparising PSI values at each genomic coordinates between groups of single cells. Parametric methods include student t-test and analysis of variance. Non-parametric methods include wilcoxon rank sum test and Kruskal-Wallis test.
#' @param multiple.testing Method for adjusting p-values for multiple comparisons.
#' @param Plots Folder to output PSI plots.
#' @param EventType Indicates the type of alternative splicing event to plot.
#' @param Groups Indicate the number of groups of single cells.
#' @param plot.width Width of outplot plots.
#' @param plot.height Height of outplot plots.
#' @export
#' @return For each alternative splicing event, a single plot consisting of three subplots arranged from top to bottom is returned. Bottom subplot is a line graph of PSI means at each genomic coordinate for the two groups of single cells. Middle subplot is a line graph of p-values corresponding to the comparison of PSI values at each genomic coordinate between the groups of single cells. Top subplot is a heatmap of PSI values at each genomic coordinate across all single cells. Location of plots as per specified in the \code{Plots} argument.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @examples
#' PSI <- readRDS(system.file("extdata/PSI", "PSI_RED_Two_Groups_small.rds", package="VALERIE"))
#' PlotPSI(PSI, SampleInfo=system.file("extdata/Sample_Info",
#'   "Sample_Info_RED_Two_Groups.txt", package="VALERIE"),
#'   ExonInfo=system.file("extdata/Exon_Info", "Exon_Info_RED_small.txt", package="VALERIE"),
#'   statistical.test="wilcox",  multiple.testing="fdr",
#'   Plots=tempdir(),
#'   plot.width=5, plot.height=8, EventType="SE", Groups=2)
#' @importFrom plyr join
#' @import ggplot2
#' @import pheatmap
#' @import ggplotify
#' @import ggpubr
#' @import scales

PlotPSI <- function(object, SampleInfo, ExonInfo, statistical.test=c("wilcox", "t.test", "KW", "ANOVA"), multiple.testing=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), Plots, plot.width, plot.height, EventType=c("SE", "MXE", "RI", "A5SS", "A3SS"), Groups) {
    
    object <- object
    SampleInfo <- SampleInfo
    ExonInfo <- ExonInfo
    statistical.test <- statistical.test
    multiple.testing <- multiple.testing
    Plots <- Plots
    plot.width <- plot.width
    plot.height <- plot.height
    EventType <- EventType
    Groups <- Groups
    
    if(EventType=="SE" & Groups==2) {
        
        PlotPSI.SE.TwoGroups(object=object, SampleInfo=SampleInfo, ExonInfo=ExonInfo, statistical.test=statistical.test, multiple.testing=multiple.testing, Plots=Plots, plot.width=plot.width, plot.height=plot.height)
        
    } else if(EventType=="SE" & Groups > 2) {
        
        PlotPSI.SE.MultiGroups(object=object, SampleInfo=SampleInfo, ExonInfo=ExonInfo, statistical.test=statistical.test, multiple.testing=multiple.testing, Plots=Plots, plot.width=plot.width, plot.height=plot.height)
        
    } else if(EventType=="MXE" & Groups==2) {
        
        PlotPSI.MXE.TwoGroups(object=object, SampleInfo=SampleInfo, ExonInfo=ExonInfo, statistical.test=statistical.test, multiple.testing=multiple.testing, Plots=Plots, plot.width=plot.width, plot.height=plot.height)
         
    } else if(EventType=="MXE" & Groups > 2) {
        
        PlotPSI.MXE.MultiGroups(object=object, SampleInfo=SampleInfo, ExonInfo=ExonInfo, statistical.test=statistical.test, multiple.testing=multiple.testing, Plots=Plots, plot.width=plot.width, plot.height=plot.height)
        
    } else if(EventType=="RI" & Groups==2) {
        
        PlotPSI.RI.TwoGroups(object=object, SampleInfo=SampleInfo, ExonInfo=ExonInfo, statistical.test=statistical.test, multiple.testing=multiple.testing, Plots=Plots, plot.width=plot.width, plot.height=plot.height)

    } else if(EventType=="RI" & Groups > 2) {
        
        PlotPSI.RI.MultiGroups(object=object, SampleInfo=SampleInfo, ExonInfo=ExonInfo, statistical.test=statistical.test, multiple.testing=multiple.testing, Plots=Plots, plot.width=plot.width, plot.height=plot.height)
        
    } else if(EventType=="A5SS" & Groups==2) {
        
        PlotPSI.A5SS.TwoGroups(object=object, SampleInfo=SampleInfo, ExonInfo=ExonInfo, statistical.test=statistical.test, multiple.testing=multiple.testing, Plots=Plots, plot.width=plot.width, plot.height=plot.height)
        
    } else if(EventType=="A5SS" & Groups > 2) {
        
        PlotPSI.A5SS.MultiGroups(object=object, SampleInfo=SampleInfo, ExonInfo=ExonInfo, statistical.test=statistical.test, multiple.testing=multiple.testing, Plots=Plots, plot.width=plot.width, plot.height=plot.height)
        
    } else if(EventType=="A3SS" & Groups==2) {
        
        PlotPSI.A3SS.TwoGroups(object=object, SampleInfo=SampleInfo, ExonInfo=ExonInfo, statistical.test=statistical.test, multiple.testing=multiple.testing, Plots=Plots, plot.width=plot.width, plot.height=plot.height)
        
    }  else {
        
        PlotPSI.A3SS.MultiGroups(object=object, SampleInfo=SampleInfo, ExonInfo=ExonInfo, statistical.test=statistical.test, multiple.testing=multiple.testing, Plots=Plots, plot.width=plot.width, plot.height=plot.height)
    
    }
    
}

