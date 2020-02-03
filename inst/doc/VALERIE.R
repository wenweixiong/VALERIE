## ----wrap-hook, echo=FALSE-----------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)

def.chunk.hook  <- knitr::knit_hooks$get("chunk")
knitr::knit_hooks$set(chunk = function(x, options) {
  x <- def.chunk.hook(x, options)
  ifelse(options$size != "normalsize", paste0("\n \\", options$size,"\n\n", x, "\n\n \\normalsize"), x)
})

## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE----------------------------------------------------------
#  library(devtools)
#  install_github("wenweixiong/VALERIE")

## ----message=FALSE, warning=FALSE----------------------------------------
library(VALERIE)

## ----message=FALSE, warning=FALSE, size="tiny"---------------------------
# Exon file (General)
exon.info <- read.table(system.file("extdata/Exon_Info", "Exon_Info_Further_Examples.txt", package="VALERIE"), sep="\t", header=FALSE, stringsAsFactors=FALSE)
print(exon.info)

# Exon file (Use for this example)
exon.info <- read.table(system.file("extdata/Exon_Info", "Exon_Info.txt", package="VALERIE"), sep="\t", header=FALSE, stringsAsFactors=FALSE)
print(exon.info)

# Sample information file
sample.info <- read.table(system.file("extdata/Sample_Info", "Sample_Info.txt", package="VALERIE"), sep="\t", header=FALSE, stringsAsFactors=FALSE)
head(sample.info) ; tail(sample.info)

# BAM files
BAM <- system.file("extdata/BAM", "", package="VALERIE")
head(list.files(BAM))

## ----eval=FALSE----------------------------------------------------------
#  # Compute PSI
#  PSI <- ComputePSI(SampleInfo=system.file("extdata/Sample_Info", "Sample_Info.txt", package="VALERIE"), ExonInfo=system.file("extdata/Exon_Info", "Exon_Info.txt", package="VALERIE"), BAM=system.file("extdata/BAM", "", package="VALERIE"), MinCoverage=10)
#  
#  # Plot PSI (Output as shown in Figure 1)
#  PlotPSI(object=PSI, SampleInfo=system.file("extdata/Sample_Info", "Sample_Info.txt", package="VALERIE"), ExonInfo=system.file("extdata/Exon_Info", "Exon_Info.txt", package="VALERIE"), statistical.test="wilcox", multiple.testing="bonferroni", Plots=tempdir(), plot.width=5, plot.height=8, EventType="SE", Groups=2)

## ----message=FALSE-------------------------------------------------------
# Check plot
output <- system.file("extdata/Plots", "1_SE_Plots_Mbp.pdf", package="VALERIE")
knitr::include_graphics(output)

