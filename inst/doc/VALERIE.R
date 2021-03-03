## ----wrap-hook, echo=FALSE----------------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)

def.chunk.hook  <- knitr::knit_hooks$get("chunk")
knitr::knit_hooks$set(chunk = function(x, options) {
  x <- def.chunk.hook(x, options)
  ifelse(options$size != "normalsize", paste0("\n \\", options$size,"\n\n", x, "\n\n \\normalsize"), x)
})

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  # From Github
#  library(devtools)
#  install_github("wenweixiong/VALERIE")
#  
#  # From CRAN
#  install.packages("VALERIE")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(VALERIE)

## ----eval=FALSE, message=FALSE, warning=FALSE---------------------------------
#  # Read sample metadata
#  path_to_file <- system.file("extdata", "BAM_PhenoData.txt", package="VALERIE")
#  BamPheno <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
#  head(BamPheno)
#  
#  # Plot
#  PlotPSI(tran_id="chr18:82554580:82554750:+@chr18:82561778:82561855:+@chr18:82572825:82572926",
#          event.type="SE",
#          strand="positive",
#          Bam=system.file("extdata/BAM", package="VALERIE"),
#          BamPheno=BamPheno,
#          cell.types=c("Ctrl", "EAE"),
#          min.coverage=10,
#          cons.exon.cutoff=100,
#          method="t.test",
#          method.adj="bonferroni",
#          cell.types.colors="ggplot.default",
#          plot.title="Mbp",
#          plot.width=5,
#          plot.height=8,
#          plot.out=system.file("extdata/Plots", "Mbp.pdf", package="VALERIE")
#          )

## ----message=FALSE, echo=FALSE------------------------------------------------
path_to_file <- system.file("extdata/Plots", "Mbp.png", package="VALERIE")
knitr::include_graphics(path_to_file)

