#' @title Percent spliced-in (PSI) visualization for skipped exons (SE)
#'
#' @description
#' \code{PlotPSI.SE} visualizes percent spliced-in (PSI) for each genomic coordinate for skipped exon (SE) event across groups of single cells.
#'
#' @details
#' This function visualizes the percent spliced-in (PSI) at each genomic coordinate encompassing the alternative exon and its flanking constitutive exons for each single cell in the form of a heatmap. The PSI mean for the respective groups are also display in the form of a line graph to summarize the PSI distributions of the respective groups. Pair-wise comparison of PSI at each genomic coordinate is performed using either the parametric ANOVA or non-parameteric Kruskal-Wallis test. The p-values can be adjusted for multiple testing using the \code{p.adjust} function.
#'
#' @param object Object of class rehab generated using \code{ComputePSI}.
#' @param SampleInfo Tab-delimited file describing the naming and grouping of the single cells. First column should contain the names of the binary alignment map (BAM) files. Second column indicates the grouping for each single cell, namely Group1, Group2, etc. Third column indicates the group names. Example file provided in extdata directory of the package.
#' @param ExonInfo Tab-delimited file describing the alternative splicing events. First column contains the alternative splicing nomenclature as per BRIE (Huang et al, Genome Biology, 2019) or MISO (Katz et al, Nature Methods, 2010). Second column indicates the type of alternative splicing event, namely SE, MXE, RI, A5SS, and A3SS. Third column contains the gene name or any personal notation. Example file provided in extdata directory of the package.
#' @param statistical.test Method for comparising PSI values at each genomic coordinates between more groups of single cells.
#' @param multiple.testing Method for adjusting p-values for multiple comparisons.
#' @param Plots Folder to output PSI plots.
#' @param plot.width Width of outplot plots.
#' @param plot.height Height of outplot plots.
#' @export
#' @return For each alternative splicing event, a single plot consisting of three subplots arranged from top to bottom is returned. Bottom subplot is a line graph of PSI means at each genomic coordinate for the groups of single cells. Middle subplot is a line graph of p-values corresponding to the comparison of PSI values at each genomic coordinate between groups of single cells. Top subplot is a heatmap of PSI values at each genomic coordinate across all single cells. Location of plots as per specified in the \code{Plots} argument.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @examples
#' PSI <- readRDS(system.file("extdata/PSI", "PSI_RED_Multi_Groups_small.rds", package="VALERIE"))
#' PlotPSI.SE.MultiGroups(PSI, SampleInfo=system.file("extdata/Sample_Info",
#'   "Sample_Info_RED_Multi_Groups.txt", package="VALERIE"),
#'   ExonInfo=system.file("extdata/Exon_Info", "Exon_Info_RED_small.txt", package="VALERIE"),
#'   statistical.test="KW",  multiple.testing="bonferroni",
#'   Plots=tempdir(), plot.width=5, plot.height=8)
#' @importFrom plyr join
#' @import ggplot2
#' @import pheatmap
#' @import ggplotify
#' @import ggpubr
#' @import scales


PlotPSI.SE.MultiGroups <- function(object, SampleInfo, ExonInfo, statistical.test=c("KW", "ANOVA"), multiple.testing=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), Plots, plot.width, plot.height) {
    
    if(!inherits(object, "rehab"))
        stop("Object must be of class 'rehab'")
 
    # Read file
        # Sample info file
        sample.info <- utils::read.table(SampleInfo, sep="\t", header=FALSE, stringsAsFactors=FALSE)
        
        # Exon file
        exon.info <- utils::read.table(ExonInfo, sep="\t", header=FALSE, stringsAsFactors=FALSE)
        
        # PSI
        psi <- unique(object)
        # psi <- unique(PSI)
        
    # Subset relevant event type
    exon.info <- exon.info[which(exon.info$V2=="SE"), ]

    # Replace header with group label
        # Subset coordinates columns
        coordinates <- psi[,c(1:2)]

        # Subset file.names columns
        file.names <- psi[,-c(1:2)]

        # Subset file names overlapping with sample info file
        file.names <- file.names[,which(names(file.names) %in% sample.info$V1)]

        # Create data frame for file names
        labels <- data.frame("V1"=names(file.names), stringsAsFactors=FALSE)

        # Annotate with group label
        labels <- join(labels, sample.info, by="V1", type="left")

        # Annotate file name data frame with group label
        names(file.names) <- labels$V2

        # Merge with coordinate columns
        psi <- cbind.data.frame(coordinates, file.names)

    # Check if chr contains chr prefix
    chr.prefix <- grep("chr", psi$Chr)
    
    # Generate plot for each event
    if(nrow(exon.info) > 5 ) {
        
        pb <- utils::txtProgressBar(1, nrow(exon.info), style=3)
        
    } else {
        
        pb <- NULL
        
    }
    
    for(k in 1:nrow(exon.info)) {
        
        # Subset event
        exons <- exon.info[k, 1]
        
        ################################# RETRIEVE PSI #################################
        
        # Determine strand
        strand <- grep("\\+", exons)
        
        if(length(strand) != 0) {
        
            # Split exons
            exons.split <- strsplit(exons, split="\\:\\+@")[[1]]
            
            # Retrieve chr
            if(length(chr.prefix)==0) {
                
                chr <- strsplit(exons.split, split="\\:")[[1]][1]
                chr <- gsub("chr", "", chr)

                } else {

                chr <- strsplit(exons.split, split="\\:")[[1]][1]

            }
                
            # Retrieve 5' constitutive exon
            exon.5.constitutive <- exons.split[1]
            exon.5.constitutive <- strsplit(exon.5.constitutive, split="\\:")
            start.5.constitutive <- as.numeric(exon.5.constitutive[[1]][2])
            end.5.constitutive <- as.numeric(exon.5.constitutive[[1]][3])
            
            # Retrieve AS exon
            exon.AS <- exons.split[2]
            exon.AS <- strsplit(exon.AS, split="\\:")
            start.AS <- as.numeric(exon.AS[[1]][2])
            end.AS <- as.numeric(exon.AS[[1]][3])
            
            # Retrieve 3' constitutive exon
            exon.3.constitutive <- exons.split[3]
            exon.3.constitutive <- strsplit(exon.3.constitutive, split="\\:")
            start.3.constitutive <- as.numeric(exon.3.constitutive[[1]][2])
            end.3.constitutive <- as.numeric(exon.3.constitutive[[1]][3])

            # Retrieve coordinates
                # Subset chr
                psi.small <- psi[which(psi$Chr==chr), ]

                # Subset 5' constitutive exon
                psi.small.exon.5.constitutive <- psi.small[which(psi.small$Position >= start.5.constitutive & psi.small$Position <= end.5.constitutive), ]
            
                # Subset AS exon
                psi.small.exon.AS <- psi.small[which(psi.small$Position >= start.AS & psi.small$Position <= end.AS), ]
            
                # Subset 3' constitutive exon
                psi.small.exon.3.constitutive <- psi.small[which(psi.small$Position >= start.3.constitutive & psi.small$Position <= end.3.constitutive), ]

                # Merge
                psi.small <- rbind.data.frame(psi.small.exon.5.constitutive, psi.small.exon.AS, psi.small.exon.3.constitutive)
                
                # Order by position
                psi.small <- psi.small[order(psi.small$Position), ]

                # Convert position to factor
                psi.small$Position <- as.factor(psi.small$Position)
                psi.small$Position <- factor(psi.small$Position, levels=psi.small$Position)
        
        } else {
            
            # Split exons
            exons.split <- strsplit(exons, split="\\:\\-@")[[1]]
            
            # Retrieve chr
            if(length(chr.prefix)==0) {
                
                chr <- strsplit(exons.split, split="\\:")[[1]][1]
                chr <- gsub("chr", "", chr)

                } else {

                chr <- strsplit(exons.split, split="\\:")[[1]][1]

            }
                
            # Retrieve 5' constitutive exon
            exon.5.constitutive <- exons.split[1]
            exon.5.constitutive <- strsplit(exon.5.constitutive, split="\\:")
            start.5.constitutive <- as.numeric(exon.5.constitutive[[1]][3])
            end.5.constitutive <- as.numeric(exon.5.constitutive[[1]][2])
            
            # Retrieve AS exon
            exon.AS <- exons.split[2]
            exon.AS <- strsplit(exon.AS, split="\\:")
            start.AS <- as.numeric(exon.AS[[1]][3])
            end.AS <- as.numeric(exon.AS[[1]][2])
            
            # Retrieve 3' constitutive exon
            exon.3.constitutive <- exons.split[3]
            exon.3.constitutive <- strsplit(exon.3.constitutive, split="\\:")
            start.3.constitutive <- as.numeric(exon.3.constitutive[[1]][3])
            end.3.constitutive <- as.numeric(exon.3.constitutive[[1]][2])

            # Retrieve coordinates
                # Subset chr
                psi.small <- psi[which(psi$Chr==chr), ]

                # Subset 5' constitutive exon
                psi.small.exon.5.constitutive <- psi.small[which(psi.small$Position <= start.5.constitutive & psi.small$Position >= end.5.constitutive), ]
            
                # Subset AS exon
                psi.small.exon.AS <- psi.small[which(psi.small$Position <= start.AS & psi.small$Position >= end.AS), ]
            
                # Subset 3' constitutive exon
                psi.small.exon.3.constitutive <- psi.small[which(psi.small$Position <= start.3.constitutive & psi.small$Position >= end.3.constitutive), ]

                # Merge
                psi.small <- rbind.data.frame(psi.small.exon.5.constitutive, psi.small.exon.AS, psi.small.exon.3.constitutive)
                
                # Order by position
                psi.small <- psi.small[order(psi.small$Position, decreasing=TRUE), ]
                
                # Convert position to factor
                psi.small$Position <- as.factor(psi.small$Position)
                psi.small$Position <- factor(psi.small$Position, levels=psi.small$Position)
            
        }
        
        ################################# COMPUTE MEAN BY GROUP #####################################

        # Compute mean PSI for each base
        groups <- c("Group1", "Group2", "Group3")
        
        mean_list <- list()

        for(j in 1:length(groups)) {

            psi.small.small <- psi.small[, which(names(psi.small) %in% groups[j])]

            mean_vector <-  NULL
            CI.lower_vector <- NULL
            CI.upper_vector <- NULL

            for(i in 1:nrow(psi.small.small)) {

                # Subset row
                row <- as.numeric(psi.small.small[i,])
                
                # Remove missing values
                row <- row[!is.na(row)]

                # Compute means
                means <- mean(row, na.rm=TRUE)
                mean_vector[i] <- means

                # Compute 95% confidence interval
                errors <- stats::qt(0.975, df=length(row)-1)*stats::sd(row)/sqrt(length(row))

                CI.lower_vector[i] <- means-errors
                CI.upper_vector[i] <- means+errors

            }
        
            mean_df <- data.frame(psi.small[,c(1:2)], "Mean"=mean_vector, "CI.lower"=CI.lower_vector, "CI.upper"=CI.upper_vector, stringsAsFactors=FALSE)

            mean_list[[j]] <- mean_df

        }

        ################################# GRAPH - GROUP #####################################

        # PSI graph
            # Define data frame
            data1 <- mean_list[[1]]
            data2 <- mean_list[[2]]
            data3 <- mean_list[[3]]

            # Merge data frame
            data1$V2 <- "Group1"
            data2$V2 <- "Group2"
            data3$V2 <- "Group3"
            data <- rbind.data.frame(data1, data2, data3, stringsAsFactors=FALSE)

            # Annotate with group name
            data <- join(data, unique(sample.info[,c("V2", "V3")]), by="V2", type="left")

            # Rename V2 and V3 columns
            names(data)[which(names(data)=="V2")] <- "Group_temp"
            names(data)[which(names(data)=="V3")] <- "Group"

            # Set factor level
            level.1 <- unique(data[which(data$Group_temp=="Group1"), "Group"])
            level.2 <- unique(data[which(data$Group_temp=="Group2"), "Group"])
            level.3 <- unique(data[which(data$Group_temp=="Group3"), "Group"])
            data$Group <- factor(data$Group, levels=c(level.1, level.2, level.3))

            # Remove intermediate column
            data$Group_temp <- NULL

            # Definitions for plot
            x <- data$Position
            y <- data$Mean
            Group <- data$Group
            ytitle <- "Mean (PSI)"

            x.1 <- c(1:length(data1$Position))
            ymin.1 <- data1$CI.lower
            ymax.1 <- data1$CI.upper

            x.2 <- c(1:length(data2$Position))
            ymin.2 <- data2$CI.lower
            ymax.2 <- data2$CI.upper
            
            x.3 <- c(1:length(data3$Position))
            ymin.3 <- data3$CI.lower
            ymax.3 <- data3$CI.upper

            accuracy <- 0.01

            # Plot
            p1 <- ggplot() +
                    geom_line(data=data, aes(x=x, y=y, color=Group, group=Group)) +
                    geom_ribbon(data=data1, aes(x=x.1, ymin=ymin.1, ymax=ymax.1), alpha=0.2, fill="red") +
                    geom_ribbon(data=data2, aes(x=x.2, ymin=ymin.2, ymax=ymax.2), alpha=0.2, fill="green") +
                    geom_ribbon(data=data3, aes(x=x.3, ymin=ymin.3, ymax=ymax.3), alpha=0.2, fill="blue") +
                    labs(y=ytitle, x=NULL) +
                    scale_y_continuous(labels=scales::number_format(accuracy=accuracy), limits=c(0, 1), position="right") +
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
                        
        ################################# GRAPH - P-VALUES #####################################

        # p-value graph
            # Perform statistical test
            positions <- psi.small$Position

            pvalues_vector <- NULL

            for(i in 1:length(positions)) {

                # Subset numeric vectors
                y.1_pvalue <- psi.small[which(psi.small$Position==positions[i]), which(names(psi.small) %in% groups[1])]
                y.2_pvalue <- psi.small[which(psi.small$Position==positions[i]), which(names(psi.small) %in% groups[2])]
                y.3_pvalue <- psi.small[which(psi.small$Position==positions[i]), which(names(psi.small) %in% groups[3])]
                y_pvalue <- c(y.1_pvalue, y.2_pvalue, y.3_pvalue)
                
                # Create factor group
                x_pvalue <- c(rep("Group1", times=length(y.1_pvalue)),
                  rep("Group2", times=length(y.2_pvalue)),
                  rep("Group3", times=length(y.3_pvalue))
                  )
                
                # Check if sufficient no. of cells with psi values
                y.1_pvalue.exp <- y.1_pvalue[which(!is.na(y.1_pvalue))]
                y.2_pvalue.exp <- y.2_pvalue[which(!is.na(y.2_pvalue))]
                y.3_pvalue.exp <- y.3_pvalue[which(!is.na(y.3_pvalue))]
                
                if(length(y.1_pvalue.exp) < 5 | length(y.2_pvalue.exp) < 5 | length(y.3_pvalue.exp) < 5) {
                    
                    pvalues_vector[i] <- 0
                    
                    } else {
                    
                
                        if(statistical.test=="KW") {
                            
                            pvalues_vector[i] <- stats::kruskal.test(as.numeric(y_pvalue) ~ as.factor(x_pvalue))$p.value
                        
                        } else {
                            
                            pvalues_vector[i] <- summary(stats::aov(as.numeric(y_pvalue) ~ as.factor(x_pvalue), alternative="two.sided"))[[1]][1,5]
                            
                        }
                        
                    }
                
            }
            
            # Replace missing values with 1.00
            pvalues_vector[is.na(pvalues_vector)] <- 1.00
            
            # Adjust for multiple testing
            pvalues_vector <- stats::p.adjust(pvalues_vector, method=multiple.testing, n=length(pvalues_vector))
          
            # Save into data frame
            pvalues_df <- data.frame("Position"=positions, "pvalue"=pvalues_vector, "pvalue_transformed"=-log10(pvalues_vector), stringsAsFactors=FALSE)
            
            # Replace infinities with arbitrary value
            inf <- pvalues_df$pvalue_transformed[!is.finite(pvalues_df$pvalue_transformed)]
            
            if(length(inf) != 0) {
                
                pvalues_df$pvalue_transformed[!is.finite(pvalues_df$pvalue_transformed)] <- 5
                
            } else {
            
                pvalues_df$pvalue_transformed <- pvalues_df$pvalue_transformed
                
            }
            
            # Replace missing values with 0
            missing <- pvalues_df$pvalue_transformed[is.nan(pvalues_df$pvalue_transformed)]
            
            if(length(missing) != 0) {
                
                pvalues_df$pvalue_transformed[is.nan(pvalues_df$pvalue_transformed)] <- 0
                
            } else {
                
                pvalues_df$pvalue_transformed <- pvalues_df$pvalue_transformed
                
            }

            # Definition
            data_pvalue <- pvalues_df
            x_pvalue <- data_pvalue$Position
            y_pvalue <- data_pvalue$pvalue_transformed
            ytitle_pvalue <- "-log10(p-value)"
            threshold <- -log10(0.05)

            ymin_pvalue <- 0
            ymax_pvalue <- max(c(y_pvalue, threshold)) + 1

            ymax_pvalue.temp <- ymax_pvalue
            ymax_pvalue.temp[is.nan(ymax_pvalue.temp)] <- 0

            if(ymax_pvalue.temp > 10) {

                accuracy_pvalue <- 0.1

                } else {

                accuracy_pvalue <- 0.01

            }
                
            # Plot
            p2 <- ggplot(data=data_pvalue, aes(x=x_pvalue, y=y_pvalue, group=1)) +
                    geom_line() +
                    geom_hline(yintercept=threshold, col="red", linetype="dashed") +
                    labs(x=NULL, y=ytitle_pvalue) +
                    scale_y_continuous(labels=scales::number_format(accuracy=accuracy_pvalue),
                                       limits=c(ymin_pvalue, ymax_pvalue), position="right") +
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

        ################################# GRAPH - SINGLE CELLS #####################################

        # PSI (individual) graph
            # Label group name
                # Save psi data frame as new object
                exp <- psi.small

                # Retrieve group label
                label <- data.frame("V2"=names(psi.small)[-c(1:2)], stringsAsFactors=FALSE)

                # Annotate with group name
                label <- join(label, unique(sample.info[,c("V2", "V3")]), by="V2", type="left")

                # Annotate psi data frame
                names(exp)[-c(1:2)] <- label$V3

            # Create expression matrix
                # Retrieve mean for each group
                mean.1 <- mean(mean_list[[1]]$Mean)
                mean.2 <- mean(mean_list[[2]]$Mean)
                mean.3 <- mean(mean_list[[3]]$Mean)

                # Define groups
                groups <- unique(sample.info[order(sample.info$V2), "V3"])

                # Retrieve groups
                group.1 <- exp[ , names(exp) %in% groups[1]]
                group.2 <- exp[ , names(exp) %in% groups[2]]
                group.3 <- exp[ , names(exp) %in% groups[3]]

                # Create data frame from each group
                exp <- cbind.data.frame(group.1, group.2, group.3)
                rownames(exp) <- psi.small[, "Position"]

            # Transpose data frame
            exp <- as.data.frame(t(exp))

            # Reorder data frame (aesthetic purpose)
                # Create genotype column
                exp$Group <- rownames(exp)

                # Split data frame by groups
                exp.1 <- exp[which(exp$Group %in% names(group.1)), ]
                exp.2 <- exp[which(exp$Group %in% names(group.2)), ]
                exp.3 <- exp[which(exp$Group %in% names(group.3)), ]

                # Remove intermediate column
                exp.1$Group <- NULL
                exp.2$Group <- NULL
                exp.3$Group <- NULL
                
                # Retrieve alternative exon
                exp.1_exon_AS <- exp.1[,which(names(exp.1) %in% c(start.AS:end.AS))]
                exp.2_exon_AS <- exp.2[,which(names(exp.2) %in% c(start.AS:end.AS))]
                exp.3_exon_AS <- exp.3[,which(names(exp.3) %in% c(start.AS:end.AS))]

                # Compute row means
                exp.1$Mean_exon_AS <- rowMeans(exp.1_exon_AS)
                exp.2$Mean_exon_AS <- rowMeans(exp.2_exon_AS)
                exp.3$Mean_exon_AS <- rowMeans(exp.3_exon_AS)

                # Reorder by means
                exp.1 <- exp.1[order(exp.1$Mean_exon_AS, decreasing=TRUE), ]
                exp.2 <- exp.2[order(exp.2$Mean_exon_AS, decreasing=TRUE), ]
                exp.3 <- exp.3[order(exp.3$Mean_exon_AS, decreasing=TRUE), ]
    
                # Merge data frames
                exp <- rbind.data.frame(exp.1, exp.2, exp.3)

                # Remove intermediate column
                exp$Mean_exon_AS <- NULL

            # Define column labels
                # Create data frame
                exon.5.constitutive <- rep("5' Cons. exon", times=(abs(end.5.constitutive-start.5.constitutive)+1))
                exon.AS <- rep("Alt. exon", times=(abs(end.AS-start.AS)+1))
                exon.3.constitutive <- rep("3' Cons. exon", times=(abs(end.3.constitutive-start.3.constitutive)+1))
                annotation.col.lab <- data.frame("Base position"=c(exon.5.constitutive, exon.AS, exon.3.constitutive))

                # Change column name
                names(annotation.col.lab) <- "Base position"

                # Create row names
                row.names(annotation.col.lab) <- colnames(exp)

            # Define row labels
                # Retrieve rownames
                annotation.row.lab <- rownames(exp)

                # Remove number suffix
                annotation.row.lab <- gsub("\\.[0-9][0-9][0-9]$", "", annotation.row.lab)
                annotation.row.lab <- gsub("\\.[0-9][0-9]$", "", annotation.row.lab)
                annotation.row.lab <- gsub("\\.[0-9]$", "", annotation.row.lab)

                # Convert to data frame
                annotation.row.lab <- data.frame("Group"=annotation.row.lab, stringsAsFactors=FALSE)

                # Create rownames
                rownames(annotation.row.lab) <- rownames(exp)

                # Set factor level
                level.1 <- unique(annotation.row.lab$Group)[1]
                level.2 <- unique(annotation.row.lab$Group)[2]
                level.3 <- unique(annotation.row.lab$Group)[3]
                annotation.row.lab$Group <- factor(annotation.row.lab$Group, levels=c(level.1, level.2, level.3))

            # Define column and colors
                # Columnm color
                annotation.col.color <- c("black", "orange", "black")
                names(annotation.col.color) <- c("5' Cons. exon", "Alt. exon", "3' Cons. exon")
                
                # Row colors
                annotation.row.color <- c("indianred1","lightgreen", "lightskyblue")
                names(annotation.row.color) <- c(level.1, level.2, level.3)

                # Create color list
                annotation.row.col.color <- list("Base position"=annotation.col.color, "Group"=annotation.row.color)

            # Color range
            color <- grDevices::colorRampPalette(c("yellow", "white", "blue"))((20))

            # Convert data frame to matrix
            exp <- as.matrix(exp)

            # Plot
            p3 <- pheatmap(exp, cluster_cols=FALSE, cluster_rows=FALSE, scale="column", show_rownames=FALSE, show_colnames=FALSE, color=color, annotation_col=annotation.col.lab, annotation_row=annotation.row.lab, annotation_colors=annotation.row.col.color, silent=TRUE, border_color=NA, fontsize=10, main=exon.info[k,3])

            # Convert pheatmap class to grob class
            p3 <- as.grob(p3)

        # Arrange plots
        plot.final <- ggarrange(p3, p1, p2, ncol=1, nrow=3, widths=0.25)

        # Save plots
        ggsave(paste(Plots, "//", k, "_SE_Plots_", exon.info[k, "V3"], ".pdf", sep=""), plot.final, device="pdf", width=plot.width, height=plot.height)
        
        # Track progress
        if(nrow(exon.info) > 5 ) {
            
            utils::setTxtProgressBar(pb, k)
            
        } else {
            
            #
            
        }

    }
    
}
