#' @title Percent spliced-in (PSI) computation
#'
#' @description
#' \code{ComputePSI} computes percent spliced-in (PSI) at each genomic coordinate for exon-level alternative splicing events, namely skipped exon (SE), mutually exclusive exons (MXE), retained intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS)
#'
#' @details
#' This function computes the percent spliced-in (PSI) at each genomic coordinate encompassing the alternative exon and its flanking constitutive exons. Formula for computing PSI is number of reads with non-N CIGAR operation divided by the total number of reads. Total number of reads is the sum of reads with non-N CIGAR operation and reads with N-CIGAR operation
#'
#' @param SampleInfo Tab-delimited file describing the naming and grouping of the single cells. First column should contain the names of the binary alignment map (BAM) files. Second column indicates the grouping for each single cell, namely Group1 and Group2. Third column indicates the group names. Example file provided in extdata directory of the package.
#' @param ExonInfo Tab-delimited file describing the alternative splicing events. First columns contains the alternative splicing nomenclature as per BRIE (Huang et al, Genome Biology, 2018) or MISO (Katz et al, Nature Methods, 2010). Second column indicates the type of alternative splicing event, namely SE, MXE, RI, A5SS, and A3SS. Third column contains the gene name or any personal notation. Example file provided in extdata directory of the package.
#' @param BAM Folder containing the BAM files sorted by genomic coordinates.
#' @param MinCoverage numeric. Coverage (Total reads) threshold below which the PSI value of the genomic coordinate is annotate as missing value, i.e. no coverage.
#' @export
#' @return A data frame of class rehab where rows are the genomic coordinates and columns are the sample names.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @examples
#' PSI <- ComputePSI(SampleInfo=system.file("extdata/Sample_Info",
#'   "Sample_Info_small.txt", package="VALERIE"),
#'   ExonInfo=system.file("extdata/Exon_Info", "Exon_Info.txt", package="VALERIE"),
#'   BAM=system.file("extdata/BAM", "", package="VALERIE"),
#'   MinCoverage=10)
#' PSI[1:5,1:4]
#' @import GenomicAlignments
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import IRanges
#' @import Rsamtools

ComputePSI <- function(SampleInfo, ExonInfo, BAM, MinCoverage) {
        
    ###################################### READ FILES ######################################

    # Retrieve BAM file list
        # Retrieve all files within directory
        files <- list.files(BAM)

        # Retrieve BAM files from directory
        files <- grep(".bam$", files, value=TRUE)
        
        # Subset BAM files present in sample info file
            # Read file
            sample.info <- utils::read.table(SampleInfo, sep="\t", header=FALSE, stringsAsFactors=FALSE)
            
            # Subset
            files <- intersect(files, sample.info$V1)

    # Read exon file
    exon.info <- utils::read.table(ExonInfo, sep="\t", header=FALSE, stringsAsFactors=FALSE)
        
    # Check if header contains chr
        # Read example file
        bamfile.GA <- readGAlignments(paste(BAM, files[1], sep="/"))
        
        # Retrieve header
        header <- names(coverage(bamfile.GA))
        
        # Check if header contains chr
        header <- grep("chr", header)
    
    ###################################### CREATE EXON LIST ######################################
    
    # Split into exons
        # SE
            # Subset relevant event
            event <- exon.info[which(exon.info$V2=="SE"), "V1"]
            
            if(length(event) != 0) {
            
                # Collapse vector
                event <- paste(event, collapse=":+@")
                
                # Split into exons
                event <- strsplit(event, split="\\:[\\+|\\-]\\@")[[1]]
                
                # Save as new object
                event.SE <- event
            
            } else {
                
                event.SE <- NULL
                
            }
            
        # MXE
            # Subset relevant event
            event <- exon.info[which(exon.info$V2=="MXE"), "V1"]
            
            if(length(event) != 0) {
            
                # Remove suffix
                event <- gsub("\\:[\\+|\\-]\\.[A|B]$", "", event)
                
                # Collapse vector
                event <- paste(event, collapse=":+@")
                
                # Split into exons
                event <- strsplit(event, split="\\:[\\+|\\-]\\@")[[1]]
                
                # Save as new object
                event.MXE <- event
                
            } else {
                
                event.MXE <- NULL
                
            }
            
        # RI
            # Subset relevant event
            event <- exon.info[which(exon.info$V2=="RI"), "V1"]
            
            if(length(event) != 0) {
            
                # Remove suffix
                event <- gsub("\\:[\\+|\\-]\\.[A|B]$", "", event)
                
                # Reformat tran_id
                event.reformatted <- character(length(event))
                
                for(i in 1:length(event)) {
                    
                    # Subset
                    event.small <- event[i]
                    
                    # Split into exons
                    event.small <- strsplit(event.small, split="\\:[\\+|\\-]\\@")[[1]]
                    
                    # Retrieve chr
                    chr <- strsplit(event.small, split="\\:")[[1]][1]
                    
                    # Retrieve smallest and largest coordinates
                    exon.1 <- as.numeric(strsplit(event.small, split="\\:")[[1]][-1])
                    exon.2 <- as.numeric(strsplit(event.small, split="\\:")[[2]][-1])
                    coord <- c(exon.1, exon.2)
                    coord.1 <- min(coord)
                    coord.2 <- max(coord)
                    
                    # Merge
                    event.reformatted[i] <- paste(chr, coord.1, coord.2, sep=":")
                
                }
                
                # Save as new object
                event.RI <- event.reformatted
                
            } else {
                
                event.RI <- NULL
                
            }
            
        # A5SS
            # Subset relevant event
            event <- exon.info[which(exon.info$V2=="A5SS"), "V1"]
            
            if(length(event) != 0) {
            
                # Remove suffix
                event <- gsub("\\:[\\+|\\-]\\.[A|B]$", "", event)
                
                # Reformat tran_id
                event.reformatted <- vector(mode="list", length=length(event))
                
                for(i in 1:length(event)) {
                    
                    # Subset
                    event.small <- event[i]
                    
                    # Split into exons
                    event.small <- strsplit(event.small, split="\\:[\\+|\\-]\\@")[[1]]
                    
                    # Retrieve chr
                    chr <- strsplit(event.small, split="\\:")[[1]][1]
                    
                    # Retrieve smallest and largest coordinates
                        # exon 1
                        coord <- as.numeric(strsplit(event.small, split="[\\:|\\|]")[[1]][-1])
                        coord.1 <- min(coord)
                        coord.2 <- max(coord)
                        exon.1 <- paste(chr, coord.1, coord.2, sep=":")
                        
                        # exon 2
                        coord <- as.numeric(strsplit(event.small, split="[\\:|\\|]")[[2]][-1])
                        coord.1 <- min(coord)
                        coord.2 <- max(coord)
                        exon.2 <- paste(chr, coord.1, coord.2, sep=":")
                        
                    # Merge
                    event.reformatted[[i]] <- c(exon.1, exon.2)
                
                }
            
                # Save as new object
                event.A5SS <- unlist(event.reformatted)
            
            } else {
                
                event.A5SS <- NULL
                
            }
            
        # A3SS
            # Subset relevant event
            event <- exon.info[which(exon.info$V2=="A3SS"), "V1"]
            
            if(length(event) != 0) {
            
                # Remove suffix
                event <- gsub("\\:[\\+|\\-]\\.[A|B]$", "", event)
                
                # Reformat tran_id
                event.reformatted <- vector(mode="list", length=length(event))
                
                for(i in 1:length(event)) {
                    
                    # Subset
                    event.small <- event[i]
                    
                    # Split into exons
                    event.small <- strsplit(event.small, split="\\:[\\+|\\-]\\@")[[1]]
                    
                    # Retrieve chr
                    chr <- strsplit(event.small, split="\\:")[[1]][1]
                    
                    # Retrieve smallest and largest coordinates
                        # exon 1
                        coord <- as.numeric(strsplit(event.small, split="[\\:|\\|]")[[1]][-1])
                        coord.1 <- min(coord)
                        coord.2 <- max(coord)
                        exon.1 <- paste(chr, coord.1, coord.2, sep=":")
                        
                        # exon 2
                        coord <- as.numeric(strsplit(event.small, split="[\\:|\\|]")[[2]][-1])
                        coord.1 <- min(coord)
                        coord.2 <- max(coord)
                        exon.2 <- paste(chr, coord.1, coord.2, sep=":")
                        
                    # Merge
                    event.reformatted[[i]] <- c(exon.1, exon.2)
                
                }
                
                # Save as new object
                event.A3SS <- unlist(event.reformatted)
                
            } else {
                
                event.A3SS <- NULL
                
            }
        
    # Merge all exons
    exons <- c(event.SE, event.MXE, event.RI, event.A5SS, event.A3SS)
    
    # Create database of chromosomes and ranges (for subsetting read-in BAM)
        # Retrieve list of chromosomes
        exons.split <- strsplit(exons, split="\\:")
        chrs <- unique(sapply(exons.split, function(x) {x[1]}))
        
        # Retrieve ranges
        ranges.list.final <- vector(mode="list", length=length(chrs))
        
        for(j in 1:length(chrs)) {
        
            # Subset exons by chrs
            exons.chr <- exons[grep(paste(chrs[j], "\\:", sep=""), exons)]
            
            # Split
            exon.split <- strsplit(exons.chr, split="\\:")
        
            # Start
            starts <- as.numeric(sapply(exon.split, function(x) {x[2]}))
        
            # End
            ends <- as.numeric(sapply(exon.split, function(x) {x[3]}))
            
            # Range
            ranges.list <- vector(mode="list", length=length(starts))
            
            for(i in 1:length(starts)) {
                
                ranges.list[[i]] <- seq(starts[i], ends[i])
                
            }
                        
            # Save into list
            ranges.list.final[[j]] <- unlist(ranges.list)
            
        }
        
        # Save chrs into list
        chrs.list.final <- as.list(chrs)
        
    # Create GRanges object of chromosomes and ranges (for subsetting BAM)
        # Retrieve list of chromosomes
        exons.split <- strsplit(exons, split="\\:")
        chrs <- sapply(exons.split, function(x) {x[1]})
        
        if(length(header) != 0) {
            
            chrs <- chrs
            
        } else {
            
            chrs <- gsub("chr", "", chrs)
            
        }
        
        # Retrieve start coordinates
        exon.split <- strsplit(exons, split="\\:")
        starts <- as.numeric(sapply(exon.split, function(x) {x[2]})) - 10000
        
        # Compute width
        exon.split <- strsplit(exons, split="\\:")
        ends <- as.numeric(sapply(exon.split, function(x) {x[3]})) + 10000
        widths <- ends - starts + 1
        
        # Create GRanges object
        gr <- GRanges(seqnames=chrs, ranges=IRanges(start=starts, width=widths))
    
    ###################################### COMPUTE PSI ######################################
    
    # Tabulate PSI for each genomic coordinate
    bamfileGA.list <- vector(mode="list", length=length(files))
    
    pb <- utils::txtProgressBar(1, length(files), style=3)
    
    for(i in 1:length(files)) {
        
        # Read file
        bamfileGA.list[[i]] <- readGAlignments(file=paste(BAM, files[i], sep="/"), index=paste(BAM, files[i], sep="/"), param=ScanBamParam(which=gr))
                
        # Track progress
        utils::setTxtProgressBar(pb, i)
        
    }
    
    # Compute PSI
    print("Reading BAM completed, computing PSI now")
    
    psi.collapsed.list <- vector(mode="list", length=length(files))
    
    pb <- utils::txtProgressBar(1, length(files), style=3)

    for(j in 1:length(files)) {

        # Retrieve GAalignment object
        bamfile.GA <- bamfileGA.list[[j]]
        
        # Retrive chrs to analyse
        chrs <- unlist(chrs.list.final)
        
        if(length(header) != 0) {
            
            chrs <- chrs
            
        } else {
            
            chrs <- gsub("chr", "", chrs)
            
        }
   
        psi.chrs.list <- vector(mode="list", length=length(chrs))
        
        for(i in 1:length(chrs)) {
            
            # Retrieve chr
            chr <- chrs[i]
            
            # Retrieve ranges for relevant chr
            ranges <- ranges.list.final[[i]]

            # Retrieve read counts
            read.counts <- as.vector(coverage(bamfile.GA)[[chrs[i]]])[ranges]

            # Retrieve read + skipped counts
            all.counts <- as.vector(coverage(granges(bamfile.GA))[[chrs[i]]])[ranges]
            
            # Set threshold for coverage
            all.counts[which(all.counts < MinCoverage)] <- NA

            # Compute PSI
            psi <- read.counts/all.counts

            # Save as data frame
            psi.df <- data.frame("Chr"=chrs[i], "Position"=ranges, "PSI"=psi, stringsAsFactors=FALSE)

            # Save PSI in list
            psi.chrs.list[[i]] <- psi.df
            
        }
        
        ###################################### COLLAPSE ##########################################

        # Collapse list into data frame
        psi.collapsed <- do.call(rbind.data.frame, psi.chrs.list)

        # Save collapsed PSI into list
        psi.collapsed.list[[j]] <- psi.collapsed

        # Remove GAlignemnts and GRanges object
        remove(bamfile.GA)

        # Track progress
        utils::setTxtProgressBar(pb, j)

    }
    
    ###################################### TABULATE PSI ######################################
    
    # Collapse list into data frame
        # Retrieve Chr column
        col.chr <- psi.collapsed.list[[1]][,1]

        # Retrieve Position column
        col.position <- psi.collapsed.list[[1]][,2]

        # Retrieve the PSI values for each sample
            # Create empty matrix
            n.col <- length(files)
            n.row <- length(col.chr)
            psi.all <- matrix(nrow=n.row, ncol=n.col)

            # Retrieve third column (PSI) of each sample data frame
            for(k in 1:length(psi.collapsed.list)) {

                psi.all[,k] <- psi.collapsed.list[[k]][,3]

            }

            # Convert matrix to data frame
            psi.all <- as.data.frame(psi.all)

            # Use file names as header
            names(psi.all) <- files

        # Merge Chr, Position, PSI columns
        psi.final <- data.frame("Chr"=as.character(col.chr), "Position"=as.character(col.position), psi.all, stringsAsFactors=FALSE)

    # Assign class
    class(psi.final) <- c("rehab", "data.frame")
    return(psi.final)
    
}

