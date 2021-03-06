#' Normalize mutational catalogues 
#' 
#' Normalizes the catalogues to a target distribution (e.g. to match the 
#' distribution of the reference signatures).
#'
#' @param mut_cat mutational catalogues (96 by n, n being the amounts of 
#' catalogues) that will be normalized. The tri-nucleotide contexts are expected 
#' to be in the default lexicographical order (see simulated data or cosmicSigs)
#' @param source_context Distribution of tri-nucleotides in the source region.
#' @param target_context Distribution of tri-nucleotides in the target region.
#' Defaults to the context frequencies of BSgenome.Hsapiens.UCSC.hg19 since that
#' corresponds to the COSMIC signatures
#'
#' @return mutational catalogues (96 by n, n being the amounts of catalogues) 
#' normalized to match the target distribution (context)
#' 
#' @examples
#' # this is a toy example:
#' #create mutational catalogue:
#' sim_data <- create_mut_catalogues(1, 500)[['catalogues']]
#' # get trinucleotide frequencies for the genome:
#' genome_context <- get_context_freq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#' #get trinucleotide frequencies for a specific region:
#' gr<-GenomicRanges::GRanges(seqnames=c("chr1"),
#'           ranges=IRanges::IRanges(start=c(100000),end=c(1000000)),
#'           strand=c("+"))
#' region_context<-get_context_freq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr)
#' #normalize data:
#' normalized_mut_cat <- normalize(sim_data, region_context, genome_context)
#' 
#' \dontrun{
#' # get the tri-nucleotide distribution of an exome region
#' exome_contexts <- get_context_freq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
#'                                    'example_exome.bed')
#' 
#' # normalize the mutational catalogue to match the COSMIC signatures
#' normalized_mut_cat <- normalize(mut_cat, exome_contexts, hg19context_freq)
#' }
#' 
#' @note The output from get_context_freq() can be used as input to this 
#' function
#'
#' @export
normalize <- function(mut_cat, source_context,
                      target_context=get(utils::data("hg19context_freq", package="SigsPack"))){
    # TODO not usable for other contexts then the default 96 for now...
    
    if (! is.numeric(mut_cat)){
        stop('The mutational catalogues areexpected to contain only numbers.',
             call. = TRUE)
    }
    m <- as.matrix(mut_cat)
    if (nrow(m) != 96){
        stop('The mutational catalogues are expected to have 96 rows.',
             call. = TRUE)
    }    
    if (! is.numeric(source_context)){
        stop('The source context is expected to contain only numbers.',
             call. = TRUE)
    }
    if (! is.numeric(target_context)){
        stop('The target_context is expected to contain only numbers.',
             call. = TRUE)
    }
    if (length(source_context) != length(target_context)){
        stop('The source and target contexts have to have the same dimensions',
             call. = TRUE)
    }
    # more checks probably needed
    
    P <- get(utils::data("cosmicSigs", package="SigsPack"))
    rownames(m) <- rownames(P)
    
    # add a triplet column to the catalogue to compare to the contexts
    for(i in seq_len(96)){
        rownames(m)[i] <- paste0(substring(rownames(P)[i],1,1),
                                 substring(rownames(P)[i],3,3),
                                 substring(rownames(P)[i],7,7))
    }
    
    #normalize
    for(triplet in rownames(as.data.frame(target_context))){
        m[which(rownames(m) == triplet),] <- m[which(rownames(m) == triplet),]/
            (source_context[triplet,])*(target_context[triplet,])
    }
    
    rownames(m) <- rownames(P)
    m <- apply(m, 2, function(x) {
        x/sum(x)
    })
    
    return(m)
}
