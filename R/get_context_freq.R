#' Extract occurence of tri-nucleotide contexts 
#' 
#' Extracts the frequencies of the tri-nucleotide contexts in a given region 
#' of the genome. These frequencies are needed to normalize a mutational 
#' catalogue. The output can be input to normalize().
#'
#' @param genome a BSgenome object
#' @param region a GRanges object, path, URL, connection or BEDFile object. 
#'
#' @return matrix containing the frequencies of the trinucleotide contexts
#' @examples 
#' gr<-GenomicRanges::GRanges(seqnames=c("chr1"),
#'           ranges=IRanges::IRanges(start=c(100000),end=c(1000000)),
#'           strand=c("+"))
#' get_context_freq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr)
#' get_context_freq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#' 
#' \dontrun{
#' get_context_freq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 'example.bed')
#' }
#' 
#' @export
get_context_freq <- function(genome, region=NULL){
    check <- lapply(c('BSgenome'), function(pack){
        if(!requireNamespace(pack, quietly = TRUE)){
            stop(paste0('Please install the library ', pack, 
                        ' to use this function.'),
                 call. = TRUE)
        }
    })
    
    if(!is(genome, 'BSgenome')){
        stop('Genome has to be a BSgenome object.',
             call. = TRUE)
    }
    #=============================================================================
    wanted <- c("ACA", "ACC", "ACG", "ACT", "ATA", "ATC", "ATG", "ATT", "CCA", 
                "CCC", "CCG", "CCT", "CTA", "CTC", "CTG", "CTT", "GCA", "GCC",
                "GCG", "GCT", "GTA", "GTC", "GTG", "GTT", "TCA", "TCC", "TCG",
                "TCT", "TTA", "TTC", "TTG", "TTT")
    #=============================================================================
    
    if (is.null(region)){
        tri_context <- matrix(ncol = 64, nrow = length(BSgenome::seqnames(genome)))
        
        # rewrite this somehow with apply?
        i <- 1
        for(chr in BSgenome::seqnames(genome)){
            tri_context[i,] <- Biostrings::trinucleotideFrequency(genome[[chr]])
            i <- i+1
        }
    }
    else{
        
        if(!is(region, 'GRanges')){
            region <- rtracklayer::import(con = region)
        }
        seq <- BSgenome::getSeq(genome, region)
        tri_context <- Biostrings::trinucleotideFrequency(seq)
    }
    
    tri_context_all <- as.data.frame(colSums(tri_context))
    rownames(tri_context_all) <- rownames(as.data.frame(
        Biostrings::trinucleotideFrequency(genome[[as.character(BSgenome::seqnames(genome)[1])]])))
    context_freq <- matrix(ncol = 1, nrow = 32)
    rownames(context_freq) <- wanted
    
    for (s in wanted){
        # get the reverse complement
        sr <- Biostrings::reverseComplement(Biostrings::DNAString(s))
        sr <- as.character(sr)
        context_freq[s,] <- tri_context_all[s,] + tri_context_all[sr,]
    }
    
    return(context_freq)
}
