#' Derive the mutational catalogue from a vcf
#'
#' Creates a matrix containing the mutational catalogue from a vcf file or 
#' object. The result can be input to the analysis functions of this package.
#'
#' @param vcf *.vcf file or a vcf object containing variant calling data for one
#'  patient
#' @param genome a BSgenome object corresponding to the genome the variants were
#' called on
#' @param name optional, a sample name 
#' @param seqs optional, a character vector containing the names of the 
#' sequences that are to be included in the mutational profile. If none is given
#' everything will we included
#'
#' @return mutational catalogue (matrix) of a patient containing SNV absolute 
#' counts (in the 96 trinucleotide context)
#' format: 1 by 96
#'
#' @note The execution can take some time, depending on the size of the vcf 
#'
#' @importFrom stats na.omit
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' vcf2mut_cat('test.vcf', BSgenome.Hsapiens.1000genomes.hs37d5)
#' }
#' @export
vcf2mut_cat <- function(vcf, genome, name=NULL, seqs=NULL){

  if(!requireNamespace('VariantAnnotation', quietly = TRUE)){
    stop(paste0('Please install the library 
                VariantAnnotation to use this function.'),
         call. = TRUE)
  }
  
   # if vcf is a file path, first import
  if(class(vcf) == 'character'){
    vcf <- VariantAnnotation::readVcf(vcf)
  }
  
  check <- lapply(c('BSgenome', 'SummarizedExperiment'), function(pack){
    if(!requireNamespace(pack, quietly = TRUE)){
      stop(paste0('Please install the library ', pack, 
                  ' to use this function.'),
           call. = TRUE)
    }
  })
  
   if (missing(vcf) || is.null(vcf) || !is(vcf, "VCF"))
     stop("Missing or illegal vcf", call. = TRUE)
   if (is.null(genome) || !is(genome, "BSgenome"))
     stop("The genome has to be a BSgenome object", call. = TRUE)
  
   # Test if genome names match
   genome_names <- sort(unique(stats::na.omit(GenomeInfoDb::genome(vcf))))
   if (length(genome_names) == 0) {
     warning("vcf has no genome name")
   } else {
     if (!all(genome_names %in% GenomeInfoDb::genome(genome)))
       warning("vcf genome \"",
               paste(genome_names, collapse=", "),
               "\" not in genome \"",
               paste(sort(unique(na.omit(GenomeInfoDb::genome(genome)))), collapse=", "),
               "\"")
   }
  
   # Remove all variants not in the genome
   vcf <- vcf[as.character(
     GenomicRanges::seqnames(vcf)) %in% BSgenome::seqnames(genome),]
   
   # Only keep variants of specified areas
   if (!is.null(seqs)){
     if(!is.character(seqs)){
       stop("seqs has to be a character vector or NULL", call. = TRUE)
     }
     vcf <- vcf[as.character(
       GenomicRanges::seqnames(vcf)) %in% seqs,]
   }
   
   if (nrow(vcf) == 0) {
     warning("no mutation overlap between vcf & genome", call. = TRUE)
     return(NULL)
   }
  
   # Remove all non-SNV variants
   vcf <- vcf[GenomicRanges::width(vcf)==1,]
   if (nrow(vcf) == 0) {
     warning("no SNV mutation", call. = TRUE)
     return(NULL)
   }
   # Set the strand to "unknown", required for promoter
   locii <- SummarizedExperiment::rowRanges(vcf)
   GenomicRanges::strand(locii) <- "*"
  
   # Extract mutation contexts
   mutation_contexts <- GenomicRanges::promoters(
     x=locii,
     upstream=1, downstream=2
   )
  
   # Remove illegal contexts (outside sequence boundaries)
   lens <- GenomeInfoDb::seqlengths(BSgenome::seqinfo(genome))
   i <- which(
     GenomicRanges::start(mutation_contexts)<1 |
     GenomicRanges::end(mutation_contexts)>lens[as.character(
       GenomicRanges::seqnames(mutation_contexts))]
   )
   if (length(i)>0) {
     warning(length(i), " mutations at sequences boundaries, ignored", 
             call. = TRUE)
     mutation_contexts <- mutation_contexts[-i,]
     if (length(i) == 0) return(NULL)
   }
   
   # Extract sequences
   seqs <- Biostrings::DNAStringSet(BSgenome::Views(genome, mutation_contexts))
   ref <- SummarizedExperiment::rowRanges(vcf)$REF
   alt <- unlist(SummarizedExperiment::rowRanges(vcf)$ALT)
   
   #format to match consensus
   mut_cat <- matrix(0L,ncol = 1, nrow = 96)
   rownames(mut_cat) <- rownames(cosmicSigs)
   for(i in seq_len(length(seqs))){
     if(as.character(ref[i]) %in% c('A', 'G')){
       seqsI <- as.character(Biostrings::reverseComplement(seqs[i]))
       altI <- Biostrings::reverseComplement(alt[i])
     }
     else{
       altI <- alt[i]
       seqsI <- as.character(seqs[i])
     }
     context <- paste0(substring(seqsI,1,1), '[',
                       substring(seqsI,2,2), '>',
                       altI, ']',
                       substring(seqsI,3,3))
     mut_cat[which(context == rownames(mut_cat))] <- mut_cat[which(context == rownames(mut_cat))] + 1
   }
   
   if (!is.null(name) && is.character(name)){
    colnames(mut_cat) <- name 
   }
   
   return(mut_cat)
}
