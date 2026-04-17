#### ====================================== ####
####    Make/define the main code options   ####
#### ====================================== ####
#wdir = getwd( ) 
source(paste(wdir,"/R/HandlePackages.R",sep = ""))

option_list = list(
  make_option(c("-f", "--file"), type="character", default="", 
              help="path to direct seq. bam file", metavar="character"),
  make_option(c("-g", "--file2"), type="character", default="", 
              help="path to control bam file", metavar="character"),
  make_option(c("-k", "--file3"), type="character", default="", 
              help="path to kmer-dependent error file", metavar="character"),
  make_option(c("-r", "--file4"), type="character", default="", 
              help="Reference genome fasta file", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default="", 
              help="path to BED file", metavar="character"),
  make_option(c("-s", "--start_position"), type="character", default=NULL, 
              help="start position", metavar="character"),
  make_option(c("-e", "--end_position"), type="character", default=NULL, 
              help="end position", metavar="character"),
  make_option(c("-c", "--chromosome"), type="character", default=NULL, 
              help="chromosome, EX: chr1", metavar="character"),
  make_option(c("-m", "--max_p"), type="character", default=NULL, 
              help="Sites with a p-value this low or smaller will be outputted", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="psi_candidates.csv", 
              help="output file name of psi candidate sites", metavar="character"),
  make_option(c("-d", "--del"), action="store_true", default=FALSE,
              help="run statistics on deletion signatures"),
  make_option(c("-a", "--all"), action="store_true", default=FALSE,
              help="examine mismatches at all nucleotides"),
  make_option(c("-x", "--read_depth"), type="character", default=NULL, 
              help="read depth cutoff, Ex: 3", metavar="character"),
  make_option(c("-y", "--percentage_mismatch"), type="character", default=NULL, 
              help="percentage mismatch delta cutoff, Ex: 10", metavar="character")
); 