#!/usr/bin/env Rscript
wdir = getwd( ) 
print(wdir)

print("Initializing the package...")
source(paste(wdir,"/R/MakeOptions_modDetect.R",sep = ""))
library(data.table)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

wdir = getwd( ) 
acc.bases = c("A","T","C","G", "-", "+")
print("------------------------------------")

print("Loading input files...")
if (opt$bed != ""){
  bedFile = read.table(as.character(opt$bed), header=FALSE, sep="\t", quote="")
  bedFile = bedFile[which(bedFile$V8 == "gene"),] # Filter to only include BED file entries that are genes
  bedFile = bedFile[, c(1,2,3,4,6,10)]
  colnames(bedFile) = c("chr", "start", "end", "gene", "strand", "metadata")
}
if (opt$file3 != ""){
  kmer.summary.df = read.csv(as.character(opt$file3))
}
fasta_file <- FaFile( as.character(opt$file4) )
Direct1.bam.dir = as.character(opt$file) 
ctrl.bam.dir    = as.character(opt$file2)
print("------------------------------------")
print("Data Preparation...")

D1b <- BamFile(Direct1.bam.dir)
ctrlb <- BamFile(ctrl.bam.dir)

finalOutput = data.frame()
totalRows = nrow(bedFile)
for (i in c(1:nrow(bedFile))){
  print(paste("Have examined this many genes:", i, "out of", totalRows))

  param <- ScanBamParam(which=GRanges(strand = bedFile$strand[i],
                                      seqnames = bedFile$chr[i],
                                      ranges = IRanges(start=bedFile$start[i], 
                                                       end=bedFile$end[i])))
  pilup_params =  Rsamtools::PileupParam(max_depth = 200000,min_mapq = 1,distinguish_nucleotides = T,
                                         ignore_query_Ns = T, min_base_quality = 1, include_insertions = T)
  
  PU = pileup(D1b, scanBamParam=param,pileupParam = pilup_params)
  # if no reads mapped to that gene, move on
  if (nrow(PU) == 0) {
    next
  }
  if (bedFile$strand[i] != "*"){
    PU = PU[which(PU$strand == bedFile$strand[i]),]
  }

  if (nrow(PU) == 0){
    next
  }
  PU$nucleotide.new = PU$nucleotide
  PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="A")] = "T"
  PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="T")] = "A"
  PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="C")] = "G"
  PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="G")] = "C"
  PU = PU[which(PU$nucleotide %in% acc.bases),]
  

  PU.ctrl = pileup(ctrlb, scanBamParam=param,pileupParam = pilup_params)
  # if no reads mapped to that gene, move on
  if (nrow(PU.ctrl) == 0) {
    next
  }
  if (bedFile$strand[i] != "*"){
    PU.ctrl = PU.ctrl[which(PU.ctrl$strand == bedFile$strand[i]),]
  }
  if (nrow(PU.ctrl) == 0){
    next
  }
  PU.ctrl$nucleotide.new = PU.ctrl$nucleotide
  PU.ctrl$nucleotide.new[which(PU.ctrl$strand == "-" & PU.ctrl$nucleotide =="A")] = "T"
  PU.ctrl$nucleotide.new[which(PU.ctrl$strand == "-" & PU.ctrl$nucleotide =="T")] = "A"
  PU.ctrl$nucleotide.new[which(PU.ctrl$strand == "-" & PU.ctrl$nucleotide =="C")] = "G"
  PU.ctrl$nucleotide.new[which(PU.ctrl$strand == "-" & PU.ctrl$nucleotide =="G")] = "C"
  PU.ctrl = PU.ctrl[which(PU.ctrl$nucleotide %in% acc.bases),]
  
  print(paste("Examining the following gene: ", bedFile$gene[i]))

  res.df = data.frame("chr" = bedFile$chr[i], "pos" = unique(PU$pos), "gene" = bedFile$gene[i], "metadata" = bedFile$metadata[i], "A.count"=0,"C.count"=0,"G.count"=0,"T.count"=0,"Deletion.count"=0, "Insertion.count"=0,
                      "A.count.ctrl"=0,"C.count.ctrl"=0,"G.count.ctrl"=0,"T.count.ctrl"=0, "Deletion.count.ctrl"=0, "Insertion.count.ctrl"=0,
                      "reference"="","kmer"="","strand"=bedFile$strand[i])
  for (j in c(1:nrow(PU))){
    ro.res.df = which(res.df$pos == PU$pos[j])
    
    if (res.df$reference[ro.res.df] != "T" & res.df$reference[ro.res.df] != "" & !opt$all){
        next
    }

    # assign reference value
    grange_ref = GRanges(bedFile$chr[i],IRanges(start = res.df$pos[ro.res.df], end = res.df$pos[ro.res.df]))
    referenceBase = getSeq(fasta_file, grange_ref)
    referenceBase = as.data.frame(referenceBase)$x
    referenceBaseNew = referenceBase
      
    if (res.df$strand[ro.res.df] == "-" & referenceBase == "A") {referenceBaseNew = "T"}
    if (res.df$strand[ro.res.df] == "-" & referenceBase == "T") {referenceBaseNew = "A"}
    if (res.df$strand[ro.res.df] == "-" & referenceBase == "G") {referenceBaseNew = "C"}
    if (res.df$strand[ro.res.df] == "-" & referenceBase == "C") {referenceBaseNew = "G"}
    res.df$reference[ro.res.df] = referenceBaseNew
    
    if (res.df$reference[ro.res.df] != "T" & !opt$all){
      next
    }

    col = 5
    if (PU$nucleotide.new[j] == "C"){col=6}
    if (PU$nucleotide.new[j] == "G"){col=7}
    if (PU$nucleotide.new[j] == "T"){col=8}
    if (PU$nucleotide.new[j] == "-"){col=9}
    if (PU$nucleotide.new[j] == "+"){col=10}

    res.df[ro.res.df, col] = res.df[ro.res.df, col] + PU$count[j]
  }
  # only extract positions with a T in the reference
  if (!opt$all){
    res.df = res.df[which(res.df$reference == "T"),]
  }
  if (nrow(res.df) == 0) {
    next
  }

  for (j in c(1:nrow(PU.ctrl))){

    rowValue = which(res.df$pos == PU.ctrl$pos[j])
    if (identical(rowValue, integer(0))){
      next
    }
    ro.res.df = rowValue
    col = 11
    if (PU.ctrl$nucleotide.new[j] == "C"){col=12}
    if (PU.ctrl$nucleotide.new[j] == "G"){col=13}
    if (PU.ctrl$nucleotide.new[j] == "T"){col=14}
    if (PU.ctrl$nucleotide.new[j] == "-"){col=15}
    if (PU.ctrl$nucleotide.new[j] == "+"){col=16}
    
    res.df[ro.res.df, col] = res.df[ro.res.df, col] + PU.ctrl$count[j]
  }
  for (j in c(1:nrow(res.df))){
    kmer = ""
    for (k in c(-2:2)){
      tryCatch({
        gr1 <- GRanges(res.df$chr[j],IRanges(start=res.df$pos[j]+k, end=res.df$pos[j]+k))
        ### Extract the kmers
        refbase <- getSeq(fasta_file, gr1)
        refbase <- as.data.frame(refbase)$x
        res.df.new = refbase
        if (refbase == "A" & res.df$strand[j] == "-"){res.df.new = "T"}
        if (refbase == "T" & res.df$strand[j] == "-"){res.df.new = "A"}
        if (refbase == "C" & res.df$strand[j] == "-"){res.df.new = "G"}
        if (refbase == "G" & res.df$strand[j] == "-"){res.df.new = "C"}
        
        kmer = paste(kmer,res.df.new,sep = "")
        if (k == 0){
          res.df$reference[j] = res.df.new
        }
      }, error = function(e) {
        
      })
    }
    if (res.df$strand[j] == "-"){
      kmer = paste(substring(kmer, 5:1, 5:1), collapse = "")
    }
    res.df$kmer[j] = kmer
  }
  if (nrow(res.df) > 0){
    finalOutput = rbind(finalOutput, res.df)
  }
  print(paste("Done examining the following gene: ", bedFile$gene[i]))
}

print("------------------------------------")
print("Calculating p-values...")

if (nrow(finalOutput) == 0) {
  stop('\r No data to display')
}

if (opt$del){
  finalOutput$NReads_test = finalOutput$Deletion.count + finalOutput$T.count + finalOutput$A.count + finalOutput$G.count + finalOutput$C.count + finalOutput$Insertion.count
  finalOutput$NReads_ctrl = finalOutput$T.count.ctrl + finalOutput$Deletion.count.ctrl + finalOutput$A.count.ctrl + finalOutput$G.count.ctrl + finalOutput$C.count.ctrl + finalOutput$Insertion.count.ctrl
  finalOutput = finalOutput[which((finalOutput$Deletion.count + finalOutput$T.count + finalOutput$A.count + finalOutput$G.count + finalOutput$C.count + finalOutput$Insertion.count) > 0),] # Filtering to only include positions with at least 1 read in TEST (avoid dividing by 0 when calculating error)
  finalOutput = finalOutput[which((finalOutput$T.count.ctrl + finalOutput$Deletion.count.ctrl + finalOutput$A.count.ctrl + finalOutput$G.count.ctrl + finalOutput$C.count.ctrl + finalOutput$Insertion.count.ctrl) > 0),] # Filtering to only include positions with at least 1 read in CTRL (avoid dividing by 0 when calculating error)
  finalOutput$mm.perc = finalOutput$Deletion.count / (finalOutput$NReads_test)*100
  finalOutput$ctrl.err = finalOutput$Deletion.count.ctrl / (finalOutput$T.count.ctrl + finalOutput$Deletion.count.ctrl + finalOutput$A.count.ctrl + finalOutput$G.count.ctrl + finalOutput$C.count.ctrl + finalOutput$Insertion.count.ctrl)*100
} else if (opt$all){
  finalOutput$NReads_test = finalOutput$Deletion.count + finalOutput$T.count + finalOutput$A.count + finalOutput$G.count + finalOutput$C.count + finalOutput$Insertion.count
  finalOutput$NReads_ctrl = finalOutput$T.count.ctrl + finalOutput$Deletion.count.ctrl + finalOutput$A.count.ctrl + finalOutput$G.count.ctrl + finalOutput$C.count.ctrl + finalOutput$Insertion.count.ctrl
  finalOutput = finalOutput[which((finalOutput$Deletion.count + finalOutput$T.count + finalOutput$A.count + finalOutput$G.count + finalOutput$C.count + finalOutput$Insertion.count) > 0),] # Filtering to only include positions with at least 1 read in TEST (avoid dividing by 0 when calculating error)
  finalOutput = finalOutput[which((finalOutput$T.count.ctrl + finalOutput$Deletion.count.ctrl + finalOutput$A.count.ctrl + finalOutput$G.count.ctrl + finalOutput$C.count.ctrl + finalOutput$Insertion.count.ctrl) > 0),] # Filtering to only include positions with at least 1 read in CTRL (avoid dividing by 0 when calculating error)
  
  finalOutput$mm.perc <- apply(finalOutput, 1, function(row) {
    count_column_name <- paste0(row["reference"], ".count")
    mismatch <- row[count_column_name]
    return(as.numeric(mismatch))
  })
  finalOutput$mm.perc <- ((finalOutput$NReads_test - finalOutput$mm.perc) / finalOutput$NReads_test) * 100

  finalOutput$ctrl.err <- apply(finalOutput, 1, function(row) {
    count_column_name <- paste0(row["reference"], ".count.ctrl")
    mismatch <- row[count_column_name]
    return(as.numeric(mismatch))
  })
  finalOutput$ctrl.err <- ((finalOutput$NReads_ctrl - finalOutput$ctrl.err) / finalOutput$NReads_ctrl) * 100
  
} else {
  finalOutput$NReads_test = finalOutput$C.count + finalOutput$T.count
  finalOutput$NReads_ctrl = finalOutput$C.count.ctrl + finalOutput$T.count.ctrl
  finalOutput = finalOutput[which(finalOutput$NReads_test > 0),] # Filtering to only include positions where at least 1 C or T was read in TEST (avoid dividing by 0 when calculating error)
  finalOutput = finalOutput[which((finalOutput$T.count.ctrl + finalOutput$C.count.ctrl) > 0),] # Filtering to only include positions where at least 1 C or T was read in CTRL (avoid dividing by 0 when calculating error)
  finalOutput$mm.perc = (finalOutput$C.count / finalOutput$NReads_test)*100
  finalOutput$ctrl.err = finalOutput$C.count.ctrl / (finalOutput$NReads_ctrl)*100
}


# Filtering to only show positions where mismatch percentage in the test is greater 0% (i.e. there is some mismatch)
finalOutput = finalOutput[which(finalOutput$mm.perc > 0),]
calc.p.val = function(n.read, mod.prob, subject.mm){
  n.err = round(subject.mm * n.read)
  n.err.lst = c(n.err:n.read)
  p.val = 0
  for (n.err in n.err.lst){
    comb = exp(lgamma(n.read+1) - lgamma(n.read-n.err+1)-lgamma(n.err+1))
    prob.analytical = comb * (mod.prob)^n.err *(1-mod.prob)^(n.read-n.err)
    p.val = p.val + prob.analytical
  }
  return (p.val)
}
finalOutput$expected.err = finalOutput$ctrl.err

if (opt$file3 != ""){
  finalOutput$kmer.err = 0
  for (i in c(1:nrow(kmer.summary.df))) {
    finalOutput$kmer.err[which(finalOutput$kmer == kmer.summary.df$kmer[i])] = kmer.summary.df$avg.mm.ivt[i]
  }
  
  finalOutput$expected.err = finalOutput$kmer.err
  # Finds all positions where control error is greater than kmer error. Then, sets the 'expected error' value for those positions to control error. 
  finalOutput$expected.err[which(finalOutput$ctrl.err>finalOutput$kmer.err)] = finalOutput$ctrl.err[which(finalOutput$ctrl.err>finalOutput$kmer.err)]
}




finalOutput$p=1
for (i in c(1:nrow(finalOutput))){
  if (finalOutput$mm.perc[i]>finalOutput$expected.err[i]){
    finalOutput$p[i] = calc.p.val(n.read = finalOutput$NReads_test[i],
                                  mod.prob = finalOutput$expected.err[i]/100, 
                                  subject.mm = finalOutput$mm.perc[i]/100)
    counter = 1
    step = 100
    while (is.na(finalOutput$p[i])){
      finalOutput$p[i] = calc.p.val(n.read = finalOutput$NReads_test[i]-counter*step,
                                    mod.prob = finalOutput$expected.err[i]/100, 
                                    subject.mm = finalOutput$mm.perc[i]/100)
      counter =  counter+1
    }
  }
}


finalOutput$NReads_test <- NULL
finalOutput$NReads_ctrl <- finalOutput$Deletion.count.ctrl + finalOutput$T.count.ctrl + finalOutput$A.count.ctrl + finalOutput$G.count.ctrl + finalOutput$C.count.ctrl + finalOutput$Insertion.count.ctrl
finalOutput$NReads_test <- finalOutput$Deletion.count + finalOutput$T.count + finalOutput$A.count + finalOutput$G.count + finalOutput$C.count + finalOutput$Insertion.count
finalOutput = finalOutput[which(finalOutput$NReads_test >= as.numeric(opt$read_depth) & finalOutput$NReads_ctrl >= as.numeric(opt$read_depth)),] # Applies read depth cutoff
finalOutput = finalOutput[which(finalOutput$mm.perc - as.numeric(opt$percentage_mismatch) >= finalOutput$expected.err),] # Applies mismatch percentage cutoff
finalOutput <- finalOutput[which(finalOutput$p <= as.numeric(opt$max_p)),] # Applies p-value filtering
finalOutput <- finalOutput[,c(ncol(finalOutput),c(1:(ncol(finalOutput)-1)))]
setDT(finalOutput)
finalOutput <- finalOutput[, lapply(.SD, function(x) paste(unique(x), collapse = ", ")), 
                     by = .(chr, pos)]
write.csv(as.data.frame(finalOutput), opt$out,row.names = F)

print("Done! Please find the output file here:")
print(opt$out)


