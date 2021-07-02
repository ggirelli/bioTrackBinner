
# x: path to file
# bins: GRanges object with genomic bins
# chromosomes: list of chromosomes to evaluate

tmp = import.bed(x,
	seqinfo = seqinfo(BSgenome.Mmusculus.UCSC.mm10),
	trackLine = FALSE)

tmp = as.data.table(tmp)
tmp[, c("start", "end", "score") := .((start+end)/2, (start+end)/2+1, 1)]

chromLevels = intersect(chromosomes, unique(tmp$seqnames))
tmp = GRanges(tmp)
seqinfo(tmp) = seqinfo(BSgenome.Mmusculus.UCSC.mm10)

tmp = keepSeqlevels(tmp, chromLevels, pruning.mode = "coarse")

tmp = coverage(tmp, weight = tmp$score)

out = rbindlist(lapply(bins, process_single_bin,
	ftype = "bed", tmp = tmp), idcol = "bins")

# out is returned to the main script
