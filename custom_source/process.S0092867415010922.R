# Pre-loaded:
# require(data.table)
# require(rtracklayer)

chromosomes = paste0("chr", c(1:7, 10:21, "9:22", "22:9", "X"))

grefDT = fread(file.path("/media/MiSo/bioRDSmaker/custom_source/",
		"hg19.chr_size.HAP1.Philadelphia_corrected.txt"),
	col.names = c("seqnames", "end"))
grefDT$start = 1
gref = GRanges(grefDT)
isCircular(gref) = rep(F, length(seqlevels(gref)))
genome(gref) = rep("hg19.HAP1.transCorr", length(seqlevels(gref)))
setkeyv(grefDT, "seqnames")
seqlengths(gref) = grefDT[seqlevels(gref), end]
gref = keepSeqlevels(gref, chromosomes, pruning.mode = "coarse")
gref = sortSeqlevels(gref)

bins = list(
	bins.1MbSize.100kbStep = import.bed(file.path(
		"/media/MiSo/bioRDSmaker/bins/",
		"hg19.bins.1MbSize.100kbStep.HAP1.transCorrected.bed"
	)),
	bins.100kbSize.100kbStep = import.bed(file.path(
		"/media/MiSo/bioRDSmaker/bins/",
		"hg19.bins.100kbSize.100kbStep.HAP1.transCorrected.bed"
	))
)
for( i in seq_along(bins) ){
	bins[[i]] = keepSeqlevels(bins[[i]], chromosomes, pruning.mode = "coarse")
	bins[[i]] = sortSeqlevels(bins[[i]])
	seqinfo(bins[[i]]) = seqinfo(gref)
}

# x: path to file
# bins: GRanges object with genomic bins
# chromosomes: list of chromosomes to evaluate

tmp = fread(x, header = T, key = c("seqnames", "start", "end"))
tmp = tmp[, .(score = rowMeans(.SD, na.rm = T)), by = key(tmp)]

tmp2 = GRanges(tmp)
tmp2 = keepSeqlevels(tmp2, chromosomes, pruning.mode = "coarse")
tmp2 = sortSeqlevels(tmp2)
seqinfo(tmp2) = seqinfo(gref)
tmp2 = trim(tmp2)
tmp2 = coverage(tmp2, weight = tmp2$score)

out = rbindlist(lapply(bins, process_single_bin,
	ftype = "custom", tmp = tmp2), idcol = "bins")
setnames(out, "seqnames", "chrom")

# out is returned to the main script
