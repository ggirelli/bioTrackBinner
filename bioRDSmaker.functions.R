
getFullExt = function(x) {
	fileName = tools::file_path_sans_ext(x, compression = TRUE)
	ext = gsub(fileName, "", x)
	ext = gsub("^\\.", "", ext)
	return(ext)
}

process_single_bin = function(b, ftype, tmp){
	b = keepSeqlevels(b, intersect(seqlevels(b), seqlevels(tmp)), pruning.mode = "coarse")
	sco = sapply(seq(1, length(b), 1e3),
		function(bi) {
			bi = as.numeric(as.character(bi))
			sub = b[bi:min(bi + 1e3 - 1, length(b)),]
			mean(tmp[sub], na.rm = T)
		})
	sco = unlist(sco)
	data.table(as.data.frame(b)[,1:3], type = ftype, score = sco)
}

processBed = function(x, bins, chromosomes) {
	tmp = import.bed(x,
		seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19),
		trackLine = FALSE)
	tmp = keepSeqlevels(tmp, chromosomes, pruning.mode = "coarse")
	if ( !"score" %in% names(tmp) ) tmp$score = 1
	tmp = coverage(tmp, weight = tmp$score)
	rbindlist(lapply(bins, process_single_bin,
		ftype = "bed", tmp = tmp), idcol = "bins")
}

processBedGraph = function(x, bins, chromosomes) {
	tmp = import.bedGraph(x, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19))
	tmp = keepSeqlevels(tmp, chromosomes, pruning.mode = "coarse")
	tmp = coverage(tmp, weight = tmp$score)
	rbindlist(lapply(bins, process_single_bin,
		ftype = "bedGraph", tmp = tmp), idcol = "bins")
}

processBigWig = function(x, bins, chromosomes){
	tmp = import.bw(x, as = "Rle")
	if ( all(!grepl("^chr", names(tmp))) ) {
		seqnames(seqinfo(tmp)) = paste0("chr", seqnames(seqinfo(tmp)))
		names(tmp) = paste0("chr", names(tmp))
	}
	chromosomes = intersect(chromosomes, seqlevels(tmp))
	tmp = keepSeqlevels(tmp, chromosomes, pruning.mode = "coarse")
	seqinfo(tmp) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)[chromosomes]
	rbindlist(lapply(bins, process_single_bin,
		ftype = "bigWig", tmp = tmp), idcol = "bins")
}

processCod = function(x, bins, chromosomes, compressed = F){
	if ( compressed ) {
		tmp = fread(cmd = paste("gunzip -c", x),
			header = TRUE, sel = c(2:4, 15),
			col.names = c("chrom", "start", "end", "score"))
	} else {
		tmp = fread(x,
			header = TRUE, sel = c(2:4, 15),
			col.names = c("chrom", "start", "end", "score"))
	}
	tmp = with(tmp, GRanges(chrom, IRanges(start, end),
		score = score, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19)))
	tmp = keepSeqlevels(tmp, chromosomes, pruning.mode = "coarse")
	tmp = coverage(tmp, weight = tmp$score)
	rbindlist(lapply(bins, process_single_bin,
		ftype = "Cod", tmp = tmp), idcol = "bins")
}

processText = function(x, bins, chromosomes){
	tmp = fread(x)
	tmp = tmp[, list(chrom = paste0("chr", Chromosome),
		start = Start, end = End, score = `Wild-type`)]
	tmp = keepSeqlevels(with(tmp, GRanges(chrom, IRanges(start, end),
			score = score, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19))),
		chromosomes, pruning.mode = "coarse")
	tmp = coverage(tmp, weight = tmp$score)
	rbindlist(lapply(bins, process_single_bin,
		ftype = "txt", tmp = tmp), idcol = "bins")
}

processXlsx = function(x, bins, chromosomes) {
	tmp = data.table(read_excel(x, sheet = 1, na = "NA",
		col_types = c("text", rep("numeric", 4))))
	tmp = tmp[, list(chrom = paste0("chr", Chromosome),
		start = Start, end = End, score = `Wild-type`)]
	tmp = keepSeqlevels(with(tmp, GRanges(chrom, IRanges(start, end),
		score = score, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19))),
		chromosomes, pruning.mode = "coarse")
	tmp = coverage(tmp, weight = tmp$score)
	rbindlist(lapply(bins, process_single_bin,
		ftype = "xlsx", tmp = tmp), idcol = "bins")
}
