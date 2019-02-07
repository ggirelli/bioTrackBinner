#!/usr/bin/env Rscript

require(argparser)
parser = arg_parser(paste0(
		'Bins biological data tracks and generates corresponding RDS files. ',
		'The reference table should be a tab-separated no-header table with ',
		'one line per track and the following columns: ',
		'cellLine, markerType, accessionID, replicateID, pathToTrack. ',
		'The bins table should be a tab-separated no-header table with one ',
		'line per bin file and the following columns: label, pathToBins. ',
		'Bin files are expected to be in bed format.'
	), name = 'raw2rds.R')
parser = add_argument(parser, arg = 'referenceTable',
	help = 'Path to reference table.')
parser = add_argument(parser, arg = 'binsTable',
	help = 'Path to bins table.')
parser = add_argument(parser, arg = 'outputFolder',
	help = 'Path to output folder, created if missing.')
parser = add_argument(parser, arg = '--elective', short = '-e', type = class(0),
	help = 'Elective argument.', default = 1, nargs = 1)
p = parse_args(parser)
attach(p['' != names(p)])

chromosomes = paste0("chr", c(1:22, "X", "Y"))

require(data.table)
require(pbapply)
require(readxl)
require(rtracklayer)

require(BSgenome.Hsapiens.UCSC.hg19)

source("/media/MiSo/GPSeq/analysis/GG/src/functions.common.R")
source("bioRDSmaker.functions.R")

# hg19 bins
bins = list(
	bin.100kbSize.100kbStep  = import.bed(file.path("bins",
		"hg19.bins.100kbSize.100kbStep.bed")),
	bin.1MbSize.100kbStep = import.bed(file.path("bins",
		"hg19.bins.1MbSize.100kbStep.bed"))
)
for( i in seq_along(bins) ){
	bins[[i]] = keepSeqlevels(bins[[i]], chromosomes,
		pruning.mode = "coarse")
	seqinfo(bins[[i]]) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
}

ref = fread(cmd = sprintf('grep -v "^#" %s', referenceTable),
	header = FALSE,
	col.names = c("cell_line", "marker", "accession", "rep", "path"),
	key = "accession")

bw.files = with(ref, setNames(path, accession))
allFilesExist = file.exists(bw.files)
if ( !all(file.exists(bw.files)) )
	print(bw.files[!allFilesExist])
stopifnot( all(allFilesExist) )

epig = rbindlist(lapply(seq_along(bw.files),
	function(i, bins) {
		x = bw.files[i]
		cat(sprintf("Processing '%s'\n", x))
		ext = getFullExt(x)
		fileList = file.path(outputFolder,
			paste0(names(bw.files)[i], ".epig.", names(bins), ".rds"))
		if ( any(!file.exists(fileList)) ) {
			processed = F
			if ( ext %in% c("bed", "bed.gz") ) {
				out = processBed(x, bins, chromosomes)
				processed = T
			}
			if ( "bdg" == ext ) {
				out = processBedGraph(x, bins, chromosomes)
				processed = T
			}
			if ( "bigWig" == ext ) {
				out = processBigWig(x, bins, chromosomes)
				processed = T
			}
			if ( "cod.gz" == ext ) {
				out = processCod(x, bins, chromosomes, T)
				processed = T
			}
			if ( "cod" == ext ) {
				out = processCod(x, bins, chromosomes)
				processed = T
			}
			if ( "txt" == ext ) {
				out = processText(x, bins, chromosomes)
				processed = T
			}
			if ( "xlsx" == ext ) {
				out = processXlsx(x, bins, chromosomes)
				processed = T
			}
			if ( processed ) {
				out$accession = names(bw.files)[i]
				return(out)
			} else {
				cat("Skipped dataset entirely...\n")
			}
		} else {
			cat("Already processed, reading...\n")
			for ( fileName in fileList )
				return(readRDS(fileName))
		}
	}, bins
))

setkeyv(epig, "accession")
epig = epig[ref[, list(accession, cell_line, marker, rep)], nomatch = 0]
setnames(epig, "seqnames", "chrom")
setcolorder(epig, c("chrom", "start", "end", "score",
	"cell_line", "marker", "type", "rep", "accession"))
epig[, `:=`(start = as.numeric(start), end = as.numeric(end))]

epig = getTranslocationCoordinates_binned(epig)

dt = split(epig, epig[, bins])
for ( i in seq_along(bins) ) {
	dt[[i]][, bins := NULL]

	cond = (dt[[i]][, end - start]) != (width(bins[[i]])[1]-1)
	dt[[i]][cond, end := end + (width(bins[[i]])[1]-1 - (end - start))]

	at = split(dt[[i]], unique(dt[[i]][, accession]))
	for ( j in seq_along(at) ) {
		fName = sprintf("%s.epig.%s.rds", names(at)[j], names(bins)[i])
		if ( !file.exists(fName) )
			saveRDS(at[[j]], file = file.path(outputFolder, fName))
	}

	saveRDS(dt[[i]], file = sprintf("epig.%s.rds", names(bins)[i]))
}