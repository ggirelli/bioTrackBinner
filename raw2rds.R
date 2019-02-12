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
	help = paste0('Path to output folder for intermediate single-track RDS',
		', created if missing.'))
p = parse_args(parser)
attach(p['' != names(p)])

chromosomes = paste0("chr", c(1:22, "X", "Y"))

require(data.table)
require(pbapply)
require(readxl)
require(rtracklayer)
require(BSgenome.Hsapiens.UCSC.hg19)
source("coordCorrect.functions.R")
source("bioRDSmaker.functions.R")

cat("Reading bins...\n")
binRef = fread(cmd = sprintf('grep -v "^#" %s', binsTable), header = F,
	col.names = c("label", "path"))
bins = lapply(1:nrow(binRef), FUN = function(binID) {
	return(import.bed(binRef[binID, path]))
})
names(bins) = binRef$label
for( i in seq_along(bins) ){
	bins[[i]] = keepSeqlevels(bins[[i]], chromosomes,
		pruning.mode = "coarse")
	seqinfo(bins[[i]]) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
}

cat("Reading references...\n")
ref = fread(cmd = sprintf('grep -v "^#" %s', referenceTable),
	header = FALSE,
	col.names = c("cell_line", "marker", "accession", "rep", "path", "src"),
	key = "accession")

bw.files = with(ref, setNames(path, accession))
allFilesExist = file.exists(bw.files)
if ( !all(file.exists(bw.files)) )
	print(bw.files[!allFilesExist])
stopifnot( all(allFilesExist) )

epig = rbindlist(lapply(seq_along(bw.files),
	function(i, bins) {
		x = bw.files[i]
		accession = names(bw.files)[i]
		cat(sprintf("Processing '%s'\n", x))
		ext = getFullExt(x)
		fileList = file.path(outputFolder,
			paste0(accession, ".epig.", names(bins), ".rds"))
		if ( any(!file.exists(fileList)) ) {
			out = NULL
			if ( !is.na(ref[accession, src]) ) {
				if ( file.exists(ref[accession, src]) )
					source(ref[accession, src], local = T)
			} else {
				if ( ext %in% c("bed", "bed.gz") )
					out = processBed(x, bins, chromosomes)
				if ( "bdg" == ext )
					out = processBedGraph(x, bins, chromosomes)
				if ( "bigWig" == ext | "bw" == ext)
					out = processBigWig(x, bins, chromosomes)
				if ( "cod.gz" == ext )
					out = processCod(x, bins, chromosomes, T)
				if ( "cod" == ext )
					out = processCod(x, bins, chromosomes)
				if ( "txt" == ext )
					out = processText(x, bins, chromosomes)
				if ( "xlsx" == ext )
					out = processXlsx(x, bins, chromosomes)
				# If it reaches here, no format recognized
			}
			if ( !is.null(out) ) {
				out$accession = accession

				setkeyv(out, "accession")
				out = out[ref[, list(accession, cell_line, marker, rep)], nomatch = 0]
				setnames(out, "seqnames", "chrom")
				setcolorder(out, c("chrom", "start", "end", "score",
					"cell_line", "marker", "type", "rep", "accession"))
				out[, `:=`(start = as.numeric(start), end = as.numeric(end))]

				dt = split(out, out[, bins])
				for ( k in seq_along(bins) ) {
					dt[[k]][(dt[[k]][, end - start]) != (width(bins[[k]])[1]-1),
						end := end + (width(bins[[k]])[1]-1 - (end - start))]

					dt[[k]] = getTranslocationCoordinates_binned(dt[[k]])

					at = split(dt[[k]], unique(dt[[k]][, accession]))
					for ( j in seq_along(at) ) {
						fName = sprintf("%s.epig.%s.rds", names(at)[j], names(bins)[k])
						if ( !file.exists(fName) )
							saveRDS(at[[j]], file = file.path(outputFolder, fName))
					}
				}
				return(rbindlist(dt))
			} else {
				cat(" Skipped dataset entirely...\n")
			}
		} else {
			cat(" Already processed, reading...\n")
			out = list()
			for ( binLabel in names(bins) ) {
				out[[binLabel]] = readRDS(file.path(outputFolder,
					paste0(names(bw.files)[i], ".epig.", binLabel, ".rds")))
				out[[binLabel]]$bins = binLabel
			}
			return(rbindlist(out))
		}
	}, bins
))

cat("Merging tracks...\n")
dt = split(epig, epig[, bins])
for ( i in seq_along(bins) )
	saveRDS(dt[[i]], file = sprintf("epig.%s.rds", names(bins)[i]))
