require(data.table)

#source("/home/agostif/Projects/GPSeq_analysis/FA/common.functions.R")
source("/media/MiSo/GPSeq/analysis/GG/src/functions.common.R")

#setwd("/media/MiSo/Epigenetic_marks")
setwd("/media/MiSo/Epigenetic_marks")

epig = list(B1000 = readRDS("/media/MiSo/Epigenetic_marks/hist_1Mb.rds"),
			B100 = readRDS("/media/MiSo/Epigenetic_marks/hist_100Kb.rds"))

epig[[1]] = rbindlist(list(epig[[1]], readRDS("/media/MiSo/Epigenetic_marks/HeLa/hela_1Mb.rds")))
epig[[2]] = rbindlist(list(epig[[2]], readRDS("/media/MiSo/Epigenetic_marks/HeLa/hela_100Kb.rds")))

epig[[1]] = rbindlist(list(epig[[1]], readRDS("/media/MiSo/Epigenetic_marks/epig_rds/atac_1Mb.rds")))
epig[[2]] = rbindlist(list(epig[[2]], readRDS("/media/MiSo/Epigenetic_marks/epig_rds/atac_100Kb.rds")))

epig[[1]] = rbindlist(list(epig[[1]], readRDS("/media/MiSo/Epigenetic_marks/epig_rds/meth_1Mb.rds")))
epig[[2]] = rbindlist(list(epig[[2]], readRDS("/media/MiSo/Epigenetic_marks/epig_rds/meth_100Kb.rds")))

epig[[1]] = rbindlist(list(epig[[1]], readRDS("/media/MiSo/Epigenetic_marks/HAP1/hap1_epig_1Mb.rds")))
epig[[2]] = rbindlist(list(epig[[2]], readRDS("/media/MiSo/Epigenetic_marks/HAP1/hap1_epig_100Kb.rds")))

for( i in seq_along(epig) )
	epig[[i]][, chrom := paste0("chr", gsub("24", "Y", gsub("23", "X", chrom)))]

epig = lapply(epig, getTranslocationCoordinates_binned)

epig = lapply(epig,
	function(x){
		x[, `:=`(rep = as.numeric(rep), score = as.numeric(score), cell_line = gsub("Hap1", "HAP1", cell_line))]
		x[accession=="ENCFF485YQE", rep := 2]
		setkeyv(x, c("accession", "chrom", "start", "end"))
		unique(x, by = key(epig))
	})

saveRDS(epig[[1]], file="/media/MiSo/Epigenetic_marks/epig_1Mb.rds")
saveRDS(epig[[2]], file="/media/MiSo/Epigenetic_marks/epig_100kb.rds")
