# Correct coordinate -----------------------------------------------------------

site2coords = function(site) {
    # Extract information from site label (e.g., chrN:MMM...).
    
    stopifnot(grepl("^chr[0-9XY]+:[0-9]+$", site))

    coords = unlist(strsplit(site, ":", fixed = T))
    coords = data.table(data.frame(
        chrom = coords[1],
        start = as.numeric(coords[2]),
        stringsAsFactors = F
    ))

    return(coords)
}

rmBinsOverlappingSite = function(coords, site, value = NA) {
    # Set coordinates of bins overlapping a site to value.
    
    stopifnot(all(c("chrom", "start", "end") %in% colnames(data)))
    
    siteCoords = site2coords(site)

    overlapping = levels(coords$chrom)[coords$chrom] == siteCoords$chrom
    overlapping = overlapping & coords$start <= siteCoords$start
    overlapping = overlapping & coords$end >= siteCoords$start

    if ( 0 == sum(overlapping, na.rm = T) ) return(coords)

    coords[overlapping, chrom := value]
    coords[overlapping, start := value]
    coords[overlapping, end := value]

    return(coords)
}

mergeChroms = function(chr1, chr2) {
    # Assemble chromosome labels for merged chromosomes (e.g., chrA:B)
    paste0(chr1, ":", substr(chr2, 4, nchar(chr2)))
}

correctTranslocation_binned = function(coords, bin_step, siteA, siteB,
    oneside = F) {
    # Correct coordinates a chromosome translocation. Use oneside to do it only
    # on one side of the translocation site, i.e., generates only the first part
    # of chrA:B and the last of chrB:A.
    
    siteDF1 = site2coords(siteA)
    siteDF2 = site2coords(siteB)

    transloced = levels(coords$chrom)[coords$chrom] == siteDF1$chrom
    transloced = transloced & coords$start >= siteDF1$start
    transloced[is.na(transloced)] = F
    coords[transloced, chrom := mergeChroms(siteDF2$chrom, siteDF1$chrom)]
    oldStart = siteDF1$start + bin_step - siteDF1$start %% bin_step
    newStart = siteDF2$start + bin_step - siteDF2$start %% bin_step
    coords$start[transloced] = coords$start[transloced] - oldStart + newStart
    coords$end[transloced] = coords$end[transloced] - oldStart + newStart

    transloced = levels(coords$chrom)[coords$chrom] == siteDF1$chrom
    transloced = transloced & coords$start <= siteDF1$start
    transloced[is.na(transloced)] = F
    coords[transloced, chrom := mergeChroms(siteDF1$chrom, siteDF2$chrom)]

    if ( !oneside ) coords = correctTranslocation_binned(
        coords, bin_step, siteB, siteA, oneside = T)

    return(coords)
}

getTranslocationCoordinates_binned = function(data,
    site1 = "chr9:133681295", site2 = "chr22:23632359",
    removeOverlapping = F, oneIndexed = T, suffix = "phil") {
    # Correct coordinates for chromosomes arm translocation.
    # 
    # Provide a table with columns chrom, start and end.
    # Translocation site should be in chrNN:MMM... format.
    # Default translocation site is for HAP1 Philadelphia chromosome.
    # 
    # Use removeOverlapping to remove (set to NA) bins overlapping the site.
    # Set oneIndexed for classical bed format (F) or GRanges (T).
    # 
    # Returns a table with chrom, start and end coordinates corrected.
    
    stopifnot(all(grepl("^chr[0-9XY]+:[0-9]+$", c(site1, site2))))
    stopifnot(all(c("chrom", "start", "end") %in% colnames(data)))

    bin_size = unique(data$end - data$start)
    if ( oneIndexed ) bin_size = bin_size + 1
    stopifnot(1 == length(bin_size))

    bin_step = unique(diff(data$end))
    bin_step = bin_step[bin_step >= 0]
    stopifnot(1 == length(bin_step))

    newCoords = data[, .(chrom, start, end)]
    if ( removeOverlapping ) {
        newCoords = rmBinsOverlappingSite(newCoords, site1)
        newCoords = rmBinsOverlappingSite(newCoords, site2)
    }

    newCoords = correctTranslocation_binned(newCoords, bin_step, site1, site2)

    setnames(newCoords, paste0(suffix, c("Chrom", "Start", "End")))
    data = cbind(data, newCoords)

    return(data)
}
