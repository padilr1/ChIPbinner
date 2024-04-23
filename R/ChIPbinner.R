#' Generate normalized bigwig.
#'
#' @param wd working directory.
#' @param chromSizes chromosome sizes.
#' @param binned_files_path directory for binned files. Files should be in the format $line.$samp.$mark.10kb.bed.
#' @param treatment treatment file name.
#' @param control control file name,
#' @param table table for sequencing depth.
#' @param table_samp_name name in sequencing depth table.
#' @param cell_line cell line.
#' @param histone_mark histone mark.
#' @param blacklist blacklist file.
#' @param addition add.
#' @param raw_count_cutoff raw cutoff.
#' @param use_input whether to include input in calculating normalized bigwigs or not.
#'
#' @return normalized bigwig files.
#' @examples generate_norm_bw(ChIP_files)
#' @export
generate_norm_bw <- function(wd,chromSizes,binned_files_path,treatment,control,table,table_samp_name,cell_line,histone_mark,blacklist,addition,raw_count_cutoff,use_input){
  #need to set wd
  setwd(wd)
  getwd()
  # params
  chrom.sizes=paste0(chromSizes)
  path <- paste0(binned_files_path)
  addition <- as.numeric(paste0(addition))
  raw_count_cutoff <- as.numeric(paste0(raw_count_cutoff))
  # sample info
  cell_line=paste0(cell_line)
  mark=paste0(histone_mark)
  samp=paste0(treatment)
  input=paste0(control) #watch spaces in names
  # RX #
  ########## make sure column name for chip library is "chip" and column name for input library is "inp", as well as "rx_ratio" for rx_ratio #########
  rx=paste0(table)
  rx_samp=paste0(table_samp_name)
  rx <- fread(rx)
  chipLibrary <- rx$CR_MappedReads[rx$samp == rx_samp]
  inputLibrary <- rx$IgG_MappedReads[rx$samp == rx_samp]
  # rxRatio <- rx$rx_ratio[rx$samp == rx_samp]
  # scalingFactor <- rxRatio
  # blacklist
  blacklist=paste0(blacklist)
  bl <- import.bed(blacklist)

  ###################################

  # generate filter for k_final
  s <- list.files(path = path, pattern = '.10kb.',full.names = FALSE,recursive = FALSE) %>%
    tibble(f = .) %>%
    separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
    mutate(f = file.path(path, f))
  # deframe
  d <- deframe(s[,c('samp', 'f')])
  # read in chrom sizes and chr list
  keep <- fread(chrom.sizes) %>%
    setNames(c("chr","seqlength"))
  gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
  # loop through each file and keep only signal for each 10kb bin
  raw <- lapply(d,function(x){
    inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
    out <- dplyr::semi_join(inp,keep,by="chr")
    final <- out$score
  })
  # read in bin
  bs_final <- import.bed("~/Documents/NSD2i/ref/10kb.windows.bed") %>% as.data.frame() %>% dplyr::select(1:3) %>%
    setNames(c('chr','start','end')) %>%
    dplyr::semi_join(keep,by="chr") %>%
    dplyr::filter(chr != "chrM") %>%
    makeGRangesFromDataFrame()
  # get max values
  mxs_final <- bind_cols(raw) %>% apply(1, max)
  # kfinal
  k_final <- mxs_final > raw_count_cutoff & !overlapsAny(bs_final,bl)
  # k_final <- !overlapsAny(bs_final,bl)

  ####################################

  # load exclusion factor
  # load("")
  # list files
  s <- list.files(path = path, pattern = '.10kb.',full.names = FALSE,recursive = FALSE) %>%
    tibble(f = .) %>%
    separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
    mutate(f = file.path(path, f)) %>%
    dplyr::filter(line==cell_line) %>%
    dplyr::filter(samp==samp | samp==input)
  d <- deframe(s[,c('samp', 'f')])
  # read in chrom sizes and chr list
  keep <- fread(chrom.sizes) %>%
    setNames(c("chr","seqlength"))
  gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
  # loop through each file
  raw <- lapply(d,function(x){
    inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
    out <- dplyr::semi_join(inp,keep,by="chr")
    out <- out[k_final]
  })
  pre_bl <- lapply(d,function(x){
    inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
    imd <- dplyr::semi_join(inp,keep,by="chr")
    out <- imd[,1:3] %>% mutate(start = start + 1)
    out <- out[k_final]
    out <- out %>%
      makeGRangesFromDataFrame(seqinfo = gn)
  })
  bw <- lapply(d,function(x){
    inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
    imd <- dplyr::semi_join(inp,keep,by="chr")
    out <- imd[,1:3] %>% mutate(start = start + 1)
    out <- out[k_final]
    out <- out %>%
      makeGRangesFromDataFrame(seqinfo = gn)
    k <- !overlapsAny(out, bl)
    out <- out[!overlapsAny(out, bl)]
  })
  # normalize by input and scale by a factor
  if (use_input == TRUE){
    bw[[samp]]$score <- (log2(((raw[[samp]]$score/chipLibrary) + addition))*1e6)[!overlapsAny(pre_bl[[samp]], bl)]
  } else {
    bw[[samp]]$score <- (log2(((raw[[samp]]$score/chipLibrary) + addition))*1e6)
  }
  # bw[[samp]]$score <- (log2(((raw[[samp]]$score/chipLibrary) + addition))*1e6)[!overlapsAny(pre_bl[[samp]], bl)]
  # bw[[samp]]$score <- (log2(((raw[[samp]]$score/chipLibrary) + addition))*1e6)
  # drop NA values
  bw[[samp]] <- bw[[samp]][!is.na(bw[[samp]]$score)]
  # Create the 'norm.bw' directory if it doesn't exist
  if (!file.exists("data/norm.bw")) {
    dir.create("data/norm.bw")
  }
  # export bw
  return(export.bw(bw[[samp]],con=sprintf("data/norm.bw/%s.%s.%s.10kb.norm.bw",cell_line,samp,mark)))
}

#' Pre-clustering.
#'
#' @param wd working directory.
#' @param path_for_norm_bw path for normalized bigiwg.
#' @param typeof_file norm bigwig.
#' @param treatment treatment sample.
#' @param control control sample.
#' @param cell_line cell line.
#' @param histone_mark histone mark.
#'
#' @return comparison matrices.
#' @export
#'
#' @examples treatment vs control sample.
pre_clustering <- function(wd,path_for_norm_bw,typeof_file,treatment,control,cell_line,histone_mark){
  # parameters
  setwd(wd)
  getwd()
  path <- paste0(path_for_norm_bw)
  pattern = paste0(typeof_file)
  # gene <- import.bed('ensembl/gene.bed')
  # igr <- import.bed('ensembl/intergenic.bed')
  control=paste0(control)
  test=paste0(treatment)
  cell_line=paste0(cell_line)
  mark=paste0(histone_mark)
  #
  s <- list.files(path = path, pattern = pattern,full.names = FALSE,recursive = FALSE) %>%
    tibble(f = .) %>%
    separate(f,c("line","samp","mark",NA,NA,NA),'\\.',F) %>%
    mutate(f = file.path(path, f)) %>%
    dplyr::filter(line == cell_line) %>%
    dplyr::filter(mark == mark) %>%
    dplyr::filter(samp == control | samp == test)
  # filter(cond != 'MT') %>%
  # filter(samp != "Cal27_KO17") %>%
  # filter(samp != "Cal27_NSD2OE")
  #
  #
  # NEED TO ORDER SAMPLE FIRST THEN PARENTAL #
  odr <- c(test,control)
  s <- s %>%
    dplyr::slice(match(odr,samp))
  #
  d <- deframe(s[,c('samp', 'f')]) %>%
    lapply(import.bw) #original ranges for KO = 283852
  #
  lapply(d,length)
  #
  #
  d[[control]] <- subsetByOverlaps(d[[control]],d[[test]])
  d[[test]]<- subsetByOverlaps(d[[test]],d[[control]])
  #
  lapply(d,length)
  #
  #
  #
  r <- lapply(d, function(y) y[y$score != 0]) %>%
    Reduce(function(a, b) a[overlapsAny(a, b)], .) %>%
    granges()
  #
  #
  cell_line <- as.character(s$line[[1]])
  mark <- as.character(s$mark[[1]])
  #
  # Create the 'mat.csv' directory if it doesn't exist
  if (!file.exists("data/mat.csv")) {
    dir.create("data/mat.csv")
  }
  # Create the 'pooled.bed' directory if it doesn't exist
  if (!file.exists("data/pooled.bed")) {
    dir.create("data/pooled.bed")
  }
  # Create the figs directory if it doesn't exist
  if (!file.exists("figs/")) {
    dir.create("figs/")
  }
  # split
  split(d, s$line) %>%
    lapply(function(x) {
      o <- lapply(x, function(y) {
        findOverlaps(r, y) %>%
          to() %>%
          {y[.]} %>%
          score()
      }) %>%
        bind_cols() %>%
        `names<-`(c('x', 'y'))
      ok <- o$x > quantile(o$x, .01) &
        o$x < quantile(o$x, .99) &
        o$y > quantile(o$y, .01) &
        o$y < quantile(o$y, .99)
      write_csv(o[ok,], sprintf('data/mat.csv/%s.%s.%s.%s.mat.csv',cell_line,control,test,mark), col_names = F)
      export.bed(r[ok], sprintf('data/pooled.bed/%s.%s.%s.%s.pooled.bed',cell_line,control,test,mark))
    })

  lfc <- fread(sprintf('data/mat.csv/%s.%s.%s.%s.mat.csv',cell_line,control,test,mark))
  # 1. Open jpeg file
  jpeg(filename = sprintf('figs/%s.%s.%s.%s.smoothScatter.jpeg',cell_line,control,test,mark))
  # 2. Create the plot
  smoothScatter(y = lfc$V1,x=lfc$V2,xlab = control,ylab=test)
  # 3. Close the file
  dev.off()
  return(list(s=s,mark=mark))
}
