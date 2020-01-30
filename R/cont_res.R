
#' results concatenation
#'
#' @param pattern file prefix of the results files
#' @param res.fld folder that contains the result files
#' @param filter.EPIC do filtering for EPIC chip data?
#' @description This function allows you to concatenate the results generated from mcLME()
#' @keywords LME
#' @export
#' @examples
#' cont_res('PLEvsCLE_','results')

cont_res = function(pattern, res.fld, filter.EPIC=F){
  require(data.table)
  require(dplyr)

  files = list.files(res.fld,full=TRUE)
  files = files[grep(paste0(res.fld,'/',pattern), files)]
  res = lapply(files,fread) %>% rbindlist()


  if(filter.EPIC){
    data(cross_reactive);data(polymorphic);data(baseext)
    cross_reactive.probe=cross_reactive$V1
    polymorphic.probe=polymorphic$PROBE
    baseext.probe=baseext$PROBE
    res = res %>%
      mutate(cross.reactive.hit = ifelse(Probe %in% cross_reactive.probe, TRUE, FALSE)) %>%
      mutate(polymorphic.hit = ifelse(Probe %in% polymorphic.probe, TRUE, FALSE)) %>%
      mutate(baseext.hit = ifelse(Probe %in% baseext.probe, TRUE, FALSE)) %>%
      dplyr::filter(!cross.reactive.hit & !polymorphic.hit & !baseext.hit)

    res = res %>% mutate(FDR = p.adjust(P, method = 'BH')) %>%
      dplyr::select(Probe, Beta, SE, stat, P, FDR, everything())

    rm(cross_reactive);rm(polymorphic);rm(baseext);
    rm(cross_reactive.probe);rm(polymorphic.probe);rm(baseext.probe);
  }else{
    res = res %>% mutate(FDR = p.adjust(P, method = 'BH')) %>%
      dplyr::select(Probe, Beta, SE, stat, P, FDR)
  }

  gc();gc()

  file.remove(files)
  return(res %>% arrange(P))
}

