#' RLM in parallel
#'
#' The features that you are working with will be called 'FeatureVal' within the function.
#' Make sure you use 'FeatureVal' in 'form' (model formula) and 'var.interest' (coefficient to extract, if applicable).
#'
#' @param pDat phenotype table (in the long format if it is longitudinal)
#' @param form model formula for rlm
#' @param var.interest coefficient to extract (append the non-reference level after in the case of a factor, eg. SexFemale)
#' @param chunk.list file names of your chunks
#' @param chunk.prefix folder that contains your chunks, make sure it ends with a '/'
#' @param output.prefix folder name and file name prefix of your result files
#' @param num.chunks number of chunks you have
#' @param mc.core cores that you want to use
#' @param Feature_Name name of the column that contains the feature names (in each chunk)
#' @param Sample_Name name of the column that contains the sample ID (in the phenotype table)
#' @description This function allows you to run RLM in parallel for high dimensional data (methylation, SNP, expression...)
#' @keywords RLM
#' @export
#' @examples
#' chunk.list=paste0('chunk_',1:32,'.csv')
#'
#' mcLM(pDat, form = 'Exprs ~ FeatureVal*Group',
#' var.interest = 'FeatureVal:GroupCLE', chunk.list=chunk.list,
#' chunk.prefix = 'processed_data/Beta_chunks/',
#' output.prefix='tmpres/PLECLEInt_',
#' num.chunks=32, mc.core=32, Feature_Name='cpg_names', Sample_Name='Sample_Name')


mcRLM = function(pDat, form, var.interest, chunk.list, chunk.prefix, output.prefix, num.chunks=32, mc.core=32, Feature_Name=NULL, Sample_Name=NULL){
  require(nlme)
  require(doMC)
  require(data.table)
  require(dplyr)
  registerDoMC(mc.core)

  foreach(i = 1:num.chunks)  %dopar% {
    pDat <- pDat %>% setnames(old=Sample_Name, new='Sample_Name') %>% as.data.table()

    dat <- fread(paste0(chunk.prefix, chunk.list[i]))
    dat <- dat %>% setnames(old=Feature_Name, new='Feature_Name') %>% as.data.table()

    RLM_res <- data.table(Probe=dat[,Feature_Name],Beta=as.numeric(NA),
                          SE=as.numeric(NA),stat=as.numeric(NA), P=as.numeric(NA))

    for(k in 1:nrow(RLM_res)){
      if(k%%500==0) print(k)
      pDat[,FeatureVal:=as.numeric(dat[k,pDat[,Sample_Name],with=FALSE])]
      form = as.formula(form)

      mod <- rlm(as.formula(form), maxit=200, data=pDat)

      cf <- try({
        result = coeftest(mod, vcov=vcovHC(mod, type='HC0'))
      },silent=T)

      if(class(cf)!="try-error"){
        RLM_res[k, c('Beta','SE','stat','P'):= as.list(cf[var.interest,c(1,2,3,4)])]
      }


      pDat[,FeatureVal:=NULL]
      rm(mod);rm(cf)
    }# end for

    fwrite(RLM_res,paste0(output.prefix,chunk.list[i]))

    cat(paste0('\n\n\n ###########  done for testing chunk ',i,'\n'))
    rm(RLM_res)
    rm(dat)
    gc();gc()
  } # end foreach
} # end function
