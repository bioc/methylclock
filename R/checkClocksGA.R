#' Check wheter input data contains the required CpGs for the implemented clocks for Gestational Age.
#' @param x data.frame or tibble (Individual in columns, CpGs in rows, CpG names in first colum - i.e. Horvath's format), ExpressionSet or GenomicRatioSet. A matrix is also possible having the CpG names in the rownames.
#'
#' @details To be supplied
#'
#' @export

checkClocksGA <- function(x,  ...){
  if (inherits(x, "data.frame") & !inherits(x, c("tbl", "tbl_df")))
    cpg.names <- x[,1]
  else if (inherits(x, "matrix"))
    cpg.names <- rownames(x)
  else if (inherits(x, c("tbl", "tbl_df")))
    cpg.names <- pull(MethylationData,1)
  else if (inherits(x, "ExpressionSet"))
    cpg.names <- Biobase::featureNames(x)
  else if (inherits(x, "GenomicRatioSet"))
    cpgs.names <- Biobase::featureNames(x)

  checkKnigth <- coefKnigthGA$CpGmarker[-1][!coefKnigthGA$CpGmarker[-1]%in%cpg.names]
  checkBohlin <- coefBohlinGA$CpGmarker[!coefBohlinGA$CpGmarker%in%cpg.names]
  checkMayne <- coefMayneGA$CpGmarker[-1][!coefMayneGA$CpGmarker[-1]%in%cpg.names]
  checkLee <- coefLeeGA$CpGmarker[-1][!coefLeeGA$CpGmarker[-1]%in%cpg.names]
  
  
  sizes <- c(length(checkKnigth), length(checkBohlin),
             length(checkMayne), length(checkLee))
  
  n <- c(nrow(coefKnigthGA[-1]), nrow(coefBohlinGA),
         nrow(coefMayneGA[-1]), nrow(coefLeeGA[-1]))
  
  df <- data.frame(clock = c("Knigth", "Bohlin", "Mayne", "Lee"),
                   Cpgs_in_clock = n,
                   missing_CpGs = sizes,
                   percentage = round((sizes/n)*100, 1))
  
  print(df)
  
  if (any(sizes!=0)){
    cat("There are some clocks that cannot be computed since your data do not contain the required CpGs. 
        These are the total number of missing CpGs for each clock : \n \n")
    print(df)
    
    out <- list(Knigth=checkKnigth, Bohlin=checkBohlin,
                Mayne=checkMayne, Lee=checkLee)
  }
  else {
    cat("Your data contain the required CpGs for all clocks")
    out <- NULL
  }
  return(out)
}