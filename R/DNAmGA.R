#' Gestational DNAm age estimation using different DNA methylation clocks.
#' @param x data.frame (Individual in columns, CpGs in rows, CpG names in 
#' first colum - i.e. Horvath's format), matrix (individuals in columns 
#' and Cpgs in rows having CpG names in the rownames), ExpressionSet or 
#' GenomicRatioSet.
#' @param toBetas Should data be transformed to beta values? Default is FALSE. 
#' If TRUE, it implies data are M values.
#' @param fastImp Is fast imputation performed if necessary? (see details). 
#' Default is FALSE
#' @param normalize Is Horvath's normalization performed? By default is FALSE
#' @param age individual's chronological age. Required to compute gestational 
#' age difference output
#' @param cell.count Are cell counts estimated? Default is TRUE.
#' @param cell.count.reference Used when 'cell.count' is TRUE. Default is 
#' "blood gse35069 complete". See 'meffil::meffil.list.cell.count.references()' 
#' for possible values.
#' @param min.perc Indicates the minimum conicidence percentage required 
#' between CpGs in or dataframee x and CpGs in clock coefficients to perform 
#' the calculation. If min.prec is too low, the estimated gestational DNAm 
#' age can be poor
#' @param ... Other arguments to be passed through impute package
#'
#' @details Imputation is performed when having missing data. 
#' Fast imputation is performed by ... 
#' what about imputing only when CpGs for the clock are missing?
#'
#' @examples
#' TestDataset <- get_TestDataset()
#' TestDataset[1:5, ]
#' ga.test <- DNAmGA(TestDataset)
#' @return the estimated gestational DNAm age
#'
#' @import impute dplyr tidyverse tibble
#' @importFrom minfi getBeta
#' @importFrom Biobase featureNames exprs
#' @importFrom minfi getBeta
#'
#' @export

DNAmGA <- function(x, toBetas = FALSE,
                    fastImp = FALSE,
                    normalize = FALSE,
                    age,
                    cell.count = TRUE,
                    cell.count.reference = "andrews and bakulski cord blood",
                    min.perc = 0.8,
                    ...) {
    #..# cpgs <- getInputCpgValues(x, toBetas)

    if (inherits(x, "data.frame")) {
        cpgs.names <- as.character(x[, 1, drop = TRUE])
        if (length(grep("cg", cpgs.names)) == 0) {
            stop("First column should contain CpG names")
        }
        cpgs <- t(as.matrix(x[, -1]))
        colnames(cpgs) <- cpgs.names
    }
    else if (inherits(x, "matrix")) {
        cpgs <- t(x)
    }
    else if (inherits(x, "ExpressionSet")) {
        cpgs <- t(Biobase::exprs(x))
    }
    else if (inherits(x, "GenomicRatioSet")) {
        cpgs <- t(minfi::getBeta(x))
    }
    else {
        stop("x must be a data.frame, matrix, 'GenomicRatioSet' 
             or an 'ExpressionSet' object")
    }
    
    if (!all(c("coefKnightGA", "coefBohlin", "coefMayneGA", 
               "coefLeeGA", "coefEPIC") %in% ls(.GlobalEnv))) {
        load_DNAmGA_Clocks_data()
    }
    
    if (toBetas) {
        toBeta <- function(m) {
            2^m / (2^m + 1)
        }
        cpgs <- toBeta(cpgs)
    }
    
    if (any(cpgs < -0.1 | cpgs > 1.1, na.rm = TRUE)) {
        stop("Data seems to do not be beta values. Check your data or set 'toBetas=TRUE'")
    }

    cpgs.all <- c( as.character(coefKnightGA$CpGmarker[-1]), 
                   as.character(coefBohlin$CpGmarker[-1]),
                   as.character(coefMayneGA$CpGmarker[-1]),
                   as.character(coefLeeGA$CpGmarker[-1]),
                   as.character(coefEPIC$CpGmarker[-1]))

    if (any(!cpgs.all %in% colnames(cpgs))) {
        warning("CpGs in all Gestational Age clocks are not present in your
        data. Try 'checkClocksGA' function to find the missing CpGs of
                each method.")
    }

    cpgs.in <- intersect(cpgs.all, colnames(cpgs))

    miss <- apply(cpgs[, cpgs.in], 2, function(x) any(is.na(x)))

    if (any(miss)) {
        if (fastImp) {
            message("Imputing missing data of", sum(miss), "CpGs .... \n")
            mm <- apply(cpgs[, cpgs.in], 2, median, na.rm = TRUE)
            cpgs.imp <- sweep(cpgs[, cpgs.in], 2, STATS = mm, 
                                FUN = function(x, s) ifelse(is.na(x), s, x) )
        }
        else {
            quiet <- function(x) {
                sink(tempfile())
                on.exit(sink())
                invisible(force(x))
            }
            message("Imputing missing data of the entire matrix .... \n")
            cpgs.imp <- quiet(t(impute.knn(t(cpgs), ...)$data))
        }
        message("Data imputed. Starting DNAm clock estimation ... \n")
    }
    else {
        cpgs.imp <- cpgs
    }

    # --------------> Knight
    Knight <- predAge(cpgs.imp, coefKnightGA, intercept = TRUE, min.perc)
    Knight <- data.frame( id = rownames(cpgs.imp), Knight = Knight )

    # --------------> Bohlin
    bohlin <- predAge(cpgs.imp, coefBohlin, intercept = TRUE, min.perc)
    Bohlin <- data.frame( id = rownames(cpgs.imp), Bohlin = bohlin / 7 )

    # --------------> Mayne
    mayne <- predAge(cpgs.imp, coefMayneGA, intercept = TRUE, min.perc)
    Mayne <- data.frame( id = rownames(cpgs.imp), Mayne = mayne )
    
    # --------------> EPIC
    epic <- predAge(cpgs.imp, coefEPIC, intercept = TRUE, min.perc)
    EPIC <- data.frame( id = rownames(cpgs.imp), EPIC = epic )

    # --------------> Lee
    if (  sum(coefLeeGA$CpGmarker[-1] %in% colnames(cpgs.imp)) / 
          length(coefLeeGA$CpGmarker) >= min.perc) {
        coefLeeSel <- data.frame(CpGmarker = coefLeeGA$CpGmarker, 
                            CoefficientTraining = coefLeeGA$Coefficient_RPC)
        Lee.RPC <- predAge(cpgs.imp, coefLeeSel, intercept = TRUE, min.perc)
        
        coefLeeSel <- data.frame(CpGmarker = coefLeeGA$CpGmarker, 
                            CoefficientTraining = coefLeeGA$Coefficient_CPC)
        Lee.CPC <- predAge(cpgs.imp, coefLeeSel, intercept = TRUE, min.perc)
        
        coefLeeSel <- data.frame(CpGmarker = coefLeeGA$CpGmarker, 
                        CoefficientTraining = coefLeeGA$Coefficient_refined_RPC)
        Lee.refRPC <- predAge(cpgs.imp, coefLeeSel, intercept = TRUE, min.perc)
        
        if (!missing(age)) {
            Lee <- data.frame(
                id = rownames(cpgs.imp),
                Lee.RPC = Lee.RPC,
                Lee.CPC = Lee.CPC,
                Lee.refRPC = Lee.refRPC,
                Lee.RPCdiff = Lee.RPC - age,
                Lee.CPCdiff = Lee.CPC - age,
                Lee.refRPCdiff = Lee.refRPC - age
            )
        } else {
            Lee <- data.frame(
                id = rownames(cpgs.imp),
                Lee.RPC = Lee.RPC,
                Lee.CPC = Lee.CPC,
                Lee.refRPC = Lee.refRPC
            )
        }
    }
    else {
        warning("The number of missing CpGs for Lee clocks exceeds 80%.
                \n  ---> This DNAm clock will be NA.")
        Lee <- data.frame(id=rownames(cpgs.imp), Lee = rep(NA, nrow(cpgs.imp)))
    }

    # --------------> output

    if (!missing(age)) {
        if (!cell.count) {
            Knight <- ageAcc1(Knight, age, lab = "Knight")
            Bohlin <- ageAcc1(Bohlin, age, lab = "Bohlin")
            Mayne <- ageAcc1(Mayne, age, lab = "Mayne")
            EPIC <- ageAcc1(EPIC, age, lab = "EPIC")
            if( dim(Lee)[2]>2){  #. 04/01/2022.#
                Lee.RPC <- ageAcc1(Lee[,c("id", "Lee.RPC")], 
                                   age, lab = "Lee.RPC") #. 04/01/2022.#
                Lee.CPC <- ageAcc1(Lee[,c("id", "Lee.CPC")], 
                                   age, lab = "Lee.CPC") #. 04/01/2022.#
                Lee.refRPC <- ageAcc1(Lee[,c("id", "Lee.refRPC")], 
                                      age, lab = "Lee.refRPC") #. 04/01/2022.#
            }
        }
        else {
            if( cell.count.reference != 'placenta') {
                # Cell Counts using meffil estimate cell counts
                cell.counts <- try(meffil.estimate.cell.counts.from.betas(
                    t(cpgs),
                    cell.count.reference
                ), TRUE)
                if (inherits(cell.counts, "try-error")) {
                    warning("cell countscannot be estimated since your data 
                            have missing CpGs for meffil")
                } else {
                    ok <- which(apply(cell.counts, 2, IQR) > 10e-6)
                    cell.counts <- cell.counts[, ok]
                    df <- data.frame(age = age, cell.counts)
                    
                    Knight <- ageAcc2(Knight, df, lab = "Knight")
                    Bohlin <- ageAcc2(Bohlin, df, lab = "Bohlin")
                    Mayne <- ageAcc2(Mayne, df, lab = "Mayne")
                    EPIC <- ageAcc2(EPIC, df, lab = "EPIC")
                    if(exists( "Lee.RPC" )){  
                        Lee.RPC <- ageAcc2(Lee[,c("id", "Lee.RPC")], 
                                           df, lab = "Lee.RPC") }
                    if(exists( "Lee.CPC" )){  
                        Lee.CPC <- ageAcc2(Lee[,c("id", "Lee.CPC")], 
                                           df, lab = "Lee.CPC") }
                    if(exists( "Lee.refRPC" )) {
                        Lee.refRPC <- ageAcc2(Lee[,c("id", "Lee.refRPC")], 
                                              df, lab = "Lee.refRPC") }
                }
                
            } else {
                
                if (!requireNamespace("planet", quietly = TRUE))
                    install.packages("planet")
                
                data("plCellCpGsThird")
                
                # Apply planet package to deconvolute data 
                cell.counts <- try(
                    
                    minfi:::projectCellType(
                        t(cpgs)[which( rownames(t(cpgs)) %in% 
                                           rownames(plCellCpGsThird)), ],
                        plCellCpGsThird,
                        lessThanOne = FALSE
                    )
                    , TRUE)
                if (inherits(cell.counts, "try-error")) {
                    warning("cell counts cannot be estimated since your data 
                            have missing CpGs for meffil")
                } else {
                    ok <- which(apply(cell.counts, 2, IQR) > 10e-6)
                    cell.counts <- cell.counts[, ok]
                    df <- data.frame(age = age, cell.counts)
                    
                    Knight <- ageAcc2(Knight, df, lab = "Knight")
                    Bohlin <- ageAcc2(Bohlin, df, lab = "Bohlin")
                    Mayne <- ageAcc2(Mayne, df, lab = "Mayne")
                    EPIC <- ageAcc2(EPIC, df, lab = "EPIC")
                    if(exists( "Lee.RPC" )) { 
                        Lee.RPC <- ageAcc2(Lee[,c("id", "Lee.RPC")], 
                                           df, lab = "Lee.RPC") }
                    if(exists( "Lee.CPC" )) { 
                        Lee.CPC <- ageAcc2(Lee[,c("id", "Lee.CPC")], 
                                           df, lab = "Lee.CPC") }
                    if(exists( "Lee.refRPC" )) { 
                        Lee.refRPC <- ageAcc2(Lee[,c("id", "Lee.refRPC")], 
                                              df, lab = "Lee.refRPC") }
                }
            }
        }
    } else {
        cell.count <- FALSE
    }
    
    out <- Knight %>%
        full_join(Bohlin, by = "id") %>%
        full_join(Mayne, by = "id") %>%
        full_join(EPIC, by = "id")
    
    if( exists( "Lee.RPC" )) {
        
        if (inherits(Lee.RPC, "matrix")) {
            Lee.RPC <- data.frame( cbind(id = rownames(Lee.RPC), 
                                         Lee.RPC = Lee.RPC[,1]))
        }
        if (inherits(Lee.CPC, "matrix")) {
            Lee.CPC <- data.frame( cbind(id = rownames(Lee.CPC), 
                                         Lee.CPC = Lee.CPC[,1]))
        }
        if (inherits(Lee.refRPC, "matrix")) {
            Lee.refRPC <- data.frame( cbind(id = rownames(Lee.refRPC), 
                                            Lee.refRPC = Lee.refRPC[,1]))
        }
        
        out <- out %>%
            full_join(Lee.RPC, by = "id") %>%
            full_join(Lee.CPC, by = "id") %>%
            full_join(Lee.refRPC, by = "id")
        
    } else {
        out <- out %>%
            full_join(Lee, by = "id")
    }
    
    out <- tibble::as_tibble(out)
    
    if (!missing(age)) { out <- add_column(out, age = age) }
    if (cell.count) { attr(out, "cell_proportion") <- cell.counts }
    
    out
}
