PROBINR_TXT_LITERAL = """#' profile_binarization suite
#'
#' parallel implementation of compute_criteria()
#' implemented by Gustavo Magaña López, as part
#' of his M1 BIBS course at Université Paris-Saclay.
#'
#' The rest of the functions found within this script are either
#' verbatim copies or adaptations of the work performed by
#' the following researchers :
#'
#' Beal, Jonas
#' Montagud, Arnau
#' Traynard, Pauline
#' Barillot, Emmanuel
#' Calzone, Laurence
#'
#' at the Computational Systems Biology of Cancer group at Institut Curie
#' website : https://sysbio.curie.fr/
#'   email :  contact-sysbio@curie.fr
#'
#' This is the repository containing the original implementation in
#' Rmarkdown notebooks :
#'
#' https://github.com/sysbio-curie/PROFILE
#'
#' I leave the following BibTex citations referencing their work :
#'
#' @article{Beal2019,
#' abstract = {Logical models of cancer pathways are typically built by mining the literature for relevant experimental observations. They are usually generic as they apply for large cohorts of individuals. As a consequence, they generally do not capture the heterogeneity of patient tumors and their therapeutic responses. We present here a novel framework, referred to as PROFILE, to tailor logical models to a particular biological sample such as a patient tumor. This methodology permits to compare the model simulations to individual clinical data, i.e., survival time. Our approach focuses on integrating mutation data, copy number alterations (CNA), and expression data (transcriptomics or proteomics) to logical models. These data need first to be either binarized or set between 0 and 1, and can then be incorporated in the logical model by modifying the activity of the node, the initial conditions or the state transition rates. The use of MaBoSS, a tool based on Monte-Carlo kinetic algorithm to perform stochastic simulations on logical models results in model state probabilities, and allows for a semi-quantitative study of the model phenotypes and perturbations. As a proof of concept, we use a published generic model of cancer signaling pathways and molecular data from METABRIC breast cancer patients. For this example, we test several combinations of data incorporation and discuss that, with these data, the most comprehensive patient-specific cancer models are obtained by modifying the nodes' activity of the model with mutations, in combination or not with CNA data, and altering the transition rates with RNA expression. We conclude that these model simulations show good correlation with clinical data such as patients' Nottingham prognostic index (NPI) subgrouping and survival time. We observe that two highly relevant cancer phenotypes derived from personalized models, Proliferation and Apoptosis, are biologically consistent prognostic factors: patients with both high proliferation and low apoptosis have the worst survival rate, and conversely. Our approach aims to combine the mechanistic insights of logical modeling with multi-omics data integration to provide patient-relevant models. This work leads to the use of logical modeling for precision medicine and will eventually facilitate the choice of patient-specific drug treatments by physicians.},
#' author = {Beal, Jonas and Montagud, Arnau and Traynard, Pauline and Barillot, Emmanuel and Calzone, Laurence},
#' doi = {10.3389/fphys.2018.01965},
#' file = {:home/gml/Documents/Master/M1_BIBS/S2/TER/biblio/personal-models-omics.pdf:pdf},
#' issn = {1664042X},
#' journal = {Frontiers in Physiology},
#' keywords = {Breast cancer,Data discretization,Logical models,Personalized mechanistic models,Personalized medicine,Stochastic simulations},
#' mendeley-groups = {TER : BNeDiction},
#' number = {JAN},
#' title = {{Personalization of logical models with multi-omics data allows clinical stratification of patients}},
#' volume = {10},
#' year = {2019}
#' }
#' @article{Beal2019a,
#' abstract = {Logical models of cancer pathways are typically built by mining the literature for relevant experimental observations. They are usually generic as they apply for large cohorts of individuals. As a consequence, they generally do not capture the heterogeneity of patient tumors and their therapeutic responses. We present here a novel framework, referred to as PROFILE, to tailor logical models to a particular biological sample such as a patient tumor. This methodology permits to compare the model simulations to individual clinical data, i.e., survival time. Our approach focuses on integrating mutation data, copy number alterations (CNA), and expression data (transcriptomics or proteomics) to logical models. These data need first to be either binarized or set between 0 and 1, and can then be incorporated in the logical model by modifying the activity of the node, the initial conditions or the state transition rates. The use of MaBoSS, a tool based on Monte-Carlo kinetic algorithm to perform stochastic simulations on logical models results in model state probabilities, and allows for a semi-quantitative study of the model phenotypes and perturbations. As a proof of concept, we use a published generic model of cancer signaling pathways and molecular data from METABRIC breast cancer patients. For this example, we test several combinations of data incorporation and discuss that, with these data, the most comprehensive patient-specific cancer models are obtained by modifying the nodes' activity of the model with mutations, in combination or not with CNA data, and altering the transition rates with RNA expression. We conclude that these model simulations show good correlation with clinical data such as patients' Nottingham prognostic index (NPI) subgrouping and survival time. We observe that two highly relevant cancer phenotypes derived from personalized models, Proliferation and Apoptosis, are biologically consistent prognostic factors: patients with both high proliferation and low apoptosis have the worst survival rate, and conversely. Our approach aims to combine the mechanistic insights of logical modeling with multi-omics data integration to provide patient-relevant models. This work leads to the use of logical modeling for precision medicine and will eventually facilitate the choice of patient-specific drug treatments by physicians.},
#' author = {Beal, Jonas and Montagud, Arnau and Traynard, Pauline and Barillot, Emmanuel and Calzone, Laurence},
#' doi = {10.3389/fphys.2018.01965},
#' file = {:home/gml/Documents/Master/M1_BIBS/S2/TER/biblio/additional-personal-models-omics.pdf:pdf},
#' issn = {1664042X},
#' journal = {Frontiers in Physiology},
#' keywords = {Breast cancer,Data discretization,Logical models,Personalized mechanistic models,Personalized medicine,Stochastic simulations},
#' mendeley-groups = {TER : BNeDiction},
#' number = {JAN},
#' pages = {1--23},
#' title = {{Personalization of logical models with multi-omics data allows clinical stratification of patients}},
#' volume = {10},
#' year = {2019}
#' }

# pacman resolves most dependencies in an easy to use way it should be replaced
# on the final release as I consider it to be a bad practice (installs to user
# library instead of a renv environment).
if (!require("pacman")) install.packages("pacman")

# Declare CRAN dependencies :
r_dependencies <- c("mclust", "diptest", "moments", "magrittr", "tidyr", "dplyr", 
    "tibble", "bigmemory", "doSNOW", "foreach", "parallel", "glue", "docstring")

# load dependencies via pacman
pacman::p_load(r_dependencies, character.only = TRUE)

tbl_to_df <- function(x, column) {
    #' Cast dplyr::tible to base::data.frame
    x <- x %>%
        as.data.frame() %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames(var = column)
    x
}

random_words <- function(n = 1, word_length = 5) {
    #' Create n random words of length 'word_length'
    a <- do.call(paste0, replicate(word_length, sample(LETTERS, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

tbl_transpose <- function(df) {
    df %>%
        tibble::rownames_to_column() %>%
        tidyr::pivot_longer(-rowname) %>%
        tidyr::pivot_wider(names_from = rowname, values_from = value)
}

chunks_of_n <- function(x, n) {
    #' Split a vector x in ~m pieces of size n
    split(x, ceiling(seq_along(x)/n))
}

split_in_n <- function(x, n) {
    #' Split a vector in n pieces
    suppressWarnings(split(x, 1:n))
}

BI <- function(dataset) {
    #' BI : Bimodality Index
    #'
    #' Function to compute the Bimodality Index described in Wang et al. (2009)
    x <- dataset
    mc <- mclust::Mclust(na.omit(x), G = 2, modelNames = "E", verbose = FALSE)
    if (is.null(mc)) {
        bi <- NA
    } else {
        sigma <- sqrt(mc$parameters$variance$sigmasq)
        delta <- abs(diff(mc$parameters$mean))/sigma
        pi <- mc$parameters$pro[1]
        bi <- delta * sqrt(pi * (1 - pi))
    }
    bi
}

OSclass <- function(exp_dataset, ref_dataset = exp_dataset) {
    #' Function to binarise the tails of the distribution
    #' Based on inter-quartile range (IQR)
    #' similar to methods described in teh outlier-sum statistic
    #' (Tibshirani and Hastie, 2007).
    #' Can be called with a reference dataset

    classif <- rep(NA, length(exp_dataset))
    q25 <- quantile(ref_dataset, 0.25, na.rm = T)
    q75 <- quantile(ref_dataset, 0.75, na.rm = T)
    IQR <- q75 - q25  # Inter-Quartile Range

    classif[exp_dataset > IQR + q75] <- 1
    classif[exp_dataset < q25 - IQR] <- 0
    return(classif)
}

BIMclass <- function(exp_dataset, ref_dataset = exp_dataset) {
    #' Function to to binarise bimodal distributions
    #' based on a 2-modes gaussian mixture model (with equal variances).
    #' Can be called with a reference dataset
    mc <- mclust::Mclust(na.omit(ref_dataset), modelNames = "E", G = 2, verbose = FALSE)
    classif <- rep(NA, length(exp_dataset))
    if (diff(mc$parameters$mean) > 0) {
        thresh_down <- max(mc$data[mc$classification == 1 & mc$uncertainty <= 0.05])
        thresh_up <- min(mc$data[mc$classification == 2 & mc$uncertainty <= 0.05])
        classif[exp_dataset <= thresh_down] <- 0
        classif[exp_dataset >= thresh_up] <- 1
    } else if (diff(mc$parameters$mean) < 0) {
        thresh_down <- max(mc$data[mc$classification == 2 & mc$uncertainty <= 0.05])
        thresh_up <- min(mc$data[mc$classification == 1 & mc$uncertainty <= 0.05])
        classif[exp_dataset <= thresh_down] <- 0
        classif[exp_dataset >= thresh_up] <- 1
    }
    return(classif)
}

norm_fun_lin <- function(xdat, reference = xdat) {
    #' Function for normalisation of zero-inflated data
    x_proc <- (xdat - quantile(reference, 0.01, na.rm = T))/quantile(xdat - quantile(reference, 
        0.01, na.rm = T), 0.99, na.rm = T)
    x_proc[x_proc < 0] <- 0
    x_proc[x_proc > 1] <- 1
    x_proc
}

norm_fun_sig <- function(xdat, reference = xdat) {
    #' Function for normalisation of unimodal data
    xdat <- xdat - median(reference, na.rm = T)
    lambda <- log(3)/mad(reference, na.rm = T)
    transformation <- function(x) {
        y <- 1/(1 + exp(-lambda * x))
        y
    }
    transformation(xdat)
}

norm_fun_bim <- function(xdat, reference = xdat) {
    #' Function for normalisation of unimodal data
    not_na_xdat <- !is.na(xdat)
    not_na_ref <- !is.na(reference)
    mc <- mclust::Mclust(reference[not_na_ref], modelNames = "E", G = 2, verbose = FALSE)
    pred <- mclust::predict.Mclust(mc, xdat[not_na_xdat])
    normalization <- rep(NA, length(xdat))
    if (diff(mc$parameters$mean) > 0) {
        normalization[not_na_xdat] <- pred$z[, 2]
    } else if (diff(mc$parameters$mean) < 0) {
        normalization[not_na_xdat] <- pred$z[, 1]
    }
    normalization
}


criteria_iter <- function(columns, data, genes) {
    #' Compute criteria for a subset of genes
    #'
    #' Data should be generated calling
    #' bigmemory::attach.big.matrix()
    #' on a big_matrix_descriptor.
    #'
    #' This is an auxiliary function intended
    #' to be called by compute_criteria
    #' and not directly by the user.

    criterix <- foreach::foreach(i = columns, .combine = rbind) %do% {
        x <- na.omit(unlist(data[, i]))
        criteria.iter <- list(Gene = genes[i], Dip = NA, BI = NA, Kurtosis = NA, 
            DropOutRate = NA, MeanNZ = NA, DenPeak = NA, Amplitude = max(x) - min(x))

        if (criteria.iter$Amplitude != 0) {
            criteria.iter$Dip <- diptest::dip.test(x)$p.value
            criteria.iter$BI <- BI(x)
            criteria.iter$Kurtosis <- moments::kurtosis(x) - 3
            criteria.iter$DropOutRate <- sum(x == 0)/length(x)
            criteria.iter$MeanNZ <- sum(x)/sum(x != 0)
            den <- density(x, na.rm = T)
            criteria.iter$DenPeak <- den$x[which.max(den$y)]
        }

        as.data.frame(criteria.iter)
    }
    criterix
}

compute_criteria <- function(exp_dataset, n_threads = parallel::detectCores() - 2, 
    descriptor_filename = NULL) {
    #' Function used to compute all statistical tools and criteria
    #' needed to perform the classification of distributions
    #' in the following categories:
    #'  * discarded
    #'  * zero-inflated
    #'  * unimodal
    #'  * bimodal
    #'

    exp_dataset <- exp_dataset %>%
        tibble::rownames_to_column("individual_id")

    .remove_descriptor <- FALSE
    if (is.null(descriptor_filename)) {
        descriptor_filename <- glue::glue("{date()} {random_words()}")
        .remove_descriptor <- TRUE
    }

    backing_file <- glue::glue("{descriptor_filename}.bin")
    descriptor_file <- glue::glue("{descriptor_filename}.desc")

    genes <- exp_dataset %>%
        dplyr::select(-individual_id) %>%
        colnames()
    n_genes <- length(genes)

    cl <- snow::makeSOCKcluster(names = rep("localhost", n_threads))
    doSNOW::registerDoSNOW(cl)

    big_exp_dataset <- exp_dataset %>%
        dplyr::select(-individual_id) %>%
        as.data.frame %>%
        bigmemory::as.big.matrix(type = "double", separated = FALSE, backingfile = backing_file, 
            descriptorfile = descriptor_file)

    big_exp_descriptor <- bigmemory::describe(big_exp_dataset)
    gene_iterator <- split_in_n(1:ncol(big_exp_dataset), n_threads)

    criteria <- foreach::foreach(i = gene_iterator, .combine = rbind, .inorder = TRUE, 
        .export = c("criteria_iter", "BI")) %dopar% {
        require(foreach)
        require(mclust)
        yy <- bigmemory::attach.big.matrix(big_exp_descriptor)
        criteria_iter(i, yy, genes)
    }

    stopCluster(cl)

    threshold <- median(criteria$Amplitude)/10
    # Added `tibble` call to enable the use of dplyr operators.
    criteria <- criteria %>%
        dplyr::tibble() %>%
        dplyr::mutate(Category = ifelse(Amplitude < threshold | DropOutRate > 0.95, 
            "Discarded", NA)) %>%
        dplyr::mutate(Category = ifelse(is.na(Category) & (BI > 1.5 & Dip < 0.05 & 
            Kurtosis < 1), "Bimodal", Category)) %>%
        dplyr::mutate(Category = ifelse(is.na(Category) & DenPeak < threshold, "ZeroInf", 
            Category)) %>%
        dplyr::mutate(Category = ifelse(is.na(Category), "Unimodal", Category))

    if (.remove_descriptor) {
        unlink(backing_file)
        unlink(descriptor_file)
    }

    criteria %>%
        tibble::column_to_rownames("Gene")
}

# TODO : refactor or remove this function (superseded by compute_criteria)
compute_criteria_sequential <- function(exp_dataset, individual_id) {
    #' Quasi exact copy of the original compute_criteria
    #' function, only that the hardcoded 'PATIENT_ID' has
    #' been replaced by a function argument in
    #' Here we compute all statistical tools and criteria
    #' needed to perform the classification of distributions
    #' in the following categories:
    #'  discarded,
    #'  zero-inflated,
    #'  unimodal,
    #'  bimodal

    # defuse (quote) individual_id :
    individual_id <- tryCatch(rlang::enquo(individual_id), error = function(e) e)

    exp_dataset <- exp_dataset %>%
        select(-{
            {
                individual_id
            }
        })

    criteria <- tibble(Gene = colnames(exp_dataset), Dip = NA, BI = NA, Kurtosis = NA, 
        DropOutRate = NA, MeanNZ = NA, DenPeak = NA, Amplitude = NA)

    # Compute
    pb = txtProgressBar(min = 1, max = ncol(exp_dataset), initial = 1)
    for (i in 1:ncol(exp_dataset)) {
        x <- na.omit(unlist(exp_dataset[, i]))
        criteria$Amplitude[i] <- max(x) - min(x)

        if (criteria$Amplitude[i] != 0) {
            criteria$Dip[i] <- dip.test(x)$p.value
            criteria$BI[i] <- BI(x)
            criteria$Kurtosis[i] <- kurtosis(x) - 3
            criteria$DropOutRate[i] <- sum(x == 0)/length(x)
            criteria$MeanNZ[i] <- sum(x)/sum(x != 0)
            den <- density(x, na.rm = T)
            criteria$DenPeak[i] <- den$x[which.max(den$y)]
        }

        setTxtProgressBar(pb, i)
    }

    threshold <- median(criteria$Amplitude)/10
    criteria <- criteria %>%
        mutate(Category = ifelse(Amplitude < threshold | DropOutRate > 0.95, "Discarded", 
            NA)) %>%
        mutate(Category = ifelse(is.na(Category) & (BI > 1.5 & Dip < 0.05 & Kurtosis < 
            1), "Bimodal", Category)) %>%
        mutate(Category = ifelse(is.na(Category) & DenPeak < threshold, "ZeroInf", 
            Category)) %>%
        mutate(Category = ifelse(is.na(Category), "Unimodal", Category))

    return(criteria)
}

binarize_exp <- function(exp_dataset, ref_dataset, ref_criteria, gene) {
    #' function to apply the proper binarization method depending
    #' on the gene expression distribution category
    #' Sorts ref_criteria rows according to exp_dataset

    .col_names <- colnames(exp_dataset)

    # Boolean flags to verify there is a reference for each and every single gene of
    # the exp_dataset
    .is_subset_of_ref_dataset <- Reduce(`&`, (.col_names %in% colnames(ref_dataset)))
    .is_subset_of_criteria <- Reduce(`&`, (.col_names %in% rownames(ref_criteria)))
    stopifnot(.is_subset_of_criteria, .is_subset_of_ref_dataset)

    # Sort the criteria according to the order of appearance in exp_dataset
    ref_criteria <- ref_criteria[.col_names, ]
    ref_dataset <- ref_dataset[, .col_names]

    if (!missing(gene)) {

        stopifnot(gene %in% rownames(ref_criteria), gene %in% colnames(ref_dataset))

        gene_cat <- ref_criteria[gene, "Category"]
        x <- unlist(dplyr::select(exp_dataset, gene))
        x_ref <- unlist(dplyr::select(ref_dataset, gene))

        if (gene_cat == "Discarded") {
            gene_bin <- rep(NA, length(x))
        } else if (gene_cat == "Bimodal") {
            gene_bin <- BIMclass(x, x_ref)

        } else {
            gene_bin <- OSclass(x, x_ref)
        }

        names(gene_bin) <- colnames(exp_dataset)

        return(gene_bin)

    } else {

        # given the prior subset checks, this condition seems unreachable.  TODO : verify
        # and remove it
        if (dim(exp_dataset)[2] != dim(ref_criteria)[1]) {
            stop("Different number of genes")
        }

        # these vectors should match the order
        logi_dis <- ref_criteria$Category == "Discarded"
        logi_OS <- ref_criteria$Category == "Unimodal" | ref_criteria$Category == 
            "ZeroInf"
        logi_bim <- ref_criteria$Category == "Bimodal"

        exp_dataset[, logi_dis] <- sapply(exp_dataset[, logi_dis], function(x) rep(NA, 
            length(x)))
        exp_dataset[, logi_OS] <- mapply(OSclass, as.data.frame(exp_dataset[, logi_OS]), 
            as.data.frame(ref_dataset[, logi_OS]))
        exp_dataset[, logi_bim] <- mapply(BIMclass, as.data.frame(exp_dataset[, logi_bim]), 
            as.data.frame(ref_dataset[, logi_bim]))

        colnames(exp_dataset) <- .col_names
        return(exp_dataset)
    }

}


normalize_exp <- function(exp_dataset, ref_dataset, ref_criteria, gene) {
    #' Normalise gene expression dataset
    #' using a reference dataset and reference
    #' criteria.

    .col_names <- colnames(exp_dataset)

    # Boolean flags to verify there is a reference for each and every single gene of
    # the exp_dataset
    .is_subset_of_ref <- Reduce(`&`, (.col_names %in% colnames(ref_dataset)))
    .is_subset_of_criteria <- Reduce(`&`, (.col_names %in% rownames(ref_criteria)))
    stopifnot(.is_subset_of_criteria, .is_subset_of_ref)

    ref_criteria <- ref_criteria[.col_names, ]
    ref_dataset <- ref_dataset[, .col_names]

    if (!missing(gene)) {

        stopifnot(gene %in% rownames(ref_criteria), gene %in% colnames(ref_dataset))

        gene_cat <- ref_criteria[gene, "Category"]
        x <- unlist(dplyr::select(exp_dataset, gene))
        x_ref <- unlist(dplyr::select(ref_dataset, gene))

        if (gene_cat == "Discarded") {
            gene_bin <- rep(NA, length(x))

        } else if (gene_cat == "Bimodal") {
            gene_bin <- norm_fun_bim(x, x_ref)

        } else if (gene_cat == "Unimodal") {
            gene_bin <- norm_fun_sig(x, x_ref)

        } else {
            gene_bin <- norm_fun_lin(x, x_ref)
        }

        names(gene_bin) <- colnames(exp_dataset)

        return(gene_bin)

    } else {

        # given the prior subset checks, this condition seems unreachable.  TODO : verify
        # and remove it
        if (dim(exp_dataset)[2] != dim(ref_criteria)[1]) {
            stop("Different number of genes")
        }

        logi_dis <- ref_criteria$Category == "Discarded"
        logi_uni <- ref_criteria$Category == "Unimodal"
        logi_zero <- ref_criteria$Category == "ZeroInf"
        logi_bim <- ref_criteria$Category == "Bimodal"

        exp_dataset[, logi_dis] <- sapply(exp_dataset[, logi_dis], function(x) rep(NA, 
            length(x)))
        exp_dataset[, logi_uni] <- mapply(norm_fun_sig, as.data.frame(exp_dataset[, 
            logi_uni]), as.data.frame(ref_dataset[, logi_uni]))
        exp_dataset[, logi_zero] <- mapply(norm_fun_lin, as.data.frame(exp_dataset[, 
            logi_zero]), as.data.frame(ref_dataset[, logi_zero]))
        exp_dataset[, logi_bim] <- mapply(norm_fun_bim, as.data.frame(exp_dataset[, 
            logi_bim]), as.data.frame(ref_dataset[, logi_bim]))

        return(exp_dataset)
    }

}
"""
