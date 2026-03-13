suppressPackageStartupMessages({
library(lucida)
library(SingleCellExperiment)
library(GenomicDataStream)
library(tidyverse)
library(nebula)
library(dreamlet)
library(reformulas)
library(muscat)
library(glmGamPoi)
})

run_nebula = function(sce, formula, cluster_id, method="LN", nthreads = 1){

  lapply( unique(sce[[cluster_id]]), function(CT){
    message(CT)
    # cells of given type
    idx = (sce[[cluster_id]] == CT)

    # grouping variable
    ran_var <- all.vars(findbars(formula)[[1]])

    # order cells by grouping variable
    sceSub = sce[,order(sce[[ran_var]])]

    data = colData(sceSub[,idx])
    data = droplevels(data)

    success = TRUE
    if( nlevels(data[[ran_var]]) < 3 ){
      success = FALSE
    }

    tryCatch(
      {design = model.matrix(nobars(formula), data)},
      error = function(e) success <<- FALSE)

    if( ! success ) return(NULL)

    fit = nebula( 
      count = counts(sceSub[,idx]), 
      id = data[[ran_var]], 
      pred = design, 
      offset = sceSub$libSize[idx], 
      method = method, 
      ncore = nthreads)
    
    tibble(cluster_id = CT, 
            ID = fit$summary$gene,
            fit$summary %>% select(-gene, -gene_id), 
            fit$overdispersion) %>%
      rename(sigSq_g = "Subject") %>%
      mutate(theta = 1 / Cell) %>%
      select(-Cell)
    }) %>%
    bind_rows
}

# modify to extract theta
.DESeq2_new <- function(x, k, design, contrast, ct, cs) {
  library(DESeq2)
  cd <- colData(x)
  y <- assay(x, k)
  mode(y) <- "integer"
  y <- DESeqDataSetFromMatrix(y, cd, design)
  y <- suppressMessages(DESeq(y))
  tbl <- lapply(cs, function(c) {
      tbl <- results(y, contrast[, c])
      tbl <- muscat:::.res_df(tbl, k, ct, c)
      tbl$theta <- 1 / DESeq2::dispersions(y)
      rename(tbl, logFC = "log2FoldChange", 
          p_val = "pvalue", p_adj.loc = "padj")
  })
  list(table = tbl, data = y)
}
assignInNamespace(".DESeq2", .DESeq2_new, "muscat")


# modify to extract theta
.edgeR_new <- function(x, k, design, coef, contrast, ct, cs, treat) {
    library(edgeR)
    y <- assay(x, k)
    y <- suppressMessages(DGEList(y, 
        group = x$group_id[colnames(y)], 
        remove.zeros = TRUE))
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    # treat: test for DE relative to logFC threshold
    # else:  genewise NB GLM with quasi-likelihood test
    .fun <- ifelse(treat, glmTreat, glmQLFTest)
    tbl <- lapply(cs, function(c) {
        fit <- .fun(fit, coef[[c]], contrast[, c])
        tbl <- topTags(fit, n = Inf, sort.by = "none")
        # combine tables & reformat
        tbl <- rename(tbl$table, p_val = "PValue", p_adj.loc = "FDR")
        tbl <- muscat:::.res_df(tbl, k, ct, c)
        tbl$theta <- 1 / y$tagwise.dispersion
        tbl
    })
    list(table = tbl, data = y, fit = fit)
}
assignInNamespace(".edgeR", .edgeR_new, "muscat")


run_analysis <- function( sce.sim, formula, coefTest, cluster_id, methods, nthreads = 1, include_metadata = TRUE){

  validMethods <- c(  
    "lucida",
    "lucida [1 step]",
    "lucida [Bayesian]",
    "lucida [pb]",
    "nebula",
    "nebula (HL)",
    "dreamlet",
    "DESeq2",
    "edgeR",
    "glmGamPoi",
    "MAST")

  stopifnot( all(methods %in% validMethods))

  df <- tibble()
  df.time <- list()

  # lucida 
  if( "lucida" %in% methods ){
    df.time[["lucida"]] <- system.time({
    fit.lucida <- lucida(sce.sim, formula, cluster_id, nthreads = nthreads)
    })

    # merge with expression magnitude information
    df_mu <- lapply( names(fit.lucida), function(x){
      tibble(cluster_id = x,
        ID = rownames(fit.lucida[[x]]$fit),
        mu = fit.lucida[[x]]$fit$mu_mean)
    }) %>%
      bind_rows

    df <- bind_rows(df,
            lucida::results(fit.lucida, coefTest, expand=TRUE) %>%
            mutate(Method = "lucida"))
  }

  # lucida one step
  if( "lucida [1 step]" %in% methods ){

    df.time[["lucida [1 step]"]] <- system.time({
    fit.lucida1 <- lucida(sce.sim, formula, cluster_id, shrinkDispersion=FALSE, nthreads = nthreads)
    })

    df <- bind_rows(df,
            lucida::results(fit.lucida1, coefTest, expand=TRUE) %>%
            mutate(Method = "lucida [1 step]"))
  }

  if( any(c("dreamlet", "lucida [pb]", "DESeq2", "edgeR") %in% methods) ){
    df.time[["pseudobulk"]] <- system.time({
    sce.tmp = SingleCellExperiment(list(
                counts = counts(sce.sim)), 
                colData = colData(sce.sim))
   
    sce.tmp$id <- lapply(all.vars(formula), function(x){
      colData(sce.tmp)[,x]
      }) %>%
      bind_cols %>%
      apply(1, function(x) paste(x, collapse="_"))

    pb <- aggregateToPseudoBulk(
          sce.tmp ,
         cluster_id = cluster_id,
         sample_id = "id")   
    })
  }

  # lucida pseudobulk
  if( "lucida [pb]" %in% methods ){
    df.time[["lucida [pb]"]] <- system.time({
    pb2 <- lapply(assayNames(pb), function(x){
      SingleCellExperiment(list(counts = assay(pb, x)), 
                colData = data.frame(colData(pb), 
                  cluster_id = x))
      })
    pb2 <- do.call(cbind, pb2)

    pb2$libSize = colSums(counts(pb2))
    pb2 = pb2[,pb2$libSize > 0]

    fit.pb = lucida(pb2, ~ Dx, cluster_id = 'cluster_id', nthreads = nthreads)
    })

    df <- bind_rows(df,
            lucida::results(fit.pb, coefTest, expand=TRUE) %>%
            mutate(Method = "lucida [pb]"))
  }

  # nebula
  if( "nebula" %in% methods ){

    df.time[["nebula"]] <- system.time({
    res.neb <- run_nebula(sce.sim, formula, cluster_id, nthreads = nthreads)
    })

    df <- bind_rows(df,
            res.neb %>%
            rename(logFC = paste0("logFC_", coefTest), 
            P.Value = paste0("p_", coefTest)) %>%
            select(cluster_id, ID, logFC, P.Value, sigSq_g, theta) %>%
            mutate(FDR = p.adjust(P.Value)) %>%
            mutate(Method = "nebula"))
  }

  if( "nebula (HL)" %in% methods ){
    df.time[["nebula (HL)"]] <- system.time({
    res.neb.HL <- run_nebula(sce.sim, formula, cluster_id, method="HL", nthreads = nthreads)
    })

    df <- bind_rows(df,
            res.neb.HL %>%
            rename(logFC = paste0("logFC_", coefTest), 
            P.Value = paste0("p_", coefTest)) %>%
            select(cluster_id, ID, logFC, P.Value, sigSq_g, theta) %>%
            mutate(FDR = p.adjust(P.Value)) %>%
            mutate(Method = "nebula (HL)"))
  }

  # dreamlet
  if( "dreamlet" %in% methods ){
    df.time[["dreamlet"]] <- system.time({
    res.proc <- processAssays(pb, ~ Dx)
    res.dl <- dreamlet(res.proc, ~ Dx)
    })

    df <- bind_rows(df,
            topTable(res.dl, coefTest, number=Inf) %>%
            as_tibble %>%
            mutate(Method = "dreamlet") %>%
            rename(FDR = "adj.P.Val", cluster_id = "assay") %>%
            select(-t, -B, -z.std, -AveExpr))
  }


  # if( "DESeq2" %in% methods ){

  #   library(DESeq2)
  #   library(edgeR)

  #   df.time[["DESeq2"]] <- system.time({

  #   res.deq = lapply(assayNames(pb), function(CT){
  #     message(CT)
  #     # get counts
  #     countMatrix = assay(pb, CT)

  #     # filter samples
  #     libSize = colSums2(countMatrix)
  #     keep = libSize > 100
  #     countMatrix <- countMatrix[,keep,drop=FALSE]
  #     info = colData(pb)[keep,]

  #     # filter genes
  #     keep <- filterByExpr(countMatrix)
  #     countMatrix <- countMatrix[keep,,drop=FALSE]

  #     dds <- DESeqDataSetFromMatrix(
  #              countMatrix,
  #             DataFrame(info), ~ Dx)
  #     dds <- DESeq(dds, quiet=TRUE)
  #     tbl <- results(dds, name = "Dx_Disease_vs_Control") 
  #     tbl$theta <- 1 / DESeq2::dispersions(dds)

  #     tbl %>%
  #       as.data.frame %>%
  #       rownames_to_column("ID") %>%
  #       mutate(cluster_id = CT) %>%
  #       as_tibble %>%
  #       filter(!is.na(log2FoldChange))
  #     }) %>%
  #     bind_rows 
  #   })

  #   df = res.deq %>%
  #     mutate(logFC = log2FoldChange, 
  #       P.Value = pvalue, 
  #       FDR = p.adjust(P.Value),
  #       Method = "DESeq2") %>%
  #     select(ID, cluster_id, logFC, P.Value, FDR, theta, Method) %>%
  #     bind_rows(df, .)
  # }

  # if( "edgeR" %in% methods ){

  #   library(edgeR)

  #   df.time[["edgeR"]] <- system.time({

  #   res.edgeR = lapply(assayNames(pb), function(CT){

  #     message(CT)
  #     # get counts
  #     countMatrix = assay(pb, CT)

  #     # filter samples
  #     libSize = colSums2(countMatrix)
  #     keep = libSize > 100
  #     countMatrix <- countMatrix[,keep,drop=FALSE]
  #     info = colData(pb)[keep,]

  #     # filter genes
  #     keep <- filterByExpr(countMatrix)
  #     countMatrix <- countMatrix[keep,,drop=FALSE]

  #     design <- model.matrix(~Dx, info)

  #     y <- suppressMessages(DGEList(
  #         countMatrix, 
  #         remove.zeros = TRUE))
  #     y <- calcNormFactors(y)
  #     y <- estimateDisp(y, design)
  #     fit <- glmQLFit(y, design)
  #     res <- glmQLFTest(fit, coef = "DxDisease")
  #     tbl <- topTags(res, n = Inf, sort.by = "none") %>%
  #               as.data.frame
  #     tbl$theta <- 1 / y$tagwise.dispersion

  #     tbl %>%
  #       as.data.frame %>%
  #       rownames_to_column("ID") %>%
  #       mutate(cluster_id = CT) %>%
  #       as_tibble 
  #     }) %>%
  #     bind_rows 
  #   })

  #   df <- res.edgeR %>%
  #     mutate( 
  #       P.Value = PValue, 
  #       FDR = p.adjust(P.Value),
  #       Method = "edgeR") %>%
  #     select(ID, cluster_id, logFC, P.Value, FDR, theta, Method) %>%
  #     bind_rows(df, .)
  # }

  # muscat: edgeR, DESeq2
  if( "DESeq2" %in% methods ){

    # hypothesis test on _LAST_ fixed effect variable
    grpVariable = all.vars(nobars(formula))
    grpVariable = grpVariable[length(grpVariable)]

    sce.tmp2 <- prepSCE(sce.tmp, 
      kid = cluster_id, 
      sid = "id",
      gid = grpVariable)

    pb <- aggregateData(sce.tmp2)

    pb[[grpVariable]] = pb$group_id

    design <- model.matrix(nobars(formula), colData(pb))

    tab.muscat <- lapply( c("edgeR", "DESeq2"), function(method){
      res.muscat = pbDS(pb, method = method, design, min_cells=2, filter="both")

      tab = res.muscat$table[[coefTest]] %>%
        bind_rows %>%
        as_tibble %>%
        mutate(Method = method) %>%
        rename(ID = "gene", 
          P.Value = "p_val", 
          FDR = 'p_adj.glb')

      if( method == "DESeq2"){
        tab = tab %>%
          select(-baseMean, -lfcSE, -stat, -p_adj.loc)
      }
      if( method == "edgeR"){
        tab = tab %>%
          select(-logCPM, -F, -p_adj.loc, -contrast)
      }
      tab
    }) %>%
      bind_rows

    df <- bind_rows(df, tab.muscat)
  }

  # glmGamPoi
  if( "glmGamPoi" %in% methods ){

    df.time[["glmGamPoi"]] <- system.time({

    x1 <- all.vars(nobars(formula))[1]
    x2 <- all.vars(findbars(formula)[[1]])
    sce.sim$id <- paste(colData(sce.sim)[,x1], 
                      colData(sce.sim)[,x2])

    sce.pb <- pseudobulk( sce.sim, group_by = vars(id, !!sym(cluster_id), Dx))

    # coefID <- paste0(x1, levels(factor(colData(sce.pb)[[x1]]))[2])

    CTs <- unique(colData(sce.pb)[[cluster_id]])
    res.gp <- lapply( CTs, function(CT){
      sceSub <- sce.pb[,colData(sce.pb)[[cluster_id]] == CT]
      fit.gp <- glm_gp(sceSub, nobars(formula))

      test_de( fit.gp, coefTest, compute_lfc_se=TRUE) %>%
        as_tibble %>%
        rename(ID = name, P.Value = pval, logFC = lfc, se = lfc_se) %>%
        mutate(cluster_id = CT) %>%
        select(cluster_id, ID, logFC, se, P.Value) %>%
        filter(!is.na(se)) %>%
        inner_join(data.frame(
          ID = names(fit.gp$overdispersions),
          theta = 1/fit.gp$overdispersions))
    }) %>%
      bind_rows %>%
      mutate( FDR = p.adjust(P.Value, "BH")) %>%
      mutate(Method = "glmGamPoi")
      })

    df <- bind_rows(df, res.gp)
  }

  if( include_metadata ){

    df <- df %>% 
      inner_join(metadata(sce.sim) %>%
        as_tibble, by=c('cluster_id', "ID"))
  } 

  df = df %>%
    mutate(Method = factor(Method, validMethods)) %>%
    droplevels

  if( "lucida" %in% methods ){
    df <- df %>%
          left_join(df_mu, by=c("cluster_id", "ID"))
  }

  # join df.time
  df.time2 = lapply(names(df.time), function(Method){
    data.frame(Method = Method, t(data.frame(df.time[[Method]])))
    }) %>%
  bind_rows %>%
  tibble

  list(df = df, df.time = df.time2)
}





