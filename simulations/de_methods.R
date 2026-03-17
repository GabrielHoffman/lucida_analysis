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
            fit$summary %>% 
      dplyr::select(-gene, -gene_id), 
            fit$overdispersion) %>%
      dplyr::rename(sigSq_g = "Subject") %>%
      mutate(theta = 1 / Cell) %>%
      dplyr::select(-Cell)
    }) %>%
    bind_rows
}

# Add dispersion to edgeR, DESeq2
# in pbDS(), remove design matrix from filterExpr()
# devtools::install_github("GabrielHoffman/muscat")
stopifnot(packageVersion("muscat") == "1.25.4")


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

    fit.pb = lucida(pb2, nobars(formula), cluster_id = 'cluster_id', nthreads = nthreads)
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

    if( nrow(res.neb) > 0){
      df <- bind_rows(df,
            res.neb %>%
            dplyr::rename(logFC = paste0("logFC_", coefTest), 
            P.Value = paste0("p_", coefTest)) %>%
            dplyr::select(cluster_id, ID, logFC, P.Value, sigSq_g, theta) %>%
            mutate(FDR = p.adjust(P.Value)) %>%
            mutate(Method = "nebula"))
    }
  }

  if( "nebula (HL)" %in% methods ){
    df.time[["nebula (HL)"]] <- system.time({
    res.neb.HL <- run_nebula(sce.sim, formula, cluster_id, method="HL", nthreads = nthreads)
    })

    if( nrow(res.neb.HL) > 0){
      df <- bind_rows(df,
            res.neb.HL %>%
            dplyr::rename(logFC = paste0("logFC_", coefTest), 
            P.Value = paste0("p_", coefTest)) %>%
            dplyr::select(cluster_id, ID, logFC, P.Value, sigSq_g, theta) %>%
            mutate(FDR = p.adjust(P.Value)) %>%
            mutate(Method = "nebula (HL)"))
    }
  }

  # dreamlet
  if( "dreamlet" %in% methods ){
    df.time[["dreamlet"]] <- system.time({
    res.proc <- processAssays(pb, nobars(formula))
    res.dl <- dreamlet(res.proc, nobars(formula))
    })

    df <- bind_rows(df,
            topTable(res.dl, coefTest, number=Inf) %>%
            as_tibble %>%
            mutate(Method = "dreamlet") %>%
            dplyr::rename(FDR = "adj.P.Val", cluster_id = "assay") %>%
            dplyr::select(-t, -B, -z.std, -AveExpr))
  }

  # muscat: edgeR, DESeq2
  if( any(c("edgeR", "DESeq2") %in% methods) ){

    df.time[["pb2"]] <- system.time({   
    # hypothesis test on _LAST_ fixed effect variable
    grpVariable = all.vars(nobars(formula))
    grpVariable = grpVariable[length(grpVariable)]

    sce.tmp2 <- prepSCE(sce.tmp, 
      kid = cluster_id, 
      sid = "id",
      gid = grpVariable)

    pb <- aggregateData(sce.tmp2)
    })

    pb[[grpVariable]] = pb$group_id
    pb$group_id = 1
    pb$group_id[seq(floor(ncol(pb)/2))] = 0

    # make sure numeric variables are numeric
    for(x in all.vars(nobars(formula))){
      if( is.numeric(sce.tmp[[x]]) ){
        pb[[x]] = as.numeric(as.character(pb[[x]]))
      }
    }

    design <- model.matrix(nobars(formula), colData(pb))



    tab.muscat <- lapply( c("edgeR", "DESeq2")[c("edgeR", "DESeq2") %in% methods],
      function(method){
      df.time[[method]] <<- system.time({

        success <- TRUE        
        tryCatch({
        res.muscat <- pbDS(pb, 
          method = method, 
          design = design, 
          coef = which(coefTest == colnames(design)), 
          min_cells = 2, 
          filter = "both")}, 
          error = function(e) success <<- FALSE)
      })

      if( ! success ) return( NULL )

      tab = res.muscat$table[[coefTest]] %>%
        bind_rows %>%
        as_tibble %>%
        mutate(Method = method) %>%
        dplyr::rename(ID = "gene", 
          P.Value = "p_val", 
          FDR = 'p_adj.glb')

      if( method == "DESeq2"){
        tab = tab %>%
          dplyr::select(-baseMean, -lfcSE, -stat, -p_adj.loc)
      }
      if( method == "edgeR"){
        tab = tab %>%
          dplyr::select(-logCPM, -F, -p_adj.loc)
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

    # include = setdiff(c("id", cluster_id, all.vars(formula)), x1)

    # info = colData(sce.pb) %>% 
    #         as.data.frame %>% 
    #         left_join( colData(sce.sim)[,include] %>%
    #           unique %>%
    #           as.data.frame, by=c("id", cluster_id))

    # i = match(paste(info$id, info[[cluster_id]]),
    #           paste(sce.pb$id, sce.pb[[cluster_id]]))

    # colData(sce.pb) = DataFrame(info)


    # coefID <- paste0(x1, levels(factor(colData(sce.pb)[[x1]]))[2])

    CTs <- unique(colData(sce.pb)[[cluster_id]])
    res.gp <- lapply( CTs, function(CT){
      sceSub <- sce.pb[,colData(sce.pb)[[cluster_id]] == CT]
      fit.gp <- glm_gp(sceSub, nobars(formula))

      test_de( fit.gp, coefTest, compute_lfc_se=TRUE) %>%
        as_tibble %>%
        dplyr::rename(ID = name, P.Value = pval, logFC = lfc, se = lfc_se) %>%
        mutate(cluster_id = CT) %>%
        dplyr::select(cluster_id, ID, logFC, se, P.Value) %>%
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

  if( nrow(df) > 1 ){
    df = df %>%
      mutate(Method = factor(Method, validMethods)) %>%
      droplevels
  }

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





