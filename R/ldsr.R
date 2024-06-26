req_cols <- c("SNP", "Z", "N")
ldsc_h2 <- function(sumstat, weights=NULL, M=NULL) {
  stopifnot("sumstat has to be a data.frame or tbl" = "data.frame" %in% class(sumstat))
  stopifnot("SNP, Z and N are required in `sumstat`" = all(req_cols %in% colnames(sumstat)))

  if(is.null(weights)) {
    weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))
    M <- 1173569
  } else {
    weights <- weights
    stopifnot(!is.null(M))
  }
  sumstat <- dplyr::select(sumstat, dplyr::all_of(req_cols))
  merged <- dplyr::inner_join(weights, sumstat, by = "SNP")

  results <- ldscore(
    y = merged$Z^2,
    x = merged$L2 |> as.matrix(),
    w = merged$L2,
    N = merged$N,
    M = M
  )

  dplyr::tibble(h2 = results$tot, h2_se = results$tot_se)


}
ldsc_rg <- function(){

}

celltype_analysis <- function(sumstat, weights, baseline_ldscores,baseline_M, celltype_ldscores, celltype_M) {
  stopifnot("SNP, Z and N are required in `sumstat`" = all(req_cols %in% colnames(sumstat)))
  # stopifnot(all(file.exists(c(weights, baseline_ldscores, celltype_ldscores))))

  sumstat <- dplyr::select(sumstat, dplyr::all_of(req_cols)) |>
    tidyr::drop_na()
  w <- dplyr::arrange(weights, CHR, BP) |>
    dplyr::select(dplyr::all_of(c("SNP", "L2")))

  merged <- dplyr::inner_join(w, sumstat, by = "SNP") |>
    dplyr::inner_join(baseline_ldscores, by = "SNP") |>
    dplyr::inner_join(celltype_ldscores, by = "SNP") |>

  non_ldscores <- c(colnames(sumstat),colnames(w))
  celltype_ldscore_names <-  colnames(dplyr::select(celltype_ldscores, -SNP))


  res <- purrr::map(celltype_ldscore_names, \(celltype)
             .celltype(
               celltype = celltype,
               merged = merged,
               baseline_M = baseline_M,
               celltype_M = celltype_M,
               non_ldscores = non_ldscores,
               celltype_ldscore_names = celltype_ldscore_names
               )) |>
    purrr::list_rbind()


}


.celltype <- function(celltype, merged, baseline_M, celltype_M, non_ldscores,celltype_ldscore_names) {

    cols_to_remove <- unique(c(non_ldscores, celltype_ldscore_names[!celltype_ldscore_names %in% celltype]))
    celltype_M <- dplyr::filter(celltype_M, annot == celltype)



    M <- dplyr::bind_rows(baseline_M, celltype_M) |> dplyr::pull(m50)
    y <- merged$Z^2
    x <- dplyr::select(merged,-dplyr::all_of(cols_to_remove)) |> as.matrix()
    w <- merged$L2
    N <- merged$N

    stopifnot("wrong dimensions between celltype ldscores and M"= length(M) == dim(x)[2])

    res <- ldscore(
      x = x,
      y = y,
      w = w,
      N = N,
      M = M
    )

    tidy_results(res) |>
      dplyr::filter(annot == celltype)

}




