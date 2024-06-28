test_that("ldsc_h2 runs", {

  skip()

  gwas <- arrow::read_parquet("~/projects/move/speed_scz2022_eur.parquet", col_select = c("SNP", "Z","N")) |>
    dplyr::mutate(SNP = stringr::str_c("rs", as.integer(SNP)))

  res <- ldsc_h2(gwas)


})



test_that("partitioned ldscores", {
  skip()
  gwas <- arrow::read_parquet("tests/testthat/fixtures/move/speed_scz2022_eur.parquet")
  weights <- arrow::read_parquet("tests/testthat/fixtures/move/weights2.parquet") |>
    dplyr::arrange(CHR, BP)
  ld <- arrow::read_parquet("tests/testthat/fixtures/move/ld.parquet")
  annot <- arrow::read_parquet("tests/testthat/fixtures/move/annot.parquet")


  merged <- dplyr::inner_join(weights, ld, by = "SNP") |>
    dplyr::inner_join(gwas, by = "SNP")
  remove <- c(colnames(weights), colnames(gwas))

  expect_no_error(res <- ldscore(
    y = merged$Z^2,
    x = dplyr::select(merged,-dplyr::all_of(remove)) |> as.matrix(),
    w = merged$L2,
    N = merged$N,
    M = annot$m50
    ))
  tidy_results(res)


})


test_that("univariate ldsc", {
  skip()

  gwas <- arrow::read_parquet("~/projects/move/speed_scz2022_eur.parquet", col_select = c("SNP", "Z","N")) |>
    dplyr::mutate(SNP = stringr::str_c("rs", as.integer(SNP)))

  weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))
  twostep <- 30
  merge <- dplyr::inner_join(weights, gwas, by = "SNP")
  y <- merge$Z^2
  w <- merge$L2
  x <- w |> as.matrix()
  N <- merge$N
  M <- 1173569

  results <- ldscore(y,x,w,N,M)

  dplyr::tibble(
    h2 = results$tot, h2_se = results$tot_se,
    int = results$int, int_se = results$int_se,
    mean_chi2 = mean(merge$Z^2),
    lambda_gc = median(merge$Z^2) / 0.4549
  )




  n_snp <- dim(x)[1]
  n_annot <- dim(x)[2]
  M_tot <- sum(M)
  stopifnot(length(M) == n_annot)
  x_tot <- rowSums(x)

  # provide a starting estimate of heritability
  hsq <- M_tot * (mean(y) - 1) / mean((x_tot * N))

  # first update of weights
  initial_w <- get_weights(ld = x_tot, w_ld = w, N = N, M = M_tot, hsq = hsq)

  # Normalise by the mean of N and add intercept
  Nbar <- mean(N)
  x <- (N*x) / Nbar
  x <- cbind(1, x)


  # -------------------------------------------------------------------------
  ll <- univariate_ldsc(
    y = y,
    x = x,
    x_tot = x_tot,
    w = w,
    N = N,
    initial_w = initial_w,
    M = M,
    Nbar = Nbar
  )

  ll



  # -------------------------------------------------------------------------



  yp <- y
  ii <- y < twostep
  x1 <- x[ii, ]
  n1 <- sum(ii)
  yp1 <- yp[ii]
  w1 <- w[ii]
  N1 <- N[ii]
  initial_w1 <- initial_w[ii]

  # run two rounds of iterated weighted linear regression to improve weight
  update_func1 <- function(x) update_func(x, ld_tot = x1, w_ld = w1, N = N1, M = M, Nbar = Nbar)
  tmp <- iterated_weights(x = x1, y = yp1, w = initial_w1, update_function = update_func1)
  step1_res <- lstq_jackknife(tmp$x, tmp$y)
  step1_int <- step1_res$full_est[1]
  step1_int_se <- step1_res$se[1]


  # estimate h2, with pre-estimated intercept --------------------------------

  updated_s <- update_separators(step1_res$separators, ii)
  yp <- y - step1_int
  x <- x[, 2]
  c <- sum(initial_w * x) / sum(initial_w * x^2)

  update_function2 <- function(x) update_func(x, x_tot, w, N, M, Nbar, step1_int)
  tmp <- iterated_weights(x = x, y = yp, w = initial_w, update_function = update_function2)
  step2_res <- lstq_jackknife(as.matrix(tmp$x), tmp$y, separators = updated_s)


  res <- combine_twostep_jknives(step1_res, step2_res, M, c, Nbar)
  res[["int"]] <- step1_int
  res[["int_se"]] <- step1_int_se
  res



})


