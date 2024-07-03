ii <- y < twostep
x1 <- x[ii,]
y1 = y[ii]
w1 = w[ii]
N1 <- N[ii]
initial_w1 <- initial_w[ii]
# run two rounds of iterated weighted linear regression to improve weight
update_func1 <- function(x) update_func(x, x1, w1,  N1, M, Nbar)
tmp <- iterated_weights(x = x1, y = y1, w = initial_w1, update_function = update_func1)



update_func <- function(wls, ld_tot, w_ld, N, M, Nbar, intercept = NULL) {
  # ld_tot should be sum of ldscores if partitioned by annotations
  # wls: first entry is intercept, second entry is coef

  if(is.null(intercept)) {

    hsq = M * wls[2] / Nbar
    intercept = max(wls[1])
    ld <- ld_tot[, 2] # drop intercept from x

  } else {

    hsq = M * wls[1] / Nbar
    ld <- ld_tot

  }
  get_weights(ld, w_ld, N, M, hsq, intercept)

}




# -------------------------------------------------------------------------

get_weights <- function(ld, w_ld, N, M, hsq, intercept=1) {

  # hsq bounded between 1, minimum weight value is 1.
  hsq = max(hsq, 0.0)
  hsq = min(hsq, 1.0)
  ld = pmax(ld, 1)
  w_ld = pmax(w_ld, 1)

  c = hsq * N / M
  het_w = 1.0 / (2*  ((intercept + (c*ld))^2) )
  oc_w = 1.0 / w_ld
  w = het_w * oc_w
  w


}

wls <- function(x,y,w) {

  # normalize weights, weight x and y by w
  w_norm <- (w / sum(w))
  x = w_norm * x
  y = w_norm * y
  XtX = crossprod(x, x)
  Xy = crossprod(x, y)
  res <- solve(XtX, Xy)

  # format the results to return a vector
  drop(t(res))

}


get_weights(ld, w_ld, N, M, hsq, intercept)

iterated_weights <- function(x,y,w, ld_tot, w_ld, N, M, Nbar, intercept = NULL) {

  w <- sqrt(w)
  for(i in 1:2) {
    ests <- wls(x,y,w)

    if(is.null(intercept)) {

      hsq = M * ests[2] / Nbar
      intercept = max(ests[1])
      ld <- ld_tot[, 2] # drop intercept from x
    
    } else {
    
      hsq = M * ests[1] / Nbar
      ld <- ld_tot
    
    }


    updated <- get_weights(ld, w_ld, N, M, hsq, intercept)
    new_w <- sqrt(updated)
    w <- new_w
  }
  w_norm <- (w / sum(w))

  list(
    x = w_norm * x,
    y = w_norm * y
  )

}