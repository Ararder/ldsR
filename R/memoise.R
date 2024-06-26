

# truth <- t(x) %*% x
#
# l2 <- x[,54]
# X <- x[,-54]
#
# XtX <- t(X) %*% X
#
#
#
# MTM <- t(M) %*% M
# MTX <- t(M) %*% X
# XTX <- t(X) %*% X
#
# MTM_star <- rbind(
#   cbind(MTM, MTX),
#   cbind(t(MTX), XTX)
# )
#
# dim(MTM_star)
#
#
# MTM_star == truth
#
#
# get_memoised_block_values <- function(X, p, XtX) {
#
#   ptX <- t(X_base) %*% p
#   XtP <- t(X) %*% p
#
#   rbind(
#     cbind(MTM, MTX),
#     cbind(t(MTX), XTX)
#   )
#
#
# }
