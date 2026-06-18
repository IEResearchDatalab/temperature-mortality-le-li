################################################################################
#
# Period life table
#
################################################################################

lifetable <- function(mx, age, ax = 0.5) {

  ord <- order(age)
  mx <- mx[ord]
  age <- age[ord]

  n <- c(diff(age), Inf)

  qx <- (n * mx) / (1 + n * (1 - ax) * mx)
  qx[length(qx)] <- 1

  px <- 1 - qx

  lx <- vector("numeric", length(qx))
  lx[1] <- 100000
  for (i in 2:length(qx))
    lx[i] <- lx[i - 1] * px[i - 1]

  dx <- lx * qx

  Lx <- vector("numeric", length(qx))
  for (i in 1:(length(qx) - 1))
    Lx[i] <- n[i] * (lx[i] - (1 - ax) * dx[i])

  Lx[length(qx)] <- ifelse(mx[length(qx)] > 0,
    lx[length(qx)] / mx[length(qx)], 0)

  Tx <- rev(cumsum(rev(Lx)))

  ex <- Tx / lx
  ex[is.nan(ex)] <- 0

  Sx <- vector("numeric", length(qx))
  Sx[1:(length(Sx) - 1)] <- Lx[2:length(Lx)] / Lx[1:(length(Lx) - 1)]

  ax_out <- rep(ax, length(mx))

  data.table(age = age, n = n, mx = mx, qx = qx, px = px, lx = lx,
    dx = dx, Lx = Lx, Sx = Sx, Tx = Tx, ex = ex, ax = ax_out)
}