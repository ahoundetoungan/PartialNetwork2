#####################################################################################
# functions to be optimized
## part0: GX and Gy are observed
### No contextual effects
qmleopt0_0 <- function(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, seed, print.proc, S = 1L){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  out      <- NULL
  if(print.proc){
    cat("lambda: ", alpha, "\n")
    out    <- fcpal0_0(alpha, R, distr, Ilist, y, Gy, X, Kx, M, Ncum, FE)^2
    cat("objective: ", out, "\n")
    cat("*********************\n")
  } else{
    Out    <- fcpal0_0(alpha, R, distr, Ilist, y, Gy, X, Kx, M, Ncum, FE)^2
  }
  out
}

### contextual effects
qmleopt0_1 <- function(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, print.proc, seed, S = 1L){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  out      <- NULL
  if(print.proc){
    cat("lambda: ", alpha, "\n")
    out    <- fcpal0_1(alpha, R, distr, Ilist, y, Gy, X, GX1, Kx, Kx1, Kx2, M, Ncum, FE)^2
    cat("objective: ", out, "\n")
    cat("*********************\n")
  } else{
    out    <- fcpal0_1(alpha, R, distr, Ilist, y, Gy, X, GX1, Kx, Kx1, Kx2, M, Ncum, FE)^2
  }
  out
}

## part1: Gy is not observed and columns of GX are observed
### No contextual effects
qmleopt1_0 <- function(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, print.proc, seed, S = 1L){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  out      <- NULL
  if(print.proc){
    cat("lambda: ", alpha, "\n")
    out    <- fcpal1_0(alpha, R, distr, Ilist, y, X, Kx, M, Ncum, FE, S)^2
    cat("objective: ", out, "\n")
    cat("*********************\n")
  } else{
    out    <- fcpal1_0(alpha, R, distr, Ilist, y, X, Kx, M, Ncum, FE, S)^2
  }
  out
}

### contextual effects
qmleopt1_1 <- function(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, print.proc, seed, S = 1L){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  out      <- NULL
  if(print.proc){
    cat("lambda: ", alpha, "\n")
    out    <- fcpal1_1(alpha, R, distr, Ilist, y, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, S)^2
    cat("objective: ", out, "\n")
    cat("*********************\n")
    out
  } else{
    out    <- fcpal1_1(alpha, R, distr, Ilist, y, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, S)^2
  }
}

## part2: Gy is observed and columns of GX are observed
### contextual effects
qmleopt2_1 <- function(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, print.proc, seed, S = 1L){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  out      <- NULL
  if(print.proc){
    cat("lambda: ", alpha, "\n")
    out    <- fcpal2_1(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, S)^2
    cat("objective: ", out, "\n")
    cat("*********************\n")
  } else{
    out    <- fcpal2_1(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, S)^2
  }
  out
}

#####################################################################################
# functions that compute beta gamma and sigma2
## part0: GX and Gy are observed
### No contextual effects
qmlebeta0_0 <- function(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, seed, print.proc, S = 1L){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  efcpal0_0(alpha, R, distr, Ilist, y, Gy, X, Kx, M, Ncum, FE)
}

### contextual effects
qmlebeta0_1 <- function(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, print.proc, seed, S = 1L){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  efcpal0_1(alpha, R, distr, Ilist, y, Gy, X, GX1, Kx, Kx1, Kx2, M, Ncum, FE)
}

## part1: Gy is not observed and columns of GX are observed
### No contextual effects
qmlebeta1_0 <- function(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, print.proc, seed, S = 1L){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  efcpal1_0(alpha, R, distr, Ilist, y, X, Kx, M, Ncum, FE, S)
}

### contextual effects
qmlebeta1_1 <- function(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, print.proc, seed, S = 1L){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  efcpal1_1(alpha, R, distr, Ilist, y, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, S)
}

## part2: Gy is observed and columns of GX are observed
### contextual effects
qmlebeta2_1 <- function(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, print.proc, seed, S = 1L){
  assign(".Random.seed", seed, envir = .GlobalEnv)
  efcpal2_1(alpha, R, distr, Ilist, y, Gy, X, X1, GX1, X2, Kx, Kx1, Kx2, M, Ncum, FE, S)
}
