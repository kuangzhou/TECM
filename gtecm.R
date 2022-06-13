#' General Transfer Evidential c-means algorithm
#'
#' \code{tecm} computes a credal partition from a matrix of attribute data using the
#' Transfer Evidential c-means (TECM) algorithm.
#' We refer the code of ECM in EVCLUST package by Thierry Denoeux
#'
#' @param x input matrix of size n x d, where n is the number of objects and d the number of
#' attributes.
#' @param c Number of  clusters.
#' @param vs input matrix of size c x d, the prototype matrix in the source domain 
#' @param g0 Initial prototypes, matrix of size c x d. If not supplied, the prototypes are
#' initialized randomly.
#' @param type Type of focal sets ("simple": empty set, singletons and Omega;
#' "full": all \eqn{2^c} subsets of \eqn{\Omega}; "pairs": \eqn{\emptyset}, singletons,
#' \eqn{\Omega}, and all
#' or selected pairs).
#' @param Omega Logical. If TRUE (default), the whole frame is included (for types 'simple' and
#' 'pairs').
#' @param pairs Set of pairs to be included in the focal sets; if NULL, all pairs are
#' included. Used only if type="pairs".
#' @param ntrials Number of runs of the optimization algorithm (set to 1 if g0 is  supplied).
#' @param alpha Exponent of the cardinality in the cost function.
#' @param mybeta Exponent of masses (beta) in the cost function.
#' @param beta3 Parameter for the source domain knowledge 2 (distrance between the prototpyes in the source and target domain)
#' @param eta Exponent of the membership rij.
#' @param delta  Distance to the empty set in the target domain.
#' @param deltas Distance to the empty set in the source domain.
#' @param epsi Minimum amount of improvement.
#' @param disp If TRUE (default), intermediate results are displayed.
#'
#' @return The credal partition (an object of class \code{"credpart"}).
#'

library(evclust)
gtecm <- function(x, c, vs, g0 = NULL, type = "full", pairs = NULL, Omega = TRUE, ntrials = 1, alpha = 1, mybeta = 2, beta3 = 1, eta = 2,  delta = 10, deltas = 100, epsi = 0.001, disp = TRUE) {
    
  ##########################################################
  cs <-  nrow(vs)
  ct <- c
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  delta2 <- delta^2
  deltas2 <- deltas^2
  
  if ((ntrials > 1) & !is.null(g0)) {
    print("WARNING: ntrials>1 and g0 provided. Parameter ntrials set to 1.")
    ntrials <- 1
  }
  if (missing(vs)){
    cs = c
    vs <- x[sample(1:n, cs), ] + 0.1 * rnorm(cs * d, cs, d) 
    beta3 = 0
  }
  
  # the nubmer of clusters in the sourc domain
  cs <-  nrow(vs)
  ct <- c
  myF <- makeF(ct, type, pairs, Omega)
  f <- nrow(myF)
  card <- rowSums(myF[2:f, ])
  
  myFs <- makeF(cs, type, pairs, Omega)
  fs <- nrow(myFs)
  cards <- rowSums(myFs[2:f, ])
  #------------------------ iterations--------------------------------
  Jbest <- Inf
  for (itrial in 1:ntrials) {
    
    
    itrial <- 1
    if (is.null(g0)){
      g <- x[sample(1:n, ct), ] + 0.1 * rnorm(ct * d, ct, d) 
    } else {g <- g0}
    
    #calculate rij for the first time
    r = matrix(1, cs, ct)
    #distance from vi in the source domain to the j cluster in the target domain
    dst = matrix(0, cs, ct)
    for (j in 1: ct) {
      dst[, j] <- rowSums((vs - matrix(g[j, ], cs, d, byrow = TRUE))^2)
    }
    for(i in 1:cs){
      vect0 = (1/dst[i, ])^(1/(eta-1))
      for(j in 1:ct){
        r[i, j] = 1/(dst[i,j]^(1/(eta-1))*sum(vect0))
      }
    }
    
    
    pasfini <- TRUE
    Jold <- Inf
    
    # gplus: the center vector in the target domain for all the clusters including the specific and imprecise classes
    gplus <- matrix(0, f - 1, d)
    # gpluss: the center vector in the source domain for all the clusters including the specific and imprecise classes
    gpluss <- matrix(0, fs - 1, d)
    
    # The Ds matrix, the distances to the centers in the source domain
    for (i in 2:fs) {
      fi <- myFs[i, ]
      truc <- matrix(fi, cs, d)
      gpluss[i - 1, ] <- colSums(vs * truc)/sum(fi)
    }
    
    Ds <- matrix(0, n, f - 1)
    for (j in 1: (f - 1)) {
      Ds[, j] <- rowSums((x - matrix(gpluss[j, ], n, d, byrow = TRUE))^2)
    }
    
    iter <- 0
    while (pasfini) {
      iter <- iter + 1
      
      for (i in 2:f) {
        fi <- myF[i, ]
        truc <- matrix(fi, c, d)
        gplus[i - 1, ] <- colSums(g * truc)/sum(fi)
      }
      
      # calculation of distances to centers in the target domain
      D <- matrix(0, n, f - 1)
      for (j in 1: (f - 1)) {
        D[, j] <- rowSums((x - matrix(gplus[j, ], n, d, byrow = TRUE))^2) 
      }
      
      D1 = D 
      deltaT2 = delta2 
      
      # Calculation of masses
      m <- matrix(1/f, n, f - 1)
      
      # mold <- m
      for (i in 1:n) {
        vect0 <- D1[i, ] 
        for (j in 1:(f - 1)) {
          vect1 <- (rep(D1[i, j], f - 1)/vect0)^(1/(mybeta - 1))
          vect2 <- rep(card[j]^(alpha/(mybeta - 1)), f - 1)/(card^(alpha/(mybeta - 1)))
          vect3 <- vect1 * vect2
          m[i, j] <- 1/(sum(vect3) + (card[j]^alpha * D1[i, j]/deltaT2)^(1/(mybeta - 1)))
        }
      }
      
      # m <- mold
      
      
      # Calculation of centers, the H matrix in the paper 
      A <- matrix(0, c, c)
      for (k in 1:c) {
        for (l in 1:c) {
          indices = which(myF[, l] == 1 & myF[, k] == 1) # indices of all Aj including wk and wl
          indices <- indices - 1
          if (length(indices) == 0) 
            A[l, k] <- 0 else {
              for (jj in 1:length(indices)) {
                j <- indices[jj]
                mj <- m[, j]^mybeta
                A[l, k] <- A[l, k] + sum(mj) * card[j]^(alpha - 2)
              }
            }
        }
      }
      
      
      
      
      
      # H + beta2 * I matrix in the paper
      rDiag = colSums(r)
      A1 = A + beta3 * diag(rDiag)
      
      # Construction of the B matrix
      B <- matrix(0, c, d)
      for (l in 1:c) {
        indices = which(myF[, l] == 1) # indices of all Aj inclduing wl
        indices <- indices - 1
        mi <- matrix(card[indices]^(alpha - 1), n, length(indices), byrow = TRUE) * m[, indices]^mybeta
        s <- rowSums(mi)
        mats <- matrix(s, n, d)
        xim <- x * mats
        B[l, ] <- colSums(xim)
      }
      
      # B + beta2*vs in the paper
      B1 =  B + beta3 * t(r^eta) %*% vs
      
      # gold <- g
      g <- solve(A1, B1)
      # g <- gold
      
      # update rij
      #distance from vi in the source domain to the j cluster in the target domain
      for (j in 1: c) {
        dst[, j] <- rowSums((vs - matrix(g[j, ], cs, d, byrow = TRUE))^2)
      }
      # rold <- r
      for(i in 1:cs){
        vect0 = (1/dst[i, ])^(1/(eta-1))
        for(j in 1:ct){
          r[i, j] = 1/(dst[i,j]^(1/(eta-1))*sum(vect0))
        }
      }
      # r <- rold
      mvide <- 1 - rowSums(m)
      J <- sum((m^mybeta) * D * matrix(card^alpha, n, f - 1, byrow = TRUE)) + delta2 * sum(mvide^mybeta)
      J3 <- beta3 * sum(r^eta * dst)
      JT = J + J3
      # print(JT)
      if (disp) 
        cat(c(iter, JT),"\n")
      pasfini <- (abs(JT - Jold) > epsi)
      Jold <- JT
      
    }  # end while loop
    if (JT < Jbest) {
      Jbest <- JT
      mbest <- m
      gbest <- g
    }
    res <- c(itrial, JT, Jbest)
    names(res) <- NULL
    if (ntrials > 1) 
      print(res)
  }  #end for loop iter
  
  m <- cbind(1 - rowSums(mbest), mbest)
  clus <- extractMass(m, myF, g = gbest, method = "tecm", crit = Jbest)
  return(clus)
}
