similarity <- function(s1, s2, mu0, k0, nu0, L0, nh, lg) {
  
  sbar1 <- sum(s1) / nh
  sbar2 <- sum(s2) / nh
  
  Vs <- matrix(0, nrow = 2, ncol = 2)
  
  for (ii in 1:nh) {
    
    s_sbar1 <- s1[ii] - sbar1
    s_sbar2 <- s2[ii] - sbar2
    
    Vs[1, 1] <- Vs[1, 1] + s_sbar1^2
    Vs[1, 2] <- Vs[1, 2] + s_sbar1 * s_sbar2
    Vs[2, 1] <- Vs[2, 1] + s_sbar2 * s_sbar1
    Vs[2, 2] <- Vs[2, 2] + s_sbar2^2
  }
  
  kn <- k0 + nh
  nun <- nu0 + nh
  
  sbar_mu01 <- sbar1 - mu0[1]
  sbar_mu02 <- sbar2 - mu0[2]
  
  Vsbarmu <- matrix(0, nrow = 2, ncol = 2)
  Vsbarmu[1, 1] <- sbar_mu01^2
  Vsbarmu[1, 2] <- sbar_mu01 * sbar_mu02
  Vsbarmu[2, 1] <- sbar_mu02 * sbar_mu01
  Vsbarmu[2, 2] <- sbar_mu02^2
  
  Ln <- matrix(0, nrow = 2, ncol = 2)
  Ln[1, 1] <- L0[1] + Vs[1, 1] + k0 * nh / (k0 + nh) * Vsbarmu[1, 1]
  Ln[1, 2] <- L0[2] + Vs[1, 2] + k0 * nh / (k0 + nh) * Vsbarmu[1, 2]
  Ln[2, 1] <- L0[3] + Vs[2, 1] + k0 * nh / (k0 + nh) * Vsbarmu[2, 1]
  Ln[2, 2] <- L0[4] + Vs[2, 2] + k0 * nh / (k0 + nh) * Vsbarmu[2, 2]
  
  dL0 <- det(matrix(L0, nrow = 2, ncol = 2))
  dLn <- det(Ln)
  
  out <- -nh * log(pi) +
    (lmvgamma((0.5 * nun),2) - lmvgamma((0.5 * nu0),2)) +
    (0.5 * nu0 * log(dL0) - 0.5 * nun * log(dLn)) +
    (log(k0) - log(kn))
  
  if (!lg) {
    out <- exp(out)
  }
  
  return(out)
  
}
