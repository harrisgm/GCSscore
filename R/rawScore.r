rawScore <- function(diff1, diff2, SDT1, SDT2) {
  m.gamma <- 0.1
  N1 <- N2 <- 1
  rawS <- (diff1/N1-diff2/N2)/
    sqrt(m.gamma^2*(diff1^2/(N1^2) + diff2^2/(N2^2))+
           sum(SDT1^2)/(N1^2)+sum(SDT2^2)/(N2^2))
  rawS[rawS > 15.0] <- 15.0
  rawS[rawS < -15.0] <- -15.0
  return(rawS)
}
