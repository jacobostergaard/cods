# Test a cpp function call

testCpp <- function(){
  alpha = matrix(c(1,1,0,1,0,1),nrow=3, byrow = FALSE)
  qr.Q(qr(alpha),complete=TRUE)
  alphaPerp = qr.Q(qr(alpha), complete = TRUE)[,2:3]
  t(alphaPerp)%*%alpha
  testFun(alpha)
}


