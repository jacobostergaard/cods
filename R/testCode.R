# Test scripts for cods, set parameters and simulate data, then run esimtation procedures
# texLib = "~/Google Drive/PhD/TeX/JournalOfMathematicalBiology/Git/fig/"


# Set parameters
setPars <- function(model){
  # Cointegration parameters:
  if(model == 0){
    # Independent
    a = b = matrix(0,nr=3,nc=3)
    A = diag(3)
    B = diag(3)
    assign("model","Independent", envir=.GlobalEnv)
  } else if(model == 1){
    # Uni-directional
    a = matrix(c(-.5, 0,0),byrow=FALSE,nrow=3)
    b = matrix(c(1,-1,0),byrow=FALSE,nrow=3)
    A = as.matrix(c(1, 0, 0))
    B = as.matrix(c( 1,-1, 0))
    assign("model","Uni-directional", envir=.GlobalEnv)
  } else if(model == 2){
    # Bi-directional
    a = matrix(c(-.5, .5,0),byrow=FALSE,nrow=3)
    b = matrix(c(1,-1,0),byrow=FALSE,nrow=3)
    A = as.matrix(c(1, -1, 0))
    B = as.matrix(c( 1,-1, 0))
    assign("model","Bi-directional", envir=.GlobalEnv)
  } else if(model == 3){
    # All-to-all
    a = matrix(c(-.5,.25,.25,.25,-.5,.25),byrow=FALSE,nrow=3)
    b = matrix(c(1,0,-1,0,1,-1),byrow=FALSE,nrow=3)
    A = matrix(c(-5, 1, 4, 1.5, -4, 2.6), nr=3,nc=2,byrow=FALSE)
    B = matrix(c(1, 0, -1, 0, 1, -1), nr=3,nc=2,byrow=FALSE)
    assign("model","All-to-all", envir=.GlobalEnv)
  } else {
    a = b = matrix(0,nr=3,nc=3)
  }
  n.sim = 1000000
  dt = 0.1
  n = 1000
  dt.sim = n*dt/n.sim

  assign("Alpha", a, envir=.GlobalEnv)
  assign("Beta", b, envir=.GlobalEnv)
  assign("A", A, envir=.GlobalEnv)
  assign("B", B, envir=.GlobalEnv)
  assign("Pi", a%*%t(b), envir=.GlobalEnv)
  assign("omega", c(0,0,0), envir=.GlobalEnv)
  assign("n.sim", n.sim, envir=.GlobalEnv)
  assign("dt", dt, envir=.GlobalEnv)
  assign("dt.sim", dt.sim, envir=.GlobalEnv)
  assign("p", 3, envir=.GlobalEnv)
  assign("lvl", c(6,4,4), envir=.GlobalEnv)
  assign("freq", c(6,5,5), envir=.GlobalEnv)
  assign("z0", c(lvl[1],0,0,lvl[2],-lvl[3],0), envir=.GlobalEnv)     # Initial coordinates
  assign("S_phi", diag(1,p)*1, envir=.GlobalEnv)      # phi noise
  assign("S_gam", diag(1,p)*1, envir=.GlobalEnv)    # gamma noise
}

runOsc <- function(model="std"){
  # set.seed(1234)
  if(model=="std"){
    assign("freq", c(6,5,5), envir=.GlobalEnv)
    assign("lvl", c(1,1,1), envir=.GlobalEnv)
    assign("S_gam", diag(1,p)*0, envir=.GlobalEnv)
    assign("z0", c(lvl[1],0,0,lvl[2],-lvl[3],0), envir=.GlobalEnv)     # Initial coordinates
  } else if(model=="win"){
    assign("lvl", c(6,4,4), envir=.GlobalEnv)
    assign("freq", c(0,0,0), envir=.GlobalEnv)
    assign("z0", c(lvl[1],0,0,lvl[2],-lvl[3],0), envir=.GlobalEnv)     # Initial coordinates
    assign("S_gam", diag(1,p)*1, envir=.GlobalEnv)

  }

  Z = oscillator(n.sim,dt.sim,z0,Alpha,Beta,omega,freq,lvl,S_phi,S_gam,model)
  idx = seq(1,n.sim,dt/dt.sim)
  Z = data.frame(t(Z[,idx]))
  names(Z) = c("t",paste0(c("x","y"),c(rbind(1:p,1:p))))
  phi = unwrapPhase(Z)
  assign("Z", Z, envir=.GlobalEnv)
  assign("phi", phi, envir=.GlobalEnv)
}

runTest <- function(m=2){

    # Simulate and bootstrap
    setPars(m)
    #checkI1conditions()
    set.seed(12345)
    runOsc("win"); #head(Z)
    M0 = bootstrap(phi,B = 1000,dt=0.1); #head(M)

    # Pick r:
    r = rankTest(M0)$r
    M1 = johansen(phi,r,dt=0.1)
    M2 = johansen(phi,r,A=A,H=B,dt=0.1)

    SE1 = getSE(M1)
    SE2 = getSE(M2)
    LRT = LRtest(M1,M2)

    print(paste(model,"system"));
    print(paste0("Estimated r=",rankTest(M0)$r));
    print(round(rankTest(M0)$pVal,3));
    print(data.frame(M1 =round(SE1,3),M2= round(SE2,3)));
    print(LRT)
}
