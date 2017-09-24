# Wrapper functions for simulating data and performing bootstrap test and inference (restrictions on matrices)
# Output here is more user friendly and useable.


# Wrapper function for the johansen cpp function.
    johansen <-function(X,r=1,A=NULL,H=NULL,dt=0.1){
      X = as.matrix(X)
      N = nrow(X)-1
      p = ncol(X)

      if(is.null(A)){
        A = as.matrix(diag(p))
      }
      if(is.null(H)){
        H = as.matrix(diag(p))
      }
      df = r*(p-ncol(A))+r*(p-ncol(H))

      out = johansenCpp(X,r,as.matrix(A),as.matrix(H),dt)

      if(r > 0){
        a.hat = matrix(out[1:r,],nr=p,byrow = TRUE)
        b.hat = matrix(out[(r+1):(2*r),],nr=p,byrow = TRUE)
        rownames(a.hat) = paste0("x",1:p)
        colnames(a.hat) = paste0("r",1:r)
        rownames(b.hat) = paste0("x",1:p)
        colnames(b.hat) = paste0("r",1:r)
      }
      P.hat = out[2*r+1,]
      O.hat = out[(2*r+2):(2*r+1+p),]
      test  = out[(2*r+2+p),]
      eigs  = out[(2*r+2+p+1),]
      res   = out[(2*r+2+p+2):(2*r+2+p+N+1),]

      names(P.hat)    = paste0("x",1:p)
      rownames(O.hat) = paste0("x",1:p)
      colnames(O.hat) = paste0("x",1:p)
      names(test)     = paste0("r=",0:(p-1))
      names(eigs)     = paste0("l",1:p)
      colnames(res)   = paste0("x",1:p)

      if(r==0){
        a.hat = b.hat = NULL
      }
      return(list(N=N, p=p, r=r ,alpha = a.hat, beta = b.hat, Psi=P.hat, Omega = O.hat, test=test, lambda=eigs, A=A, H=H, df=df, dt=dt, res=res, data=X))
    }


# Wrapper function for the boostrap cpp function.
    bootstrap <-function(X,B=1000, dt =0.1){

      X = as.matrix(X)
      N = nrow(X)-1
      p = ncol(X)
      r = 1
      tmp = johansenCpp(X,r,diag(p),diag(p),dt)
      out = list()

      out$test = tmp[(2*r+2+p),]
      out$boot = bootstrapCpp(X,B=B,dt)
      return(out)
    }

# Wrapper function for the NEW boostrap cpp function.
    bootstrapOld <-function(X,B=1000, dt =0.1){

      X = as.matrix(X)
      N = nrow(X)-1
      p = ncol(X)
      r = 1
      tmp = johansenCpp(X,r,diag(p),diag(p),dt)
      out = list()

      out$test = tmp[(2*r+2+p),]
      out$boot = bootstrapCppOld(X,B=B,dt)
      return(out)
    }
# bootstrap <-function(X,r=1,B=1000, dt =0.1){
#   X = as.matrix(X)
#   N = nrow(X)-1
#   p = ncol(X)
#
#   out = bootstrapCpp(X,r,B,dt)
#
#   if(r > 0){
#     a.hat = matrix(out[1:r,],nr=p,byrow = TRUE)
#     b.hat = matrix(out[(r+1):(2*r),],nr=p,byrow = TRUE)
#     rownames(a.hat) = paste0("x",1:p)
#     colnames(a.hat) = paste0("r",1:r)
#     rownames(b.hat) = paste0("x",1:p)
#     colnames(b.hat) = paste0("r",1:r)
#   }
#
#   P.hat = out[2*r+1,]
#   O.hat = out[(2*r+2):(2*r+1+p),]
#   test  = out[(2*r+2+p),]
#   eigs  = out[(2*r+2+p+1),]
#   res   = out[(2*r+2+p+2):(2*r+2+p+N+1),]
#
#   names(P.hat)    = paste0("x",1:p)
#   rownames(O.hat) = paste0("x",1:p)
#   colnames(O.hat) = paste0("x",1:p)
#   names(test)     = paste0("r=",0:(p-1))
#   names(eigs)     = paste0("l",1:p)
#   colnames(res)   = paste0("x",1:p)
#
#   if(B > 0){
#     boot =  out[(2*r+2+p+N+2):(2*r+2+p+N+B+1-1),]
#     colnames(boot) = paste("r=",0:(p-1), sep="")
#   }
#
#   if(r>0){
#     if(B== 0) return(list(alpha = a.hat, beta = b.hat, Psi=P.hat, Omega = O.hat, test=test, lambda=eigs, res=res))
#     if(B > 0) return(list(alpha = a.hat, beta = b.hat, Psi=P.hat, Omega = O.hat, test=test, lambda=eigs, res=res, boot=boot))
#   } else{
#     if(B== 0) return(list(Psi=P.hat, Omega = O.hat, test=test, lambda = eigs, res=res))
#     if(B > 0) return(list(Psi=P.hat, Omega = O.hat, test=test, lambda = eigs, res=res, boot=boot))
#   }
# }


# Wrapper function for the oscillator function
# oscillator <- function(alpha,beta,phi.b,phi.sig,x0,y0,nSim,dt,r.a,r.b,r.sig){
#   nOsc = nrow(alpha)
#   alpha = as.matrix(alpha)
#   beta = as.matrix(beta)
#   phi.b = as.matrix(phi.b)
#   phi.sig = as.matrix(phi.sig)
#   x0 = as.matrix(x0)
#   y0 = as.matrix(y0)
#   r.a = as.matrix(r.a)
#   r.b = as.matrix(r.b)
#   r.sig = as.matrix(r.sig)
#   Z = linearOscillator(alpha,beta,phi.b,phi.sig,x0,y0,nSim,dt,r.a,r.b,r.sig)
#   Z = as.data.frame(Z)
#   colnames(Z) = c("t",paste(c("x","y"),c(rbind(1:nOsc,1:nOsc)),sep=""),paste("phi",1:nOsc,sep=""),paste("r",1:nOsc,sep=""))
#   #colnames(Z) = c("dt",paste(c("x","y"),c(rbind(1:nOsc,1:nOsc)),sep=""),paste("phi",1:nOsc,sep=""),paste("r",1:nOsc,sep=""),paste(c("x","y"),c(rbind(1:nOsc,1:nOsc)),"_1",sep=""),paste("phi",1:nOsc,"_1",sep=""),paste("coInt",1:nOsc,sep=""),paste("coInt",1:nOsc,"_1",sep=""),paste("coInt",1:nOsc,"_2",sep=""),paste("rot",1:nOsc,sep=""))
#   invisible(Z)
# }

