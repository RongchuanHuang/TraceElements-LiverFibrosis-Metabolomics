library(MASS)


####################################################################################################################################
###### First Order PDM
####################################################################################################################################

PDM_1 =function(x,y,m,imax, tol, theta=theta, w1=w1, interval = interval, step = step)  {
  
  I=diag(length(w1))
  
  f=function(lambda,psi,phi){
    crossprod(solve(lambda*I + psi) %*% phi) - 1  
  }
  
  #negative log
  g11=function(x,y,m,theta,w1,s1,s2){
    crossprod(y-theta[3]-x*theta[5]-m%*%w1*theta[4]) / (1/s1^2)
  }
  
  g12=function(x,y,m,theta,w1,s1,s2){
    crossprod(m%*%w1-theta[1]-x*theta[2]) / (1/s2^2)
  }
  
  count=0
  
  for (i in 1:imax) {
    
    mw1=m%*%w1
    J=rep(1,nrow(m))
    X1=cbind(J,x,mw1)
    X2=cbind(J,x)
    
    s1=sqrt( crossprod( y-X1%*%solve( crossprod(X1) )%*% crossprod(X1, y) )  / ( nrow(m)-3 ) )
    
    s2=sqrt(  crossprod(mw1-X2%*%solve(  crossprod(X2) )%*% crossprod(X2, mw1)) / (nrow(m)-2) )
    
    psi=as.numeric(1/s1^2)*crossprod(m)*(theta[4])^2
    +as.numeric(1/s2^2)*t(m)%*%m
    
    phi=as.numeric(1/s1^2)*crossprod(m, y-theta[3]-x*theta[5])*theta[4]
    +as.numeric(1/s2^2)*crossprod(m, theta[1] + x*theta[2])
    
    
    ################################################################
    ################################################################
    
    
    grid = round(2*interval/step)
    
    
    temp=numeric( step + 1 )
    
    
    for (i in 1 : (step + 1) ) {
      
      temp[i]=f(- (interval + grid ) + grid*i, psi = psi, phi = phi)
    }
    
    
    ################################################################
    ################################################################
    
    
    for (i in 1: step ){
      
      
      if ((temp[i]*temp[i+1])<0){
        l= - (interval + grid ) + grid*i
        h= - (interval + grid ) + grid*(i+1)
        
      }
    }
    
    
    tryCatch({ 
      lambda = uniroot(f, c(l, h), phi = phi, psi = psi) $ root 
    }, error=function(e){}) 
    
    
    
    
    
    
    
    w1_new=solve(lambda*I+psi)%*%phi
    
    theta_1=optim(theta,g11,x=x,y=y,m=m,s1=s1,s2=s2,w1=w1_new, method="BFGS")$par
    theta_2=optim(theta,g12,x=x,y=y,m=m,s1=s1,s2=s2,w1=w1_new, method="BFGS")$par
    theta_new=c(theta_2[1],theta_2[2],theta_1[3],theta_1[4],theta_1[5])
    
    #L = vector(length=imax)
    L = - crossprod(y-theta_new[3]-x*theta_new[5]-m%*%w1*theta_new[4]) / (2*s1^2)
    - crossprod(m%*%w1-theta_new[1]-x*theta_new[2]) / (2*s2^2)
    
    #     print(L[i])
    #     
    
    print( paste( "Calculating the First PDM:", "log-likelihood is", L, sep = " "))
    print( paste( "lambda", "=", lambda, sep = " "))
    print( paste( "norm", "=", t(w1_new)%*%w1_new, sep = " "))
    
    if (abs( crossprod(theta_new-theta) )<tol && abs( crossprod(w1_new - w1) )<tol) break
    
    else {theta=theta_new
          w1=w1_new}
    
    count=count + 1
  }
  
  
  out1=list(w1,theta,lambda)
  names(out1) <- c("w1","theta","lambda")
  out1
}


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
###### End of the First Order PDM
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

PDM_2 =function(x, y, m, imax, tol_1, tol_2, tol_3, tol_4, theta_1 = theta_1, theta_2 = theta_2, w1 = w1, w2=w2, 
                interval_1 = interval_1, step_1 = step_1, interval_2 = interval_2, step_2 = step_2
)  {
  
  PDM_1 = PDM_1(x, y, m, imax, tol = tol_1, theta = theta_1, w1 = w1, interval = interval_1, step= step_1)
  
  w1_new = PDM_1$w1
  theta_1_new = PDM_1$theta
  
  I=diag(length(w2))
  
  H = I - tcrossprod(w1_new)
  
  
  f=function(lambda, psi_2, phi_2){
    t( ginv( lambda* crossprod(H) + psi_2 ) %*% phi_2  ) %*% crossprod(H) %*% (  ginv( lambda* crossprod(H) +  psi_2 ) %*% phi_2  ) - 1  
  }
  
  
  
  #negative log
  g11 = function(x, y, m, theta_2, w1, w2, s1, s2 ){
    crossprod( y - theta_2[3] - x*theta_2[6] - m%*%w1*theta_2[4] -  m%*%w2*theta_2[5]) / (1/s1^2)
  }
  
  g12 = function(x, y, m, theta_2, w2, s1, s2){
    crossprod(m%*%w2 - theta_2[1] - x*theta_2[2]) / (1/s2^2)
  }
  
  
  
  count=0
  
  for (i in 1:imax) {
    
    
    mw1 = m%*%w1_new
    mw2 = m%*%w2
    J = rep(1,nrow(m))
    X1 = cbind(J, x, mw1, mw2)
    X2 = cbind(J, x)
    
    s1 = sqrt( crossprod(y - X1%*%ginv(crossprod(X1))%*% crossprod(X1, y)) / (nrow(m) - 4) )
    s2 = sqrt( crossprod(mw2-X2%*%ginv( crossprod(X2) )%*% crossprod(X2, mw2) ) / (nrow(m) - 3) )
    
    
    
    
    psi_2 =
      as.numeric(1/s1^2)* H%*% crossprod(m) %*% H * (theta_2[5])^2
    + as.numeric(1/s2^2)* H%*% crossprod(m) %*% H
    
    
    
    phi_2 =
      as.numeric(1/s1^2)* H%*% crossprod(m, y - theta_2[3] - x*theta_2[6] - mw1*theta_2[4]) *theta_2[5]
    + as.numeric(1/s2^2)*H %*% crossprod(m, theta_2[1] + x*theta_2[2] )
    
    
    
    ################################################################
    ################################################################
    
    
    grid = round(2*interval_2/step_2)
    
    
    temp=numeric( step_2 + 1 )
    
    
    for (i in 1 : (step_2 + 1) ) {
      
      temp[i]=f(- (interval_2 + grid ) + grid*i, psi_2 = psi_2, phi_2 = phi_2)
    }
    
    
    ################################################################
    ################################################################
    
    
    for (i in 1: step_2 ){
      
      
      if ((temp[i]*temp[i+1])<0){
        l= - (interval_2 + grid ) + grid*i
        h= - (interval_2 + grid ) + grid*(i+1)
        
      }
    }
    
    
    
    tryCatch({ 
      lambda = uniroot(f, c(l, h), phi_2 = phi_2, psi_2 = psi_2) $ root
    }, error=function(e){}) 
    
    
    w2_new=ginv( lambda* crossprod(H) + psi_2 ) %*% phi_2 
    
    theta_1_2 = optim(theta_2, g11, x=x, y=y, m=m, s1=s1, s2=s2, w1=w1_new, w2=w2_new, method="BFGS")$par
    theta_2_2 = optim(theta_2, g12, x=x, y=y, m=m, s1=s1, s2=s2, w2=w2_new, method="BFGS")$par
    theta_2_new = c(theta_2_2[1], theta_2_2[2], theta_1_2[3], theta_1_2[4], theta_1_2[5], theta_1_2[6])
    

    L2= - crossprod( y - theta_2[3] - x*theta_2[6] - m%*%w1_new*theta_2[4] -  m%*%w2_new*theta_2[5]) / (1/s1^2)
    
    - crossprod(m%*%w2_new - theta_2[1] - x*theta_2[2]) / (1/s2^2)
    
    
    print (paste( "Calculating the Second PDM:", "log-likelihood is", L2, sep = " "))
    print( paste( "lambda", "=", lambda, sep = " "))
    print( paste( "orthogonality", "=", t(w2_new)%*%w1_new, sep = " "))
    print( paste( "theta_2_new", "=", theta_2_new, sep = " "))    
    
    if ( abs( crossprod(theta_2_new[3:6] - theta_2 [3:6]) ) < tol_2 
         && abs( crossprod(theta_2_new[1] - theta_2_new[1]) )<tol_4
         && abs( crossprod(theta_2_new[2] - theta_2_new[2]) )<tol_3
         && abs( crossprod(w2_new - w2) )<tol_2 ) break
    
    else {theta_2 = theta_2_new
          w2=w2_new}
    
    count=count + 1
  }
  
  w2 = w2/as.numeric(sqrt(crossprod(w2)))
  
  out1=list(w1_new, w2, theta_1_new, theta_2, lambda)
  names(out1) <- c("w1","w2", "theta_1", "theta_2", "lambda")
  out1
}


