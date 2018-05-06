# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

function_CCD <- function(n,mu,sigma,m,mm,division,alpha){


  Generate_data <- function(n,mu,sigma){
    data <- list()
    for(ii in 1:length(n)){
      data[[ii]] <- rnorm(n[ii],mu[ii],sigma[ii])
    }

    return(data)

  }

  Generate_sample_value <- function(n,data_value) {

    s_mean <- NULL
    s_variance <- NULL
    for(ii in 1:length(n)){
      s_mean[ii] <- mean(data_value[[ii]])
      s_variance[ii] <- var(data_value[[ii]])
    }
    return(list(s_mean,s_variance))
  }


  Generate_R_theta <- function(n,xx,ss,m){

    u <- rchisq(m,n-1)
    z <- rnorm(m,0,1)
    v <- rchisq(m,n-1)
    R_theta <- xx - (z/ (sqrt(u)/ sqrt(n-1))) * sqrt(ss)/ sqrt(n) + ss /( 2* v /(n-1) )

    return(R_theta)
  }

  Calculate_R_all <- function(n,xx,ss,m){
    k <- length(n)
    R_all <- array(dim = c(k,m))
    for(ii in 1:k){
      R_all[ii,] <- Generate_R_theta(n=n[ii],xx=xx[ii],ss=ss[ii],m)
    }
    return(R_all)
  }

  Calculate_H_i <- function(theta, R_theta){
    H_i <- sum (R_theta <= theta)
    b <- length(R_theta)
    return(H_i/b)
  }



  Calculate_Hc_value <- function(n,R_all,theta){
    H_prob <- array(dim=c(length(n)))
    for(ii in 1:length(n)){
      H_prob[ii] <- Calculate_H_i(theta=theta,R_theta = R_all[ii,])
    }
    Hc_value <- (-2)*sum(log (1-H_prob))
    return(Hc_value)
  }

  Calculate_theta_all <- function(R_all,division ){

    theta_all <- array(dim = c(division))
    a0 <-  max(apply(R_all,1,min))
    leftp <- 0              #min(-5,a0)   #gai
    b0 <-  min(apply(R_all,1,max))
    rightp <- 10             #min(50,b0)
    gaps <- (rightp-leftp)/division
    for(jj in 1:division){
      theta_all[jj] <- leftp + jj * gaps
    }
    return(theta_all)
  }

  Calculate_Hc_all <- function(n,theta_all,R_all){
    Hc_all <- array(dim=c(length(theta_all)))
    for(ii in 1:length(theta_all)){
      Hc_all[ii] <- Calculate_Hc_value(n,R_all,theta = theta_all[ii])
    }
    return(Hc_all)
  }


  Calculate_interval <- function(n,Hc_all,theta_all,alpha,mu,sigma){
    FC <- pchisq(Hc_all,2*length(n))
    FC1 <- sort(FC)
    FC2 <- order(FC)

    number1 <- sum(FC1 <= alpha/2)
    if(abs(FC1[number1] - alpha/2) >= abs(FC1[number1+1]) - alpha/2){
      number1 <- number1 + 1
    }

    number2 <- sum(FC1 <= 1-alpha/2 )
    if(abs(FC1[number2]-(1-alpha/2)) >= abs(FC1[number2+1])){
      number2 <- number2 + 1
    }

    left <- theta_all[FC2[number1]]
    right <- theta_all[FC2[number2]]

    res <- c(mu[1]+(sigma[1]^2)/2 , left , right)
    return(res)
  }


  count <- 0
  Changdu <- 0
  proof <- array(dim = c(mm,3))




  for(rep in 1:mm){
    data_value <- Generate_data(n,mu,sigma)
    sample_value <- Generate_sample_value(n,data_value)
    xx <- sample_value[[1]]
    ss <- sample_value[[2]]

    R_all <- Calculate_R_all(n,xx,ss,m)
    theta_all <- Calculate_theta_all(R_all,division)
    Hc_all <- Calculate_Hc_all(n,theta_all,R_all)
    interval <- Calculate_interval(n,Hc_all,theta_all,alpha,mu,sigma)
    proof[rep,] <- interval
    Changdu <- Changdu + (interval[3]-interval[2])
    if(interval[1]>=interval[2] && interval[1]<=interval[3]){
      count <- count +1
    }
  }
  return(list(Empirical_coverage=count/mm,Interval_length=Changdu/rep))

}
