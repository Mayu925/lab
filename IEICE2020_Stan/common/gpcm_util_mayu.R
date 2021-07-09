convert_category_estimates <- function(category_prm, N, K){
  for(n in 1:N){
    prm = category_prm[((n-1)*(K-2)+1):((n-1)*(K-2)+(K-2))]
    prm = append(prm, -1*sum(prm))
    if(n == 1){
      mat = t(data.frame(prm))
    } else {
      mat = rbind(mat, t(data.frame(prm)))
    }
  }  
  return(mat)
}

convert_alpha_rt <- function(alpha_rt, T, R){
  for(r in 1:R){
    alrt = alpha_rt[((r-1)*T+1):((r-1)*T+T)]
    if(r == 1){
      mat = t(data.frame(alrt))
    }else{
      mat = rbind(mat,t(data.frame(alrt)))
    }
  }
  return(mat)
}

prob <-function(param, k, x){
  all_sum <- 0
  tmp <- 0
 
    tmp <- tmp + 1.7 * logit(param, param$category_prm, x)
   
      target <- exp(tmp)
      all_sum <- all_sum + exp(tmp)

    
  
  return(target/all_sum)
}

gen_category_param <- function(N, K){
  category <- matrix(0, nrow=N, ncol=(K-1))
  for (k in 1:(K-1)){
    category[,k] <- rnorm(N, 0, 1)
  }
  for (i in 1:N){
#    category[i,] = sort(category[i,])
    category[i,] = category[i,] - mean(category[i,])
  } 
  return(category)
}

fisher_information_no_alpha_multiplication <- function(param, K, theta) {
  for (k in 0:(K-1)) {
    z <- prob(param, k+1, theta)
    if(k == 0){
      f1 <- k^2 * z
      f2 <- k * z
    } else {
      f1 = f1 + k^2 * z
      f2 = f2 + k * z
    }
  }
  return(f1 -  f2 ^ 2)
}

draw_icc <- function(param, K, def){
  x<- seq(-def$xlim, def$xlim, length=(10))
  curve(prob(param, 1, x), 
        from=-def$xlim, to=def$xlim, ylim=c(0,def$ylim), xlab=def$xlab, ylab=def$ylab, main=def$title, 
        col=def$color[1],lty=def$style[1], cex.lab=def$lcex, lwd=def$llwd, cex.axis=def$axcx, cex.main=def$maincx, yaxt= "n")
  par(new=T)
  plot(x, prob(param, 1, x), xlim=c(-def$xlim, def$xlim), ylim=c(0,def$ylim), xlab="", ylab="",
       col=def$color[1],lty=def$style[1], cex.lab=def$lcex, lwd=def$llwd, axes=FALSE, pch=def$pchs[1], cex=1.5,  yaxt= "n")
  for(k in 2:K){
    par(new=T)
    curve(prob(param, k, x), 
          from=-def$xlim, to=def$xlim, ylim=c(0,def$ylim), xlab="", ylab="",  
          col=def$color[k],lty=def$style[k], cex.lab=def$lcex, lwd=def$llwd, axes=FALSE,  yaxt= "n")
    par(new=T)
    plot(x, prob(param, k, x), xlim=c(-def$xlim, def$xlim), ylim=c(0,def$ylim), xlab="", ylab="",
         col=def$color[k],lty=def$style[k], cex.lab=def$lcex, lwd=def$llwd, axes=FALSE, pch=def$pchs[k], cex=1.5,  yaxt= "n")
  }
  axis(2, cex.axis = def$caxcx)
}
