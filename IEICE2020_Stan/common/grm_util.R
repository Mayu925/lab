convert_category_estimates <- function(category_prm, N, K){
  for(n in 1:N){
    prm = category_prm[((n-1)*(K-1)+1):((n-1)*(K-1)+(K-1))]
    if(n == 1){
      mat = t(data.frame(prm))
    } else {
      mat = rbind(mat, t(data.frame(prm)))
    }
  }  
  return(mat)
}

prob <-function(param, k, x){
  category_prm <- param$category_prm
  if (k == 1) {
    return( 1.0 - bcc(param, category_prm[1], x) )
  } else if(k == length(category_prm)+1){
    return( bcc(param, category_prm[k-1], x) )
  } else {
    return( bcc(param, category_prm[k-1], x) - bcc(param, category_prm[k], x) )
  }
}

fisher_information_no_alpha_multiplication <-function(param, theta) {
  category_prm <- param$category_prm
  p <- bcc(param, category_prm[1], theta)
  Z <- (- p * (1-p))^2 / (1.0 - p)
  for (k in 2:(length(category_prm))) {
    p <- bcc(param, category_prm[k], theta)
    p_1 <- bcc(param, category_prm[k-1], theta)
    Z <- Z +  (p_1 * (1-p_1) - p * (1-p))^2 / (p_1 - p)
  }
  p <- bcc(param, category_prm[k-1], theta)
  Z <- Z + (p * (1-p))^2 / p
  return(Z)
}

gen_category_param <- function(N, K){
  prior_sd <- 1.0
  prior_mean <- seq(-0.5, 0.5, length=(K-1))
  category <- matrix(0, nrow=N, ncol=(K-1))
  for (k in 1:(K-1)){
    category[,k] <- rnorm(N, prior_mean[k], prior_sd)
  }
  for (i in 1:N){
    category[i,] = sort(category[i,])
  } 
  return(category)
}

draw_grm_rubric <- function(param, K, def){
  draw_icc(param,K, def)
#  par(new=T)
#  draw_se(param, K, def)
  draw_legend(def)
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

draw_se <- function(param, K, def){
  curve(standard_error(param, x),
        from=-def$xlim, to=def$xlim, ylim=c(0,def$SeYlim), xlab="", ylab="",  
        col=def$color[K+1],lty=def$style[K+1], cex.lab=def$lcex, lwd=def$llwd, axes=FALSE)
  axis(4, cex.axis = def$caxcx)
}

draw_legend <- function(def){
  legend(def$legendPosition, cex = def$axcx, lwd=def$llwd, legend=def$legends, col=def$color, lty=def$style, bg = "white", ncol = def$legendNcol, pch=def$pchs)
}
