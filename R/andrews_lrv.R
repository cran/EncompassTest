#' Long-run covariance estimation using Andrews quadratic spectral kernel.
#' 
#' Given a vector of residuals, it generates the Heteroskedastic Long run variance.
#'
#' @param e a vector of residual series, for which we recommend to use the recursive residuals from larger model.
#' @return  a vector of Long run variance using Andrews quadratic spectral kernel HAC.
#' @references Andrews, D. W. (1991). Heteroskedasticity and autocorrelation consistent covariance matrix estimation. Econometrica: Journal of the Econometric Society, 817-858.
#' @import pracma
#' @examples
#' set.seed(2014)
#' x<- rnorm(15);
#' #Andrews kernel HAC
#' andrews_lrcov = andrews_lrv(x)
#' @export



andrews_lrv=function(e){
  e=as.vector(e)
  t=length(e)
  
  a  =  pracma::mldivide(e[2:t],e[1:(t-1)])
  k = 0.8
  c_left = 1.1447*(((4*(a^2)*t)/(((1+a)^2)*((1-a)^2)))^(1/3))
  c_right = 1.1447*(((4*(k^2)*t)/(((1+k)^2)*((1-k)^2)))^(1/3))
  veco = c(c_left,c_right)
  l  =  min(veco)
  l = floor(l)
  lrv = t(e)%*%e/t
  
  #fix the empty element error
  if (l >= t){
    l = t
    
    for (i in 1:(l-1)){
      w = (1-(i/(1+l)))
      lrv = lrv + 2*t(e[1:(t-i)])%*%e[(1+i):t]*w/t
    }
    
  }else if (l > 0){
    for (i in 1:l){
      w = (1-(i/(1+l)))
      lrv = lrv + 2*t(e[1:(t-i)])%*%e[(1+i):t]*w/t
    }
  }
  
  
  
  return(lrv)
}
