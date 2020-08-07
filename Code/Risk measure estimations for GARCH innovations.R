# Risk measure estimations for GARCH innovations

# Our application requires that we estimate risk measures in two stages: first, 
# for the innovations of the fitted GARCH and, later, for the time series studied. 
# Therefore, two routines are presented to estimate these risk measures, 
# one for each mentioned stage.
# 
# The routine created to estimate risk measures for the innovations of a 
# fitted GARCH model consists of the following three functions.

# Expected Shortfall for Normal model

# The ES.norm function estimates the Expected Shortfall of the Normal distribution. 
# Its parameters are:
#   alpha: the confidence level at which you want to estimate the Expected
#          Shortfall.
#      mu: the mean of the Normal model for which the Expected Shortfall is 
#          being estimated.
#   sigma: the standard deviation of the Normal model for which the Expected
#          Shortfall is being estimated.
# 
# The output of the ES.norm function is the expected shortfall for a Normal
# model, estimated at an alpha confidence level. The function can compute
# the Expected Shortfall for different values of alpha in one call.


ES.norm<-function(alpha, mu, sigma) {
  n<-length(alpha)
  ES<-c()
  for(i in 1:n) {
    ES[i]<-mu+(sigma*(dnorm(qnorm(alpha[i],mu,sigma),
                            mu,sigma)/(1-alpha[i]))) }
  return(ES) }

# Expected Shortfall for Student-t model

# The ES.t function estimates the Expected Shortfall of the Student-t distribution. 
# Its parameters are:
#   alpha: The confidence level at which you want to estimate the 
#          Expected Shortfall.
#      df: The degrees of freedom of the Student-t model for which 
#          the Expected Shortfall is being estimated.
# 
# The output of the ES.t function is the expected shortfall for a Student-t model, 
# estimated at an alpha confidence level. The function can compute the Expected 
# Shortfall for different values of alpha in one call.

ES.t<-function(alpha, df) {
  n<-length(alpha)
  ES<-c()
  for(i in 1:n) {
    ES[i]<-(dt(qt(alpha[i],df),df)/(1-alpha[i])) *
                         ((df+(qt(alpha[i],df))^2)/(df-1)) }
  return(ES) }

# Innovations risk measure estimations - Synthesis

# To estimate the VaR Expected Shortfall of a fitted GPD model we use the function riskmeasures
# of the  evir library. The quantiles of the Normal and Student-t distributions can be obtained
# using R's functions qnorm and qt, respectively.

# The above functions already allow us to calculate our two risk measures, for an innovation 
# series for any of the three innovation models considered in our work. However, in order to
# synthesize these estimations in a single procedure, we created the function risk_innov, 
# which is the last function of our routine for innovations.


# The risk_innov function has four parameters: 
#       res: the innovations of a GARCH model; 
#   measure: character string indicating which risk measure you want to estimate, "VaR" or "ES"; 
#     alpha: confidence level at which the risk measures are to be estimated. This can be a
#            numerical vector containing values between 0 and 1
#    umbral: the threshold on which the fitting of the GPD model is based. This can be expressed
#            numerically or based on the amount of data that exceeds the threshold.
# 
# The function risk_innov returns a matrix The rows of this matrix correspond to the models from 
# which the VaR and Expected Shortfall are estimated. Its columns correspond to the confidence 
# levels at which it is requested to estimate these risk measures. 


riesgo_innov<-function(res,medida,alfa,umbral=NA,extremos) {
  res<-as.numeric(res)
  par.norm<-fitdistr(res,"normal")$estimate 
  par.t<-fitdistr(res,"t")$estimate 
  if(is.na(umbral)) 
    gpd<-gpd(res,nextremes = extremos)
  else 
    gpd<-gpd(res,threshold = umbral)
  if(medida == "VaR") {           
      return(matrix(c(qnorm(alfa, par.norm["mean"],
                            par.norm["sd"]),
                      qt(alfa, par.t["df"]),
                      riskmeasures(gpd,alfa)[,2]),
                    nrow=3,col=length(alfa),
                    dimnames=list(c("VaR, modelo Normal",
                                    "VaR, modelo t-Student",
                                    "VaR, modelo DPG"),
                    as.character(alfa)),byrow=TRUE)) }
  if(medida == "ES") {
      return(matrix(c(ES.norm(alfa, par.norm["mean"],
                              par.norm["sd"]),
                      ES.t(alfa, par.t["df"]),
                      riskmeasures(gpd,alfa)[,3]),
                      nrow=3, ncol=length(alfa),
                      dimnames = list(c("ES, modelo Normal",
                                      "ES, modelo t-Student",
                                      "ES, modelo DPG"),
                     as.character(alfa)),byrow=TRUE)) } }