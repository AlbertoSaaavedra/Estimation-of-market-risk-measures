# Modelling GARCH innovations

# The function compara_modelos fits a GPD model to the innovations of a GARCH
# model and compares the quality of the obtained fit against that given by
# using Normal and Student-t models. 
# The compara_modelos function has the parameters:
#      res: the innovations of the adjusted GARCH; and
#  umbral: either expressed as a number (in the threshold variable) or
#          as a function of the amount of excesses to it (in the 
#          extreme variable).
# 
# There are two options displayed by the function compare_models: 
#   i) GPD fit based on the given threshold; and
#  ii) Comparison GPD vs Normal vs Student-t.
# 
# Selecting the first option displays an additional menu, which offers 
# four plots to evaluate the goodness of fit of the GPD model to the data. 
# These options are given by default by the plot.gpd function of the evir 
# library. In our work we use only the first two options of this menu.
# 
# On the other hand, the second option of the function compara_modelos 
# allows us to make a plot that compares the quality of fit of the GPD 
# model to the (right) tail of the analyzed innovations against that
# given by a the Normal and a Student-t models fited via maximum likelihood.

compara_modelos<-function(res,umbral=NA,extremos) {
  library(evir) 
  res <- as.numeric(res)
  #Verificando si el umbral se indica numéricamente
  if(is.na(umbral)) 
    gpd<-gpd(res,nextremes = extremos) #Ajustando la DPG
  #Verificando si el umbral se determina en función de la
  #cantidad de observaciones que lo exceden
  else
    gpd<-gpd(res,threshold = umbral) #Ajustando la DPG
  #Enumerando las opciones de la función
  choices <- c("Ajuste de DPG a con base en el umbral dado", 
               "Comparación DPG vs Normal vs t-Student") 
  tmenu <- paste("Opción:", choices)
  pick <- 1
  lastcurve <- NULL
  while (pick > 0) {
    pick <- menu(tmenu, title = "\nSelecciona una opción (o 0 
                                para salir):")
    if (pick == 1) {
      out1<-gpd(res[res>0],threshold=umbral,
                nextremes=extremos)
      plot.gpd(out1) }
    if (pick == 2) {
      #Estimación MV de parámetros del modelo normal
      par.norm<-fitdistr(res,"normal")$estimate 
      #Estimación MV de parámetros del modelo t-Student
      par.t<-fitdistr(res,"t")$estimate 
      par(mfrow=c(1,1))
      tailplot(gpd,main="Innovaciones GARCH: 
                        comparación de modelos")
      curve(1-pnorm(x,par.norm['mean'],par.norm['sd']),
            from=0.6089704,to=30,lty="dashed",
            col="red",add=TRUE)
      curve(1-pt(x,par.t["df"]),from=0.6089704,
            to=30,lty="dotted",col="blue",add=TRUE)
      legend("topright", 
        c("Modelo GPD","Modelo Normal","Modelo t-Student"),
        col=c("black","red","blue"),
        lty=c("solid","dashed","dotted"),bty="n",cex=0.8) } }
  invisible(lastcurve) }