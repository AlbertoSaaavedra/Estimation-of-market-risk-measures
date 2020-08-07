# Exploring GARCH innovations with an extreme value theory approach
# 
# The explora_evt function explores the innovations of a GARCH model with an extreme 
# value theory approach. The explora_evt function has the parameters:
#    res: the standardized residuals of the adjusted GARCH model; and
#  alpha: variable that indicates the level of the quantile to be plotted
#         in the exploratory analysis.
# 
# When the explora_evt function is run, a menu is displayed to the user, 
# their options are:
#   i) graph the sample mean excess function;
#  ii) plot the variability of the parameter xi when varying the threshold; and
# iii) graph the variablility, before changes in the threshold, of the 
#      alpha-quantile of a DPG model fitted to innovations.

explora_evt<-function(res,alfa) {
  # Library with functions for extreme value theory methods
  library(evir) 
  res<-as.numeric(res)
  # Numbering the choices given by the function
  choices<-c("Gráfica de la función muestral media de 
               exceso", 
              "Gráfica del parámetro de forma de la DPG al 
               variar el umbral", 
              "Gráfica del cuantil alfa al variar el umbral")
  tmenu<-paste("Opción:", choices) 
  pick<-1 
  lastcurve<-NULL
  while(pick > 0) {
    pick<-menu(tmenu, title = "\nSelecciona una opción (o 0 
                                 para salir):")
    if(pick == 1)
        meplot(res[res>0])
    if(pick == 2)
        shape(res[res>0])
    if(pick == 3) 
        quant(res[res>0], p=alfa) }
  invisible(lastcurve) }