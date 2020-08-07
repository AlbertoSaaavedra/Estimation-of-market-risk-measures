# Estimation of market risk measures

### Contents

This repository contains some elements of the work I carried out during the development of my BSc Thesis at the National Autonomous University of Mexico (UNAM). Specifically, in this repository you can fin:
- Sample codes of the programming I did to estimate and backtest conditional Value-at-Risk and Expected Shortfall on Mexican financial time series at very high confidence levels.
- Sample of the data I used for this research, which is freely available at the Central Bank of Mexico's webpage:
  - For the USD/MXN exchange rate, [here](https://www.banxico.org.mx/SieInternet/consultarDirectorioInternetAction.do?sector=6&accion=consultarCuadro&idCuadro=CF102&locale=es).
  - For the index of the Mexican Stock Exchange, [here](https://www.banxico.org.mx/SieInternet/consultarDirectorioInternetAction.do?sector=7&accion=consultarCuadroAnalitico&idCuadro=CA54&locale=es).
- A paper I published at the Mexican Journal of Economics and Finance with main findings on my research on this topic.

### Abstract of my BSc Thesis

The objectives of this work are to investigate whether: i) a GARCH model with Generalized Pareto Distribution (GPD) innovations, complemented with an EWMA volatility forecast in order to consider practical problems that might arise in GARCH applications that comprise long periods of time, appropriately estimate risk measures (VaR and Expected Shortfall) for Mexican financial series, at high confidence levels; ii) the estimates yielded by such model are better than those given by a GARCH with Gaussian or Student-t innovations. Our quality assessment and comparison between models consist of backtests of the risk measures estimates yielded by each method used in this paper. Our results show that: i) the methodology used in this paper appropriately estimates our two risk measures; ii) the GARCH-GPD model yields better results than the GARCH-Gaussian and GARCH-t-Student models. Our results are limited to one-day risk measures estimates. As far as we know, our results on the Expected Shortfall are the first of its kind for Mexican series. We conclude that the study achieved its objectives and there are important areas of opportunity for further studies
