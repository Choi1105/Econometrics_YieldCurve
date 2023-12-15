## DGP and Application of Space State Model(Code with Data) ##
## ※ Notice ※
___
+ **Textbook**
  + [State-Space Models with Regime Switching -Classical and Gibbs-Sampling Approaches with Applications](https://mitpress.mit.edu/9780262535502/state-space-models-with-regime-switching)
  + By Chang-Jin Kim and Daniel C R. Halbert 
  + Publisher: The MIT Press

+ **Code**
  + This code is based on the work of Professor [Young-min Kim](https://kimymecon.weebly.com/matlab.html), an assistant professor at Jeonbuk National University.

+ **Data**
  + Some data were referred to Professor [Chang-jin Kim's data](http://econ.korea.ac.kr/~cjkim/), and others were obtained from [KOSIS](https://kosis.kr/statHtml/statHtml.do?orgId=214&tblId=DT_214N_Z01900&vw_cd=MT_ZTITLE&list_id=214_21404&seqNo=&lang_mode=ko&language=kor&obj_var_id=&itm_id=&conn_path=MT_ZTITLE).

+ **Reference**
  + [Diebold, F. X., and C. L. Li, 2006, Forecasting the term structure of government bond yields,Journal of Econometrics130, pp. 337–364.](https://www.sas.upenn.edu/~fdiebold/papers/paper49/Diebold-Li.pdf)
  + [Kang, K. H., 2012a, Forecasting the term structure of Korean Government bond tields usingthe dynamic Nelson-Siegel class models,Asia-Pacific Journal of Financial Studies41, pp.765–787.](https://onlinelibrary.wiley.com/doi/full/10.1111/ajfs.12000)
  + [Clark, P., 1987,The cyclical component of U.S. economic activity, Quarterly Journal of Economics 102, 797-814](https://www.jstor.org/stable/1884282)
  + [Nelson, C. R., & Plosser, C. R. (1982). Trends and random walks in macroeconomic time series: Some evidence and implications. Journal of Monetary Economics, 10(2), 139–162.](https://www.sciencedirect.com/science/article/abs/pii/0304393282900125)
  + [Burmeister, Edwin, and Kent Wall, "Kalman Filtering Estimation of Unobserved Rational Expectations with an Application to the German Hyperinflation," Journal of Econometrics, November 1982, 20, 255–84.](https://www.sciencedirect.com/science/article/abs/pii/0304407682900215)
  + [Stock, J.H., Watson, M.W., 1991. A probability model of the coincident economic indicators. In: Moore, G., Lahiri, K. (Eds.), The Leading Economic Indicators: New Approaches and Forecasting Records. Cambridge
University Press, Cambridge, pp. 63–89](https://www.nber.org/papers/w2772)
  + [Kim, C.J., Nelson, C.R., 1989. The time-varying parameter model for modeling
changing conditional variance: the case of the Lucas hypothesis. J. Bus. Econ.
Stat. 7 (4), 433–440](https://www.tandfonline.com/doi/abs/10.1080/07350015.1989.10509755)

+ **The meaning of each DGP, DGP_Practice, Practice is as follows.**
  + DGP : Data Generating Process (Basic Code)
  + DGP_Practice : Data Generating Process for models to be run in Practice (Modified code for Application)
  + Practice : Modified code and Real Data(Include result_plot)
  
___

## ※ Index ※
___

## Regime Switching Models
+ **Serially Uncorrelated Data and Switching.**
  + Uncorrelated Data and Switching DGP
    + Beta_1_Sig_2
    + Beta_2_Sig_1
    + Beta_2_Sig_2
    + Beta_2_Sig_2_OLS
   
  + Independent Switching DGP
    + Independent Switching Non Exogenous Variable
    + Independent Switching With Exogenous Variable
    + Independent Switching Non Exogenous Variable With Paramconst
      
  + Markov Switching DGP
    + The case of Markov Switching
   
      
+ **Serially Correlated Data and Switching.**
  + DGP
    + Dealing with the Problem of Unobserved $S_t$

+ **Structure Break.**
  + DGP
    + 2Structure Break
    + 3Structure Break
    + Expected Duration of a Regime in a Markov-Switching Model
      
+ **Kim's Smoothing Algorithm.**
  + DGP
    + Markov-Switching Dependent Model with Kim's Smoothing
    + Markov-Switching Serially Correlated Model with Kim's Smoothing
___
## ※ Summary ※

## **Hamilton's Markov Switching Replication.**
### Model 1.
___
