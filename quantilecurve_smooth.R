quantilecurve_smooth<- function(qc, x, tau, x_out, tau_out = c(0.05, 1:9*0.1, 0.95), xmin = 4, xmax = 16, ymin = 0, ymax = 32){
  # -------------------------------------------------------------------
  #  @Function Function Name
  # 
  #  @Descriptions:  
  #
  # 
  #  @Parameters:
  #     qc    = matrix of fitted quantile curves -row - age; column - quantile 
  #     x     = age corresponding to row
  #     tau   = quantile levels corresponding to columns
  #     x_out = Age forfitted quantiles to output 
  #
  #  @Returns:  
  #     qc_long       = quantile curves in long format
  #     curve_plt     = ggplot objective of quantile curves
  #     q_v_out       = for the grid points of interest
  # --------------------------------------------------------------------
  
  
  qc= matrix(qc, nrow = dim(qc)[1], ncol = dim(qc)[2], 
             dimnames = list(x, tau))
  q_curve =reshape2::melt(qc)
  names(q_curve)= c('AGE','Quantile', 'AVAL')
  q_curve$q = as.numeric( q_curve$Quantile)
  q_curve$Quantile = paste0(q_curve$q*100%>%round, '%')
  
  q_curve$SmoothVAL = rep(0, length(q_curve$AGE))
  for(ii in unique(q_curve$Quantile)){
    id1 = which(q_curve$Quantile == ii)
    temp.set = q_curve[id1,]
    model1 = loess(temp.set$AVAL~temp.set$AGE)
    temp.vec = model1$fitted
    id2 = which(temp.vec < 0)
    if(length(id2) > 0){
      temp.vec [id2[1]:length(temp.vec)] = 0
    }
    q_curve$SmoothVAL[id1] = temp.vec
  }
  
  #NSAA quantile curves
  curve_plt = ggplot(data = q_curve, aes_string(x = 'AGE', y = 'SmoothVAL')) + 
    geom_path(aes_string(group = 'Quantile', color = 'Quantile')) +
    #geom_smooth(method = "loess",aes_string(group = 'Quantile', color = 'Quantile'))+
    ylab(label ='Smoothed Estimated Score') + xlab(label = 'Age')  + ylim(ymin, ymax) + xlim(xmin,xmax)+
    theme(
      plot.title = element_text(hjust = 0.5,size=30, face = "bold"),
      axis.title=element_text(size=30),
      axis.text=element_text(size=26),
      legend.text=element_text(size=15))
  # # age quantile plot by age
  # 
  # quantile_plt = ggplot(data = filter(q_curve, AGE%in%x_out)%>%mutate(AGE = factor(AGE)), 
  #                    aes_string(x = 'q', y = 'AVAL')) + 
  #   geom_path(aes_string(group = 'AGE', color = 'AGE')) 
  # 
  # output the fitted quantile at grid age points
  qc[qc<0] = 0
  q_v_out = qc[as.character(x_out), ]%>%apply(2, round, digit = 1)
  q_v_out_col = colnames(q_v_out)
  q_v_out = data.frame(q_v_out,  check.names = F )%>%tibble::rownames_to_column( "AGE")
  
  
  
  return(list(qc_long = q_curve, curve_plt = curve_plt, q_v_out= q_v_out, 
              qc = qc))
  
}
