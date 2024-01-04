
# -------------------------------------------------------------------
#  @FileName: EC_chart_fun.R
# 
#  @Descriptions:  Define the necessary functions in this section
#
# 
#  @Function List:
#       Function1    = the array of earthquake magnitudes
#       Function2    = the array of earthquake magnitudes
#       Function3    = the array of earthquake magnitudes
#       Function4    = the array of earthquake magnitudes
#
# ------------------Functions start here-----------------------------


# quantilecurve --------------------------------------------------------------
quantilecurve<- function(qc, x, tau, x_out, tau_out = c(0.05, 1:9*0.1, 0.95),ylim_max = 34){
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
  
  #NSAA quantile curves
  curve_plt = ggplot(data = q_curve, aes_string(x = 'AGE', y = 'AVAL')) + 
    geom_path(aes_string(group = 'Quantile', color = 'Quantile')) +
    labs(title = 'NSAA') + ylab(label ='') + xlab(label = 'Years')  + ylim(0, ylim_max)
  
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

# Function 2 --------------------------------------------------------------
# quadric functions
c_fun<-function(a, b, c, x){a+b*x+c*x^2}


# LOA immputation ---------------------------------------------------------
LOA_impute<-function(dat, loa, seed = 12346){
  # -------------------------------------------------------------------
  #  @Function Function Name
  #   qunatile smoothing spline from qc curves
  #  @Descriptions:  
  #
  # 
  #  @Parameters:
  #     dat    = include AGEV, AVAL
  #     loa     = age corresponding to row
  #     tau   = quantile levels corresponding to columns
  #     x_out = Age forfitted quantiles to output 
  #
  #  @Returns:  
  #     qc_long       = quantile curves in long format
  #     curve_plt     = ggplot objective of quantile curves
  #     q_v_out       = for the grid points of interest
  # --------------------------------------------------------------------
  
  set.seed(seed = seed)
  
  # break interval baseline on survial probability
  
  loa = arrange(loa, age)
  dat = filter(dat, (!is.na(AGEV))&(!is.na(AVAL)))
  
  dat$AGE_interval =findInterval(dat$AGEV, loa$age)
  
  loa$count =  table(dat$AGE_interval)[as.character(1:dim(loa)[1])]
  loa$count[is.na(loa$count)]= 0
  loa$age_up = c(na.omit(lead(loa$age))%>%as.vector, 
                         max(max(dat$AGEV, na.rm= T),max(loa$age)) )
  loa = filter(loa, loa>0)%>%mutate(impute =round( count/(1-loa)*loa))
  
  age_impute = mapply(runif, n= loa$impute, min = loa$age, max = loa$age_up)%>%unlist
  dat_impute = data.frame(AGEV= age_impute, AVAL = runif(length(age_impute), 0, 0.01))
  dat_impute = rbind(select(dat, AGEV, AVAL), dat_impute)
  return(list(dat_impute= dat_impute, loa = loa))
}

# Quantile interpolating function -----------------------------------------

Quantile_estimate<-function(qc, AVAL, AGEV, up = 34, low = 0){
  require(SiZer)
  #cat(AVAL, AGEV,'\n')
  
  x = rownames(qc)%>%as.numeric
  qv = colnames(qc)%>%as.numeric
  
  # predict quantile value at AGEV
  pq<-function(y, x, AGEV){
    sp = smooth.spline(x = x, y = y)
    yout = predict(sp,x= AGEV )$y
    yout = min(yout, up)
    yout = max(yout, low)
    return(yout)
  }
  
  q_fit = apply(qc, 2, pq, x= x, AGEV= AGEV)
  # out = predict(piecewise.linear(x = q_fit,y = qv, CI=FALSE),
  #               AVAL)
  out = approx(q_fit, qv, xout = AVAL)$y
  
  out = max(out, 0)
  out = min(out, 1)
  return(out)
  
}

# estimate quantile value at a given age
Quantile_value<-function(qc,  AGEV, ql, up = 34, low = 0){
  
  #cat(AVAL, AGEV,'\n')
  require(SiZer)
  x = rownames(qc)%>%as.numeric
  qv = colnames(qc)%>%as.numeric  
  # predict quantile value at AGEV
  pq<-function(y, x, AGEV){
    sp = smooth.spline(x = x, y = y)
    yout = predict(sp,x= AGEV )$y
    yout = min(yout, up)
    yout = max(yout, low)
    return(yout)
  }
  
  if(AGEV%in%x){
    q_fit = qc[as.character(AGEV),]
    
  }else{ 
    q_fit = apply(qc, 2, pq, x= x, AGEV= AGEV)
    }
  

  # out = predict(piecewise.linear(y = q_fit,x = qv, CI=FALSE),
  #               ql)
  out =  approx(qv, q_fit, xout = ql)$y
  
  out = max(out, low)
  out = min(out, up)
  return(out)
  
}



## QC curve 



# ------------------Functions end here-------------------------------
