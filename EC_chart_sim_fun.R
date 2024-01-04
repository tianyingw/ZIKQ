# -------------------------------------------------------------------
#  @FileName: EC_chart_sim_fun.R
# 
#  @Descriptions:  simulation function for EC chart analysis 
#
# 
#  @Function List:
#       Function1    = the array of earthquake magnitudes
#       Function2    = the array of earthquake magnitudes
#       Function3    = the array of earthquake magnitudes
#       Function4    = the array of earthquake magnitudes
#
# ------------------Functions start here-----------------------------

# 1- H(t)

loa_s = read.csv(paste0(data_path,'fig2.csv'))%>%filter(skip =='Others'&Y<1)%>%
  transmute(age = round(X), loa =  1-Y)%>%arrange(age, loa )
loa_s = loa_s[!duplicated(loa_s$age),]

#ggplot(data = loa_s, aes(x = age, y = 1-loa))+geom_step()+geom_point()

range_data = data.frame(age = 4:15, max = c(30, 32, 33, 34, 34, 34,33, 33, 33, 33, 33, 32 ),
                        min = c(15, 10, 5,rep(0, 9) ))




LOA_f <-function(age, loa = loa_s){
  require(dplyr)
  loa = arrange(loa, age)
  # age_int = findInterval(age, loa$age)
  # loa_out = c(0, loa$loa)[age_int+1]
  if(age < min(loa$age)){
    loa_out = 0
  }else{
    if(age >max(loa$age)){loa_out = loa$loa[loa$age ==max(loa$age)]}
    else{
      loa_out = approx(loa$age, loa$loa, xout = age)$y
    }
    
  }
  return(loa_out)
}


vector_bound <-function(x, l=0.01, u = 0.99){
  return(sapply(x, function(y) min(max(y, l), u)))
}


# simulate age data from the curve 
pq<-function(y, x, AGEV, up =34, low =0){
  # y - grid points value 
  # x - grid age
  
  sp = smooth.spline(x = x, y = y)
  yout = predict(sp,x= AGEV )$y
  yout = min(yout, up)
  yout = max(yout, low)
  return(yout)
}


score_estimate<- function(age, q, loa = loa_s, range_dat=range_data){
  #age - a given age 
  #q - a given quantile 
  #loa - LOA function 
  # range_dat - 

  loa_rate_here = LOA_f(age = age, loa = loa)
  
  if(q<= loa_rate_here){
    out = 0
  }else{
    wt = c(1-q, q-loa_rate_here)/(1-loa_rate_here)
    # if(age%in%range_dat$age){
    #   range_h = range_dat[range_dat$age== age, c('min','max')]
    #   
    # }else{
      range_h = c(pq(y = range_dat$min, x = range_dat$age, AGEV = age), 
                  pq(y = range_dat$max, x = range_dat$age, AGEV = age)
                  )
    # }
    out =sum(wt*range_h)
    
  }
  return(out)
  
}

## validate the growth data
# q = seq(from = 0, to = 1, by = .05)
# q = q[q>0&q<1]
# 
# val_par = expand.grid(age = seq(from =4, to =14, by =0.1), q = q)
# val_par$v= mapply(score_estimate, age =val_par$age, q = val_par$q)
# val_par$q = factor(val_par$q) 
# plt =ggplot(data = val_par, aes_string(x = 'age', y = 'v')) + 
#   geom_path(aes_string(group = 'q'))+
#   geom_point(data =range_data, 
#              aes_string(x = 'age', y = 'min'))+ 
#   geom_point(data =range_data,
#              aes_string(x = 'age', y = 'max'),colour ='red')
#   
# ggplotly(plt)

