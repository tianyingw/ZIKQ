# Main functions
#' Function to preprocess the data
#' @param dat dataset, contains "pid": patient ID, "index": observation index, "age": observed age, "val": observed value
#' @param delta.age a number indicating age interval length. The default is 0.25.
#' @import dplyr
#' @export
dat_pre_process <- function(dat, delta.age = 0.25){

  # order, clean excessive zeros
  dat = filter(dat, val ==0)%>%group_by(pid)%>%mutate(ind =first(index))%>%dplyr::select(pid, ind)%>%
    distinct%>%full_join(dat)%>%filter(is.na(ind)|(index<=ind))%>%arrange(pid, age)%>%data.frame
  
  #dat = filter(dat, val ==0)%>%group_by(pid)%>%mutate(ind =first(index))%>%select(pid, ind)%>%
  #  distinct%>%full_join(dat)%>%filter(is.na(ind)|(index<=ind))%>%arrange(pid, age)%>%data.frame

  # add rounded age
  dat$age.int = round(dat$age/delta.age)*delta.age

  #dat1 = dat %>% group_by(pid) %>% summarise(
  #  pid = pid,
  #  ind = ind,
  #  index = index,
  #  age = age,
  #  val = val,
  #  age.int = age.int,
  #  dup.age = duplicated(age.int)
  #)
  # reduce the repeated observations
  #dat1 = dat1 %>% group_by(pid, age.int)  %>% summarise(
  #  pid = pid[1],
  #  ind = ind[1],
  #  age.int = mean(age.int),
  #  val = mean(val)
  #)
  dat1 = dat %>% group_by(pid) %>%
    filter(!duplicated(age.int))
  return(list(dat1 = dat1, dat = dat))
}

#' Function to estimate the survival rate
#' @param dat data frame, contains "pid": patient ID, "index": observation index, "age": observed age, "val": observed value.
#' @param thres the threshold for manually imputation: if one's last observation is less than the threshold, its next observation is imputed to be zero.
#' @export
St_est <- function(dat, thres = 5){
  age.int = sort(unique(dat$age.int))

  t0 = age.int
  di = rep(0, length(age.int)) # having event at time t
  ni = rep(0, length(age.int)) # survived at time t (does not include the samples joining beyond time t)

# impute missing observations
  dat.short = dat %>% group_by(pid) %>% summarise(
    max_age = max(age.int),
    min_age = min(age.int),
    event_ind = ind[1],
    last_val = val[length(val)]
  )


  # calculate the individuals whose observations were observed between this age range
  for(ia in 1:length(age.int)){

    # use longitudinal format
    #ni[ia] = length(which(dat1$age.int == age.int[ia] & dat1$val > 0))
    #di[ia] = length(which(dat1$age.int == age.int[ia] & dat1$val == 0 & dat1$ind > 1))

    # use summarize format
    ni[ia] = length(which(dat.short$min_age <= age.int[ia] & dat.short$max_age > age.int[ia]))+ length(which(dat.short$max_age == age.int[ia] & is.na(dat.short$event_ind)))
    di[ia] = length(which(dat.short$max_age == age.int[ia] & dat.short$event_ind > 1))
    #
    if(ia > 1){
     di[ia] = di[ia] + length(which(dat.short$max_age == age.int[ia-1] & is.na(dat.short$event_ind) & dat.short$last_val < thres))
    }
  }

  hi = di/(ni+di)
  p1 =  1-hi # survival rate
  S1 = cumprod(p1)


   # S1 is the estimated ht, then we need to interpolate
  St.fun = approxfun(age.int, S1)

  return(St.fun)
}

#' Estimate the values at quantile sequence q
#' @param dat0 data frame, includes "pid": patient ID, "index": observation index, "age": observed age, "val": observed value.
#' @param h_mean bandwidth for regression mean estimation. If not specified, the technique of Ruppert et al (1995) will be implemented.
#' @param h_tau bandwidth for quantile regression. If not specified, h_tau will be calculated based on h_mean.
#' @param q a sequence of quantile levels you want to estimate.
#' @param delta.age a number indicating age interval length. The default is 0.25.
#' @param thres a threshold for imputation. If one's last observation is positive but less than the threshold, we manually impute 0 as the next observation. The default value is 0.
#' @param St.fun a function which provides the survival rate at certain age range. The range should cover all ages observed in the data.
#' @import quantreg KernSmooth
#' @export

ZIKQ <- function(dat0, h_mean = NULL,  h_tau = NULL, q=seq(from = 0, to = 1, by = .05), age.range = c(4,13), delta.age = 0.25, thres = 0, St.fun = NULL, Smooth.St = FALSE, df.smooth.st = 5){
  
  
  q = q[q>0&q<1]
  
  dat.list  = dat_pre_process(dat = dat0, delta.age = delta.age)
  dat = dat.list$dat # did not delete duplicated age
  dat1 = dat.list$dat1 # delete duplicated age
  
  if(is.null(St.fun)){
    St.fun = St_est(dat = dat1,thres = thres)
  }
  if(is.null(h_mean)){
    h_mean = dpill(dat$age, dat$val)
    
  }
  
  
  
  ind.temp = which((dat$val)==0)
  simu_temp = dat[-ind.temp, ] # only consider positive Y
  #t0 = sort(unique(dat$age.int))
  grid.x = seq(age.range[1], age.range[2], delta.age)
  # calculate one row at a time (one age for all tau at a time)
  fitqr = lprq2(x = simu_temp$age, y = simu_temp$val, h_mean_org = h_mean, h_tau = h_tau, grid.x = grid.x, tau = q, St.fun = St.fun, Smooth.St = Smooth.St, df.smooth.st = df.smooth.st)
  
  est.Q = fitqr$fv
  rownames(est.Q) = fitqr$xx
  colnames(est.Q) = q
  
  
  
  return(est.Q)
}



#' Supprot function for kernel esetimation
#' @export
lprq2 <-function(x, y, h_mean_org, grid.x , h_tau = NULL, tau=.5, St.fun, Smooth.St = FALSE, df.smooth.st = 5){
  xx <- grid.x
  ht0_org = St.fun(xx)
  if(length(which(is.na(ht0_org))) > 0){
    xx = xx[-which(is.na(ht0_org))]
    ht0_org = ht0_org[-which(is.na(ht0_org))]
  }
  
  if(Smooth.St){
    
    #temp.approx = approx(xx, ht0_org, xout = grid.x)
    #ht0 = temp.approx$y
    ht0 = predict(lm(ht0_org~bs(xx, df = df.smooth.st)), data.frame(xx = xx))
    # check if h0 is monotone decreasing
    ht0[which(ht0 > 1)] = 1
    ht0[which(ht0 < 0)] = 0
    
    for(ii in 2:length(xx)){
      if(ht0[ii] > ht0[ii-1]){
        ht0[ii] = ht0[ii-1]
      }
    }
    
  }else{
    ht0 = ht0_org
  }
  
  #ht0 = rep(1, length(ht0))
  #if(is.null(h_tau)){
  #  h_tau = h_mean*(tau*(1-tau)/dnorm(qnorm(tau))^2)^(1/5)
  #}
  
  
  #print(h_mean)
  fv<-matrix(0, nrow = length(xx), ncol = length(tau))
  dv<-fv
  
  
#  for(i in 1:length(xx)) {
#    new.tau = (tau-(1-ht0[i]))/ht0[i]
#    pos.tau.id = which(new.tau > 0)
    
#    if(length(pos.tau.id) > 0){
#      z <- x - xx[i]
      
#      for(k in 1:length(pos.tau.id)){
#        wx <- dnorm(z/h_tau[pos.tau.id[k]])
#        r <- rq(y~z, weights = wx, tau = new.tau[pos.tau.id[k]], ci = FALSE)
#        fv[i,pos.tau.id[k]] = r$coef[1]
#        dv[i,pos.tau.id[k]] = r$coef[2]
#      }
      
      
      
 #   }
    
 # }
  
  
  #par(mar=c(4.5,4.5,2,2))
  
  for(i in 1:length(tau)) {
    new.tau = (tau[i]-(1-ht0))/ht0
    pos.tau.id = which(new.tau > 0)
    
    if(length(pos.tau.id) > 0){
      
      # each column is an age point, each row is a subject
      z.set = apply(as.matrix(xx[1:length(pos.tau.id)]),1,function(yy){x - yy})
      tau.set = new.tau[pos.tau.id]
      
      #wx.set = apply(z.set, 2, function(yy){dnorm(yy/h_tau[i])})
      
      
      for(k in 1:length(pos.tau.id)){
        if(is.null(h_mean_org)){
          h_mean = dpill(y, z.set[,k])
        }else{
          h_mean = h_mean_org
        }
        
        h_tau = h_mean*(tau.set[k]*(1-tau.set[k])/dnorm(qnorm(tau.set[k]))^2)^(1/5)
        wx.set = 1*(abs(z.set[,k])<3*h_tau)#dnorm(z.set[,k]/(h_tau))
        r <- rq(y~z.set[,k], weights = wx.set, tau = tau.set[k], ci = FALSE)
        #print(summary(r))
        fv[pos.tau.id[k],i] = r$coef[1]
        dv[pos.tau.id[k],i] = sum(abs(wx.set)>1e-2)#r$coef[2]
      }
      
      
      
    }
    #plot(xx,dv[,i], xlab = "Age", ylab = "number of effective weight (abs>1e-2)", main = paste0("tau = ", tau[i]), cex.main = 2, cex.lab = 2, cex.axis = 2)
    
  }
  list(xx = xx, fv = fv, dv = dv) }




#' Return estimated value at given age and quantile level
#' @param est.Q estimated quantile curve
#' @param age estimated value at given age
#' @param q estimated value at given quantile level
#' @param up upper bound, the default is 34.
#' @param low lower bound, the default is 0.
#' @import SiZer dplyr
#' @export
EstQuantile_value<-function(est.Q, age, q, up = 34, low = 0){

  x = rownames(est.Q)%>%as.numeric
  qv = colnames(est.Q)%>%as.numeric
  # predict quantile value at AGEV
  pq<-function(y, x, age){
    sp = smooth.spline(x = x, y = y)
    yout = predict(sp,x= age )$y
    yout = min(yout, up)
    yout = max(yout, low)
    return(yout)
  }

  if(age%in%x){
    q_fit = est.Q[as.character(age),]

  }else{
    q_fit = apply(est.Q, 2, pq, x= x, age= age)
  }


  # out = predict(piecewise.linear(y = q_fit,x = qv, CI=FALSE),
  #               ql)
  out =  approx(qv, q_fit, xout = q)$y

  out = max(out, low)
  out = min(out, up)
  return(out)

}


lprq <-function(x, y, h, m=50 , tau=.5){
      xx <- seq(min(x),max(x),length=m)
      if(length(tau) > 1){
        fv<-matrix(0, nrow = m, ncol = length(tau))
        dv<-fv
      }else{
        fv<-xx
        dv<-xx
      }
     
     for(i in 1:length(xx)) {
       z <- x - xx[i]
       wx <- dnorm(z/h)
       r <- rq(y~z, weights = wx, tau = tau, ci = FALSE)
       if(length(tau) > 1){
         fv[i,] = r$coef[1,]
         dv[i,] = r$coef[2,]
       }else{
         fv[i] = r$coef[1]
         dv[i] = r$coef[2]
       }
       
      }
     list(xx = xx, fv = fv, dv = dv) }

