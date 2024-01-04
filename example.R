rm(list = ls())
library.list <-c('ggplot2', 'plotly',  'dplyr','reshape2','statVisual', 'GGally', 'gridExtra'
                 , 'ggcorrplot','haven', 'SiZer','quantregGrowth','quantreg','splines','KernSmooth')
a= sapply(library.list, require, character.only = T, quietly = T)

data_path = "~/Dropbox/CensoredGC/Submission_JASA/Code/"
code_path = "~/Dropbox/CensoredGC/Submission_JASA/Code/"
source(paste0(code_path, 'quantilecurve_smooth.R'))
source(paste0(code_path,'EC_chart_sim_fun.R'))
source(paste0(code_path,'EC_chart_fun.R'))
source(paste0(code_path,'mainfunctions.R'))



# get the original curve value
delta.age = 0.5
x = seq(from = 4, to = 13, by =delta.age)
t0 = seq(4,13,delta.age)
x_out = round(unique(t0))
q_out = seq(from = 0.1, to = 0.9, by = .1)

val_par = expand.grid(age = as.numeric(x_out), q = round(q_out,2))
val_par$v= mapply(score_estimate, age =val_par$age, q = val_par$q)


# major parameters
N = 1000
thres = 0
age_lower = 4
age_upper = 13
q_width = 0.05
set.seed(12346)
curve_est_all = NULL

#######################################
# simulate data
#######################################
censor.rate =0
adjust = FALSE
n_lower = 1
n_upper = 6
repN = 10
ht_est = matrix(0, nrow = repN, ncol= 21)
curve_est_all =curve_est_all2 =curve_est_all3 =NULL
avg_rate = rep(0, repN)
si = 1
#simulate measurement #, start age, quantile level
set.seed(192223)
# set the first observational age for each individual
age_start = runif(N, min = age_lower, max = 11)

# decide how many observations each individual has

ni_all = round((age_upper - age_start)/ delta.age)+1

if(censor.rate > 0){
  
  id.censor = sample(c(1:N), censor.rate*N)
  ni_censor = round(runif(length(id.censor), min = n_lower, max = ni_all[id.censor])) # ni observations per individual
  
  ni_all[id.censor] = ni_censor
}
ni_all = ni_all+1 # so every individual has at least one observation

obs_N = sum(ni_all) # total number of observations
ni_seq = unlist(sapply(ni_all, seq))


qi = runif(N, 0.01, 0.99)


# generate a dataset
sim_dat = data.frame(s = rep(1:N, times = ni_all),
                     i = ni_seq,
                     q = vector_bound( rep(qi, times = ni_all) + runif(obs_N, -q_width, q_width)), # quantile level
                     AGEV = rep(age_start, times = ni_all) + (ni_seq-1)*delta.age)

# simulate function scores
sim_dat$AVAL= mapply(score_estimate, age =sim_dat$AGEV, q = sim_dat$q)


# create a new dataset with patient_ID, observed_age and observed_value
sim_dat1 = data.frame(pid = sim_dat$s, index = sim_dat$i, age = sim_dat$AGEV, val = sim_dat$AVAL)

#########################################
# fit ZIKQ
#########################################
est.Q = ZIKQ(dat0 = sim_dat1, h_mean = NULL, h_tau = NULL, q = seq(0.01,0.99, 0.05), delta.age = 0.1, thres = thres)

est.age = as.numeric(row.names(est.Q))
est.age.int = intersect(est.age, 1:20)
est.tau = as.numeric(colnames(est.Q))
qs_loa_out<-quantilecurve_smooth(est.Q, x = est.age, tau =est.tau, x_out = est.age, xmin = 4, xmax = 13 )

qs_loa_out
