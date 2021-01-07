######################
### Assessing the impact of preventive mass vaccination campaigns on yellow fever outbreaks in Africa : a population-level self-controlled case-series study
### Jean et al., 2020 
### contact: Kévin JEAN, kevin.jean@lecnam.net

library(lubridate)
library(survival)
library(gnm)
library(geepack)

rm(list=ls(all=TRUE)) 


source("functions_PMVC_impact.R")





# time period considered
yy0 = 2004 # 2004 means 2005 is the 1st year considered
yy1 = 2019 # 2019 means 2018 is the last year considered

#############
## Load data sets
############
dn1 = readRDS("formatted_data/dn1.RDS"); length(dn1) # vector of identifiers for provinces considered in analyses

### vaccination coverage - retrieved from https://shiny.dide.imperial.ac.uk/polici/ (Hamlet et al, Vaccine 2019)
vc = read.csv("formatted_data/vc.pop.level_2d.csv", h=T, stringsAsFactors = F)
rownames(vc) = dn1
dim(vc)
vc = vc[,paste0("X", yy0:yy1)];dim(vc)
rownames(vc) = dn1



### environmental data - retrieved from Hamlet et al, PLOS NTDs 2018: https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0006284
env = read.csv("formatted_data/YF_seasonality_env.csv", h=T, stringsAsFactors = F);dim(env)
colnames(env) = c("adm1", "report", "logpop", "survqual", "temp.suit",
                  "rainfall", "interaction", "EVI")
env$adm0_adm1 = paste0(substr(env$adm1, 1,3),"_",substr(env$adm1, 4,6))
env = env[env$adm0_adm1 %in% dn1,];dim(env)
env_meca = env

env_meca_var = c("logpop", "survqual", "temp.suit", "rainfall", "EVI")
all_var = env_meca_var


### outbreak data
out = read.csv("formatted_data/outbreaks_1980s-2018.csv", stringsAsFactors = F); dim(out)
colnames(out)
out = out[out$year>yy0,];dim(out)
table(out$year)


### PMVC data
pmvc = read.csv("formatted_data/camp_cleaned.csv", stringsAsFactors = F); dim(pmvc)
colnames(pmvc)
pmvc = pmvc[pmvc$year>yy0,];dim(pmvc)



### outbreak data
out = read.csv("formatted_data/outbreaks_1980s-2018.csv", stringsAsFactors = F); dim(out)
colnames(out)
out = out[out$year>yy0,];dim(out)
table(out$year)
# do a sensitivity analysis excluding AGO, that weights strongly in terms of outbreaks


### PMVC data
pmvc = read.csv("formatted_data/camp_cleaned.csv", stringsAsFactors = F); dim(pmvc)
colnames(pmvc)
pmvc = pmvc[pmvc$year>yy0,];dim(pmvc)



#############
## formate data set for analyses
############
first.out.only = F
missing.out.month = 1
out.day = 2
PMVC.month = 12
PMVC.day = 30
camp = pmvc

dat = load_out_and_pmvc_data(first.out.only, missing.out.month, out.day, PMVC.month, PMVC.day,
                             out, camp = pmvc, env = env,
                             dn1 = dn1)
head(dat)
table(dat$nb_outbreak)
table(dat$PMVC01)



#########

pf = format_cohort_calc_pf(dat)
head(pf)


hist(pf$di)
hist(pf$Ni)
hist(pf$Mi)
table(pf$Mi)
###### ###### ###### ###### ###### ###### 
###### load IRR value ########
###### ###### ###### ###### ###### ###### 
alpha = .05
zq <- qnorm(1-alpha/2)
### irr_value = readRDS("IRR_value_1st_out_only.RDS");irr_value # here, import an estimate of IRR (on the log-scale) with standard error

# results from the main SCCS analysis: IRR = 0.14 (95% CI, CI: 0.06 - 0.34)
irr_value = c("coef" = -1.9530889, "se"=0.4451579); irr_value
IRR = exp(irr_value[1]);IRR
IRR_low = exp(irr_value[1] - zq*irr_value[2]);IRR_low
IRR_sup = exp(irr_value[1] + zq*irr_value[2]);IRR_sup


###### ###### ###### ###### ###### ###### 
###### Calculation of outbreaks prevented ########
###### ###### ###### ###### ###### ###### 

#calculate the number of outbreak prevented
calc_Nav = function(data, irr_value= irr_value){
  IRR = exp(rnorm(1, mean = irr_value[1], sd = irr_value[2]));IRR
  Nav =  data$Ni*(data$Ti-data$di)*(1-IRR)/data$di
  return(sum(Nav))
}


N = 10000
vec_resNav = NULL
for(i in 1:N){
  i_resample = sample(1:nrow(pf), nrow(pf), replace=T)
  res = calc_Nav(data=pf[i_resample,], irr_value= irr_value)
  vec_resNav = c(vec_resNav, res)
}
hist(vec_resNav)
summary(vec_resNav)
quantile(vec_resNav, probs = c(0.025, 0.5, 0.975))


# prevented fraction
N_obs = sum(pf$Ni)+sum(pf$Mi);N_obs

PF = 100*(1- N_obs/(N_obs+median(vec_resNav)));PF
PF_low =  100*(1- N_obs/(N_obs+quantile(vec_resNav, probs = c(0.025))));PF_low
PF_sup = 100*(1- N_obs/(N_obs+quantile(vec_resNav, probs = c(0.975))));PF_sup

