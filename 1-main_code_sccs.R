######################
### Assessing the impact of preventive mass vaccination campaigns on yellow fever outbreaks in Africa : a population-level self-controlled case-series study
### Jean et al., 2020 : https://doi.org/10.1101/2020.07.09.20147355 
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


### outbreak data
out = read.csv("formatted_data/outbreaks_1980s-2018.csv", stringsAsFactors = F); dim(out)
colnames(out)
out = out[out$year>yy0,];dim(out)
table(out$year)


### PMVC data
pmvc = read.csv("formatted_data/camp_cleaned.csv", stringsAsFactors = F); dim(pmvc)
colnames(pmvc)
pmvc = pmvc[pmvc$year>yy0,];dim(pmvc)



#############
## formate data set for analyses
############
first.out.only = T
missing.out.month = 1
out.day = 2
PMVC.month = 12
PMVC.day = 30
camp = pmvc

dat = load_out_and_pmvc_data(first.out.only, missing.out.month,out.day,PMVC.month, PMVC.day,
                             out, camp = pmvc, env = env,
                             dn1 = dn1)
head(dat)

cas = dat[dat$nb_outbreak>0,]  # select only provinces with outbreaks
cas_exp = dat[dat$nb_outbreak>0 & dat$PMVC01 ==1,]; dim(cas_exp)  #select only provinces with outbreaks and PMVC, ie SCCS sample



#############################################
#############################################
#############################################
#############################################


########
# descriptive analysis 

table(dat$nb_outbreak)
table(dat$nb_PMVC )
table(dat$nb_outbreak>0, dat$PMVC01)

sum(cas$nb_outbreak)
sum(cas_exp$nb_outbreak) # 43 for all outbreaks, 33 for first outbreak only


table(cas_exp$nb_outbreak) # 33 outbreaks considered
str(cas_exp)
dmy(cas_exp$outbreak1_date)
table(dmy(cas_exp$outbreak1_date) < dmy(cas_exp$PMVC_date1 ))


#############
## create pseudo-observation data sets
############

# SCCS
pseudo.cas = pseudo.for.sccs(cas_exp) # build up the pseudo-observations data
table(pseudo.cas$even)
sum(pseudo.cas$even)

table(pseudo.cas$expo)
table(pseudo.cas$expo, pseudo.cas$even)
100*prop.table(table(pseudo.cas$expo, pseudo.cas$even),2)
mod = gnm(even ~ expo  + offset(loginterval), eliminate = indiv, family = "poisson",
          data = pseudo.cas)
SummaryGNM(mod)

irr_vec = c(coef = summary(mod)$coefficients[,1], se =summary(mod)$coefficients[,2])



# if (first.out.only == T){
#   saveRDS(irr_vec, file = "IRR_value_1st_out_only.RDS")
# } else {saveRDS(irr_vec, file ="IRR_value.RDS") }



# SCCS with pre-expo
pseudo.cas.preexpo = pseudo.for.sccs.preexp(cas_exp, Nyears = 3)
table(pseudo.cas.preexpo$indiv)
table(pseudo.cas.preexpo$even)
sum(pseudo.cas.preexpo$even)

table(pseudo.cas.preexpo$expo)
table(pseudo.cas.preexpo$expo,pseudo.cas.preexpo$even)

mod_preexpo = gnm(even ~ expo  + offset(loginterval), eliminate = indiv, family = "poisson",
          data = pseudo.cas.preexpo)
SummaryGNM(mod_preexpo)





# SCCS with vac cov
pseudo.cas.vc = pseudo.for.sccs.vc(cas_exp); dim(pseudo.cas.vc)
table(pseudo.cas.vc$even)
hist(pseudo.cas.vc$vc)
pseudo.cas.vc$vc10 = pseudo.cas.vc$vc * 10 # for a 10% increase in coverage
mod <- clogit(even ~ vc10  + strata(indiv) + offset(loginterval), data = pseudo.cas.vc)
Summaryclogit(mod)


cut_vc = seq(0,1, by = .2)
pseudo.cas.vc$vc_cat = findInterval(pseudo.cas.vc$vc, cut_vc)
pseudo.cas.vc$vc_cat = as.factor(pseudo.cas.vc$vc_cat )
table(pseudo.cas.vc$vc_cat,pseudo.cas.vc$even)
sum(pseudo.cas.vc$even)


pseudo.cas.vc$vc_cat = relevel(pseudo.cas.vc$vc_cat, ref = 3)
mod_vc_cat <- clogit(even ~ vc_cat  + strata(indiv) + offset(loginterval), data = pseudo.cas.vc)
Summaryclogit(mod_vc_cat)




#############
## re-sample to assess the effect of spatial autocorrelation
## re-sample only 1 province per country
############
dat = load_out_and_pmvc_data(first.out.only= T, missing.out.month= 7,PMVC.month=12, PMVC.day=30,
                             out, camp = pmvc, env = env,
                             dn1 = dn1)
head(dat);dim(dat)
cas = dat[dat$nb_outbreak>0,];dim(cas)
cas_exp = dat[dat$nb_outbreak>0 & dat$PMVC01 == 1,];dim(cas_exp)


alpha =.05
zq <- qnorm(1-alpha/2)

resampled_sccs = function(dat){
  # dat is a data set generated by the load_out_and_pmvc_data function
  cas = dat[dat$nb_outbreak>0,];dim(cas)
  cou = unique(substr(cas$adm1,0,3)); length(cou)
  re_cas = NULL
  for (i in 1:length(cou)){
    pro = cas$adm1[grep(cou[i], cas$adm1)]
    line = cas[cas$adm1 == sample(pro, 1),]
    re_cas = rbind(re_cas, line)
  }
  re_pseudo.cas = pseudo.for.sccs(re_cas)
  re_mod = gnm(even ~ expo  + offset(loginterval), eliminate = indiv, family = "poisson",
               data = re_pseudo.cas)
 # print(SummaryGNM(re_mod))
  return(re_mod)
}




N = 100
res_sample = NULL
for (i in 1:N){
  sa = resampled_sccs(dat)
  coe = coef(summary(sa))[,"Estimate"]
  ste = coef(summary(sa))[,"Std. Error"]
  irr = exp(coe)
  low_b = exp(coe - zq*ste)
  up_b = exp(coe + zq*ste)
  line = c(coe, ste, irr, low_b, up_b)
  res_sample = rbind(res_sample,line)
}
res_sample = data.frame(res_sample)
colnames(res_sample)= c("coef", "se", "IRR", "IRR_lower", "IRR_upper")
dim(res_sample)
res_sample =res_sample[!is.infinite(res_sample$IRR_upper),];dim(res_sample) # retrieve resampling with random zero

theta_b = mean(res_sample$se); theta_b
pooled_IRR = exp(mean(res_sample[,1])) 
mean_low_b = exp(mean(res_sample[,1])- zq*theta_b )
mean_upper_b = exp(mean(res_sample[,1])+ zq*theta_b )
print( paste0(pooled_IRR, "[", mean_low_b, " - ", mean_upper_b, "]"))

