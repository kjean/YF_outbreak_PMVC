######################
### Assessing the impact of preventive mass vaccination campaigns on yellow fever outbreaks in Africa : a population-level self-controlled case-series study
### Jean et al., 2020 : https://www.medrxiv.org/content/10.1101/2020.07.09.20147355v2
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
missing.out.month = 7
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




#### test expo model
res_expo = NULL

all_var
for (i in 1:length(all_var)){
  print(i)
  modn = geeglm(dat$PMVC01 ~dat[,all_var[i]],family = "poisson", data = dat,
                id=1:nrow(dat), corstr = "exchangeable")
  SummaryGEEGLM(modn)
  line = SummaryGEEGLM(modn)[2,]
  rownames(line) = all_var[i]
  res_expo = rbind(res_expo, line)
}
res_expo




# cohort
pseudo.cohort = pseudo.for.cohort(dat)
head(pseudo.cohort)
table(pseudo.cohort$even)
table(pseudo.cohort$expo)
table(pseudo.cohort$expo, pseudo.cohort$even)

#############
# test expo to PMVC - NS
mod = geeglm(event ~ expo  + offset(loginterval),family = "poisson", data = pseudo.cohort,
             id=1:nrow(pseudo.cohort), corstr = "exchangeable")
SummaryGEEGLM(mod)



# explore env covariates

model_choice = "stat"
#model_choice = "mecha"

if(model_choice == "stat") {
  var= env_stat_var
}else{var = env_meca_var}


res_uni = NULL
for (i in 1:length(var)){
  modn = geeglm(event ~ pseudo.cohort[,var[i]] +  offset(loginterval),family = "poisson", data = pseudo.cohort,
                id=1:nrow(pseudo.cohort), corstr = "exchangeable")
  SummaryGEEGLM(modn)
  line = SummaryGEEGLM(modn)[2,]
  rownames(line) = var[i]
  res_uni = rbind(res_uni, line)
}
res_uni




formula = as.formula(paste0("event ~ expo +", paste(var, collapse = "+"), "+ offset(loginterval)"))
mod_multi = geeglm(formula,family = "poisson", data = pseudo.cohort,
              id=1:nrow(pseudo.cohort), corstr = "exchangeable")
SummaryGEEGLM(mod_multi)





################# with vc
pseudo.cohort.vc = pseudo.for.cohort.vc(dat)
head(pseudo.cohort.vc)
summary (pseudo.cohort.vc)
str(pseudo.cohort.vc)
mod = geeglm(event ~ vc  + offset(loginterval),family = "poisson", data = pseudo.cohort.vc,
             id=1:nrow(pseudo.cohort.vc), corstr = "exchangeable")
SummaryGEEGLM(mod)

formula = as.formula(paste0("event ~ vc +", paste(var, collapse = "+"), "+ offset(loginterval)"))
mod.vc = geeglm(formula, family = "poisson", data = pseudo.cohort.vc,
                id=1:nrow(pseudo.cohort.vc), corstr = "exchangeable")
SummaryGEEGLM(mod.vc)


cut_vc = seq(0,1, by = .2)
pseudo.cohort.vc$vc_cat = findInterval(pseudo.cohort.vc$vc, cut_vc)
pseudo.cohort.vc$vc_cat = as.factor(pseudo.cohort.vc$vc_cat )
table(pseudo.cohort.vc$vc_cat,pseudo.cohort.vc$even)

pseudo.cohort.vc$vc_cat = relevel(pseudo.cohort.vc$vc_cat, ref = 3)
formula = as.formula(paste0("event ~ vc_cat +", paste(var, collapse = "+"), "+ offset(loginterval)"))
mod.vc.cat = geeglm(formula,family = "poisson", data = pseudo.cohort.vc,
                    id=1:nrow(pseudo.cohort.vc), corstr = "exchangeable")
SummaryGEEGLM(mod.vc.cat)

