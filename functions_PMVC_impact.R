######################
### Assessing the impact of preventive mass vaccination campaigns on yellow fever outbreaks in Africa : a population-level self-controlled case-series study
### Jean et al., 2020 : https://doi.org/10.1101/2020.07.09.20147355 
### contact: Kévin JEAN, kevin.jean@lecnam.net




load_out_and_pmvc_data = function(first.out.only = F, missing.out.month = 7,out.day = 5 ,PMVC.month = 12, PMVC.day=30,
                                  out = out, camp = pmvc, env = env,
                                  dn1 = dn1){
  
  # first.out.only : if T, then censore after the first outbreak
  # missing.out.month: value for outrbeak month when missing
  # PMVC.month = assumed month of PMVC implementation 
  # out is outbreak data
  # camp = vaccination campaigns data
  # env = environmental dataset
  # dn1 = vector of provices considered
  
  out$month[is.na(out$month)] = missing.out.month
  out = out[order(out$adm0_adm1),]
  out$year_month = paste0(out$year, out$month)
  out = out[order(out$year_month),]
  out = out[order(out$adm0_adm1),]
  
  # retrive duplicated outbreaks
  uniq_out = paste(out$adm0_adm1, out$year_month, sep="_")
  out = out[!duplicated(uniq_out),]
  
  # if required: retrieve duplicates = keep only the 1st outbreak
  if(first.out.only){
    cbind(out$adm0_adm1, out$year, out$year_month, duplicated(out$adm0_adm1))
    out = out[!duplicated(out$adm0_adm1),]; dim(out)
  }
  
  # retrive duplicated campaigns
  
  camp$year_place = paste0(camp$adm0_adm1, "_", camp$year)
  camp = camp[!duplicated(camp$year_place),]; dim(camp)
  camp = camp[order(camp$year),]
  camp = camp[order(camp$adm0_adm1),]
  # table(table(camp$adm0_adm1)) # max 3 PMVC per province
  # camp = camp[!duplicated(camp$adm0_adm1),]; dim(camp)
  
  
  
  #### create the merged dataset
  nr = length(dn1)
  data = data.frame(adm1 = dn1, nb_outbreak = rep(0, nr), 
                    outbreak1_date = rep(NA, nr), 
                    outbreak2_date = rep(NA, nr),
                    outbreak3_date = rep(NA, nr),stringsAsFactors = F)
  data = data[order(data$adm1),]
  for (i in 1:nrow(data)){
    if(data$adm1[i] %in% out$adm0_adm1){
      data$nb_outbreak[i] = table(out$adm0_adm1)[data$adm1[i]]
      adm = data$adm1[i]
      data$outbreak1_date[i] = paste(out.day,out$month[out$adm0_adm1 ==adm][1],
                                     out$year[out$adm0_adm1 ==adm][1], sep ="/") # we put the outbreak at the date of 5th to avoid issue based on bisextile years...
      if(data$nb_outbreak[i] >1){
        data$outbreak2_date[i] = paste(out.day,out$month[out$adm0_adm1 ==adm][2],
                                       out$year[out$adm0_adm1 ==adm][2], sep ="/")
      } 
      if(data$nb_outbreak[i] >2){ #max 3 outbreaks per province
        data$outbreak3_date[i] = paste(out.day,out$month[out$adm0_adm1 ==adm][3],
                                       out$year[out$adm0_adm1 ==adm][3], sep ="/")
      } 
    }
  }
  
  ## add date of PMVC
  data$PMVC01 = NA
  data$PMVC_date3 = data$PMVC_date2 = data$PMVC_date1 = NA
  data$nb_PMVC= 0
  for (i in 1:nrow(data)){
    date_month1 = date_month2 =date_month3 =NA
    if(data$adm1[i] %in% camp$adm0_adm1){
      data$nb_PMVC[i] = table(camp$adm0_adm1)[data$adm1[i]]
      data$PMVC01[i] = 1
      date_month1 = paste(PMVC.day,PMVC.month , camp$year[camp$adm0_adm1 == data$adm1[i]][1], sep= "/")
      data$PMVC_date1[i] = date_month1
      
      if(data$nb_PMVC[i] >1){
        date_month2 = paste(PMVC.day,PMVC.month , camp$year[camp$adm0_adm1 == data$adm1[i]][2], sep= "/")
        data$PMVC_date2[i] = date_month2
      }
      if(data$nb_PMVC[i] >2){
        date_month3 = paste(PMVC.day,PMVC.month , camp$year[camp$adm0_adm1 == data$adm1[i]][3], sep= "/")
        data$PMVC_date3[i] = date_month3
      }
    }
  }
  data$PMVC01[is.na(data$PMVC01)] =0
  
  # merge with environmental dataset
  data = merge(data, env, by.x = "adm1", by.y= "adm0_adm1"); dim(data)
  
  return(data)
}



######### 
## formate for SCCS analysis, pseudo obs
#########
pseudo.for.sccs = function(dat){
  
  # keep only cases
  cas = dat[dat$nb_outbreak>0,]
  
  # add start and stop dates
  t0 = "31/12/2004"
  t1 = "01/01/2019"
  cas$pre_start = 0
  cas$end = interval(dmy(t0),dmy(t1))/days(1)
  
  cas$PMVC01[is.na(cas$PMVC01)] =0
  cas$PMVC_date1[cas$PMVC01 == 0 ] = "01/01/2019" # for no PMVC, put PMVC date to date end
  table(cas$PMVC01)
  
  cas$out1_day = interval(dmy(t0),dmy(cas$outbreak1_date))/days(1) 
  cas$out2_day = interval(dmy(t0),dmy(cas$outbreak2_date))/days(1) # if first.out.only = T, NA for all
  cas$out3_day = interval(dmy(t0),dmy(cas$outbreak3_date))/days(1) # if first.out.only = T, NA for all
  
  
  ## 
  nevents <- nrow(cas)
  ncuts = 3 # start, stop, start_expo
  ind = rep(cas$adm1, times = ncuts)
  
  #create an ordered list of individual events and 
  #cut points for start, end, exposure and age groups
  expo = interval(dmy(t0),dmy(cas$PMVC_date1))/days(1) 
  
  cutp <- c(as.matrix(cas$pre_start), as.matrix(cas$end), expo)
  o <- order(ind, cutp)
  ind = as.factor(ind[o])
  cutp = cutp[o]
  
  
  #calculate interval lengths, set to 0 if before start or after end
  interval <- c(0, cutp[2:length(ind)]-cutp[1:length(ind)-1])
  interval <- ifelse(cutp<=cas$pre_start[ind], 0, interval)
  interval <- ifelse(cutp>cas$end[ind], 0, interval)
  
  #event = 1 if event occurred in interval, otherwise 0
  event1 <- ifelse(cas$out1_day[ind]>cutp-interval, 1, 0)
  event1 <- ifelse(cas$out1_day[ind]<=cutp, event1, 0)
  
  event2 <- ifelse(!is.na(cas$out2_day[ind]) & cas$out2_day[ind]>cutp-interval,1, 0)
  event2 <- ifelse(!is.na(cas$out2_day[ind]) & cas$out2_day[ind]<=cutp, event2, 0)
  
  event3 <- ifelse(!is.na(cas$out3_day[ind]) & cas$out3_day[ind]>cutp-interval,1,0) 
  event3 <- ifelse(!is.na(cas$out3_day[ind]) & cas$out3_day[ind]<=cutp, event3, 0)
  
  tot_even = event1 + event2 + event3
  
  #exposure groups
  exgr <- rep(0, nevents*ncuts)
  exgr <- ifelse(cutp > expo[ind], 1, exgr)
  exgr <- as.factor(exgr)
  
  
  #put all data in a data frame, take out data with 0 interval lengths
  pseudo <- data.frame(indiv = ind[interval!=0], even = tot_even[interval!=0], interval = interval[interval!=0],  expo = exgr[interval!=0], loginterval = log(interval[interval!=0]))
  
  return(pseudo)
}


## consider exposure as continuous and changing each year
pseudo.for.sccs.vc = function(cas){
  #same function but add vac cov for each year
  
  cas$pre_start = 0
  t0 = "31/12/2004"
  t1 = "01/01/2019"
  cas$end = interval(dmy(t0),dmy(t1))/days(1)
  cas$PMVC_date1[cas$PMVC01 == 0 ] = "01/01/2019" # for no PMVC, put PMVC date to date end
  
  cas$out1_day = interval(dmy(t0),dmy(cas$outbreak1_date))/days(1)
  cas$out2_day = interval(dmy(t0),dmy(cas$outbreak2_date))/days(1) # if first.out.only = T, NA for all
  cas$out3_day = interval(dmy(t0),dmy(cas$outbreak3_date))/days(1) # if first.out.only = T, NA for all
  
  t_year = paste0("01/01/", 2005:2018)# add the day of each 1st january
  year_change = interval(dmy(t0),dmy(t_year))/days(1)
  
  ncuts = 4 + ncol(vc)# start, stop, start_expo + nb_years
  
  pseudo = NULL
  for (adm in cas$adm1){
    r = cas[cas$adm1==adm,] # r stands for row
    pre_start = 0
    end = interval(dmy(t0),dmy(t1))/days(1)
    expo = interval(dmy(t0),dmy(r$PMVC_date1))/days(1) 
    cutp = c(pre_start, end, expo, year_change)
    names(cutp) = c("pre_start", "end", "expo", paste0("X", 2005:2018))
    o = order(c(pre_start, end, expo, year_change))
    cutp = cutp[o]
    year_cutp = year(dmy(t0)+days(cutp)) # keep track of the year of each event to match with vc
    
    interval = c(0,365,cutp[3:(length(cutp))] - cutp[2:(length(cutp)-1)]) #interval is the lenght since last event
    names(interval)[2] = "X2005"
    #event = 1 if event occurred in interval, otherwise 0
    event1 = event2 = event3 =rep(0, length(cutp))
    names(event1)=names(event2)=names(event3)=names(cutp)
    
    index1 = max(which(r$out1_day>cutp))
    event1[index1] = 1
    
    index2 = ifelse( !is.na(r$out2_day), max(which(r$out2_day>cutp)), NA)
    event2[index2] = 1
    index3 = ifelse( !is.na(r$out3_day), max(which(r$out3_day>cutp)), NA)
    event3[index3] = 1
    
    tot_even = event1 + event2 + event3
    # cbind(interval, event1, event2, event3)
    
    exp = ifelse(cutp>expo,1,0)
    
    ind = rep(adm, times = ncuts)
    datr = data.frame(indiv = rep(adm, length(cutp)), even=tot_even, expo =exp, interval=interval,
                      years = year_cutp)
    datr$vc = t(vc[rownames(vc) == adm, match( paste0("X", datr$years),colnames(vc) ) ] )
    datr = datr[datr$interval>1,]
    
    pseudo = rbind(pseudo, datr)
  }
  pseudo$loginterval = log(pseudo$interval)
  
  return(pseudo)
}


###### include a pre-exposure period of N years
pseudo.for.sccs.preexp = function(cas, Nyears){
  
  # keep only cases
  #cas = dat[dat$nb_outbreak>0,]
  
  # add start and stop dates
  t0 = "31/12/2004"
  t1 = "01/01/2019"
  cas$pre_start = 0
  cas$end = interval(dmy(t0),dmy(t1))/days(1)
  
  cas$PMVC01[is.na(cas$PMVC01)] =0
  cas$PMVC_date1[cas$PMVC01 == 0 ] = "01/01/2019" # for no PMVC, put PMVC date to date end
  table(cas$PMVC01)
  
  # create the date of pre_expo start
  cas$PMVC_datepreexp = interval(dmy(t0),dmy(cas$PMVC_date1) - years(Nyears))/days(1)
  cas$PMVC_datepreexp[cas$PMVC01 == 0 ] =  interval(dmy(t0),dmy(t1))/days(1) # for no PMVC, put PMVC pre-expo date to date end
  
  cas$out1_day = interval(dmy(t0),dmy(cas$outbreak1_date))/days(1) 
  cas$out2_day = interval(dmy(t0),dmy(cas$outbreak2_date))/days(1) # if first.out.only = T, NA for all
  cas$out3_day = interval(dmy(t0),dmy(cas$outbreak3_date))/days(1) # if first.out.only = T, NA for all
  
  ## 
  nevents <- nrow(cas)
  ncuts = 4 # start, stop, start_pre_expo, start_expo
  ind = rep(cas$adm1, times = ncuts)
  
  #create an ordered list of individual events and 
  #cut points for start, end, exposure and age groups
  expo = interval(dmy(t0),dmy(cas$PMVC_date1))/days(1) 
  
  cutp <- c(as.matrix(cas$pre_start), as.matrix(cas$end), as.matrix(cas$PMVC_datepreexp), expo)
  o <- order(ind, cutp)
  ind = as.factor(ind[o])
  cutp = cutp[o]
  
  #calculate interval lengths, set to 0 if before start or after end
  interval <- c(0, cutp[2:length(ind)]-cutp[1:length(ind)-1])
  interval <- ifelse(cutp<=cas$pre_start[ind], 0, interval)
  interval <- ifelse(cutp>cas$end[ind], 0, interval)
  
  #event = 1 if event occurred in interval, otherwise 0
  event1 <- ifelse(cas$out1_day[ind]>cutp-interval, 1, 0)
  event1 <- ifelse(cas$out1_day[ind]<=cutp, event1, 0)
  
  event2 <- ifelse(!is.na(cas$out2_day[ind]) & cas$out2_day[ind]>cutp-interval,1, 0)
  event2 <- ifelse(!is.na(cas$out2_day[ind]) & cas$out2_day[ind]<=cutp, event2, 0)
  
  event3 <- ifelse(!is.na(cas$out3_day[ind]) & cas$out3_day[ind]>cutp-interval,1,0) 
  event3 <- ifelse(!is.na(cas$out3_day[ind]) & cas$out3_day[ind]<=cutp, event3, 0)
  
  tot_even = event1 + event2 + event3
  
  #exposure groups
  exgr <- rep(0, nevents*ncuts)
  exgr <- ifelse(cutp > cas$PMVC_datepreexp[ind], 1, exgr)
  exgr <- ifelse(cutp > expo[ind], 2, exgr)
  exgr <- as.factor(exgr)
  
  #put all data in a data frame, take out data with 0 interval lengths
  pseudo <- data.frame(indiv = ind[interval!=0], even = tot_even[interval!=0], interval = interval[interval!=0],  expo = exgr[interval!=0], loginterval = log(interval[interval!=0]))
  
  return(pseudo)
  
}



######### 
## formate for cohort analysis, pseudo obs
#########
######### 
## formate for cohort analysis, pseudo obs
#########
pseudo.for.cohort = function(dat){
  
  t0 = "31/12/2004"
  t1 = "01/01/2019"
  dat$pre_start = 0
  dat$end = interval(dmy(t0),dmy(t1))/days(1)
  dat$PMVC_date1[dat$PMVC01 == 0 ] = "01/01/2019" # for no PMVC, put PMVC date to date end
  dat$outbreak1_date[is.na(dat$outbreak1_date) ] = "02/01/2019" # for no outbreak, put outreak date to date end
  
  dat$out_day = interval(dmy(t0),dmy(dat$outbreak1_date))/days(1) 
  
  ncuts = 4 # start, stop, start_expo, event
  pseudo = NULL
  for (adm in dat$adm1){
    r = dat[dat$adm1==adm,] # r stands for row
    pre_start = 0
    end = interval(dmy(t0),dmy(t1))/days(1)
    ev = r$out_day
    expo = interval(dmy(t0),dmy(r$PMVC_date1))/days(1) 
    
    cutp = c(pre_start, ev, end, expo)
    o = order(c(pre_start, ev, end, expo))
    cutp = cutp[o]
    interval = c(0,cutp[2:length(cutp)] - cutp[1:(length(cutp)-1)])
    event = ifelse(ev==cutp,1,0)
    exp = ifelse(cutp>expo,1,0)
    to_censor = ifelse(cutp>ev,1,0) # indicative variable flagging lines to be censored, ie after outbreaks
    
    ind = rep(adm, times = ncuts)
    datr = data.frame(indiv = ind, event=event, expo =exp, interval=interval)
    # need to sensor and retrieve interval = 0
    datr = datr[to_censor==0 & interval>1,]
    
    pseudo = rbind(pseudo, datr)
    
    
  }
  pseudo$loginterval = log(pseudo$interval)
  
  pseudo = merge(pseudo, env, by.x="indiv", by.y = "adm0_adm1")
  
  return(pseudo)
}


#### same for all outbreaks
format_cohort_calc_pf = function(dat){
  
  t0 = "31/12/2004"
  t1 = "01/01/2019"
  dat$pre_start = 0
  dat$end = interval(dmy(t0),dmy(t1))/days(1)
  dat$PMVC_date1[dat$PMVC01 == 0 ] = "01/01/2019" # for no PMVC, put PMVC date to date end
  dat$outbreak1_date[is.na(dat$outbreak1_date) ] = "02/01/2019" # for no outbreak, put outreak date to date end
  dat$outbreak2_date[is.na(dat$outbreak2_date) ] = "02/01/2019"
  dat$outbreak3_date[is.na(dat$outbreak3_date) ] = "02/01/2019"
  
  dat$out1_day = interval(dmy(t0),dmy(dat$outbreak1_date))/days(1) 
  dat$out2_day = interval(dmy(t0),dmy(dat$outbreak2_date))/days(1) # if first.out.only = T, NA for all
  dat$out3_day = interval(dmy(t0),dmy(dat$outbreak3_date))/days(1) # if first.out.only = T, NA for all
  
  di =  interval(dmy(t0),dmy(dat$PMVC_date1))/days(1) 
  Ti = interval(dmy(t0),dmy(t1))/days(1)
  
  Ni = rep(0, nrow(dat))
  Ni = ifelse(dat$out1_day< di, Ni+1, Ni)
  Ni = ifelse(dat$out2_day< di, Ni+1, Ni)
  Ni = ifelse(dat$out3_day< di, Ni+1, Ni)
  
  Mi = rep(0, nrow(dat))
  Mi = ifelse(dat$out1_day>= di & dat$out1_day<dat$end, Mi+1, Mi)
  Mi = ifelse(dat$out2_day>= di & dat$out2_day<dat$end, Mi+1, Mi)
  Mi = ifelse(dat$out3_day>= di & dat$out3_day<dat$end, Mi+1, Mi)
  
  vec = data.frame(adm1 = dat$adm1, di, Ti, Ni, Mi)
  return(vec)
}







#### same with vaccination coverage as continuous
pseudo.for.cohort.vc = function(dat){
  
  t0 = "31/12/2004"
  t1 = "01/01/2019"
  dat$pre_start = 0
  dat$end = interval(dmy(t0),dmy(t1))/days(1)
  dat$PMVC_date1[dat$PMVC01 == 0 ] = "01/01/2019" # for no PMVC, put PMVC date to date end
  dat$outbreak1_date[is.na(dat$outbreak1_date) ] = "02/01/2019" # for no outbreak, put outreak date to date end
  
  dat$out_day = interval(dmy(t0),dmy(dat$outbreak1_date))/days(1) 
  t_year = paste0("01/01/", 2005:2018)# add the day of each 1st january
  year_change = interval(dmy(t0),dmy(t_year))/days(1)
  
  ncuts = 4 + length(year_change)# start, stop, start_expo + nb_years
  
  pseudo = NULL
  for (adm in dat$adm1){
    r = dat[dat$adm1==adm,] # r stands for row
    pre_start = 0
    end = interval(dmy(t0),dmy(t1))/days(1)
    ev = r$out_day
    expo = interval(dmy(t0),dmy(r$PMVC_date1))/days(1) 
    
    cutp = c(pre_start, ev, end, expo, year_change)
    names(cutp) = c("pre_start", "ev", "end", "expo", paste0("X", 2005:2018))
    o = order(c(pre_start, ev, end, expo, year_change))
    cutp = cutp[o]
    year_cutp = year(dmy(t0)+days(cutp)) # keep track of the year of each event to match with vc
    interval = c(0,cutp[2:length(cutp)] - cutp[1:(length(cutp)-1)])
    event = ifelse(ev==cutp,1,0)
    exp = ifelse(cutp>expo,1,0)
    to_censor = ifelse(cutp>ev,1,0) # indicative variable flagging lines to be censored, ie after outbreaks
    
    ind = rep(adm, times = ncuts)
    datr = data.frame(indiv = ind, event=event, expo =exp, interval=interval,  years = year_cutp)
    datr$vc = t(vc[rownames(vc) == adm, match( paste0("X", datr$years),colnames(vc) ) ] )
    datr = datr[to_censor==0 & interval>1,]
    pseudo = rbind(pseudo, datr)
    
    
  }
  pseudo$loginterval = log(pseudo$interval)
  
  pseudo = merge(pseudo, env, by.x="indiv", by.y = "adm0_adm1")
  
  return(pseudo)
}






##################################
##################################
##################################
SummaryGEEGLM <- function(fit, alpha=.05, dig=2, p.dig=4){
  # output a summary data frame of GEE results
  #fit is a fitted geese() object. dig is number of digits to report.
  zq <- qnorm(1-alpha/2)
  estimate <- fit$coefficients
  lower    <- fit$coefficients - zq*coef(summary(fit))[,"Std.err"]
  upper    <- fit$coefficients + zq*coef(summary(fit))[,"Std.err"]
  robust.se <- round(coef(summary(fit))[,"Std.err"], dig)
  p <- coef(summary(fit))[,"Pr(>|W|)"]
  p <- round(p, digits=p.dig)
  Sig <- ifelse(p<0.05, '*', ifelse(p<0.01,'**', ifelse(p<0.001, '***', '')))
  RR       <- round(exp(estimate), dig) #incidence rate ratio
  RR.lower <- round(exp(lower), dig)
  RR.upper <- round(exp(upper), dig)
  estimate  <- round(estimate, dig)
  return(data.frame(cbind(estimate,robust.se,RR,RR.lower,RR.upper,p, Sig)))	
}


SummaryGEEGLM_v2 <- function(fit, alpha=.05, dig=2, p.dig=4){
  # output a summary data frame of GEE results
  #fit is a fitted geese() object. dig is number of digits to report.
  zq <- qnorm(1-alpha/2)
  estimate <- fit$coefficients[-1]
  lower    <- fit$coefficients[-1] - zq*coef(summary(fit))[,"Std.err"][-1]
  upper    <- fit$coefficients[-1] + zq*coef(summary(fit))[,"Std.err"][-1]
  
  p <- coef(summary(fit))[,"Pr(>|W|)"][-1]
  p <- round(p, digits=p.dig)
  #Sig <- ifelse(p<0.05, '*', ifelse(p<0.01,'**', ifelse(p<0.001, '***', '')))
  RR       <- round(exp(estimate), dig) #incidence rate ratio
  RR.lower <- round(exp(lower), dig)
  RR.upper <- round(exp(upper), dig)
  
  return(data.frame(cbind(RR, paste0(RR.lower, "-", RR.upper) , p)))	
}

SummaryGNM <- function(fit, alpha=.05, dig=2, p.dig=4){
  # output a summary data frame of GNM results
  #f it is a fitted gnm() object. dig is number of digits to report.
  zq <- qnorm(1-alpha/2)
  estimate <- summary(fit)$coefficients[,1]
  se <- summary(fit)$coefficients[,2]
  lower <- estimate -zq*se
  upper <- estimate +zq*se
  
  RR <- round(exp(estimate), dig) #incidence rate ratio
  RR.lower <- round(exp(lower), dig)
  RR.upper <- round(exp(upper), dig)
  return(data.frame(cbind(RR, paste0(RR.lower, "-", RR.upper) )))	
  
}

Summaryclogit <- function(fit, alpha=.05, dig=2, p.dig=4){
  # output a summary data frame of GNM results
  #f it is a fitted gnm() object. dig is number of digits to report.
  zq <- qnorm(1-alpha/2)
  estimate <- summary(fit)$coefficients[,1]
  se <- summary(fit)$coefficients[,3]
  lower <- estimate -zq*se
  upper <- estimate +zq*se
  
  RR <- round(exp(estimate), dig) #incidence rate ratio
  RR.lower <- round(exp(lower), dig)
  RR.upper <- round(exp(upper), dig)
  return(data.frame(cbind(RR, paste0(RR.lower, "-", RR.upper) )))	
  
}



##################################
  ##################################
##################################



