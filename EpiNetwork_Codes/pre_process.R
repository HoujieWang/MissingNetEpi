# 05/14/2019
# pre-processing complete data
# to make recovery event times missing
# and produce "status report" regularly

# 1. extract all infection/recovery events
# 2. process them in chunks: delete recovery records, 
#    and update everybody's status until day of report
# 3. output event log w/ missingness & status reports

# 05/22/2019
# add the true parameter values to the miss_recovery output!

miss_recovery <- function(dat, interval = 7, miss_prop = 1, model = "SIR"){
  # `interval`: the regular period of status report
  # `miss_prop`: the proportion/probability of missing recovery times
  
  events = dat$events
  G0 = dat$G0; I0 = dat$I0
  
  N = nrow(G0)
  epid = rep(0,N); epid[I0] = 1
  tmax = max(events$time)
  report.times = seq(from = 0, to = (tmax %/% interval + 1) * interval, by = interval)
  report = epid
  
  # select all epi events
  # events.epi = events[events$event %in% c(1,2),]
  
  recov.ind = NULL
  for(ix in c(2:length(report.times))){
    lb = report.times[ix-1]; ub = report.times[ix]
    st = min(which(events$time > lb)); en = max(which(events$time <= ub))
    if(st == Inf){ break }
    for(r in st:en){
      z = events$event[r]
      p = events$per1[r]
      if(z==1){
        epid[p] = 1
      }
      if(z==2){
        if(model=="SIR"){
          epid[p] = -1
        }else{
          epid[p] = 0
        }
        recov.ind = c(recov.ind, r)
      }
    }
    # status report
    report = rbind(report, epid)
  }
  
  # records to remove
  m = length(recov.ind)
  remove = as.logical(rbinom(m,1,miss_prop))
  events = events[-recov.ind[remove],]
  
  # return new dataset and status report
  row.names(report) = as.character(report.times)
  
  if("truth" %in% names(dat)){
    return(list(G0=G0, I0=I0, report.times=report.times, 
                events=events, report = report, 
                truth = dat$truth))
  }else{
    cat("True parameters unavailable for this dataset!\n")
    return(list(G0=G0, I0=I0, report.times=report.times, 
                       events=events, report = report))
  }
}

# # try it out
# dats = readRDS("~/Documents/EpiNet_coupled_1.rds")
# miss1 = miss_recovery(dats)
# events.m = miss1$events
# events.o = dats$events
# sum(events.m$event==2)
# sum(events.o$event==2)
# nrow(events.o)
# nrow(events.m)
# ## the row numbers checked out
# miss1$report.times
# miss1$report
# sum(miss1$report["14",]==1) # 17 sick up to 14 days
# events.o[max(which(events.o$time < 14)),] # preval=.17
# ## prevalence checked out
# 
# # save it
# # 100% missing in recovery times
# saveRDS(miss1, "~/Documents/miss_recov.rds")

# # try out partially missing recovery times
# # eg. 50%
# dats = readRDS("~/Documents/EpiNet_coupled_4.rds")
# miss4 = miss_recovery(dats, miss_prop = 0.5)
# events.m = miss4$events
# events.o = dats$events
# sum(events.m$event==2)
# sum(events.o$event==2)
# nrow(events.o)
# nrow(events.m)
# miss4$report.times
# miss4$report
# sum(miss4$report["28",]==1) # 7 sick up to 28 days
# events.o[max(which(events.o$time < 28)),] # preval=.07
# # seems to be working well!
# # save it
# saveRDS(miss4, "~/Documents/miss_recov50_coupled4.rds")
