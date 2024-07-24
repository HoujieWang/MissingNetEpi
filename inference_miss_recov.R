# 05/16/2019
# Inference on partially observed data
# --Missing (some) recovery times 

# 05/22/2019
# add another function that does the entire process:
# 1) read in a complete dataset
# 2) generate missingness in recovery times
# 3) do inference (and make plots) in 3 different imputation methods
# 4) calculate ESS and conduct Geweke diagnostics

# preparations: loading stuff

source("./inference_util.R")
source("./pre_process.R")

# try some other prior
# just to shake things up
pr = data.frame(count = rep(1,4), avg = c(0.02, 0.1, 0.004, 0.06))


# 1. function to do Bayesian inference on data w/ missing recovery times
# output/plot results every `output.sams` recorded samples
# priors: data frame of vars `count` & `avg`
# assume "quarantine" case!

# 05/21/2019
# add a third method for proposing recovery times: ``..._MH'': 
# keep the previous sample if new proposal is the new proposal is no good

# 05/22/2019
# pull true params from the dataset, if such info is available

# 06/04/2019
# use `parse_augment2` function 
infer_miss_recov <- function(dats, priors, init.params = NULL,
                             verbose = T, plot = T, output.sams = 100, 
                             samples = 1000, burn = 100, thin = 1,
                             impute = "filter", model="SIR", 
                             timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dats$G0; I0 = dats$I0; events = dats$events
  reports = dats$report; report.times = dats$report.times
  
  if("truth" %in% names(dats)){
    true.params = dats$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, reports, report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # a "use new proposal" flag for the MH method
  # use.new=T -> need to recompute suff.stats.
  if(impute=="MH"){
    use.new = T
  }
  
  # run iterations
  S = samples * thin + burn
  for(it in 1:S){
    
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    if(impute == "MH"){
      if(it > 1){
        prev.times = times
      }else{
        prev.times = NULL
      }
    }
    times = NULL
    for(ix in 1:length(MR$recover)){
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      if (impute == "filter"){
        imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      }else if (impute == "MH"){
        # use filter to get a viable imputation in the 1st interation
        if(it == 1){
          imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
        }else{
          imputed = propose_recov_MH(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
        }
      }else{
        imputed = propose_recov_rej(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      }
      times = c(times, imputed)
    }
    if(impute == "MH" & it > 1){
      # update the flag
      use.new = !all(is.na(times))
      # use the previous proposals if new proposals are not good
      # for each particular interval!
      times[is.na(times)] = prev.times[is.na(times)]
    }
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    # (2) compute event counts and summations
    # PA = parse_augment2(G0, I0, events, recover.dat)
    PA = parse_augment_cpp(G0, I0, events, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  # traceplot
  # modification: add in 95% posterior credible intervals
  for(ix in 1:8){
    v = vars[ix]
    sams = params[,ix]
    if(!is.null(true.params)){
      truth = true.params[ix]
      yl = range(sams, truth)
      plot(sams~c(1:samples), main=paste("Traceplot for",v), 
           xlab = "sample", ylab = v, type="l", ylim = yl)
      abline(h=truth,col="red",lwd=2)
    }else{
      yl = range(sams)
      plot(sams~c(1:samples), main=paste("Traceplot for",v), 
           xlab = "sample", ylab = v, type="l", ylim = yl)
    }
    bounds = quantile(sams, c(.025, .975))
    abline(h=bounds[1],col="gray",lty=2,lwd=2)
    abline(h=bounds[2],col="gray",lty=2,lwd=2)
  }
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  return(params)
}

# 11/20/2022 Added by Houjie To accomodate the observation pattern of 
# only observing the network at discrete times and the people the patients meet before
# they get sick
# 11/21/2022 Added one more step to make network events consistent with imputed recovery times
infer_miss_recov2 = function(dat, priors, obs_time,
                             init.params = NULL,
                             verbose = T, plot = T, output.sams = 100, 
                             samples = 1000, burn = 100, thin = 1,
                             impute = "filter", model="SIR", 
                             timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  reports = dat$report; report.times = dat$report.times
  
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, report = reports, times = report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # a "use new proposal" flag for the MH method
  # use.new=T -> need to recompute suff.stats.
  if(impute=="MH"){
    use.new = T
  }
  
  # run iterations
  S = samples * thin + burn
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    if(impute == "MH"){
      if(it > 1){
        prev.times = times
      }else{
        prev.times = NULL
      }
    }
    times = NULL
    for(ix in 1:length(MR$recover)){
      # ix = 1
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      if (impute == "filter"){
        imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      }else if (impute == "MH"){
        # use filter to get a viable imputation in the 1st interation
        if(it == 1){
          imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
        }else{
          imputed = propose_recov_MH(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
        }
      }else{
        imputed = propose_recov_rej(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      }
      times = c(times, imputed)
    }
    if(impute == "MH" & it > 1){
      # update the flag
      use.new = !all(is.na(times))
      # use the previous proposals if new proposals are not good
      # for each particular interval!
      times[is.na(times)] = prev.times[is.na(times)]
    }
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    # (2) compute event counts and summations
    #PA = parse_augment(G0, I0, events, recover.dat)
    
    # PA = parse_augment2(G0, I0, events, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    
    events = events[,-c(5:6)]
    recover.dat$per2 = NA
    recover.dat$event = 2
    recover.dat = recover.dat[,names(events)]
    
    events_imputed = rbind(events, recover.dat)
    events_imputed = events_imputed[order(events_imputed$time),]
    events_imputed = unique(events_imputed)
    
    epi_events = events_imputed[events_imputed$event %in% 1: 2, ]
    epi_events = rbind(c(0, 1, I0, NA), epi_events)
    sick_status = matrix(0, nrow = (length(obs_time)-1), ncol = nrow(miss_dats2$G0))
    rownames(sick_status) = obs_time[-1]
    colnames(sick_status) = 1: ncol(sick_status)
    for (i in 1: nrow(recover.dat)){
      # i = 3
      per = recover.dat$per1[i]
      epi_time_i = epi_events$time[epi_events$per1 == per]
      sick_status[, per] = epi_time_i[1] <= obs_time[-1] & obs_time[-1] <= epi_time_i[2]
    }
    
    for (i in 1: nrow(recover.dat)){
      # i = 2
      per = recover.dat$per1[i]
      selected_idx = which((events_imputed$per1 == per | events_imputed$per2 == per) &
                             !events_imputed$event %in% 1: 2 &
                             events_imputed$time %% 1 == 0)
      net_events_temp = events_imputed[selected_idx, ]
      net_events_temp[net_events_temp$per1 < per, c("per1", "per2")] = net_events_temp[net_events_temp$per1 < per, c("per2", "per1")]
      events_imputed$event[selected_idx] = sick_status[as.character(net_events_temp$time), per] + 
        sick_status[cbind(as.character(net_events_temp$time), as.character(net_events_temp$per2))] + 
        ifelse(net_events_temp$event %in% 3: 5, 3, 6)
    }
    
    
    event.counts = table(factor(events_imputed$event, levels = 1:8))
    event.counts[1] = event.counts[1] - 1
    big.sum2 = big.sumer(G0, I0, events_imputed)
    big.sums = colSums(big.sum2[, 1: 8] * big.sum2$time_diff)
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  # traceplot
  # modification: add in 95% posterior credible intervals
  # for(ix in 1:8){
  #   v = vars[ix]
  #   sams = params[,ix]
  #   if(!is.null(true.params)){
  #     truth = true.params[ix]
  #     yl = range(sams, truth)
  #     plot(sams~c(1:samples), main=paste("Traceplot for",v), 
  #          xlab = "sample", ylab = v, type="l", ylim = yl)
  #     abline(h=truth,col="red",lwd=2)
  #   }else{
  #     yl = range(sams)
  #     plot(sams~c(1:samples), main=paste("Traceplot for",v), 
  #          xlab = "sample", ylab = v, type="l", ylim = yl)
  #   }
  #   bounds = quantile(sams, c(.025, .975))
  #   abline(h=bounds[1],col="gray",lty=2,lwd=2)
  #   abline(h=bounds[2],col="gray",lty=2,lwd=2)
  # }
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  return(params)
}

infer_miss_recov3 = function(dat, priors, obs_time,
                             init.params = NULL,
                             verbose = T, plot = T, output.sams = 100, 
                             samples = 1000, burn = 100, thin = 1,
                             impute = "filter", model="SIR", 
                             timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  reports = dat$report; report.times = dat$report.times
  
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, report = reports, times = report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # a "use new proposal" flag for the MH method
  # use.new=T -> need to recompute suff.stats.
  if(impute=="MH"){
    use.new = T
  }
  
  # run iterations
  S = samples * thin + burn
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    if(impute == "MH"){
      if(it > 1){
        prev.times = times
      }else{
        prev.times = NULL
      }
    }
    times = NULL
    for(ix in 1:length(MR$recover)){
      # ix = 1
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      if (impute == "filter"){
        imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      }else if (impute == "MH"){
        # use filter to get a viable imputation in the 1st interation
        if(it == 1){
          imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
        }else{
          imputed = propose_recov_MH(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
        }
      }else{
        imputed = propose_recov_rej(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      }
      times = c(times, imputed)
    }
    if(impute == "MH" & it > 1){
      # update the flag
      use.new = !all(is.na(times))
      # use the previous proposals if new proposals are not good
      # for each particular interval!
      times[is.na(times)] = prev.times[is.na(times)]
    }
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    # (2) compute event counts and summations
    #PA = parse_augment(G0, I0, events, recover.dat)
    
    # # Make sure the event types of all events are correctly labeled
    infect_times = rep(max(c(events$time, recover.dat$time))+1, length = nrow(G0))
    infect_times[I0] = 0
    infect_times[events$per1[events$event == 1]] = events$time[events$event == 1]
    recov_times = rep(max(c(events$time, recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    net_event_idx = which(events$event %in% 3: 8)
    full_net_events = events[net_event_idx, ]
    full_CD_type = ifelse(full_net_events$event %in% 3: 5, 3, 6)
    full_net_events$event = (infect_times[full_net_events$per1] < full_net_events$time & full_net_events$time <= recov_times[full_net_events$per1]) +
      (infect_times[full_net_events$per2] < full_net_events$time & full_net_events$time <= recov_times[full_net_events$per2]) +
      ifelse(full_net_events$event %in% 3: 5, 3, 6)
    events[net_event_idx, ] = full_net_events
    # full_events_w_recov = rbind.data.frame(events, cbind.data.frame(time = recover.dat$time, event = 2, per1 = recover.dat$per1, per2 = NA))
    # full_events_w_recov = full_events_w_recov[order(full_events_w_recov$time), ]
    # parse_augment3(G0, I0, full_events_w_recov[full_events_w_recov$event != 2, ], full_events_w_recov[full_events_w_recov$event == 2, ])[1:2]
    # a = big.sumer(G0, I0, full_events_w_recov)
    # colSums(a[, 1: 8] * a$time_diff)
    PA = parse_augment_cpp(G0, I0, events, recover.dat)
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    # event.counts = PA$event.counts
    
    big.sums = PA$big.sums
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  # traceplot
  # modification: add in 95% posterior credible intervals
  # for(ix in 1:8){
  #   v = vars[ix]
  #   sams = params[,ix]
  #   if(!is.null(true.params)){
  #     truth = true.params[ix]
  #     yl = range(sams, truth)
  #     plot(sams~c(1:samples), main=paste("Traceplot for",v), 
  #          xlab = "sample", ylab = v, type="l", ylim = yl)
  #     abline(h=truth,col="red",lwd=2)
  #   }else{
  #     yl = range(sams)
  #     plot(sams~c(1:samples), main=paste("Traceplot for",v), 
  #          xlab = "sample", ylab = v, type="l", ylim = yl)
  #   }
  #   bounds = quantile(sams, c(.025, .975))
  #   abline(h=bounds[1],col="gray",lty=2,lwd=2)
  #   abline(h=bounds[2],col="gray",lty=2,lwd=2)
  # }
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  return(params)
}

infer_miss_recov4_3_test = function(dat, priors, obs_time,
                                    init.params = NULL,
                                    verbose = T, plot = T, output.sams = 100, 
                                    samples = 1000, burn = 100, thin = 1,
                                    impute = "filter", model="SIR", 
                                    timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  reports = dat$report; report.times = dat$report.times
  window_length = obs_time[2] - obs_time[1]
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, report = reports, times = report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  net_events = events[events$event %in% 3: 8, ]
  net_events[net_events$per2 < net_events$per1, c("per1", "per2")] = 
    net_events[net_events$per2 < net_events$per1, c("per2", "per1")]
  
  # record the events whose time of appearance will be sampled
  net_events_to_sample = net_events[net_events$time %% 1 == 0, ] 
  net_events_to_sample$obs_window = net_events_to_sample$time
  infect_times = rep(-1, length = nrow(G0))
  infect_times[I0] = 0
  infect_times[events$per1[events$event == 1]] = events$time[events$event == 1]
  CD_type = ifelse(net_events_to_sample$event %in% 3: 5, 3, 6)
  # run iterations
  S = samples * thin + burn
  for(it in 1:S){
    # it = 1
    net_events_to_sample_temp = net_events_to_sample
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    
    times = NULL
    for(ix in 1:length(MR$recover)){
      # ix = 1
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      
      imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      
      times = c(times, imputed)
    }
    
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    # (2) compute event counts and summations
    #PA_test = parse_augment3_test(G0, I0, events, recover.dat)
    
    
    recov_times = rep(max(c(events$time, recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    # Create a matrix to store the bounds for each network times
    event_bounds = cbind.data.frame(lb = net_events_to_sample_temp$obs_window-window_length, 
                                    recov1 = recov_times[net_events_to_sample_temp$per1], 
                                    recov2 = recov_times[net_events_to_sample_temp$per2],
                                    ub = net_events_to_sample_temp$obs_window)
    rownames(event_bounds) = rownames(net_events_to_sample_temp)
    # Set the lower bounds of the network events to be the infection times of the first patient involved 
    # if the infection time is in the obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per1] & 
                       infect_times[net_events_to_sample_temp$per1] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per1[temp_idx]] 
    
    # Set the lower bounds of the network events to be the infection times of the second patient involved
    # if the second patient is infected later than the first patient in the same obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per2] & 
                       infect_times[net_events_to_sample_temp$per2] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per2[temp_idx]] 
    
    # Row sort the recovery times
    event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov1", "recov2")] = 
      event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov2", "recov1")]
    
    # Remove the recovery times if the recovery times are not in the corresponding obs window
    event_bounds$recov1[event_bounds$recov1 <= event_bounds$lb | event_bounds$recov1 >= event_bounds$ub] = NA
    event_bounds$recov2[event_bounds$recov2 <= event_bounds$lb | event_bounds$recov2 >= event_bounds$ub] = NA
    
    infect_times[infect_times < 0] = max(c(events$time, recover.dat$time))+1
    
    initial_link_type = rowSums(cbind(infect_times[net_events_to_sample_temp$per1] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per1], 
                                      infect_times[net_events_to_sample_temp$per2] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per2])) + CD_type
    
    # Find the number of NAs in each row. If there are two NAs, 
    # then the health status of the two persons are not changed within the interval. 
    # If there is one NA, the health status of one of the person changes within the interval.
    # If there are no NAs, the health status of both person changes within the interval.
    n_texp_dist = rowSums(is.na(event_bounds))
    
    idx1 = which(n_texp_dist == 2) # indices of the events with only one interval
    idx2 = which(n_texp_dist == 1) # indices of the events with two intervals
    idx3 = which(n_texp_dist == 0) # indices of the events with three intervals
    
    
    # net_events_to_sample_temp$time[idx1] = sam_texp(length(idx1),
    #                                                 params.cur[initial_link_type[idx1]],
    #                                                 a = event_bounds$lb[idx1],
    #                                                 b = event_bounds$ub[idx1])
    event_bounds1 = cbind(event_bounds$lb[idx1], event_bounds$ub[idx1])
    sampled_time1 = sam_texp(length(idx1),
                             params.cur[initial_link_type[idx1]],
                             a = event_bounds1[, 1],
                             b = event_bounds1[, 2])
    sampled_time1[initial_link_type[idx1] == 7] = (sampled_time1[initial_link_type[idx1] == 7] - event_bounds1[initial_link_type[idx1] == 7, 1]) + event_bounds1[initial_link_type[idx1] == 7, 1]
    net_events_to_sample_temp$time[idx1] = sampled_time1
    
    # idx1_sub = idx1[which(runif(length(idx1)) <= dpois(1, exp(-params.cur[initial_link_type[idx1]] * (event_bounds$ub[idx1] - net_events_to_sample_temp$time[idx1]))))]
    # (1 - exp(-params.cur[initial_link_type[idx1]] * (event_bounds$ub[idx1] - net_events_to_sample$time[idx1])))
    # new_events1 = net_event_imputer(input_event = net_events_to_sample_temp[idx1_sub, c("time", "event", "per1", "per2")],
    #                                 ub = event_bounds[idx1_sub, "ub"],
    #                                 params = params.cur)
    new_events1 = c()
    
    middle_bound = event_bounds[idx2, c("recov1", "recov2")]
    middle_bound[is.na(middle_bound)] = 0
    middle_bound = rowSums(middle_bound)
    event_bounds2 = cbind(event_bounds$lb[idx2], middle_bound, event_bounds$ub[idx2])
    net_events_to_sample_temp$time[idx2] = sam_texp2(n = length(idx2),
                                                     rate = cbind(params.cur[initial_link_type[idx2]],
                                                                  params.cur[initial_link_type[idx2]-1]),
                                                     a = event_bounds2[, 1],
                                                     b = event_bounds2[, 2],
                                                     c = event_bounds2[, 3])
    # Find the indices of the events happening after the middle bound (i.e. the connection happen after one of the patients is recovered). 
    # The event type is "minus 1" compared to the events happening before the middle bound
    idx2_A = idx2[which(middle_bound < net_events_to_sample_temp$time[idx2] & net_events_to_sample_temp$time[idx2] <= event_bounds$ub[idx2])]
    net_events_to_sample_temp$event[idx2_A] = net_events_to_sample_temp$event[idx2_A] - 1
    
    
    net_events_to_sample_temp$time[idx3] = sam_texp3(n = length(idx3), 
                                                     rate = cbind(params.cur[initial_link_type[idx3]], 
                                                                  params.cur[initial_link_type[idx3]-1], 
                                                                  params.cur[initial_link_type[idx3]-2]), 
                                                     a = event_bounds$lb[idx3], 
                                                     b = event_bounds$recov1[idx3], 
                                                     c = event_bounds$recov2[idx3], 
                                                     d = event_bounds$ub[idx3])
    idx3_A = idx3[which(event_bounds$recov1[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$recov2[idx3])]
    net_events_to_sample_temp$event[idx3_A] = net_events_to_sample_temp$event[idx3_A] - 1
    idx3_B = idx3[which(event_bounds$recov2[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$ub[idx3])]
    net_events_to_sample_temp$event[idx3_B] = net_events_to_sample_temp$event[idx3_B] - 2
    
    events[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    events[rownames(net_events_to_sample), "event"] = net_events_to_sample$event
    
    full_events = rbind(events,new_events1)
    full_events = full_events[order(full_events$time), ]
    
    # Make sure the event types of all events are correctly labeled
    net_event_idx = which(full_events$event %in% 3: 8)
    full_net_events = full_events[net_event_idx, ]
    full_CD_type = ifelse(full_net_events$event %in% 3: 5, 3, 6)
    
    full_net_events$event = (infect_times[full_net_events$per1] < full_net_events$time & full_net_events$time <= recov_times[full_net_events$per1]) + 
      (infect_times[full_net_events$per2] < full_net_events$time & full_net_events$time <= recov_times[full_net_events$per2]) + 
      ifelse(full_net_events$event %in% 3: 5, 3, 6)
    full_events[net_event_idx, ] = full_net_events
    
    PA2 = parse_augment2_test(G0, I0, full_events, recover.dat)
    event.counts = PA2$event.counts
    big.sums = PA2$big.sums
    
    # PA = parse_augment2_test(dats$G0, dats$I0, dats$events[dats$events$event != 2, ], dats$events[dats$events$event == 2, c(1, 3)])
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  # traceplot
  # modification: add in 95% posterior credible intervals
  # for(ix in 1:8){
  #   v = vars[ix]
  #   sams = params[,ix]
  #   if(!is.null(true.params)){
  #     truth = true.params[ix]
  #     yl = range(sams, truth)
  #     plot(sams~c(1:samples), main=paste("Traceplot for",v), 
  #          xlab = "sample", ylab = v, type="l", ylim = yl)
  #     abline(h=truth,col="red",lwd=2)
  #   }else{
  #     yl = range(sams)
  #     plot(sams~c(1:samples), main=paste("Traceplot for",v), 
  #          xlab = "sample", ylab = v, type="l", ylim = yl)
  #   }
  #   bounds = quantile(sams, c(.025, .975))
  #   abline(h=bounds[1],col="gray",lty=2,lwd=2)
  #   abline(h=bounds[2],col="gray",lty=2,lwd=2)
  # }
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  complete_events = rbind.data.frame(full_events, 
                                     cbind.data.frame(time = recover.dat$time, 
                                                      event = 2, 
                                                      per1 = recover.dat$per1, 
                                                      per2 = NA))
  complete_events = complete_events[order(complete_events$time), ]
  return(list(parm = params,
              imputed_events = complete_events))
}
infer_miss_recov4 = function(dat, priors, obs_time,
                             init.params = NULL,
                             verbose = T, plot = T, output.sams = 100, 
                             samples = 1000, burn = 100, thin = 1,
                             impute = "filter", model="SIR", 
                             timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  reports = dat$report; report.times = dat$report.times
  window_length = obs_time[2] - obs_time[1]
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, report = reports, times = report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  net_events = events[events$event %in% 3: 8, ]
  net_events[net_events$per2 < net_events$per1, c("per1", "per2")] = 
    net_events[net_events$per2 < net_events$per1, c("per2", "per1")]
  
  # record the events whose time of appearance will be sampled
  net_events_to_sample = net_events[net_events$time %% 1 == 0, ] 
  net_events_to_sample$obs_window = net_events_to_sample$time
  infect_times = rep(-1, length = nrow(G0))
  infect_times[I0] = 0
  infect_times[events$per1[events$event == 1]] = events$time[events$event == 1]
  
  net_snapshots = net_clipper(net_events, G0, obs_time)
  
  multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
  
  sick_pers = c(1: nrow(G0))[which(infect_times > 0)]
  
  
  # run iterations
  S = samples * thin + burn
  all_big_sums = matrix(0, nrow = S, ncol = 5)
  for(it in 1:S){
    # it = 1
    net_events_to_sample_temp = net_events_to_sample
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    
    times = NULL
    for(ix in 1:length(MR$recover)){
      # ix = 1
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      
      imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      
      times = c(times, imputed)
    }
    
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    # PA = parse_augment2_test(G0, I0, events, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # all_big_sums[it, 1] = big.sums[1]
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    recov_times = rep(max(c(events$time, recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    extra_events = c()
    multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
    for (i in 1: length(multiple_nei)) {
      # i = 1
      sick_per = as.integer(names(multiple_nei[i]))
      
      sick_nei_i = multiple_nei[[i]]
      for (j in 1: length(sick_nei_i)) {
        # j = 2
        sick_link = sort(c(sick_nei_i[j], sick_per))
        
        related_net_events = net_events[net_events$per1 == sick_link[1] & net_events$per2 == sick_link[2], ]
        
        t0 = obs_time[which(obs_time > infect_times[sick_per])[1] - 1]
        t1 = obs_time[which(obs_time > infect_times[sick_per])[1]]
        
        if (sum(related_net_events$time %% 1 != 0) == 1 & nrow(related_net_events) != 0){
          # Find the lastes time we observe this link being connected
          time_of_connection = ifelse(net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > infect_times[sick_per])[1] - 1],
                                      t0, related_net_events$time[related_net_events$time %% 1 != 0])
        }else{
          time_of_connection = t0
        }
        
        
        # If a disconnection is sampled before the infection
        if (runif(1) < exp(-(infect_times[sick_per] - time_of_connection)*params.cur[7])*(infect_times[sick_per] - time_of_connection)*params.cur[7]){
          extra_events = rbind.data.frame(extra_events,
                                          data.frame(time = sam_texp(1, params.cur[7], time_of_connection, infect_times[sick_per]),
                                                     event = 7,
                                                     per1 = sick_link[1],
                                                     per2 = sick_link[2]))
          
          if (net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > infect_times[sick_per])[1]]){ # If this link is still connected at t1, we will sample a connection back
            # Prepare the sampling bounds determined by the infection and recovery times and t1
            temp_bounds = c(infect_times[sick_per], sort(recov_times[sick_link]), t1) #
            temp_bounds = temp_bounds[temp_bounds[1] <= temp_bounds & temp_bounds <= temp_bounds[length(temp_bounds)]]
            
            
            
            # Prepare the parameters used for each interval
            parms_temp = params.cur[5: (5 - (length(temp_bounds)-2))]
            probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
            probs = probs/sum(probs)
            
            # Sample the interval the connection will take place
            interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
            
            # Generate the connection event
            extra_events = rbind.data.frame(extra_events,
                                            data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                       event = interval_idx + 3,
                                                       per1 = sick_link[1],
                                                       per2 = sick_link[2]))
          } else{ # If this link is not connected at t1, we remove it from net_events_to_sample since the disconnection has been sampled already
            net_events_to_sample_temp = net_events_to_sample_temp[-which(net_events_to_sample_temp$per1 == sick_link[1] & net_events_to_sample_temp$per2 == sick_link[2]), ]
          }
        }
      }
    }
    
    # Create a matrix to store the bounds for each network times
    event_bounds = cbind.data.frame(lb = net_events_to_sample_temp$obs_window-window_length, 
                                    recov1 = recov_times[net_events_to_sample_temp$per1], 
                                    recov2 = recov_times[net_events_to_sample_temp$per2],
                                    ub = net_events_to_sample_temp$obs_window)
    rownames(event_bounds) = rownames(net_events_to_sample_temp)
    
    # Set the lower bounds of the network events to be the infection times of the first patient involved 
    # if the infection time is in the obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per1] & 
                       infect_times[net_events_to_sample_temp$per1] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per1[temp_idx]] 
    
    # Set the lower bounds of the network events to be the infection times of the second patient involved
    # if the second patient is infected later than the first patient in the same obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per2] & 
                       infect_times[net_events_to_sample_temp$per2] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per2[temp_idx]] 
    
    # Row sort the recovery times
    event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov1", "recov2")] = 
      event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov2", "recov1")]
    
    # Remove the recovery times if the recovery times are not in the corresponding obs window
    event_bounds$recov1[event_bounds$recov1 <= event_bounds$lb | event_bounds$recov1 >= event_bounds$ub] = NA
    event_bounds$recov2[event_bounds$recov2 <= event_bounds$lb | event_bounds$recov2 >= event_bounds$ub] = NA
    
    infect_times[infect_times < 0] = max(c(events$time, recover.dat$time))+1
    
    # Find the type of links at the lower bounds
    initial_link_type = rowSums(cbind(infect_times[net_events_to_sample_temp$per1] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per1], 
                                      infect_times[net_events_to_sample_temp$per2] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per2])) + 
      ifelse(net_events_to_sample_temp$event %in% 3: 5, 3, 6)
    
    # Find the number of NAs in each row. If there are two NAs, 
    # then the health status of the two persons are not changed within the interval. 
    # If there is one NA, the health status of one of the person changes within the interval.
    # If there are no NAs, the health status of both person changes within the interval.
    n_texp_dist = rowSums(is.na(event_bounds))
    
    idx1 = which(n_texp_dist == 2) # indices of the events with only one interval
    idx2 = which(n_texp_dist == 1) # indices of the events with two intervals
    idx3 = which(n_texp_dist == 0) # indices of the events with three intervals
    
    
    event_bounds1 = cbind(event_bounds$lb[idx1], event_bounds$ub[idx1])
    sampled_time1 = sam_texp(length(idx1),
                             params.cur[initial_link_type[idx1]],
                             a = event_bounds1[, 1],
                             b = event_bounds1[, 2])
    sampled_time1[initial_link_type[idx1] == 7] = (sampled_time1[initial_link_type[idx1] == 7] - event_bounds1[initial_link_type[idx1] == 7, 1]) + event_bounds1[initial_link_type[idx1] == 7, 1]
    net_events_to_sample_temp$time[idx1] = sampled_time1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 2] = big.sums[1]
    
    
    middle_bound = event_bounds[idx2, c("recov1", "recov2")]
    middle_bound[is.na(middle_bound)] = 0
    middle_bound = rowSums(middle_bound)
    event_bounds2 = cbind(event_bounds$lb[idx2], middle_bound, event_bounds$ub[idx2])
    net_events_to_sample_temp$time[idx2] = sam_texp2(n = length(idx2),
                                                     rate = cbind(params.cur[initial_link_type[idx2]],
                                                                  params.cur[initial_link_type[idx2]-1]),
                                                     a = event_bounds2[, 1],
                                                     b = event_bounds2[, 2],
                                                     c = event_bounds2[, 3])
    # Find the indices of the events happening after the middle bound (i.e. the connection happen after one of the patients is recovered). 
    # The event type is "minus 1" compared to the events happening before the middle bound
    idx2_A = idx2[which(middle_bound < net_events_to_sample_temp$time[idx2] & net_events_to_sample_temp$time[idx2] <= event_bounds$ub[idx2])]
    net_events_to_sample_temp$event[idx2_A] = net_events_to_sample_temp$event[idx2_A] - 1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 3] = big.sums[1]
    
    net_events_to_sample_temp$time[idx3] = sam_texp3(n = length(idx3), 
                                                     rate = cbind(params.cur[initial_link_type[idx3]], 
                                                                  params.cur[initial_link_type[idx3]-1], 
                                                                  params.cur[initial_link_type[idx3]-2]), 
                                                     a = event_bounds$lb[idx3], 
                                                     b = event_bounds$recov1[idx3], 
                                                     c = event_bounds$recov2[idx3], 
                                                     d = event_bounds$ub[idx3])
    idx3_A = idx3[which(event_bounds$recov1[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$recov2[idx3])]
    net_events_to_sample_temp$event[idx3_A] = net_events_to_sample_temp$event[idx3_A] - 1
    idx3_B = idx3[which(event_bounds$recov2[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$ub[idx3])]
    net_events_to_sample_temp$event[idx3_B] = net_events_to_sample_temp$event[idx3_B] - 2
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 4] = big.sums[1]
    
    events_temp = events
    events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    events_temp = rbind.data.frame(events_temp, extra_events)
    events_temp = events_temp[order(events_temp$time), ]
    
    PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    # all_big_sums[it, 5] = big.sums[1]
    # print(all_big_sums[it, ])
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  # complete_events = rbind.data.frame(full_events, 
  #                                    cbind.data.frame(time = recover.dat$time, 
  #                                                     event = 2, 
  #                                                     per1 = recover.dat$per1, 
  #                                                     per2 = NA))
  # complete_events = complete_events[order(complete_events$time), ]
  # return(list(params, all_big_sums))
  return(params)
}
infer_miss_recov5 = function(dat, priors, obs_time,
                             init.params = NULL,
                             verbose = T, plot = T, output.sams = 100, 
                             samples = 1000, burn = 100, thin = 1,
                             impute = "filter", model="SIR", 
                             timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  reports = dat$report; report.times = dat$report.times
  window_length = obs_time[2] - obs_time[1]
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, report = reports, times = report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  net_events = events[events$event %in% 3: 8, ]
  net_events[net_events$per2 < net_events$per1, c("per1", "per2")] = 
    net_events[net_events$per2 < net_events$per1, c("per2", "per1")]
  
  # record the events whose time of appearance will be sampled
  net_events_to_sample = net_events[net_events$time %% 1 == 0, ] 
  net_events_to_sample$obs_window = net_events_to_sample$time
  infect_times = rep(-1, length = nrow(G0))
  infect_times[I0] = 0
  infect_times[events$per1[events$event == 1]] = events$time[events$event == 1]
  
  net_snapshots = net_clipper(net_events, G0, obs_time)
  
  multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
  
  sick_pers = c(1: nrow(G0))[which(infect_times > 0)]
  
  G = as.matrix(G0)
  event_temp = events
  event_temp = rbind.data.frame(data.frame(time = 0, 
                                           event = 1,
                                           per1 = I0, 
                                           per2 = NA), event_temp)
  
  infect_times[infect_times < 0] = max(events$time)+1
  
  SI_links_to_break = list()
  for (i in 1: nrow(event_temp)) {
    # i = 1
    events_i = event_temp[i, ]
    sick_link = as.integer(events_i[3: 4])
    if (events_i$event == 1){
      sick_per = events_i$per1
      t0 = obs_time[which(obs_time > infect_times[sick_per])[1] - 1]
      t1 = obs_time[which(obs_time > infect_times[sick_per])[1]]
      victims = which(G[events_i$per1, ] == 1)
      victims = victims[infect_times[victims] >= t1]
      if (sum(infect_times[victims] >= t1) > 0){
        for (j in 1: length(victims)) {
          # j = 1
          sus_per = victims[j]
          # Find the link status at each observation time point
          link_status = net_snapshots[sick_per, sus_per, ]
          
          # Find the observation intervals where the link remains connected
          sampling_intl = matrix(cbind(obs_time[-length(obs_time)], obs_time[-1])[link_status[-length(link_status)] == 1 & link_status[-1] == 1, ], ncol = 2)
          
          # Remove the observation intervals earlier than the sick person getting sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl <= events_i$time) < 2, ], ncol = 2)
          # Change the lower bound of the interval where the sick person becomes sick to be the sick time
          sampling_intl[rowSums(sampling_intl <= events_i$time) == 1, 1] = events_i$time
          # Remove the intervals when and after the susceptible person becomes sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl < infect_times[sus_per]) == 2, ], ncol = 2)
          
          if (nrow(sampling_intl) >= 1){
            sampling_intl = list(sampling_intl)
            names(sampling_intl) = paste0(sick_per, "-", sus_per)
            # strsplit(names(break_points), split = "-") %>% unlist() %>% as.integer()
            SI_links_to_break = append(SI_links_to_break, sampling_intl)
          }
        }
      }
    }
    if (events_i$event %in% 3: 5){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 1
    }
    if (events_i$event %in% 6: 8){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 0
    }
  }
  
  infect_times[infect_times == max(events$time)+1] = -1
  # run iterations
  S = samples * thin + burn
  # all_big_sums = matrix(0, nrow = S, ncol = 5)
  for(it in 1:S){
    # it = 1
    net_events_to_sample_temp = net_events_to_sample
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    
    times = NULL
    for(ix in 1:length(MR$recover)){
      # ix = 1
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      
      imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      
      times = c(times, imputed)
    }
    
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    # PA = parse_augment2_test(G0, I0, events, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # all_big_sums[it, 1] = big.sums[1]
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    recov_times = rep(max(c(events$time, recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    extra_events = c()
    multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
    if (length(multiple_nei) > 0){
      for (i in 1: length(multiple_nei)) {
        # i = 2
        sick_per = as.integer(names(multiple_nei[i]))
        
        sick_nei_i = multiple_nei[[i]]
        for (j in 1: length(sick_nei_i)) {
          # j = 2
          sick_link = sort(c(sick_nei_i[j], sick_per))
          
          related_net_events = net_events[net_events$per1 == sick_link[1] & net_events$per2 == sick_link[2], ]
          
          t0 = obs_time[which(obs_time > infect_times[sick_per])[1] - 1]
          t1 = obs_time[which(obs_time > infect_times[sick_per])[1]]
          
          if (sum(related_net_events$time %% 1 != 0) == 1 & nrow(related_net_events) != 0){
            # Find the lastes time we observe this link being connected
            time_of_connection = ifelse(net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > infect_times[sick_per])[1] - 1],
                                        t0, related_net_events$time[related_net_events$time %% 1 != 0])
          }else{
            time_of_connection = t0
          }
          
          
          # If a disconnection is sampled before the infection
          if (runif(1) < exp(-(infect_times[sick_per] - time_of_connection)*params.cur[7])*(infect_times[sick_per] - time_of_connection)*params.cur[7]){
            extra_events = rbind.data.frame(extra_events,
                                            data.frame(time = sam_texp(1, params.cur[7], time_of_connection, infect_times[sick_per]),
                                                       event = 7,
                                                       per1 = sick_link[1],
                                                       per2 = sick_link[2]))
            
            if (net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > infect_times[sick_per])[1]]){ # If this link is still connected at t1, we will sample a connection back
              # Prepare the sampling bounds determined by the infection and recovery times and t1
              temp_bounds = c(infect_times[sick_per], sort(recov_times[sick_link]), t1) #
              temp_bounds = temp_bounds[temp_bounds[1] <= temp_bounds & temp_bounds <= temp_bounds[length(temp_bounds)]]
              
              
              
              # Prepare the parameters used for each interval
              parms_temp = params.cur[5: (5 - (length(temp_bounds)-2))]
              probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
              probs = probs/sum(probs)
              
              # Sample the interval the connection will take place
              interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
              
              # Generate the connection event
              extra_events = rbind.data.frame(extra_events,
                                              data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                         event = interval_idx + 3,
                                                         per1 = sick_link[1],
                                                         per2 = sick_link[2]))
            } else{ # If this link is not connected at t1, we remove it from net_events_to_sample since the disconnection has been sampled already
              net_events_to_sample_temp = net_events_to_sample_temp[-which(net_events_to_sample_temp$per1 == sick_link[1] & net_events_to_sample_temp$per2 == sick_link[2]), ]
            }
          }
        }
      }
    }
    
    extra_events2 = c()
    if (length(SI_links_to_break) > 0){
      for (i in 1: length(SI_links_to_break)) {
        # i = 92
        sampling_intl = SI_links_to_break[[i]]
        linked_pers = strsplit(names(SI_links_to_break[i]), split = "-") %>% unlist() %>% as.integer()
        sick_per = linked_pers[1]
        sus_per = linked_pers[2]
        
        
        # Remove the observation intervals later than the sick person recovered
        sampling_intl = matrix(sampling_intl[rowSums(sampling_intl >= recov_times[sick_per]) < 2, ], ncol = 2)
        # Change the lower bound of the interval where the sick person becomes sick to be the sick time
        sampling_intl[rowSums(sampling_intl >= recov_times[sick_per]) == 1, 2] = recov_times[sick_per]
        if (nrow(sampling_intl) > 0){
          for (j in 1: nrow(sampling_intl)) {
            # j = 1
            if (runif(1) < exp(-params.cur[7] * (sampling_intl[j, 2] - sampling_intl[j, 1])) * params.cur[7]*(sampling_intl[j, 2] - sampling_intl[j, 1])){
              new_disconn = data.frame(time = sam_texp(1, params.cur[7], sampling_intl[j, 1], sampling_intl[j, 2]),
                                       event = 7,
                                       per1 = sort(linked_pers)[1],
                                       per2 = sort(linked_pers)[2])
              extra_events2 = rbind.data.frame(extra_events2, new_disconn)
              
              if (!sampling_intl[j, 2] %in% obs_time){ # If the sick person recovers before t1
                # Prepare the sampling bounds determined by the disconnection and t1
                t1 = obs_time[which(obs_time > sampling_intl[j, 2])[1]]
                temp_bounds = c(new_disconn$time, sampling_intl[j, 2], t1) 
                
                # Prepare the parameters used for each interval
                parms_temp = params.cur[c(4, 3)]
                probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
                probs = probs/sum(probs)
                
                # Sample the interval the connection will take place
                interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
                
                # Generate the connection event
                extra_events2 = rbind.data.frame(extra_events2,
                                                 data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                            event = c(4, 3)[interval_idx],
                                                            per1 = sort(linked_pers)[1],
                                                            per2 = sort(linked_pers)[2]))
              } else{ # If the sick person does not recover before t1 
                extra_events2 = rbind.data.frame(extra_events2,
                                                 data.frame(time = sam_texp(1, params.cur[4], new_disconn$time, sampling_intl[j, 2]),
                                                            event = 4,
                                                            per1 = sort(linked_pers)[1],
                                                            per2 = sort(linked_pers)[2]))
              }
              
            }
          }
        }
      }
    }
    
    
    # Create a matrix to store the bounds for each network times
    event_bounds = cbind.data.frame(lb = net_events_to_sample_temp$obs_window-window_length, 
                                    recov1 = recov_times[net_events_to_sample_temp$per1], 
                                    recov2 = recov_times[net_events_to_sample_temp$per2],
                                    ub = net_events_to_sample_temp$obs_window)
    rownames(event_bounds) = rownames(net_events_to_sample_temp)
    
    # Set the lower bounds of the network events to be the infection times of the first patient involved 
    # if the infection time is in the obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per1] & 
                       infect_times[net_events_to_sample_temp$per1] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per1[temp_idx]] 
    
    # Set the lower bounds of the network events to be the infection times of the second patient involved
    # if the second patient is infected later than the first patient in the same obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per2] & 
                       infect_times[net_events_to_sample_temp$per2] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per2[temp_idx]] 
    
    # Row sort the recovery times
    event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov1", "recov2")] = 
      event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov2", "recov1")]
    
    # Remove the recovery times if the recovery times are not in the corresponding obs window
    event_bounds$recov1[event_bounds$recov1 <= event_bounds$lb | event_bounds$recov1 >= event_bounds$ub] = NA
    event_bounds$recov2[event_bounds$recov2 <= event_bounds$lb | event_bounds$recov2 >= event_bounds$ub] = NA
    
    infect_times[infect_times < 0] = max(c(events$time, recover.dat$time))+1
    
    # Find the type of links at the lower bounds
    initial_link_type = rowSums(cbind(infect_times[net_events_to_sample_temp$per1] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per1], 
                                      infect_times[net_events_to_sample_temp$per2] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per2])) + 
      ifelse(net_events_to_sample_temp$event %in% 3: 5, 3, 6)
    
    # Find the number of NAs in each row. If there are two NAs, 
    # then the health status of the two persons are not changed within the interval. 
    # If there is one NA, the health status of one of the person changes within the interval.
    # If there are no NAs, the health status of both person changes within the interval.
    n_texp_dist = rowSums(is.na(event_bounds))
    
    idx1 = which(n_texp_dist == 2) # indices of the events with only one interval
    idx2 = which(n_texp_dist == 1) # indices of the events with two intervals
    idx3 = which(n_texp_dist == 0) # indices of the events with three intervals
    
    
    event_bounds1 = cbind(event_bounds$lb[idx1], event_bounds$ub[idx1])
    sampled_time1 = sam_texp(length(idx1),
                             params.cur[initial_link_type[idx1]],
                             a = event_bounds1[, 1],
                             b = event_bounds1[, 2])
    sampled_time1[initial_link_type[idx1] == 7] = (sampled_time1[initial_link_type[idx1] == 7] - event_bounds1[initial_link_type[idx1] == 7, 1]) + event_bounds1[initial_link_type[idx1] == 7, 1]
    net_events_to_sample_temp$time[idx1] = sampled_time1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 2] = big.sums[1]
    
    
    middle_bound = event_bounds[idx2, c("recov1", "recov2")]
    middle_bound[is.na(middle_bound)] = 0
    middle_bound = rowSums(middle_bound)
    event_bounds2 = cbind(event_bounds$lb[idx2], middle_bound, event_bounds$ub[idx2])
    net_events_to_sample_temp$time[idx2] = sam_texp2(n = length(idx2),
                                                     rate = cbind(params.cur[initial_link_type[idx2]],
                                                                  params.cur[initial_link_type[idx2]-1]),
                                                     a = event_bounds2[, 1],
                                                     b = event_bounds2[, 2],
                                                     c = event_bounds2[, 3])
    # Find the indices of the events happening after the middle bound (i.e. the connection happen after one of the patients is recovered). 
    # The event type is "minus 1" compared to the events happening before the middle bound
    idx2_A = idx2[which(middle_bound < net_events_to_sample_temp$time[idx2] & net_events_to_sample_temp$time[idx2] <= event_bounds$ub[idx2])]
    net_events_to_sample_temp$event[idx2_A] = net_events_to_sample_temp$event[idx2_A] - 1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 3] = big.sums[1]
    
    net_events_to_sample_temp$time[idx3] = sam_texp3(n = length(idx3), 
                                                     rate = cbind(params.cur[initial_link_type[idx3]], 
                                                                  params.cur[initial_link_type[idx3]-1], 
                                                                  params.cur[initial_link_type[idx3]-2]), 
                                                     a = event_bounds$lb[idx3], 
                                                     b = event_bounds$recov1[idx3], 
                                                     c = event_bounds$recov2[idx3], 
                                                     d = event_bounds$ub[idx3])
    idx3_A = idx3[which(event_bounds$recov1[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$recov2[idx3])]
    net_events_to_sample_temp$event[idx3_A] = net_events_to_sample_temp$event[idx3_A] - 1
    idx3_B = idx3[which(event_bounds$recov2[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$ub[idx3])]
    net_events_to_sample_temp$event[idx3_B] = net_events_to_sample_temp$event[idx3_B] - 2
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 4] = big.sums[1]
    # 
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # events_temp = rbind.data.frame(events_temp, extra_events)
    # events_temp = events_temp[order(events_temp$time), ]
    # events_temp = events_temp[events_temp$time %% 1 != 0, ]
    # 
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    
    # 
    # 
    # all_big_sums[it, 4] = big.sums[1]
    
    events_temp = events
    events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    events_temp = rbind.data.frame(events_temp, extra_events, extra_events2)
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[events_temp$time %% 1 != 0, ]
    
    PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    # all_big_sums[it, 5] = big.sums[1]
    # print(all_big_sums[it, ])
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  # complete_events = rbind.data.frame(full_events, 
  #                                    cbind.data.frame(time = recover.dat$time, 
  #                                                     event = 2, 
  #                                                     per1 = recover.dat$per1, 
  #                                                     per2 = NA))
  # complete_events = complete_events[order(complete_events$time), ]
  # return(list(params, all_big_sums))
  return(params)
}
infer_miss_recov5_test = function(dat, priors, obs_time,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  impute = "filter", model="SIR", 
                                  timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  reports = dat$report; report.times = dat$report.times
  window_length = obs_time[2] - obs_time[1]
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, report = reports, times = report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  net_events = events[events$event %in% 3: 8, ]
  net_events[net_events$per2 < net_events$per1, c("per1", "per2")] = 
    net_events[net_events$per2 < net_events$per1, c("per2", "per1")]
  
  # record the events whose time of appearance will be sampled
  net_events_to_sample = net_events[net_events$time %% 1 == 0, ] 
  net_events_to_sample$obs_window = net_events_to_sample$time
  infect_times = rep(-1, length = nrow(G0))
  infect_times[I0] = 0
  infect_times[events$per1[events$event == 1]] = events$time[events$event == 1]
  
  net_snapshots = net_clipper(net_events, G0, obs_time)
  
  multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
  
  sick_pers = c(1: nrow(G0))[which(infect_times > 0)]
  
  G = as.matrix(G0)
  event_temp = events
  event_temp = rbind.data.frame(data.frame(time = 0, 
                                           event = 1,
                                           per1 = I0, 
                                           per2 = NA), event_temp)
  
  infect_times[infect_times < 0] = max(events$time)+1
  
  SI_links_to_break = list()
  for (i in 1: nrow(event_temp)) {
    # i = 1
    events_i = event_temp[i, ]
    sick_link = as.integer(events_i[3: 4])
    if (events_i$event == 1){
      sick_per = events_i$per1
      t0 = obs_time[which(obs_time > infect_times[sick_per])[1] - 1]
      t1 = obs_time[which(obs_time > infect_times[sick_per])[1]]
      victims = which(G[events_i$per1, ] == 1)
      victims = victims[infect_times[victims] >= t1]
      if (sum(infect_times[victims] >= t1) > 0){
        for (j in 1: length(victims)) {
          # j = 1
          sus_per = victims[j]
          # Find the link status at each observation time point
          link_status = net_snapshots[sick_per, sus_per, ]
          
          # Find the observation intervals where the link remains connected
          sampling_intl = matrix(cbind(obs_time[-length(obs_time)], obs_time[-1])[link_status[-length(link_status)] == 1 & link_status[-1] == 1, ], ncol = 2)
          
          # Remove the observation intervals earlier than the sick person getting sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl <= events_i$time) < 2, ], ncol = 2)
          # Change the lower bound of the interval where the sick person becomes sick to be the sick time
          sampling_intl[rowSums(sampling_intl <= events_i$time) == 1, 1] = events_i$time
          # Remove the intervals when and after the susceptible person becomes sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl < infect_times[sus_per]) == 2, ], ncol = 2)
          
          if (nrow(sampling_intl) >= 1){
            sampling_intl = list(sampling_intl)
            names(sampling_intl) = paste0(sick_per, "-", sus_per)
            # strsplit(names(break_points), split = "-") %>% unlist() %>% as.integer()
            SI_links_to_break = append(SI_links_to_break, sampling_intl)
          }
        }
      }
    }
    if (events_i$event %in% 3: 5){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 1
    }
    if (events_i$event %in% 6: 8){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 0
    }
  }
  
  infect_times[infect_times == max(events$time)+1] = -1
  # run iterations
  S = samples * thin + burn
  # all_big_sums = matrix(0, nrow = S, ncol = 5)
  for(it in 1:S){
    # it = 1
    net_events_to_sample_temp = net_events_to_sample
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    
    times = NULL
    for(ix in 1:length(MR$recover)){
      # ix = 1
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      
      imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      
      times = c(times, imputed)
    }
    
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    # PA = parse_augment2_test(G0, I0, events, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # all_big_sums[it, 1] = big.sums[1]
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    recov_times = rep(max(c(events$time, recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    extra_events = c()
    multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
    if (length(multiple_nei) > 0){
      for (i in 1: length(multiple_nei)) {
        # i = 2
        sick_per = as.integer(names(multiple_nei[i]))
        
        sick_nei_i = multiple_nei[[i]]
        for (j in 1: length(sick_nei_i)) {
          # j = 2
          sick_link = sort(c(sick_nei_i[j], sick_per))
          
          related_net_events = net_events[net_events$per1 == sick_link[1] & net_events$per2 == sick_link[2], ]
          
          t0 = obs_time[which(obs_time > infect_times[sick_per])[1] - 1]
          t1 = obs_time[which(obs_time > infect_times[sick_per])[1]]
          
          if (sum(related_net_events$time %% 1 != 0) == 1 & nrow(related_net_events) != 0){
            # Find the lastes time we observe this link being connected
            time_of_connection = ifelse(net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > infect_times[sick_per])[1] - 1],
                                        t0, related_net_events$time[related_net_events$time %% 1 != 0])
          }else{
            time_of_connection = t0
          }
          
          
          # If a disconnection is sampled before the infection
          if (runif(1) < 1 - exp(-(infect_times[sick_per] - time_of_connection)*params.cur[7])){
            extra_events = rbind.data.frame(extra_events,
                                            data.frame(time = sam_texp(1, params.cur[7], time_of_connection, infect_times[sick_per]),
                                                       event = 7,
                                                       per1 = sick_link[1],
                                                       per2 = sick_link[2]))
            
            if (net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > infect_times[sick_per])[1]]){ # If this link is still connected at t1, we will sample a connection back
              # Prepare the sampling bounds determined by the infection and recovery times and t1
              temp_bounds = c(infect_times[sick_per], sort(recov_times[sick_link]), t1) #
              temp_bounds = temp_bounds[temp_bounds[1] <= temp_bounds & temp_bounds <= temp_bounds[length(temp_bounds)]]
              
              
              
              # Prepare the parameters used for each interval
              parms_temp = params.cur[5: (5 - (length(temp_bounds)-2))]
              probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
              probs = probs/sum(probs)
              
              # Sample the interval the connection will take place
              interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
              
              # Generate the connection event
              extra_events = rbind.data.frame(extra_events,
                                              data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                         event = interval_idx + 3,
                                                         per1 = sick_link[1],
                                                         per2 = sick_link[2]))
            } else{ # If this link is not connected at t1, we remove it from net_events_to_sample since the disconnection has been sampled already
              net_events_to_sample_temp = net_events_to_sample_temp[-which(net_events_to_sample_temp$per1 == sick_link[1] & net_events_to_sample_temp$per2 == sick_link[2]), ]
            }
          }
        }
      }
    }
    
    extra_events2 = c()
    if (length(SI_links_to_break) > 0){
      for (i in 1: length(SI_links_to_break)) {
        # i = 92
        sampling_intl = SI_links_to_break[[i]]
        linked_pers = strsplit(names(SI_links_to_break[i]), split = "-") %>% unlist() %>% as.integer()
        sick_per = linked_pers[1]
        sus_per = linked_pers[2]
        
        
        # Remove the observation intervals later than the sick person recovered
        sampling_intl = matrix(sampling_intl[rowSums(sampling_intl >= recov_times[sick_per]) < 2, ], ncol = 2)
        # Change the lower bound of the interval where the sick person becomes sick to be the sick time
        sampling_intl[rowSums(sampling_intl >= recov_times[sick_per]) == 1, 2] = recov_times[sick_per]
        if (nrow(sampling_intl) > 0){
          for (j in 1: nrow(sampling_intl)) {
            # j = 1
            if (runif(1) < 1 - exp(-params.cur[7] * (sampling_intl[j, 2] - sampling_intl[j, 1])) ){
              new_disconn = data.frame(time = sam_texp(1, params.cur[7], sampling_intl[j, 1], sampling_intl[j, 2]),
                                       event = 7,
                                       per1 = sort(linked_pers)[1],
                                       per2 = sort(linked_pers)[2])
              extra_events2 = rbind.data.frame(extra_events2, new_disconn)
              
              if (!sampling_intl[j, 2] %in% obs_time){ # If the sick person recovers before t1
                # Prepare the sampling bounds determined by the disconnection and t1
                t1 = obs_time[which(obs_time > sampling_intl[j, 2])[1]]
                temp_bounds = c(new_disconn$time, sampling_intl[j, 2], t1) 
                
                # Prepare the parameters used for each interval
                parms_temp = params.cur[c(4, 3)]
                probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
                probs = probs/sum(probs)
                
                # Sample the interval the connection will take place
                interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
                
                # Generate the connection event
                extra_events2 = rbind.data.frame(extra_events2,
                                                 data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                            event = c(4, 3)[interval_idx],
                                                            per1 = sort(linked_pers)[1],
                                                            per2 = sort(linked_pers)[2]))
              } else{ # If the sick person does not recover before t1 
                extra_events2 = rbind.data.frame(extra_events2,
                                                 data.frame(time = sam_texp(1, params.cur[4], new_disconn$time, sampling_intl[j, 2]),
                                                            event = 4,
                                                            per1 = sort(linked_pers)[1],
                                                            per2 = sort(linked_pers)[2]))
              }
              
            }
          }
        }
      }
    }
    
    
    # Create a matrix to store the bounds for each network times
    event_bounds = cbind.data.frame(lb = net_events_to_sample_temp$obs_window-window_length, 
                                    recov1 = recov_times[net_events_to_sample_temp$per1], 
                                    recov2 = recov_times[net_events_to_sample_temp$per2],
                                    ub = net_events_to_sample_temp$obs_window)
    rownames(event_bounds) = rownames(net_events_to_sample_temp)
    
    # Set the lower bounds of the network events to be the infection times of the first patient involved 
    # if the infection time is in the obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per1] & 
                       infect_times[net_events_to_sample_temp$per1] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per1[temp_idx]] 
    
    # Set the lower bounds of the network events to be the infection times of the second patient involved
    # if the second patient is infected later than the first patient in the same obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per2] & 
                       infect_times[net_events_to_sample_temp$per2] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per2[temp_idx]] 
    
    # Row sort the recovery times
    event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov1", "recov2")] = 
      event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov2", "recov1")]
    
    # Remove the recovery times if the recovery times are not in the corresponding obs window
    event_bounds$recov1[event_bounds$recov1 <= event_bounds$lb | event_bounds$recov1 >= event_bounds$ub] = NA
    event_bounds$recov2[event_bounds$recov2 <= event_bounds$lb | event_bounds$recov2 >= event_bounds$ub] = NA
    
    infect_times[infect_times < 0] = max(c(events$time, recover.dat$time))+1
    
    # Find the type of links at the lower bounds
    initial_link_type = rowSums(cbind(infect_times[net_events_to_sample_temp$per1] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per1], 
                                      infect_times[net_events_to_sample_temp$per2] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per2])) + 
      ifelse(net_events_to_sample_temp$event %in% 3: 5, 3, 6)
    
    # Find the number of NAs in each row. If there are two NAs, 
    # then the health status of the two persons are not changed within the interval. 
    # If there is one NA, the health status of one of the person changes within the interval.
    # If there are no NAs, the health status of both person changes within the interval.
    n_texp_dist = rowSums(is.na(event_bounds))
    
    idx1 = which(n_texp_dist == 2) # indices of the events with only one interval
    idx2 = which(n_texp_dist == 1) # indices of the events with two intervals
    idx3 = which(n_texp_dist == 0) # indices of the events with three intervals
    
    
    event_bounds1 = cbind(event_bounds$lb[idx1], event_bounds$ub[idx1])
    sampled_time1 = sam_texp(length(idx1),
                             params.cur[initial_link_type[idx1]],
                             a = event_bounds1[, 1],
                             b = event_bounds1[, 2])
    sampled_time1[initial_link_type[idx1] == 7] = (sampled_time1[initial_link_type[idx1] == 7] - event_bounds1[initial_link_type[idx1] == 7, 1]) + event_bounds1[initial_link_type[idx1] == 7, 1]
    net_events_to_sample_temp$time[idx1] = sampled_time1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 2] = big.sums[1]
    
    
    middle_bound = event_bounds[idx2, c("recov1", "recov2")]
    middle_bound[is.na(middle_bound)] = 0
    middle_bound = rowSums(middle_bound)
    event_bounds2 = cbind(event_bounds$lb[idx2], middle_bound, event_bounds$ub[idx2])
    net_events_to_sample_temp$time[idx2] = sam_texp2(n = length(idx2),
                                                     rate = cbind(params.cur[initial_link_type[idx2]],
                                                                  params.cur[initial_link_type[idx2]-1]),
                                                     a = event_bounds2[, 1],
                                                     b = event_bounds2[, 2],
                                                     c = event_bounds2[, 3])
    # Find the indices of the events happening after the middle bound (i.e. the connection happen after one of the patients is recovered). 
    # The event type is "minus 1" compared to the events happening before the middle bound
    idx2_A = idx2[which(middle_bound < net_events_to_sample_temp$time[idx2] & net_events_to_sample_temp$time[idx2] <= event_bounds$ub[idx2])]
    net_events_to_sample_temp$event[idx2_A] = net_events_to_sample_temp$event[idx2_A] - 1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 3] = big.sums[1]
    
    net_events_to_sample_temp$time[idx3] = sam_texp3(n = length(idx3), 
                                                     rate = cbind(params.cur[initial_link_type[idx3]], 
                                                                  params.cur[initial_link_type[idx3]-1], 
                                                                  params.cur[initial_link_type[idx3]-2]), 
                                                     a = event_bounds$lb[idx3], 
                                                     b = event_bounds$recov1[idx3], 
                                                     c = event_bounds$recov2[idx3], 
                                                     d = event_bounds$ub[idx3])
    idx3_A = idx3[which(event_bounds$recov1[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$recov2[idx3])]
    net_events_to_sample_temp$event[idx3_A] = net_events_to_sample_temp$event[idx3_A] - 1
    idx3_B = idx3[which(event_bounds$recov2[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$ub[idx3])]
    net_events_to_sample_temp$event[idx3_B] = net_events_to_sample_temp$event[idx3_B] - 2
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 4] = big.sums[1]
    # 
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # events_temp = rbind.data.frame(events_temp, extra_events)
    # events_temp = events_temp[order(events_temp$time), ]
    # events_temp = events_temp[events_temp$time %% 1 != 0, ]
    # 
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    
    # 
    # 
    # all_big_sums[it, 4] = big.sums[1]
    
    events_temp = events
    events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    events_temp = rbind.data.frame(events_temp, extra_events, extra_events2)
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[events_temp$time %% 1 != 0, ]
    
    PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    # all_big_sums[it, 5] = big.sums[1]
    # print(all_big_sums[it, ])
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  # complete_events = rbind.data.frame(full_events, 
  #                                    cbind.data.frame(time = recover.dat$time, 
  #                                                     event = 2, 
  #                                                     per1 = recover.dat$per1, 
  #                                                     per2 = NA))
  # complete_events = complete_events[order(complete_events$time), ]
  # return(list(params, all_big_sums))
  return(params)
}
infer_miss_recov6 = function(dat, priors, obs_time,
                             init.params = NULL,
                             verbose = T, plot = T, output.sams = 100, 
                             samples = 1000, burn = 100, thin = 1,
                             impute = "filter", model="SIR", 
                             timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  reports = dat$report; report.times = dat$report.times
  window_length = obs_time[2] - obs_time[1]
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, report = reports, times = report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  net_events = events[events$event %in% 3: 8, ]
  net_events[net_events$per2 < net_events$per1, c("per1", "per2")] = 
    net_events[net_events$per2 < net_events$per1, c("per2", "per1")]
  
  # record the events whose time of appearance will be sampled
  net_events_to_sample = net_events[net_events$time %% 1 == 0, ] 
  net_events_to_sample$obs_window = net_events_to_sample$time
  infect_times = rep(-1, length = nrow(G0))
  infect_times[I0] = 0
  infect_times[events$per1[events$event == 1]] = events$time[events$event == 1]
  
  net_snapshots = net_clipper(net_events, G0, obs_time)
  
  multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
  
  sick_pers = c(1: nrow(G0))[which(infect_times > 0)]
  
  G = as.matrix(G0)
  event_temp = events
  event_temp = rbind.data.frame(data.frame(time = 0, 
                                           event = 1,
                                           per1 = I0, 
                                           per2 = NA), event_temp)
  
  infect_times[infect_times < 0] = max(events$time)+1
  
  SI_links_to_break = list()
  for (i in 1: nrow(event_temp)) {
    # i = 1
    events_i = event_temp[i, ]
    sick_link = as.integer(events_i[3: 4])
    if (events_i$event == 1){
      sick_per = events_i$per1
      t0 = obs_time[which(obs_time > infect_times[sick_per])[1] - 1]
      t1 = obs_time[which(obs_time > infect_times[sick_per])[1]]
      victims = which(G[events_i$per1, ] == 1)
      victims = victims[infect_times[victims] >= t1]
      if (sum(infect_times[victims] >= t1) > 0){
        for (j in 1: length(victims)) {
          # j = 1
          sus_per = victims[j]
          # Find the link status at each observation time point
          link_status = net_snapshots[sick_per, sus_per, ]
          
          # Find the observation intervals where the link remains connected
          sampling_intl = matrix(cbind(obs_time[-length(obs_time)], obs_time[-1])[link_status[-length(link_status)] == 1 & link_status[-1] == 1, ], ncol = 2)
          
          # Remove the observation intervals earlier than the sick person getting sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl <= events_i$time) < 2, ], ncol = 2)
          # Change the lower bound of the interval where the sick person becomes sick to be the sick time
          sampling_intl[rowSums(sampling_intl <= events_i$time) == 1, 1] = events_i$time
          # Remove the intervals when and after the susceptible person becomes sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl < infect_times[sus_per]) == 2, ], ncol = 2)
          
          if (nrow(sampling_intl) >= 1){
            sampling_intl = list(sampling_intl)
            names(sampling_intl) = paste0(sick_per, "-", sus_per)
            # strsplit(names(break_points), split = "-") %>% unlist() %>% as.integer()
            SI_links_to_break = append(SI_links_to_break, sampling_intl)
          }
        }
      }
    }
    if (events_i$event %in% 3: 5){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 1
    }
    if (events_i$event %in% 6: 8){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 0
    }
  }
  
  infect_times[infect_times == max(events$time)+1] = -1
  # run iterations
  S = samples * thin + burn
  # all_big_sums = matrix(0, nrow = S, ncol = 5)
  for(it in 1:S){
    # it = 1
    net_events_to_sample_temp = net_events_to_sample
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    
    times = NULL
    for(ix in 1:length(MR$recover)){
      # ix = 1
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      
      imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      
      times = c(times, imputed)
    }
    
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    # PA = parse_augment2_test(G0, I0, events, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # all_big_sums[it, 1] = big.sums[1]
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    recov_times = rep(max(c(events$time, recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    extra_events = c()
    multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
    if (length(multiple_nei) > 0){
      for (i in 1: length(multiple_nei)) {
        # i = 2
        sick_per = as.integer(names(multiple_nei[i]))
        
        sick_nei_i = multiple_nei[[i]]
        for (j in 1: length(sick_nei_i)) {
          # j = 2
          sick_link = sort(c(sick_nei_i[j], sick_per))
          
          related_net_events = net_events[net_events$per1 == sick_link[1] & net_events$per2 == sick_link[2], ]
          
          t0 = obs_time[which(obs_time > infect_times[sick_per])[1] - 1]
          t1 = obs_time[which(obs_time > infect_times[sick_per])[1]]
          
          if (sum(related_net_events$time %% 1 != 0) == 1 & nrow(related_net_events) != 0){
            # Find the lastes time we observe this link being connected
            time_of_connection = ifelse(net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > infect_times[sick_per])[1] - 1],
                                        t0, related_net_events$time[related_net_events$time %% 1 != 0])
          }else{
            time_of_connection = t0
          }
          
          
          # If a disconnection is sampled before the infection
          if (runif(1) < 1 - exp(-(infect_times[sick_per] - time_of_connection)*params.cur[7])){
            extra_events = rbind.data.frame(extra_events,
                                            data.frame(time = sam_texp(1, params.cur[7], time_of_connection, infect_times[sick_per]),
                                                       event = 7,
                                                       per1 = sick_link[1],
                                                       per2 = sick_link[2]))
            
            if (net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > infect_times[sick_per])[1]]){ # If this link is still connected at t1, we will sample a connection back
              # Prepare the sampling bounds determined by the infection and recovery times and t1
              temp_bounds = c(infect_times[sick_per], sort(recov_times[sick_link]), t1) #
              temp_bounds = temp_bounds[temp_bounds[1] <= temp_bounds & temp_bounds <= temp_bounds[length(temp_bounds)]]
              
              
              
              # Prepare the parameters used for each interval
              parms_temp = params.cur[5: (5 - (length(temp_bounds)-2))]
              probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
              probs = probs/sum(probs)
              
              # Sample the interval the connection will take place
              interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
              
              # Generate the connection event
              extra_events = rbind.data.frame(extra_events,
                                              data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                         event = interval_idx + 3,
                                                         per1 = sick_link[1],
                                                         per2 = sick_link[2]))
            } else{ # If this link is not connected at t1, we remove it from net_events_to_sample since the disconnection has been sampled already
              net_events_to_sample_temp = net_events_to_sample_temp[-which(net_events_to_sample_temp$per1 == sick_link[1] & net_events_to_sample_temp$per2 == sick_link[2]), ]
            }
          }
        }
      }
    }
    
    extra_events2 = c()
    if (length(SI_links_to_break) > 0){
      for (i in 1: length(SI_links_to_break)) {
        # i = 92
        sampling_intl = SI_links_to_break[[i]]
        linked_pers = strsplit(names(SI_links_to_break[i]), split = "-") %>% unlist() %>% as.integer()
        sick_per = linked_pers[1]
        sus_per = linked_pers[2]
        
        
        # Remove the observation intervals later than the sick person recovered
        sampling_intl = matrix(sampling_intl[rowSums(sampling_intl >= recov_times[sick_per]) < 2, ], ncol = 2)
        # Change the lower bound of the interval where the sick person becomes sick to be the sick time
        sampling_intl[rowSums(sampling_intl >= recov_times[sick_per]) == 1, 2] = recov_times[sick_per]
        if (nrow(sampling_intl) > 0){
          for (j in 1: nrow(sampling_intl)) {
            # j = 1
            if (runif(1) < 1 - exp(-(infect_times[sick_per] - time_of_connection)*params.cur[7])){
              new_disconn = data.frame(time = sam_texp(1, params.cur[7], sampling_intl[j, 1], sampling_intl[j, 2]),
                                       event = 7,
                                       per1 = sort(linked_pers)[1],
                                       per2 = sort(linked_pers)[2])
              extra_events2 = rbind.data.frame(extra_events2, new_disconn)
              
              if (!sampling_intl[j, 2] %in% obs_time){ # If the sick person recovers before t1
                # Prepare the sampling bounds determined by the disconnection and t1
                t1 = obs_time[which(obs_time > sampling_intl[j, 2])[1]]
                temp_bounds = c(new_disconn$time, sampling_intl[j, 2], t1) 
                
                # Prepare the parameters used for each interval
                parms_temp = params.cur[c(4, 3)]
                probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
                probs = probs/sum(probs)
                
                # Sample the interval the connection will take place
                interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
                
                # Generate the connection event
                extra_events2 = rbind.data.frame(extra_events2,
                                                 data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                            event = c(4, 3)[interval_idx],
                                                            per1 = sort(linked_pers)[1],
                                                            per2 = sort(linked_pers)[2]))
              } else{ # If the sick person does not recover before t1 
                extra_events2 = rbind.data.frame(extra_events2,
                                                 data.frame(time = sam_texp(1, params.cur[4], new_disconn$time, sampling_intl[j, 2]),
                                                            event = 4,
                                                            per1 = sort(linked_pers)[1],
                                                            per2 = sort(linked_pers)[2]))
              }
              
            }
          }
        }
      }
    }
    
    
    # Create a matrix to store the bounds for each network times
    event_bounds = cbind.data.frame(lb = net_events_to_sample_temp$obs_window-window_length, 
                                    recov1 = recov_times[net_events_to_sample_temp$per1], 
                                    recov2 = recov_times[net_events_to_sample_temp$per2],
                                    ub = net_events_to_sample_temp$obs_window)
    rownames(event_bounds) = rownames(net_events_to_sample_temp)
    
    # Set the lower bounds of the network events to be the infection times of the first patient involved 
    # if the infection time is in the obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per1] & 
                       infect_times[net_events_to_sample_temp$per1] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per1[temp_idx]] 
    
    # Set the lower bounds of the network events to be the infection times of the second patient involved
    # if the second patient is infected later than the first patient in the same obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per2] & 
                       infect_times[net_events_to_sample_temp$per2] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per2[temp_idx]] 
    
    # Row sort the recovery times
    event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov1", "recov2")] = 
      event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov2", "recov1")]
    
    # Remove the recovery times if the recovery times are not in the corresponding obs window
    event_bounds$recov1[event_bounds$recov1 <= event_bounds$lb | event_bounds$recov1 >= event_bounds$ub] = NA
    event_bounds$recov2[event_bounds$recov2 <= event_bounds$lb | event_bounds$recov2 >= event_bounds$ub] = NA
    
    infect_times[infect_times < 0] = max(c(events$time, recover.dat$time))+1
    
    # Find the type of links at the lower bounds
    initial_link_type = rowSums(cbind(infect_times[net_events_to_sample_temp$per1] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per1], 
                                      infect_times[net_events_to_sample_temp$per2] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per2])) + 
      ifelse(net_events_to_sample_temp$event %in% 3: 5, 3, 6)
    
    # Find the number of NAs in each row. If there are two NAs, 
    # then the health status of the two persons are not changed within the interval. 
    # If there is one NA, the health status of one of the person changes within the interval.
    # If there are no NAs, the health status of both person changes within the interval.
    n_texp_dist = rowSums(is.na(event_bounds))
    
    idx1 = which(n_texp_dist == 2) # indices of the events with only one interval
    idx2 = which(n_texp_dist == 1) # indices of the events with two intervals
    idx3 = which(n_texp_dist == 0) # indices of the events with three intervals
    
    
    event_bounds1 = cbind(event_bounds$lb[idx1], event_bounds$ub[idx1])
    sampled_time1 = sam_texp(length(idx1),
                             params.cur[initial_link_type[idx1]],
                             a = event_bounds1[, 1],
                             b = event_bounds1[, 2])
    sampled_time1[initial_link_type[idx1] == 7] = (sampled_time1[initial_link_type[idx1] == 7] - event_bounds1[initial_link_type[idx1] == 7, 1]) + event_bounds1[initial_link_type[idx1] == 7, 1]
    net_events_to_sample_temp$time[idx1] = sampled_time1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 2] = big.sums[1]
    
    
    middle_bound = event_bounds[idx2, c("recov1", "recov2")]
    middle_bound[is.na(middle_bound)] = 0
    middle_bound = rowSums(middle_bound)
    event_bounds2 = cbind(event_bounds$lb[idx2], middle_bound, event_bounds$ub[idx2])
    net_events_to_sample_temp$time[idx2] = sam_texp2(n = length(idx2),
                                                     rate = cbind(params.cur[initial_link_type[idx2]],
                                                                  params.cur[initial_link_type[idx2]-1]),
                                                     a = event_bounds2[, 1],
                                                     b = event_bounds2[, 2],
                                                     c = event_bounds2[, 3])
    # Find the indices of the events happening after the middle bound (i.e. the connection happen after one of the patients is recovered). 
    # The event type is "minus 1" compared to the events happening before the middle bound
    idx2_A = idx2[which(middle_bound < net_events_to_sample_temp$time[idx2] & net_events_to_sample_temp$time[idx2] <= event_bounds$ub[idx2])]
    net_events_to_sample_temp$event[idx2_A] = net_events_to_sample_temp$event[idx2_A] - 1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 3] = big.sums[1]
    
    net_events_to_sample_temp$time[idx3] = sam_texp3(n = length(idx3), 
                                                     rate = cbind(params.cur[initial_link_type[idx3]], 
                                                                  params.cur[initial_link_type[idx3]-1], 
                                                                  params.cur[initial_link_type[idx3]-2]), 
                                                     a = event_bounds$lb[idx3], 
                                                     b = event_bounds$recov1[idx3], 
                                                     c = event_bounds$recov2[idx3], 
                                                     d = event_bounds$ub[idx3])
    idx3_A = idx3[which(event_bounds$recov1[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$recov2[idx3])]
    net_events_to_sample_temp$event[idx3_A] = net_events_to_sample_temp$event[idx3_A] - 1
    idx3_B = idx3[which(event_bounds$recov2[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$ub[idx3])]
    net_events_to_sample_temp$event[idx3_B] = net_events_to_sample_temp$event[idx3_B] - 2
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 4] = big.sums[1]
    # 
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # events_temp = rbind.data.frame(events_temp, extra_events)
    # events_temp = events_temp[order(events_temp$time), ]
    # events_temp = events_temp[events_temp$time %% 1 != 0, ]
    # 
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    
    # 
    # 
    # all_big_sums[it, 4] = big.sums[1]
    
    events_temp = events
    events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    events_temp = rbind.data.frame(events_temp, extra_events, extra_events2)
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[events_temp$time %% 1 != 0, ]
    
    PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    # all_big_sums[it, 5] = big.sums[1]
    # print(all_big_sums[it, ])
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  # complete_events = rbind.data.frame(full_events, 
  #                                    cbind.data.frame(time = recover.dat$time, 
  #                                                     event = 2, 
  #                                                     per1 = recover.dat$per1, 
  #                                                     per2 = NA))
  # complete_events = complete_events[order(complete_events$time), ]
  # return(list(params, all_big_sums))
  return(params)
}
infer_miss_recov6_test = function(dat, priors, obs_time,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  impute = "filter", model="SIR", 
                                  timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  reports = dat$report; report.times = dat$report.times
  window_length = obs_time[2] - obs_time[1]
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, report = reports, times = report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  net_events = events[events$event %in% 3: 8, ]
  net_events[net_events$per2 < net_events$per1, c("per1", "per2")] = 
    net_events[net_events$per2 < net_events$per1, c("per2", "per1")]
  
  # record the events whose time of appearance will be sampled
  net_events_to_sample = net_events[net_events$time %% 1 == 0, ] 
  net_events_to_sample$obs_window = net_events_to_sample$time
  infect_times = rep(-1, length = nrow(G0))
  infect_times[I0] = 0
  infect_times[events$per1[events$event == 1]] = events$time[events$event == 1]
  
  net_snapshots = net_clipper(net_events, G0, obs_time)
  
  multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
  
  sick_pers = c(1: nrow(G0))[which(infect_times > 0)]
  
  G = as.matrix(G0)
  event_temp = events
  event_temp = rbind.data.frame(data.frame(time = 0, 
                                           event = 1,
                                           per1 = I0, 
                                           per2 = NA), event_temp)
  
  infect_times[infect_times < 0] = max(events$time)+1
  
  SI_links_to_break = list()
  for (i in 1: nrow(event_temp)) {
    # i = 1
    events_i = event_temp[i, ]
    sick_link = as.integer(events_i[3: 4])
    if (events_i$event == 1){
      sick_per = events_i$per1
      t0 = obs_time[which(obs_time > infect_times[sick_per])[1] - 1]
      t1 = obs_time[which(obs_time > infect_times[sick_per])[1]]
      victims = which(G[events_i$per1, ] == 1)
      victims = victims[infect_times[victims] >= t1]
      if (sum(infect_times[victims] >= t1) > 0){
        for (j in 1: length(victims)) {
          # j = 1
          sus_per = victims[j]
          # Find the link status at each observation time point
          link_status = net_snapshots[sick_per, sus_per, ]
          
          # Find the observation intervals where the link remains connected
          sampling_intl = matrix(cbind(obs_time[-length(obs_time)], obs_time[-1])[link_status[-length(link_status)] == 1 & link_status[-1] == 1, ], ncol = 2)
          
          # Remove the observation intervals earlier than the sick person getting sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl <= events_i$time) < 2, ], ncol = 2)
          # Change the lower bound of the interval where the sick person becomes sick to be the sick time
          sampling_intl[rowSums(sampling_intl <= events_i$time) == 1, 1] = events_i$time
          # Remove the intervals when and after the susceptible person becomes sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl < infect_times[sus_per]) == 2, ], ncol = 2)
          
          if (nrow(sampling_intl) >= 1){
            sampling_intl = list(sampling_intl)
            names(sampling_intl) = paste0(sick_per, "-", sus_per)
            # strsplit(names(break_points), split = "-") %>% unlist() %>% as.integer()
            SI_links_to_break = append(SI_links_to_break, sampling_intl)
          }
        }
      }
    }
    if (events_i$event %in% 3: 5){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 1
    }
    if (events_i$event %in% 6: 8){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 0
    }
  }
  
  infect_times[infect_times == max(events$time)+1] = -1
  # run iterations
  S = samples * thin + burn
  # all_big_sums = matrix(0, nrow = S, ncol = 5)
  for(it in 1:S){
    # it = 1
    net_events_to_sample_temp = net_events_to_sample
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    
    times = NULL
    for(ix in 1:length(MR$recover)){
      # ix = 2
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      
      imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      
      times = c(times, imputed)
    }
    
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    # PA = parse_augment2_test(G0, I0, events, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # all_big_sums[it, 1] = big.sums[1]
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    recov_times = rep(max(c(events$time, recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    extra_events = c()
    multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
    if (length(multiple_nei) > 0){
      for (i in 1: length(multiple_nei)) {
        # i = 1
        sick_per = as.integer(names(multiple_nei[i]))
        
        sick_nei_i = multiple_nei[[i]]
        for (j in 1: length(sick_nei_i)) {
          # j = 3
          sick_link = sort(c(sick_nei_i[j], sick_per))
          
          related_net_events = net_events[net_events$per1 == sick_link[1] & net_events$per2 == sick_link[2], ]
          
          t0 = obs_time[which(obs_time > infect_times[sick_per])[1] - 1]
          t1 = obs_time[which(obs_time > infect_times[sick_per])[1]]
          
          if (sum(related_net_events$time %% 1 != 0) == 1 & nrow(related_net_events) != 0){
            # Find the lastes time we observe this link being connected
            time_of_connection = ifelse(net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > infect_times[sick_per])[1] - 1],
                                        t0, related_net_events$time[related_net_events$time %% 1 != 0])
          }else{
            time_of_connection = t0
          }
          time_of_SI = max(time_of_connection, infect_times[sick_nei_i[j]])
          
          if (recov_times[sick_nei_i[j]] > infect_times[sick_per]){
            # If a disconnection is sampled before the infection
            if (runif(1) < 1 - exp(-(infect_times[sick_per] - time_of_SI)*params.cur[7])){
              
              new_disconn = data.frame(time = sam_texp(1, params.cur[7], time_of_SI, infect_times[sick_per]),
                                       event = 7,
                                       per1 = sick_link[1],
                                       per2 = sick_link[2])
              extra_events = rbind.data.frame(extra_events, new_disconn)
              
              if (net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > infect_times[sick_per])[1]]){ # If this link is still connected at t1, we will sample a connection back
                # Prepare the sampling bounds determined by the time of disconnection, infection and recovery times and t1
                temp_bounds = c(new_disconn$time, infect_times[sick_per], sort(recov_times[sick_link]), t1) #
                temp_bounds = temp_bounds[temp_bounds[1] <= temp_bounds & temp_bounds <= temp_bounds[length(temp_bounds)]]
                
                
                
                # Prepare the parameters used for each interval
                
                parms_temp = params.cur[c(4, 5: (5 - (length(temp_bounds[-1])-2)))]
                probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
                probs = probs/sum(probs)
                
                # Sample the interval the connection will take place
                interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
                
                # Generate the connection event
                extra_events = rbind.data.frame(extra_events,
                                                data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                           event = interval_idx + 3,
                                                           per1 = sick_link[1],
                                                           per2 = sick_link[2]))
              } else{ # If this link is not connected at t1, we remove it from net_events_to_sample since the disconnection has been sampled already
                net_events_to_sample_temp = net_events_to_sample_temp[-which(net_events_to_sample_temp$time == t1 &
                                                                               net_events_to_sample_temp$per1 == sick_link[1] & 
                                                                               net_events_to_sample_temp$per2 == sick_link[2]), ]
              }
            }
          }
        }
      }
    }
    
    extra_events2 = c()
    if (length(SI_links_to_break) > 0){
      for (i in 1: length(SI_links_to_break)) {
        # i = 1
        sampling_intl = SI_links_to_break[[i]]
        linked_pers = strsplit(names(SI_links_to_break[i]), split = "-") %>% unlist() %>% as.integer()
        sick_per = linked_pers[1]
        sus_per = linked_pers[2]
        
        
        # Remove the observation intervals later than the sick person recovered
        sampling_intl = matrix(sampling_intl[rowSums(sampling_intl >= recov_times[sick_per]) < 2, ], ncol = 2)
        # Change the lower bound of the interval where the sick person becomes sick to be the sick time
        sampling_intl[rowSums(sampling_intl >= recov_times[sick_per]) == 1, 2] = recov_times[sick_per]
        if (nrow(sampling_intl) > 0){
          for (j in 1: nrow(sampling_intl)) {
            # j = 1
            
            if (runif(1) < 1 - exp(-(sampling_intl[j, 2] - sampling_intl[j, 1])*params.cur[7])){
              new_disconn = data.frame(time = sam_texp(1, params.cur[7], sampling_intl[j, 1], sampling_intl[j, 2]),
                                       event = 7,
                                       per1 = sort(linked_pers)[1],
                                       per2 = sort(linked_pers)[2])
              extra_events2 = rbind.data.frame(extra_events2, new_disconn)
              
              if (!sampling_intl[j, 2] %in% obs_time){ # If the sick person recovers before t1
                # Prepare the sampling bounds determined by the disconnection and t1
                t1 = obs_time[which(obs_time > sampling_intl[j, 2])[1]]
                temp_bounds = c(new_disconn$time, sampling_intl[j, 2], t1) 
                
                # Prepare the parameters used for each interval
                parms_temp = params.cur[c(4, 3)]
                probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
                probs = probs/sum(probs)
                
                # Sample the interval the connection will take place
                interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
                
                # Generate the connection event
                extra_events2 = rbind.data.frame(extra_events2,
                                                 data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                            event = c(4, 3)[interval_idx],
                                                            per1 = sort(linked_pers)[1],
                                                            per2 = sort(linked_pers)[2]))
              } else{ # If the sick person does not recover before t1 
                extra_events2 = rbind.data.frame(extra_events2,
                                                 data.frame(time = sam_texp(1, params.cur[4], new_disconn$time, sampling_intl[j, 2]),
                                                            event = 4,
                                                            per1 = sort(linked_pers)[1],
                                                            per2 = sort(linked_pers)[2]))
              }
              
            }
          }
        }
      }
    }
    
    
    # Create a matrix to store the bounds for each network times
    event_bounds = cbind.data.frame(lb = net_events_to_sample_temp$obs_window-window_length, 
                                    recov1 = recov_times[net_events_to_sample_temp$per1], 
                                    recov2 = recov_times[net_events_to_sample_temp$per2],
                                    ub = net_events_to_sample_temp$obs_window)
    rownames(event_bounds) = rownames(net_events_to_sample_temp)
    
    # Set the lower bounds of the network events to be the infection times of the first patient involved 
    # if the infection time is in the obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per1] & 
                       infect_times[net_events_to_sample_temp$per1] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per1[temp_idx]] 
    
    # Set the lower bounds of the network events to be the infection times of the second patient involved
    # if the second patient is infected later than the first patient in the same obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per2] & 
                       infect_times[net_events_to_sample_temp$per2] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per2[temp_idx]] 
    
    # Row sort the recovery times
    event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov1", "recov2")] = 
      event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov2", "recov1")]
    
    # Remove the recovery times if the recovery times are not in the corresponding obs window
    event_bounds$recov1[event_bounds$recov1 <= event_bounds$lb | event_bounds$recov1 >= event_bounds$ub] = NA
    event_bounds$recov2[event_bounds$recov2 <= event_bounds$lb | event_bounds$recov2 >= event_bounds$ub] = NA
    
    infect_times[infect_times < 0] = max(c(events$time, recover.dat$time))+1
    
    # Find the type of links at the lower bounds
    initial_link_type = rowSums(cbind(infect_times[net_events_to_sample_temp$per1] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per1], 
                                      infect_times[net_events_to_sample_temp$per2] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per2])) + 
      ifelse(net_events_to_sample_temp$event %in% 3: 5, 3, 6)
    
    # Find the number of NAs in each row. If there are two NAs, 
    # then the health status of the two persons are not changed within the interval. 
    # If there is one NA, the health status of one of the person changes within the interval.
    # If there are no NAs, the health status of both person changes within the interval.
    n_texp_dist = rowSums(is.na(event_bounds))
    
    idx1 = which(n_texp_dist == 2) # indices of the events with only one interval
    idx2 = which(n_texp_dist == 1) # indices of the events with two intervals
    idx3 = which(n_texp_dist == 0) # indices of the events with three intervals
    
    
    event_bounds1 = cbind(event_bounds$lb[idx1], event_bounds$ub[idx1])
    sampled_time1 = sam_texp(length(idx1),
                             params.cur[initial_link_type[idx1]],
                             a = event_bounds1[, 1],
                             b = event_bounds1[, 2])
    sampled_time1[initial_link_type[idx1] == 7] = (sampled_time1[initial_link_type[idx1] == 7] - event_bounds1[initial_link_type[idx1] == 7, 1]) + event_bounds1[initial_link_type[idx1] == 7, 1]
    net_events_to_sample_temp$time[idx1] = sampled_time1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 2] = big.sums[1]
    
    
    middle_bound = event_bounds[idx2, c("recov1", "recov2")]
    middle_bound[is.na(middle_bound)] = 0
    middle_bound = rowSums(middle_bound)
    event_bounds2 = cbind(event_bounds$lb[idx2], middle_bound, event_bounds$ub[idx2])
    net_events_to_sample_temp$time[idx2] = sam_texp2(n = length(idx2),
                                                     rate = cbind(params.cur[initial_link_type[idx2]],
                                                                  params.cur[initial_link_type[idx2]-1]),
                                                     a = event_bounds2[, 1],
                                                     b = event_bounds2[, 2],
                                                     c = event_bounds2[, 3])
    # Find the indices of the events happening after the middle bound (i.e. the connection happen after one of the patients is recovered). 
    # The event type is "minus 1" compared to the events happening before the middle bound
    idx2_A = idx2[which(middle_bound < net_events_to_sample_temp$time[idx2] & net_events_to_sample_temp$time[idx2] <= event_bounds$ub[idx2])]
    net_events_to_sample_temp$event[idx2_A] = net_events_to_sample_temp$event[idx2_A] - 1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 3] = big.sums[1]
    
    net_events_to_sample_temp$time[idx3] = sam_texp3(n = length(idx3), 
                                                     rate = cbind(params.cur[initial_link_type[idx3]], 
                                                                  params.cur[initial_link_type[idx3]-1], 
                                                                  params.cur[initial_link_type[idx3]-2]), 
                                                     a = event_bounds$lb[idx3], 
                                                     b = event_bounds$recov1[idx3], 
                                                     c = event_bounds$recov2[idx3], 
                                                     d = event_bounds$ub[idx3])
    idx3_A = idx3[which(event_bounds$recov1[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$recov2[idx3])]
    net_events_to_sample_temp$event[idx3_A] = net_events_to_sample_temp$event[idx3_A] - 1
    idx3_B = idx3[which(event_bounds$recov2[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$ub[idx3])]
    net_events_to_sample_temp$event[idx3_B] = net_events_to_sample_temp$event[idx3_B] - 2
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 4] = big.sums[1]
    # 
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # events_temp = rbind.data.frame(events_temp, extra_events)
    # events_temp = events_temp[order(events_temp$time), ]
    # events_temp = events_temp[events_temp$time %% 1 != 0, ]
    # 
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    
    # 
    # 
    # all_big_sums[it, 4] = big.sums[1]
    
    events_temp = events
    events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    events_temp = rbind.data.frame(events_temp, extra_events, extra_events2)
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[events_temp$time %% 1 != 0, ]
    
    PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    # all_big_sums[it, 5] = big.sums[1]
    # print(all_big_sums[it, ])
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  # complete_events = rbind.data.frame(full_events, 
  #                                    cbind.data.frame(time = recover.dat$time, 
  #                                                     event = 2, 
  #                                                     per1 = recover.dat$per1, 
  #                                                     per2 = NA))
  # complete_events = complete_events[order(complete_events$time), ]
  # return(list(params, all_big_sums))
  return(params)
}
# 02/06/23 Added by Houjie to impute the network once
infer_miss_recov7_test = function(dat, priors, obs_time,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  impute = "filter", model="SIR", 
                                  timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  reports = dat$report; report.times = dat$report.times
  window_length = obs_time[2] - obs_time[1]
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, report = reports, times = report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  net_events = events[events$event %in% 3: 8, ]
  net_events[net_events$per2 < net_events$per1, c("per1", "per2")] = 
    net_events[net_events$per2 < net_events$per1, c("per2", "per1")]
  
  # record the events whose time of appearance will be sampled
  net_events_to_sample = net_events[net_events$time %% 1 == 0, ] 
  net_events_to_sample$obs_window = net_events_to_sample$time
  infect_times = rep(-1, length = nrow(G0))
  infect_times[I0] = 0
  infect_times[events$per1[events$event == 1]] = events$time[events$event == 1]
  
  net_snapshots = net_clipper(net_events, G0, obs_time)
  
  multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
  
  sick_pers = c(1: nrow(G0))[which(infect_times > 0)]
  
  G = as.matrix(G0)
  event_temp = events
  event_temp = rbind.data.frame(data.frame(time = 0, 
                                           event = 1,
                                           per1 = I0, 
                                           per2 = NA), event_temp)
  
  infect_times[infect_times < 0] = max(events$time)+1
  
  SI_links_to_break = list()
  for (i in 1: nrow(event_temp)) {
    # i = 1
    events_i = event_temp[i, ]
    sick_link = as.integer(events_i[3: 4])
    if (events_i$event == 1){
      sick_per = events_i$per1
      i1 = infect_times[sick_per]
      t0 = obs_time[which(obs_time > i1)[1] - 1]
      t1 = obs_time[which(obs_time > i1)[1]]
      victims = which(G[events_i$per1, ] == 1)
      victims = victims[infect_times[victims] >= t1]
      if (sum(infect_times[victims] >= t1) > 0){
        for (j in 1: length(victims)) {
          # j = 1
          sus_per = victims[j]
          # Find the link status at each observation time point
          link_status = net_snapshots[sick_per, sus_per, ]
          
          # Find the observation intervals where the link remains connected
          sampling_intl = matrix(cbind(obs_time[-length(obs_time)], obs_time[-1])[link_status[-length(link_status)] == 1 & link_status[-1] == 1, ], ncol = 2)
          
          # Remove the observation intervals earlier than the sick person getting sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl <= events_i$time) < 2, ], ncol = 2)
          # Change the lower bound of the interval where the sick person becomes sick to be the sick time
          sampling_intl[rowSums(sampling_intl <= events_i$time) == 1, 1] = events_i$time
          # Remove the intervals when and after the susceptible person becomes sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl < infect_times[sus_per]) == 2, ], ncol = 2)
          
          if (nrow(sampling_intl) >= 1){
            sampling_intl = list(sampling_intl)
            names(sampling_intl) = paste0(sick_per, "-", sus_per)
            # strsplit(names(break_points), split = "-") %>% unlist() %>% as.integer()
            SI_links_to_break = append(SI_links_to_break, sampling_intl)
          }
        }
      }
    }
    if (events_i$event %in% 3: 5){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 1
    }
    if (events_i$event %in% 6: 8){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 0
    }
  }
  
  infect_times[infect_times == max(events$time)+1] = -1
  # run iterations
  S = samples * thin + burn
  # all_big_sums = matrix(0, nrow = S, ncol = 5)
  for(it in 1:S){
    # it = 1
    net_events_to_sample_temp = net_events_to_sample
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    
    times = NULL
    for(ix in 1:length(MR$recover)){
      # ix = 2
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      
      imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      
      times = c(times, imputed)
    }
    
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    # PA = parse_augment2_test(G0, I0, events, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # all_big_sums[it, 1] = big.sums[1]
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    recov_times = rep(max(c(events$time, recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    extra_events = c()
    multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
    if (length(multiple_nei) > 0){
      for (i in 1: length(multiple_nei)) {
        # i = 1
        sick_per = as.integer(names(multiple_nei[i]))
        i1 = infect_times[sick_per] # the sick time of this person
        r1 = recov_times[sick_per]
        sick_nei_i = multiple_nei[[i]]
        for (j in 1: length(sick_nei_i)) {
          
          # j = 1
          sick_link = sort(c(sick_nei_i[j], sick_per))
          i2 = infect_times[sick_nei_i[j]]; r2 = recov_times[sick_nei_i[j]] # the sick and reovery time of the neighbour
          related_net_events = net_events[net_events$per1 == sick_link[1] & net_events$per2 == sick_link[2], ]
          
          # End points of the observation interval
          t0 = obs_time[which(obs_time > i1)[1] - 1]
          t1 = obs_time[which(obs_time > i1)[1]]
          
          if (sum(related_net_events$time %% 1 != 0) == 1 & nrow(related_net_events) != 0){
            # Find the latest time we observe this link being connected
            tc = ifelse(net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > i1)[1] - 1],
                        t0, related_net_events$time[related_net_events$time %% 1 != 0])
          }else{
            tc = t0
          }
          
          # Sampling disconnection of SI links
          new_event = FALSE
          SI_end = min(i1, r2)
          if ((i2 <= tc & r2 >= i1) | (i2 <= tc & tc <= r2 & r2 <= i1)){
            if (runif(1) < pexp((SI_end - tc), rate = params.cur[7])){ # disconnection happening in (tc, min(i1, r2)]
              new_event = TRUE
              new_disconn = data.frame(time = sam_texp(1,  params.cur[7], tc, SI_end),
                                       event = 7,
                                       per1 = sick_link[1],
                                       per2 = sick_link[2])
              extra_events = rbind.data.frame(extra_events, new_disconn)
            }
          }
          
          if ((tc < i2 & i1 <= r2) | (tc < i2 & r2 <= i1)){
            if (runif(1) < pexp((i2 - tc), rate = params.cur[6])){ # disconnection happening before i2 in (tc, i2] \cup (i2, min(i1, r2)]
              new_event = TRUE
              new_disconn = data.frame(time = sam_texp(1, params.cur[6], tc, i2),
                                       event = 6,
                                       per1 = sick_link[1],
                                       per2 = sick_link[2])
              extra_events = rbind.data.frame(extra_events, new_disconn)
            } else if (runif(1) < pexp((SI_end - i2), rate = params.cur[7])){ # disconnection happening after i2 in (tc, i2] \cup (i2, min(i1, r2)]
              new_event = TRUE
              new_disconn = data.frame(time = sam_texp(1,  params.cur[7], i2, SI_end),
                                       event = 7,
                                       per1 = sick_link[1],
                                       per2 = sick_link[2])
              extra_events = rbind.data.frame(extra_events, new_disconn)
            }
            
          }
          
          if(new_event){
            if(net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > i1)[1]]){# If this link is still connected at t1, we will sample a connection back
              # Prepare the sampling bounds determined by the time of disconnection, i1, r1, i2, r2 and t1
              conn_intl = rbind(c(new_disconn$time, i2, r2, i1, r1, t1), c(0, 1, -1, 1, -1, 0))
              colnames(conn_intl) = c("td", "i2", "r2", "i1", "r1", "t1")
              conn_intl = conn_intl[, conn_intl[1, 1] <= conn_intl[1, ] & conn_intl[1, ] <= conn_intl[1, ncol(conn_intl)]] # Filter our the times not in (td, t1]
              conn_intl = conn_intl[, order(conn_intl[1, ])]
              conn_intl[2, 1] = new_disconn$event - 6 + 3 # Find the type of the disconnection, which is the initial type of the connection
              conn_type = cumsum(conn_intl[2, -ncol(conn_intl)]) # Find the connection types w.r.t each interval
              
              # Prepare the parameters used for each interval
              temp_bounds = conn_intl[1, ]
              parms_temp = params.cur[conn_type]
              probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
              probs = probs/sum(probs)
              
              # Sample the interval the connection will take place
              interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
              
              # Generate the connection event
              extra_events = rbind.data.frame(extra_events,
                                              data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                         event = conn_type[interval_idx],
                                                         per1 = sick_link[1],
                                                         per2 = sick_link[2]))
            }else{# If this link is not connected at t1, we remove it from net_events_to_sample since the disconnection has been sampled already
              net_events_to_sample_temp = net_events_to_sample_temp[-which(net_events_to_sample_temp$time == t1 &
                                                                             net_events_to_sample_temp$per1 == sick_link[1] & 
                                                                             net_events_to_sample_temp$per2 == sick_link[2]), ]
              
            }
          }
          
          
        }
      }
    }
    
    extra_events2 = c()
    if (length(SI_links_to_break) > 0){
      for (i in 1: length(SI_links_to_break)) {
        # i = 1
        sampling_intl = SI_links_to_break[[i]]
        linked_pers = strsplit(names(SI_links_to_break[i]), split = "-") %>% unlist() %>% as.integer()
        sick_per = linked_pers[1]
        sus_per = linked_pers[2]
        
        
        # Remove the observation intervals later than the sick person recovered
        sampling_intl = matrix(sampling_intl[rowSums(sampling_intl >= recov_times[sick_per]) < 2, ], ncol = 2)
        # Change the lower bound of the interval where the sick person becomes sick to be the sick time
        sampling_intl[rowSums(sampling_intl >= recov_times[sick_per]) == 1, 2] = recov_times[sick_per]
        if (nrow(sampling_intl) > 0){
          for (j in 1: nrow(sampling_intl)) {
            # j = 1
            
            if (runif(1) < 1 - exp(-(sampling_intl[j, 2] - sampling_intl[j, 1])*params.cur[7])){
              new_disconn = data.frame(time = sam_texp(1, params.cur[7], sampling_intl[j, 1], sampling_intl[j, 2]),
                                       event = 7,
                                       per1 = sort(linked_pers)[1],
                                       per2 = sort(linked_pers)[2])
              extra_events2 = rbind.data.frame(extra_events2, new_disconn)
              
              if (!sampling_intl[j, 2] %in% obs_time){ # If the sick person recovers before t1
                # Prepare the sampling bounds determined by the disconnection and t1
                t1 = obs_time[which(obs_time > sampling_intl[j, 2])[1]]
                temp_bounds = c(new_disconn$time, sampling_intl[j, 2], t1) 
                
                # Prepare the parameters used for each interval
                parms_temp = params.cur[c(4, 3)]
                probs = (1 - exp(-parms_temp*diff(temp_bounds)))*exp(-parms_temp*temp_bounds[-length(temp_bounds)])
                probs = probs/sum(probs)
                
                # Sample the interval the connection will take place
                interval_idx = sample.int(length(temp_bounds)-1, 1, prob = probs)
                
                # Generate the connection event
                extra_events2 = rbind.data.frame(extra_events2,
                                                 data.frame(time = sam_texp(1, parms_temp[interval_idx], temp_bounds[interval_idx], temp_bounds[interval_idx+1]),
                                                            event = c(4, 3)[interval_idx],
                                                            per1 = sort(linked_pers)[1],
                                                            per2 = sort(linked_pers)[2]))
              } else{ # If the sick person does not recover before t1 
                extra_events2 = rbind.data.frame(extra_events2,
                                                 data.frame(time = sam_texp(1, params.cur[4], new_disconn$time, sampling_intl[j, 2]),
                                                            event = 4,
                                                            per1 = sort(linked_pers)[1],
                                                            per2 = sort(linked_pers)[2]))
              }
              
            }
          }
        }
      }
    }
    
    
    # Create a matrix to store the bounds for each network times
    event_bounds = cbind.data.frame(lb = net_events_to_sample_temp$obs_window-window_length, 
                                    recov1 = recov_times[net_events_to_sample_temp$per1], 
                                    recov2 = recov_times[net_events_to_sample_temp$per2],
                                    ub = net_events_to_sample_temp$obs_window)
    rownames(event_bounds) = rownames(net_events_to_sample_temp)
    
    # Set the lower bounds of the network events to be the infection times of the first patient involved 
    # if the infection time is in the obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per1] & 
                       infect_times[net_events_to_sample_temp$per1] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per1[temp_idx]] 
    
    # Set the lower bounds of the network events to be the infection times of the second patient involved
    # if the second patient is infected later than the first patient in the same obs window
    temp_idx = which(event_bounds$lb <= infect_times[net_events_to_sample_temp$per2] & 
                       infect_times[net_events_to_sample_temp$per2] <= event_bounds$ub)
    event_bounds$lb[temp_idx] = infect_times[net_events_to_sample_temp$per2[temp_idx]] 
    
    # Row sort the recovery times
    event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov1", "recov2")] = 
      event_bounds[which(event_bounds$recov1 >= event_bounds$recov2), c("recov2", "recov1")]
    
    # Remove the recovery times if the recovery times are not in the corresponding obs window
    event_bounds$recov1[event_bounds$recov1 <= event_bounds$lb | event_bounds$recov1 >= event_bounds$ub] = NA
    event_bounds$recov2[event_bounds$recov2 <= event_bounds$lb | event_bounds$recov2 >= event_bounds$ub] = NA
    
    infect_times[infect_times < 0] = max(c(events$time, recover.dat$time))+1
    
    # Find the type of links at the lower bounds
    initial_link_type = rowSums(cbind(infect_times[net_events_to_sample_temp$per1] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per1], 
                                      infect_times[net_events_to_sample_temp$per2] <= event_bounds$lb & event_bounds$lb < recov_times[net_events_to_sample_temp$per2])) + 
      ifelse(net_events_to_sample_temp$event %in% 3: 5, 3, 6)
    
    # Find the number of NAs in each row. If there are two NAs, 
    # then the health status of the two persons are not changed within the interval. 
    # If there is one NA, the health status of one of the person changes within the interval.
    # If there are no NAs, the health status of both person changes within the interval.
    n_texp_dist = rowSums(is.na(event_bounds))
    
    idx1 = which(n_texp_dist == 2) # indices of the events with only one interval
    idx2 = which(n_texp_dist == 1) # indices of the events with two intervals
    idx3 = which(n_texp_dist == 0) # indices of the events with three intervals
    
    
    event_bounds1 = cbind(event_bounds$lb[idx1], event_bounds$ub[idx1])
    sampled_time1 = sam_texp(length(idx1),
                             params.cur[initial_link_type[idx1]],
                             a = event_bounds1[, 1],
                             b = event_bounds1[, 2])
    sampled_time1[initial_link_type[idx1] == 7] = (sampled_time1[initial_link_type[idx1] == 7] - event_bounds1[initial_link_type[idx1] == 7, 1]) + event_bounds1[initial_link_type[idx1] == 7, 1]
    net_events_to_sample_temp$time[idx1] = sampled_time1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 2] = big.sums[1]
    
    
    middle_bound = event_bounds[idx2, c("recov1", "recov2")]
    middle_bound[is.na(middle_bound)] = 0
    middle_bound = rowSums(middle_bound)
    event_bounds2 = cbind(event_bounds$lb[idx2], middle_bound, event_bounds$ub[idx2])
    net_events_to_sample_temp$time[idx2] = sam_texp2(n = length(idx2),
                                                     rate = cbind(params.cur[initial_link_type[idx2]],
                                                                  params.cur[initial_link_type[idx2]-1]),
                                                     a = event_bounds2[, 1],
                                                     b = event_bounds2[, 2],
                                                     c = event_bounds2[, 3])
    # Find the indices of the events happening after the middle bound (i.e. the connection happen after one of the patients is recovered). 
    # The event type is "minus 1" compared to the events happening before the middle bound
    idx2_A = idx2[which(middle_bound < net_events_to_sample_temp$time[idx2] & net_events_to_sample_temp$time[idx2] <= event_bounds$ub[idx2])]
    net_events_to_sample_temp$event[idx2_A] = net_events_to_sample_temp$event[idx2_A] - 1
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 3] = big.sums[1]
    
    net_events_to_sample_temp$time[idx3] = sam_texp3(n = length(idx3), 
                                                     rate = cbind(params.cur[initial_link_type[idx3]], 
                                                                  params.cur[initial_link_type[idx3]-1], 
                                                                  params.cur[initial_link_type[idx3]-2]), 
                                                     a = event_bounds$lb[idx3], 
                                                     b = event_bounds$recov1[idx3], 
                                                     c = event_bounds$recov2[idx3], 
                                                     d = event_bounds$ub[idx3])
    idx3_A = idx3[which(event_bounds$recov1[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$recov2[idx3])]
    net_events_to_sample_temp$event[idx3_A] = net_events_to_sample_temp$event[idx3_A] - 1
    idx3_B = idx3[which(event_bounds$recov2[idx3] < net_events_to_sample_temp$time[idx3] & net_events_to_sample_temp$time[idx3] <= event_bounds$ub[idx3])]
    net_events_to_sample_temp$event[idx3_B] = net_events_to_sample_temp$event[idx3_B] - 2
    
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # 
    # all_big_sums[it, 4] = big.sums[1]
    # 
    # events_temp = events
    # events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    # events_temp = rbind.data.frame(events_temp, extra_events)
    # events_temp = events_temp[order(events_temp$time), ]
    # events_temp = events_temp[events_temp$time %% 1 != 0, ]
    # 
    # PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    # event.counts = PA$event.counts
    # big.sums = PA$big.sums
    # big.sums
    # 
    # 
    # all_big_sums[it, 4] = big.sums[1]
    
    events_temp = events
    events_temp[rownames(net_events_to_sample_temp), "time"] = net_events_to_sample_temp$time
    events_temp = rbind.data.frame(events_temp, extra_events, extra_events2)
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[events_temp$time %% 1 != 0, ]
    
    PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    # all_big_sums[it, 5] = big.sums[1]
    # print(all_big_sums[it, ])
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  # complete_events = rbind.data.frame(full_events, 
  #                                    cbind.data.frame(time = recover.dat$time, 
  #                                                     event = 2, 
  #                                                     per1 = recover.dat$per1, 
  #                                                     per2 = NA))
  # complete_events = complete_events[order(complete_events$time), ]
  # return(list(params, all_big_sums))
  return(params)
}

# 02/06/23 Added by Houjie to impute the network as much as it allows
infer_miss_recov8_test = function(dat, priors, obs_time,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  impute = "filter", model="SIR", 
                                  timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  reports = dat$report; report.times = dat$report.times
  window_length = obs_time[2] - obs_time[1]
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, report = reports, times = report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  net_events = events[events$event %in% 3: 8, ]
  net_events[net_events$per2 < net_events$per1, c("per1", "per2")] = 
    net_events[net_events$per2 < net_events$per1, c("per2", "per1")]
  
  infect_times = rep(max(events$time)+1, length = nrow(G0))
  infect_times[I0] = 0
  infect_times[events$per1[events$event == 1]] = events$time[events$event == 1]
  
  # record the discrete events whose time of appearance will be sampled
  net_events_to_sample = net_events[net_events$time %% 1 == 0, ] 
  coarsened_events_to_sample = apply(net_events_to_sample, 1, function(event_i){
    # event_i = net_events_to_sample[1, ] %>% as.numeric()
    t1 = event_i[1]
    t0 = obs_time[which(obs_time == t1)-1]
    link = event_i[3: 4] %>% as.integer()
    # if the exact time of connection is observed
    t0 = max(net_events[net_events$per1 == link[1] & 
                          net_events$per2 == link[2] & 
                          net_events$time %% 1 != 0 &
                          net_events$time < t1 & 
                          net_events$time > t0, "time"], t0)
    c(link, t0, t1, ifelse(event_i[2] %in% 3: 5, 0, 1), 1-ifelse(event_i[2] %in% 3: 5, 0, 1))
  }) %>% t() %>% as.data.frame()
  colnames(coarsened_events_to_sample) = c("p1", "p2", "t0", "t1", "net0", "net1")
  
  # coarsened_events_to_sampleA does not contain the events where in (t0, t1) is before or after both people are sick
  # coarsened_events_to_sampleA = coarsened_events_to_sample[(coarsened_events_to_sample$i1 > coarsened_events_to_sample$t1 &
  #                                                          coarsened_events_to_sample$i2 > coarsened_events_to_sample$t1) |
  #                                                          (coarsened_events_to_sample$i1 < coarsened_events_to_sample$t0 &
  #                                                          coarsened_events_to_sample$i2 < coarsened_events_to_sample$t0), ]
  # coarsened_events_to_sampleB = coarsened_events_to_sample[!((coarsened_events_to_sample$i1 > coarsened_events_to_sample$t1 &
  #                                                             coarsened_events_to_sample$i2 > coarsened_events_to_sample$t1) |
  #                                                            (coarsened_events_to_sample$i1 < coarsened_events_to_sample$t0 &
  #                                                               coarsened_events_to_sample$i2 < coarsened_events_to_sample$t0)), ]
  
  net_snapshots = net_clipper(net_events, G0, obs_time)
  multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
  # sick_pers = c(1: nrow(G0))[which(infect_times > 0)]
  
  SI_links_to_break1 = do.call("rbind", lapply(1: length(multiple_nei), function(i){
    cbind.data.frame(p1 = as.integer(names(multiple_nei[i])), p2 = multiple_nei[[i]])
  }))
  SI_links_to_break1 = apply(SI_links_to_break1, 1, function(sick_link){
    # sick_link = SI_links_to_break1[56, ] %>% as.integer()
    sick_per = sick_link[1] # the susceptible person
    sick_nei = sick_link[2] # the infectious patient
    i1 = infect_times[sick_per] 
    i2 = infect_times[sick_nei] 
    # We only care about the interval where the susceptible person is infected
    t0 = obs_time[which(obs_time > i1)[1] - 1]
    t1 = obs_time[which(obs_time > i1)[1]]
    related_net_events = net_events[net_events$per1 == min(sick_link) & net_events$per2 == max(sick_link), ]
    if (sum(related_net_events$time %% 1 != 0) == 1 & nrow(related_net_events) != 0){
      # Find the latest time we observe this link being connected
      t0 = ifelse(net_snapshots[sick_link[1], sick_link[2], ][which(obs_time > i1)[1] - 1],
                  t0, related_net_events$time[related_net_events$time %% 1 != 0])
    }
    c(sick_per, sick_nei, t0, t1, 1, net_snapshots[sick_link[1], sick_link[2], which(obs_time > i1)[1]])
  })  %>% t()  %>% as.data.frame()
  colnames(SI_links_to_break1) = c("p1", "p2", "t0", "t1", "net0", "net1")
  SI_links_to_break1[SI_links_to_break1$p2 < SI_links_to_break1$p1, 1: 2] = SI_links_to_break1[SI_links_to_break1$p2 < SI_links_to_break1$p1, 2: 1]
  
  G = as.matrix(G0)
  event_temp = events
  event_temp = rbind.data.frame(data.frame(time = 0, event = 1, per1 = I0, per2 = NA), event_temp)
  SI_links_to_break2 = c()
  for (i in 1: nrow(event_temp)) {
    # i = 1
    events_i = event_temp[i, ]
    sick_link = as.integer(events_i[3: 4])
    if (events_i$event == 1){
      sick_per = events_i$per1
      i1 = infect_times[sick_per]
      t0 = obs_time[which(obs_time > i1)[1] - 1]
      t1 = obs_time[which(obs_time > i1)[1]]
      victims = which(G[events_i$per1, ] == 1)
      victims = victims[infect_times[victims] >= t1]
      if (sum(infect_times[victims] >= t1) > 0){
        for (j in 1: length(victims)) {
          # j = 1
          sus_per = victims[j]
          # Find the link status at each observation time point
          link_status = net_snapshots[sick_per, sus_per, ]
          
          # Find the observation intervals where the link remains connected
          sampling_intl = matrix(cbind(obs_time[-length(obs_time)], obs_time[-1])[link_status[-length(link_status)] == 1 & link_status[-1] == 1, ], ncol = 2)
          
          # Remove the observation intervals earlier than the sick person getting sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl <= events_i$time) < 2, ], ncol = 2)
          # Change the lower bound of the interval where the sick person becomes sick to be the sick time
          sampling_intl[rowSums(sampling_intl <= events_i$time) == 1, 1] = events_i$time
          # Remove the intervals when and after the susceptible person becomes sick
          sampling_intl = matrix(sampling_intl[rowSums(sampling_intl < infect_times[sus_per]) == 2, ], ncol = 2)
          
          if (nrow(sampling_intl) >= 1){
            # sampling_intl = list(sampling_intl)
            # strsplit(names(break_points), split = "-") %>% unlist() %>% as.integer()
            
            SI_links_to_break2 = rbind.data.frame(SI_links_to_break2,
                                                  cbind.data.frame(matrix(rep(c(sick_per, sus_per), each = nrow(sampling_intl)), ncol = 2), sampling_intl, 1, 1))
          }
        }
      }
    }
    if (events_i$event %in% 3: 5){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 1
    }
    if (events_i$event %in% 6: 8){
      G[sick_link[3], sick_link[4]] <- G[sick_link[4], sick_link[3]] <- 0
    }
  }
  colnames(SI_links_to_break2) = c("p1", "p2", "t0", "t1", "net0", "net1")
  SI_links_to_break2[SI_links_to_break2$p1 > SI_links_to_break2$p2, 1: 2] = 
    SI_links_to_break2[SI_links_to_break2$p1 > SI_links_to_break2$p2, 2: 1]
  SI_links_to_break2 = unique(SI_links_to_break2)
  
  # SI_links_to_break2 = cbind.data.frame(SI_links_to_break2,
  #                                      "i1" = infect_times[SI_links_to_break2$p1], 
  #                                      "i2" = infect_times[SI_links_to_break2$p2])
  
  # infect_times[infect_times == max(events$time)+1] = -1
  
  multi_events_to_sample = rbind.data.frame(SI_links_to_break1, SI_links_to_break2, coarsened_events_to_sample)
  
  multi_events_to_sample = multi_events_to_sample[!duplicated(multi_events_to_sample), ]
  multi_events_to_sample$i1 = infect_times[multi_events_to_sample$p1]
  multi_events_to_sample$i2 = infect_times[multi_events_to_sample$p2]
  # single_events_to_sample = coarsened_events_to_sampleA
  # run iterations
  S = samples * thin + burn
  for(it in 1:S){
    # it = 1
    net_events_to_sample_temp = net_events_to_sample
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    
    times = NULL
    for(ix in 1:length(MR$recover)){
      # ix = 2
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      
      imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      
      times = c(times, imputed)
    }
    
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    recov_times = rep(max(c(events$time, recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    
    extra_events1 = do.call("rbind", mclapply(1: nrow(multi_events_to_sample), function(i){
      # i = 19
      temp_edge = multi_events_to_sample[i, ] %>% as.numeric()
      
      interval_imputer(sick_link = temp_edge[1: 2],
                       sampling_intvl = temp_edge[3: 4],
                       epi_times = c(temp_edge[7], recov_times[temp_edge[1]], temp_edge[8], recov_times[temp_edge[2]]),
                       net_status = temp_edge[5: 6],
                       params.cur)
    }, mc.cores = detectCores()))
    # extra_events1 = c()
    # for(i in 1: nrow(multi_events_to_sample)){
    # # i = 19
    # temp_edge = multi_events_to_sample[i, ] %>% as.numeric()
    # 
    # extra_events1 = rbind.data.frame(extra_events1, 
    #                                  interval_imputer(sick_link = temp_edge[1: 2], 
    #                                                                  sampling_intvl = temp_edge[3: 4], 
    #                                                                  epi_times = c(temp_edge[7], recov_times[temp_edge[1]], temp_edge[8], recov_times[temp_edge[2]]), 
    #                                                                  net_status = temp_edge[5: 6],
    #                                                                  params.cur))
    # }
    
    # extra_events2 = do.call("rbind", mclapply(1: nrow(single_events_to_sample), function(i){
    #   # set.seed(101)
    #   temp_edge = single_events_to_sample[i, ] %>% as.numeric()
    #   
    #   interval_imputer(sick_link = temp_edge[1: 2], 
    #                    sampling_intvl = temp_edge[3: 4], 
    #                    epi_times = c(temp_edge[7], recov_times[temp_edge[1]], temp_edge[8], recov_times[temp_edge[2]]), 
    #                    net_status = temp_edge[5: 6])
    # }, mc.cores = 12))
    
    events_temp = rbind.data.frame(events[(events$time %% 1) != 0, ],
                                   extra_events1)
    events_temp = events_temp[order(events_temp$time), ]
    
    # G_temp = G0
    # for (i in 1: nrow(events_temp)) {
    #   event_i = events_temp[i, ]
    #   if (event_i$event %in% 3: 5){
    #     G_temp[event_i$per1, event_i$per2] <- G_temp[event_i$per2, event_i$per1] <- 1 + G_temp[event_i$per2, event_i$per1]
    #   }
    #   if (event_i$event %in% 6: 8){
    #     G_temp[event_i$per1, event_i$per2] <- G_temp[event_i$per2, event_i$per1] <- -1 + G_temp[event_i$per2, event_i$per1]
    #   }
    # }
    # which(as.matrix(G_temp)  == -1, arr.ind = TRUE)
    
    PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    # all_big_sums[it, 5] = big.sums[1]
    # print(all_big_sums[it, ])
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  # complete_events = rbind.data.frame(full_events, 
  #                                    cbind.data.frame(time = recover.dat$time, 
  #                                                     event = 2, 
  #                                                     per1 = recover.dat$per1, 
  #                                                     per2 = NA))
  # complete_events = complete_events[order(complete_events$time), ]
  # return(list(params, all_big_sums))
  return(params)
}
infer_miss_recov9_test = function(dat, priors, obs_time,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  impute = "filter", model="SIR", 
                                  timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  # dat = survey_dat
  # preparations
  G0 = dat$G0; I0 = dat$I0; 
  # events = dat$events
  reports = dat$health_report; report.times = as.integer(rownames(reports))
  window_length = obs_time[2] - obs_time[1]
  sick_events = dat$sick_events
  net_snapshots = dat$net_snapshots
  interact_report = dat$interaction_report
  
  infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
  infect_times[sick_events$per1] = sick_events$time
  
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  # if(plot){ par(mfrow=c(2,2)) }
  
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
  for (i in 1: nrow(sick_events)) {
    # i = 2
    per = sick_events$per1[i]
    
    ub = report.times[which(reports[, per] == -1)[1]]
    lb = report.times[which(reports[, per] == -1)[1]-1]
    intervals[i, ] = c(lb, ub)
  }
  
  rownames(intervals) = sick_events$per1
  colnames(intervals) = c("lb", "ub")
  # MR = get_miss_recov(reports, report.times, dats_no_recov_dis_net$events)
  # recov.persons = unlist(MR$recover)
  # intervals = MR$intervals
  
  
  
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = list()
  for (i in 2: nrow(sick_events)) {
    # i = 5
    sick_per = sick_events$per1[i]
    sick_time =  sick_events$time[i]
    time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
    
    # Infectious people at time of sick_per's infection
    other_patients = sick_events$per1[sick_events$time < sick_time]
    other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
    
    # list of persons reporting their interaction with sick_per
    interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
    
    # Identify the set of people who might be the virus transmitter
    # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
    # and those who ever connected with sick_per and infectious
    nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
                             interactd_per[interactd_per %in% other_patients]) 
  }
  names(nei_infec) = sick_events$per1[-1]
  
  multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
  SI_links_to_break1 = do.call("rbind", lapply(1: length(multiple_nei), function(i){
    cbind.data.frame(p1 = as.integer(names(multiple_nei[i])), p2 = multiple_nei[[i]])
  }))
  SI_links_to_break1 = apply(SI_links_to_break1, 1, function(sick_link){
    # sick_link = SI_links_to_break1[56, ] %>% as.integer()
    sick_per = sick_link[1] # the susceptible person
    sick_nei = sick_link[2] # the infectious patient
    i1 = infect_times[sick_per]
    i2 = infect_times[sick_nei]
    # We only care about the interval where the susceptible person is infected
    t0 = obs_time[which(obs_time > i1)[1] - 1]
    t1 = obs_time[which(obs_time > i1)[1]]
    
    interval_idx = cut(i1, obs_time, labels = 1: (length(obs_time)-1)) %>% as.integer()
    c(sick_per, sick_nei, t0, t1, net_snapshots[sick_link[1], sick_link[2], c(interval_idx, interval_idx+1)])
  })  %>% t()  %>% as.data.frame()
  colnames(SI_links_to_break1) = c("per1", "per2", "t0", "t1", "net0", "net1")
  
  edges_to_sample = do.call("rbind", lapply(1: (length(obs_time)-1), function(i){
    cbind(interact_report[[i]], t0 = i, t1 = i+1)
  }))
  edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
                                                 net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
  edges_to_sample$t0 = obs_time[edges_to_sample$t0]
  edges_to_sample$t1 = obs_time[edges_to_sample$t1]
  
  # edges_to_sample = rbind.data.frame(edges_to_sample, SI_links_to_break1)
  edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
  # nei_infec2 = get_nei_infection(dats$G0, dats_no_recov_dis_net$events, report = reports, times = report.times)
  # options(warn = 0) 
  # res_list = list()
  # for (i in 1: length(nei_infec)) {
  #   # i = 24
  #   tryCatch({
  #     res_list[[i]] = sort(nei_infec[[i]]) == sort(nei_infec2[[i]])
  #   }, error = function(e){print(i)})
  # }
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # Identify network events and intervals to sample
  # run iterations
  S = samples * thin + burn
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # Select the infectious neighbours
    
    nei_infec_temp = lapply(nei_infec, function(x){
      if (length(x) > 1){
        # sample(x, length(x), replace = FALSE)[1: sample(1: length(x), 1)]
        sample(x, length(x), replace = FALSE)[1: sample(1: 1, 1)]
      } else{
        x
      }
    })
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    
    recover.dat = c()
    for(ix in 1: (length(report.times)-1)){
      # ix = 1
      
      recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
      lb = report.times[ix]; ub = report.times[ix+1]
      
      imputed = propose_recov_filter2(lb, ub, 
                                      recovers = recovs, 
                                      events.infec = sick_events[-1, ][lb <= sick_events[-1, ]$time & sick_events[-1, ]$time <= ub, c("time","per1"), ], 
                                      nei_infec_temp, 
                                      gam = gam.cur)
      
      # times = c(times, imputed)
      recover.dat = rbind.data.frame(recover.dat, cbind.data.frame(imputed, recovs))
    }
    colnames(recover.dat) = c("time", "per1")
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    edges_to_sample_tmp = edges_to_sample
    extra_edges = c()
    for (i in 1: length(nei_infec_temp)){
      # i = 1
      sick_per = as.integer(names(nei_infec_temp)[i])
      interval_idx = cut(infect_times[sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
      
      SI_links = apply(cbind(as.integer(names(nei_infec_temp)[i]),nei_infec_temp[[i]]), 1, sort) %>% t()
      SI_links = cbind.data.frame(SI_links, t0 = obs_time[interval_idx], t1 = obs_time[interval_idx+1])
      colnames(SI_links)[1: 2] = c("per1", "per2")
      
      idx_in_report = which(rbind(SI_links, edges_to_sample_tmp[, 1: 4]) %>% duplicated()) - nrow(SI_links)
      if (length(idx_in_report) > 0){
        SI_links = SI_links[which(rbind(SI_links, edges_to_sample_tmp[, 1: 4]) %>% duplicated(fromLast = TRUE)), ]
        edges_to_sample_tmp = edges_to_sample_tmp[-idx_in_report, ]
        
        extra_edges = rbind.data.frame(extra_edges, 
                                       rbind(cbind(SI_links[, 1: 3], t1 = infect_times[sick_per]+1e-5, 
                                                   net0 = net_snapshots[as.matrix(cbind(SI_links[, 1: 2], interval_idx))], net1 = 1),
                                             cbind(SI_links[, 1: 2], t0 = infect_times[sick_per]+1e-5, t1 = obs_time[interval_idx+1],
                                                   net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links[, 1: 2], interval_idx+1))]))
        )
      }
    }
    edges_to_sample_tmp = rbind.data.frame(edges_to_sample_tmp, extra_edges)
    
    
    imputed_net_events = do.call("rbind", mclapply(1: nrow(edges_to_sample_tmp), function(i){
      # i = 1186
      temp_edge = edges_to_sample_tmp[i, ] %>% as.numeric()
      
      interval_imputer(sick_link = temp_edge[1: 2],
                       sampling_intvl = temp_edge[3: 4],
                       epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
                       net_status = temp_edge[5: 6],
                       params.cur)
    }, mc.cores = detectCores()))
    
    
    events_temp = rbind.data.frame(sick_events, 
                                   imputed_net_events)
    events_temp = events_temp[order(events_temp$time), ]
    
    
    
    PA = parse_augment2_test(G0, I0, events_temp, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    print(big.sums)
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  return(params)
}





# 2. the entire process function

# 2.1 a diagnostic function for inference results 
# (a matrix with each column being a chain for a particular parameter)
infer.diag <- function(res, method, plot=F){
  vars = colnames(res)
  
  ess = coda::effectiveSize(res)#/nrow(res)
  zscore = coda::geweke.diag(res)$z
  pvals = sapply(zscore, function(z) ifelse(z<0, pnorm(z), 1-pnorm(z))) * 2
  if(plot){
    par(mfrow=c(2,2))
    res = coda::as.mcmc(res)
    coda::geweke.plot(res, auto.layout = F)
  }
  res.dat = as.data.frame(rbind(ess, zscore, pvals))
  names(res.dat) = vars
  res.dat$method = as.character(method)
  row.names(res.dat) = c("ESS","z.score","p.values")
  return(res.dat)
}

# 06/03/2019
# modified: able to set pdfnames
# (but all plots are saved to a pdf file)
pipeline_miss_recov <- function(datname, fpath = "~/Documents/",
                                interval = 7, miss_prop = 1, miss_model = "SIR",
                                doMH = T, save_miss = T, 
                                pdfname = NULL, ...){
  dats = readRDS(paste0(fpath,datname,".rds"))
  miss_dats = miss_recovery(dats, interval, miss_prop, miss_model)
  if(save_miss){
    misspath = paste0(fpath,datname,"_miss",miss_prop*100,".rds")
    saveRDS(miss_dats, misspath)
  }
  
  diag.dat = NULL
  
  par(mfrow=c(2,2))
  if(is.null(pdfname)){
    pdf(paste0(fpath,datname,".pdf"), width = 9, height = 6)
  }else{
    pdf(paste0(fpath,pdfname,".pdf"), width = 9, height = 6)
  }
  
  res.fil = infer_miss_recov(dats = miss_dats, impute = "filter", ...)
  diag.dat = rbind(diag.dat, infer.diag(res.fil, "filter"))
  res.rej = infer_miss_recov(dats = miss_dats, impute = "rej", ...)
  diag.dat = rbind(diag.dat, infer.diag(res.rej, "reject"))
  if(doMH){
    res.MH = infer_miss_recov(dats = miss_dats, impute = "MH", ...)
    diag.dat = rbind(diag.dat, infer.diag(res.MH, "MH"))
  }
  
  dev.off()
  
  return(diag.dat)
}
infer_miss_recov9_test2 = function(dat, priors, obs_time,
                                   init.params = NULL,
                                   verbose = T, plot = T, output.sams = 100, 
                                   samples = 1000, burn = 100, thin = 1,
                                   impute = "filter", model="SIR", 
                                   timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  # dat = survey_dat
  # preparations
  G0 = dat$G0; I0 = dat$I0; 
  # events = dat$events
  reports = dat$health_report; report.times = as.integer(rownames(reports))
  window_length = obs_time[2] - obs_time[1]
  sick_events = dat$sick_events
  net_snapshots = dat$net_snapshots
  interact_report = dat$interaction_report
  
  infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
  infect_times[sick_events$per1] = sick_events$time
  
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  # if(plot){ par(mfrow=c(2,2)) }
  
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
  for (i in 1: nrow(sick_events)) {
    # i = 2
    per = sick_events$per1[i]
    
    ub = report.times[which(reports[, per] == -1)[1]]
    lb = report.times[which(reports[, per] == -1)[1]-1]
    intervals[i, ] = c(lb, ub)
  }
  
  rownames(intervals) = sick_events$per1
  colnames(intervals) = c("lb", "ub")
  # MR = get_miss_recov(reports, report.times, dats_no_recov_dis_net$events)
  # recov.persons = unlist(MR$recover)
  # intervals = MR$intervals
  
  
  
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = list()
  for (i in 2: nrow(sick_events)) {
    # i = 2
    sick_per = sick_events$per1[i]
    sick_time =  sick_events$time[i]
    time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
    
    # Infectious people at time of sick_per's infection
    other_patients = sick_events$per1[sick_events$time < sick_time]
    other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
    
    # list of persons reporting their interaction with sick_per
    interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
    
    # Identify the set of people who might be the virus transmitter
    # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
    # and those who ever connected with sick_per and infectious
    nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
                             interactd_per[interactd_per %in% other_patients]) 
  }
  names(nei_infec) = sick_events$per1[-1]
  
  multiple_nei = nei_infec[sapply(nei_infec, function(x){length(x) > 1})]
  SI_links_to_break1 = do.call("rbind", lapply(1: length(multiple_nei), function(i){
    cbind.data.frame(p1 = as.integer(names(multiple_nei[i])), p2 = multiple_nei[[i]])
  }))
  SI_links_to_break1 = apply(SI_links_to_break1, 1, function(sick_link){
    # sick_link = SI_links_to_break1[56, ] %>% as.integer()
    sick_per = sick_link[1] # the susceptible person
    sick_nei = sick_link[2] # the infectious patient
    i1 = infect_times[sick_per]
    i2 = infect_times[sick_nei]
    # We only care about the interval where the susceptible person is infected
    t0 = obs_time[which(obs_time > i1)[1] - 1]
    t1 = obs_time[which(obs_time > i1)[1]]
    
    interval_idx = cut(i1, obs_time, labels = 1: (length(obs_time)-1)) %>% as.integer()
    c(sick_per, sick_nei, t0, t1, net_snapshots[sick_link[1], sick_link[2], c(interval_idx, interval_idx+1)])
  })  %>% t()  %>% as.data.frame()
  colnames(SI_links_to_break1) = c("per1", "per2", "t0", "t1", "net0", "net1")
  
  edges_to_sample = do.call("rbind", lapply(1: (length(obs_time)-1), function(i){
    cbind(interact_report[[i]], t0 = i, t1 = i+1)
  }))
  edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
                                                 net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
  edges_to_sample$t0 = obs_time[edges_to_sample$t0]
  edges_to_sample$t1 = obs_time[edges_to_sample$t1]
  
  # edges_to_sample = rbind.data.frame(edges_to_sample, SI_links_to_break1)
  edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
  # nei_infec2 = get_nei_infection(dats$G0, dats_no_recov_dis_net$events, report = reports, times = report.times)
  # options(warn = 0) 
  # res_list = list()
  # for (i in 1: length(nei_infec)) {
  #   # i = 24
  #   tryCatch({
  #     res_list[[i]] = sort(nei_infec[[i]]) == sort(nei_infec2[[i]])
  #   }, error = function(e){print(i)})
  # }
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # Identify network events and intervals to sample
  # run iterations
  S = samples * thin + burn
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # Select the infectious neighbours
    
    nei_infec_temp = lapply(nei_infec, function(x){
      # x = nei_infec[[40]]
      if (length(x) > 1){
        # sample(x, sample(1: length(x), 1), replace = FALSE)
        sample(x, length(x), replace = FALSE)[1]
      } else{
        x
      }
    })
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    
    recover.dat = c()
    for(ix in 1: (length(report.times)-1)){
      # ix = 1
      
      recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
      lb = report.times[ix]; ub = report.times[ix+1]
      
      imputed = propose_recov_filter2(lb, ub, 
                                      recovers = recovs, 
                                      events.infec = sick_events[-1, ][lb <= sick_events[-1, ]$time & sick_events[-1, ]$time <= ub, c("time","per1"), ], 
                                      nei_infec_temp, 
                                      gam = gam.cur)
      
      # times = c(times, imputed)
      recover.dat = rbind.data.frame(recover.dat, cbind.data.frame(imputed, recovs))
    }
    colnames(recover.dat) = c("time", "per1")
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    edges_to_sample_tmp = edges_to_sample
    extra_edges = c()
    for (i in 1: length(nei_infec_temp)){
      # i = 1
      sick_per = as.integer(names(nei_infec_temp)[i])
      interval_idx = cut(infect_times[sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
      
      SI_links = apply(cbind(as.integer(names(nei_infec_temp)[i]),nei_infec_temp[[i]]), 1, sort) %>% t()
      SI_links = cbind.data.frame(SI_links, t0 = obs_time[interval_idx], t1 = obs_time[interval_idx+1])
      colnames(SI_links)[1: 2] = c("per1", "per2")
      
      idx_in_report = which(rbind(SI_links, edges_to_sample_tmp[, 1: 4]) %>% duplicated()) - nrow(SI_links)
      if (length(idx_in_report) > 0){
        SI_links = SI_links[which(rbind(SI_links, edges_to_sample_tmp[, 1: 4]) %>% duplicated(fromLast = TRUE)), ]
        edges_to_sample_tmp = edges_to_sample_tmp[-idx_in_report, ]
        
        extra_edges = rbind.data.frame(extra_edges, 
                                       rbind(cbind(SI_links[, 1: 3], t1 = infect_times[sick_per]+1e-5, 
                                                   net0 = net_snapshots[as.matrix(cbind(SI_links[, 1: 2], interval_idx))], net1 = 1),
                                             cbind(SI_links[, 1: 2], t0 = infect_times[sick_per]+1e-5, t1 = obs_time[interval_idx+1],
                                                   net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links[, 1: 2], interval_idx+1))]))
        )
      }
    }
    edges_to_sample_tmp = rbind.data.frame(edges_to_sample_tmp, extra_edges)
    
    # edges_to_sample_tmp = rbind.data.frame(edges_to_sample_tmp)
    
    imputed_net_events = do.call("rbind", mclapply(1: nrow(edges_to_sample_tmp), function(i){
      # i = 1186
      temp_edge = edges_to_sample_tmp[i, ] %>% as.numeric()
      
      interval_imputer(sick_link = temp_edge[1: 2],
                       sampling_intvl = temp_edge[3: 4],
                       epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
                       net_status = temp_edge[5: 6],
                       params.cur)
    }, mc.cores = detectCores()))
    
    events_temp = rbind.data.frame(sick_events, 
                                   imputed_net_events)
    events_temp = events_temp[order(events_temp$time), ]
    
    
    PA = parse_augment2_test(G0, I0, events_temp[-1, ], recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    # big.sums
    # table(dats$events$event)
    # dats$infer2[nrow(dats$infer2), 10: 17] / dats$infer2[nrow(dats$infer2), 2: 9]
    # big.sums / (dats$infer2[nrow(dats$infer2), 10: 17] / dats$infer2[nrow(dats$infer2), 2: 9])
    # print(big.sums)
    
    # event_coarsed = dats_no_recov_dis_net$events
    # PA_coarsed = parse_augment2_test(G0, I0, event_coarsed, recover.dat)
    # PA_coarsed$event.counts
    # PA_coarsed$big.sums
    
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  return(params)
}
infer_miss_recov10_test2 = function(dat, priors, obs_time,
                                    init.params = NULL,
                                    verbose = T, plot = T, output.sams = 100, 
                                    samples = 1000, burn = 100, thin = 1,
                                    timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; 
  # events = dat$events
  reports = dat$health_report; report.times = as.integer(rownames(reports))
  window_length = obs_time[2] - obs_time[1]
  sick_events = dat$sick_events
  net_snapshots = dat$net_snapshots
  interact_report = dat$interaction_report
  
  infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
  infect_times[sick_events$per1] = sick_events$time
  
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
  for (i in 1: nrow(sick_events)) {
    # i = 2
    per = sick_events$per1[i]
    
    ub = report.times[which(reports[, per] == -1)[1]]
    lb = report.times[which(reports[, per] == -1)[1]-1]
    intervals[i, ] = c(lb, ub)
  }
  rownames(intervals) = sick_events$per1
  colnames(intervals) = c("lb", "ub")
  
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = list()
  for (i in 2: nrow(sick_events)) {
    # i = 2
    sick_per = sick_events$per1[i]
    sick_time =  sick_events$time[i]
    time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
    
    # Infectious people at time of sick_per's infection
    other_patients = sick_events$per1[sick_events$time < sick_time]
    other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
    
    # list of persons reporting their interaction with sick_per
    interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
    
    # Identify the set of people who might be the virus transmitter
    # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
    # and those who ever connected with sick_per and infectious
    nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
                             interactd_per[interactd_per %in% other_patients]) 
  }
  names(nei_infec) = sick_events$per1[-1]
  
  
  edges_to_sample = do.call("rbind", lapply(1: (length(obs_time)-1), function(i){
    cbind(interact_report[[i]], t0 = i, t1 = i+1)
  }))
  edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
                                                 net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
  edges_to_sample$t0 = obs_time[edges_to_sample$t0]
  edges_to_sample$t1 = obs_time[edges_to_sample$t1]
  
  edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
  
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # Identify network events and intervals to sample
  # run iterations
  S = samples * thin + burn
  set.seed(1)
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # Select the infectious neighbours
    nei_infec_temp = lapply(nei_infec, function(x){
      # x = nei_infec[[40]]
      if (length(x) > 1){
        # sample(x, sample(1: length(x), 1), replace = FALSE)
        sample(x, length(x), replace = FALSE)[1]
      } else{
        x
      }
    })
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    recover.dat = c()
    for(ix in 1: (length(report.times)-1)){
      # ix = 1
      
      recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
      lb = report.times[ix]; ub = report.times[ix+1]
      
      imputed = propose_recov_filter2(lb, ub, 
                                      recovers = recovs, 
                                      events.infec = sick_events[-1, ][lb <= sick_events[-1, ]$time & sick_events[-1, ]$time <= ub, c("time","per1"), ], 
                                      nei_infec_temp, 
                                      gam = gam.cur)
      
      # times = c(times, imputed)
      recover.dat = rbind.data.frame(recover.dat, cbind.data.frame(imputed, recovs))
    }
    colnames(recover.dat) = c("time", "per1")
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    
    
    edges_to_sample_tmp2 = edges_to_sample
    
    nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
                                       sick_nei = unlist(nei_infec_temp))
    
    interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
    # SI_links2 = cbind.data.frame(rowSort(as.matrix(nei_infec_temp2)), t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
    SI_links2 = cbind.data.frame(t(apply(as.matrix(nei_infec_temp2), 1, sort)), t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
    rownames(SI_links2) = 1: nrow(SI_links2)
    colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
    temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
    SI_links2 = SI_links2[temp_idx, ]
    interval_idx2 = interval_idx2[temp_idx]
    edges_to_sample_tmp2 = edges_to_sample_tmp2[-(which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2)), ]
    
    extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
    extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                         net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
    extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                         t1 = obs_time[interval_idx2+1],
                                                         net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
    extra_edges2 = as.data.frame(extra_edges2)
    colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
    edges_to_sample_tmp2 = rbind.data.frame(edges_to_sample_tmp2, extra_edges2)
    
    
    # edges_to_sample_tmp$i1 = infect_times[edges_to_sample_tmp$per1]
    # edges_to_sample_tmp$i2 = infect_times[edges_to_sample_tmp$per2]
    # edges_to_sample_tmp$r1 = recov_times[edges_to_sample_tmp$per1]
    # edges_to_sample_tmp$r2 = recov_times[edges_to_sample_tmp$per2]
    # 
    # edges_to_sample_tmp$i1[edges_to_sample_tmp$i1 < edges_to_sample_tmp$t0 | 
    #                          edges_to_sample_tmp$i1 > edges_to_sample_tmp$t1] = NA
    # edges_to_sample_tmp$i2[edges_to_sample_tmp$i2 < edges_to_sample_tmp$t0 | 
    #                          edges_to_sample_tmp$i2 > edges_to_sample_tmp$t1] = NA
    # edges_to_sample_tmp$r1[edges_to_sample_tmp$r1 < edges_to_sample_tmp$t0 | 
    #                          edges_to_sample_tmp$r1 > edges_to_sample_tmp$t1] = NA
    # edges_to_sample_tmp$r2[edges_to_sample_tmp$r2 < edges_to_sample_tmp$t0 | 
    #                          edges_to_sample_tmp$r2 > edges_to_sample_tmp$t1] = NA
    
    imputed_net_events = c()
    
    imputed_net_events = do.call("rbind", mclapply(1: nrow(edges_to_sample_tmp2), function(i){
      temp_edge = edges_to_sample_tmp2[i, ] %>% as.numeric()
      interval_imputer(sick_link = temp_edge[1: 2],
                       sampling_intvl = temp_edge[3: 4],
                       epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
                       net_status = temp_edge[5: 6],
                       params.cur)
    }, mc.cores = detectCores()))
    events_temp = rbind.data.frame(sick_events, 
                                   imputed_net_events)
    events_temp = events_temp[order(events_temp$time), ]
    
    
    PA = parse_augment2_test(G0, I0, events_temp[-1, ], recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    if(sum(big.sums<0) > 0) break
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  return(params)
}

infer_miss_recov11_test = function(dat, priors, obs_time,
                                   init.params = NULL,
                                   verbose = T, plot = T, output.sams = 100, 
                                   samples = 1000, burn = 100, thin = 1,
                                   timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; 
  # events = dat$events
  reports = dat$health_report; report.times = as.integer(rownames(reports))
  window_length = obs_time[2] - obs_time[1]
  sick_events = dat$sick_events
  net_snapshots = dat$net_snapshots
  interact_report = dat$interaction_report
  
  infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
  infect_times[sick_events$per1] = sick_events$time
  
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
  for (i in 1: nrow(sick_events)) {
    # i = 2
    per = sick_events$per1[i]
    
    ub = report.times[which(reports[, per] == -1)[1]]
    lb = report.times[which(reports[, per] == -1)[1]-1]
    intervals[i, ] = c(lb, ub)
  }
  rownames(intervals) = sick_events$per1
  colnames(intervals) = c("lb", "ub")
  
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = list()
  for (i in 2: nrow(sick_events)) {
    # i = 2
    sick_per = sick_events$per1[i]
    sick_time =  sick_events$time[i]
    time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
    
    # Infectious people at time of sick_per's infection
    other_patients = sick_events$per1[sick_events$time < sick_time]
    other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
    
    # list of persons reporting their interaction with sick_per
    interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
    
    # Identify the set of people who might be the virus transmitter
    # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
    # and those who ever connected with sick_per and infectious
    nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
                             interactd_per[interactd_per %in% other_patients]) 
  }
  names(nei_infec) = sick_events$per1[-1]
  
  
  edges_to_sample = do.call("rbind", lapply(1: (length(obs_time)-1), function(i){
    cbind(interact_report[[i]], t0 = i, t1 = i+1)
  }))
  edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
                                                 net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
  edges_to_sample$t0 = obs_time[edges_to_sample$t0]
  edges_to_sample$t1 = obs_time[edges_to_sample$t1]
  
  edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
  
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # Identify network events and intervals to sample
  # run iterations
  S = samples * thin + burn
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # Select the infectious neighbours
    nei_infec_temp = lapply(nei_infec, function(x){
      # x = nei_infec[[40]]
      if (length(x) > 1){
        # sample(x, sample(1: length(x), 1), replace = FALSE)
        sample(x, length(x), replace = FALSE)[1]
      } else{
        x
      }
    })
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    recover.dat = c()
    for(ix in 1: (length(report.times)-1)){
      # ix = 1
      
      recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
      lb = report.times[ix]; ub = report.times[ix+1]
      
      imputed = propose_recov_filter2(lb, ub, 
                                      recovers = recovs, 
                                      events.infec = sick_events[-1, ][lb <= sick_events[-1, ]$time & sick_events[-1, ]$time <= ub, c("time","per1"), ], 
                                      nei_infec_temp, 
                                      gam = gam.cur)
      
      # times = c(times, imputed)
      recover.dat = rbind.data.frame(recover.dat, cbind.data.frame(imputed, recovs))
    }
    colnames(recover.dat) = c("time", "per1")
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    
    # Prepare the edges and their time intervals we would like to sample
    edges_to_sample_tmp2 = edges_to_sample
    
    nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
                                       sick_nei = unlist(nei_infec_temp))
    
    interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
    # SI_links2 = cbind.data.frame(rowSort(as.matrix(nei_infec_temp2)), t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
    temp_mat = t(apply(as.matrix(nei_infec_temp2), 1, sort)); rownames(temp_mat) = NULL
    SI_links2= cbind.data.frame(temp_mat, t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
    colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
    temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
    SI_links2 = SI_links2[temp_idx, ]
    interval_idx2 = interval_idx2[temp_idx]
    edges_to_sample_tmp2 = edges_to_sample_tmp2[-(which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2)), ]
    
    extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
    extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                         net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
    extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                         t1 = obs_time[interval_idx2+1],
                                                         net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
    extra_edges2 = as.data.frame(extra_edges2)
    colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
    edges_to_sample_tmp2 = rbind.data.frame(edges_to_sample_tmp2, extra_edges2)
    
    
    edges_to_sample_tmp2$i1 = infect_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$i2 = infect_times[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    
    edges_to_sample_tmp2$i1[edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$i2[edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$r1[edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$r2[edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1] = NA
    homo_time_edges = which(rowSums(is.na(edges_to_sample_tmp2[, 7:10])) == 4)
    edges_to_sample_tmp2_A = edges_to_sample_tmp2[homo_time_edges, ]
    edges_to_sample_tmp2_B = edges_to_sample_tmp2[-homo_time_edges, ]
    
    # We sample the edges where both persons don't change their health status (time-homogeneous sampling)
    # This is easier to vectorize so we do it explicitly
    imputed_net_events_A = c()
    # We first sample some events to make sure the start and end of each edge is the same
    health_status = (infect_times[edges_to_sample_tmp2_A$per1] < edges_to_sample_tmp2_A$t0 & 
                       edges_to_sample_tmp2_A$t0 < recov_times[edges_to_sample_tmp2_A$per1]) + 
      (infect_times[edges_to_sample_tmp2_A$per2] < edges_to_sample_tmp2_A$t0 & 
         edges_to_sample_tmp2_A$t0 < recov_times[edges_to_sample_tmp2_A$per2])
    # Find the events there the start and the end are not the same
    tmp_idx = which(edges_to_sample_tmp2_A$net0 != edges_to_sample_tmp2_A$net1) 
    first_type = ifelse(edges_to_sample_tmp2_A$net0[tmp_idx] == 0, 3, 6) + health_status[tmp_idx]
    imputed_net_events_A = rbind.data.frame(imputed_net_events_A,
                                            data.frame(time = sam_texp(length(tmp_idx), 
                                                                       rate = params.cur[first_type], 
                                                                       a = edges_to_sample_tmp2_A$t0[tmp_idx],
                                                                       b = edges_to_sample_tmp2_A$t1[tmp_idx]),
                                                       event = first_type,
                                                       per1 = edges_to_sample_tmp2_A$per1[tmp_idx],
                                                       per2 = edges_to_sample_tmp2_A$per2[tmp_idx])
    )
    
    edges_to_sample_tmp2_A$t0[tmp_idx] = imputed_net_events_A$time
    edges_to_sample_tmp2_A$net0[tmp_idx] = edges_to_sample_tmp2_A$net1[tmp_idx]
    while (nrow(edges_to_sample_tmp2_A) > 0) {
      next_event = ifelse(as.logical(edges_to_sample_tmp2_A$net0), 6, 3) + health_status
      proposed_time = rexp(nrow(edges_to_sample_tmp2_A), rate = params.cur[next_event]) + edges_to_sample_tmp2_A$t0
      accept_idx = which(proposed_time <= edges_to_sample_tmp2_A$t1)
      
      imputed_net_events_A = rbind.data.frame(imputed_net_events_A,
                                              data.frame(time = proposed_time[accept_idx],
                                                         event = next_event[accept_idx],
                                                         per1 = edges_to_sample_tmp2_A$per1[accept_idx],
                                                         per2 = edges_to_sample_tmp2_A$per2[accept_idx]))
      next_event = ifelse(next_event[accept_idx] %in% 3: 5, 3, -3) + next_event[accept_idx]
      
      proposed_time2 = sam_texp(length(accept_idx), 
                                rate = params.cur[next_event], 
                                a = proposed_time[accept_idx],
                                b = edges_to_sample_tmp2_A$t1[accept_idx])
      imputed_net_events_A = rbind.data.frame(imputed_net_events_A,
                                              data.frame(time = sam_texp(length(accept_idx), 
                                                                         rate = params.cur[next_event], 
                                                                         a = proposed_time2,
                                                                         b = edges_to_sample_tmp2_A$t1[accept_idx]),
                                                         event = next_event,
                                                         per1 = edges_to_sample_tmp2_A$per1[accept_idx],
                                                         per2 = edges_to_sample_tmp2_A$per2[accept_idx]))
      edges_to_sample_tmp2_A$t0[accept_idx] = proposed_time2
      edges_to_sample_tmp2_A = edges_to_sample_tmp2_A[accept_idx, ]
      health_status = health_status[accept_idx]
    }
    # Now we sample the edges and time intervals where the health status of both persons change
    imputed_net_events_B = do.call("rbind", mclapply(1: nrow(edges_to_sample_tmp2_B), function(i){
      temp_edge = edges_to_sample_tmp2_B[i, ] %>% as.numeric()
      interval_imputer(sick_link = temp_edge[1: 2],
                       sampling_intvl = temp_edge[3: 4],
                       epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
                       net_status = temp_edge[5: 6],
                       params.cur)
    }, mc.cores = 1))
    
    
    events_temp = rbind.data.frame(sick_events, 
                                   imputed_net_events_A,
                                   imputed_net_events_B)
    events_temp = events_temp[order(events_temp$time), ]
    
    
    PA = parse_augment2_test(G0, I0, events_temp[-1, ], recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    if(sum(big.sums<0) > 0) break
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  return(params)
}
infer_miss_recov12_cpp = function(dat, priors, obs_time,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42){
  
  # dat = survey_dat; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T;
  # samples = nsample; burn = nburn; thin = 1; seed = 1; init.params = NULL; timing = T
  
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; 
  # events = dat$events
  reports = dat$health_report; report.times = as.integer(rownames(reports))
  window_length = obs_time[2] - obs_time[1]
  sick_events = dat$sick_events
  net_snapshots = dat$net_snapshots
  interact_report = dat$interaction_report
  
  infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
  infect_times[sick_events$per1] = sick_events$time
  
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
  for (i in 1: nrow(sick_events)) {
    # i = 2
    per = sick_events$per1[i]
    
    ub = report.times[which(reports[, per] == -1)[1]]
    lb = report.times[which(reports[, per] == -1)[1]-1]
    intervals[i, ] = c(lb, ub)
  }
  rownames(intervals) = sick_events$per1
  colnames(intervals) = c("lb", "ub")
  
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = list()
  for (i in 2: nrow(sick_events)) {
    # i = 2
    sick_per = sick_events$per1[i]
    sick_time =  sick_events$time[i]
    time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
    
    # Infectious people at time of sick_per's infection
    other_patients = sick_events$per1[sick_events$time < sick_time]
    other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
    
    # list of persons reporting their interaction with sick_per
    interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
    
    # Identify the set of people who might be the virus transmitter
    # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
    # and those who ever connected with sick_per and infectious
    nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
                             interactd_per[interactd_per %in% other_patients]) 
  }
  names(nei_infec) = sick_events$per1[-1]
  
  
  edges_to_sample = do.call("rbind", lapply(1: (length(obs_time)-1), function(i){
    cbind.data.frame(interact_report[[i]], t0 = i, t1 = i+1)
  }))
  edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
                                                 net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
  edges_to_sample$t0 = obs_time[edges_to_sample$t0]
  edges_to_sample$t1 = obs_time[edges_to_sample$t1]
  
  edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
  
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # Identify network events and intervals to sample
  # run iterations
  S = samples * thin + burn
  imp_event_log = matrix(NA, nrow = S, ncol = 8)
  for(it in 1:S){
    # it = 34
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # Select the infectious neighbours
    nei_infec_temp = lapply(nei_infec, function(x){
      # x = nei_infec[[40]]
      if (length(x) > 1){
        # sample(x, sample(1: length(x), 1), replace = FALSE)
        sample(x, length(x), replace = FALSE)[1]
      } else{
        x
      }
    })
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    recover.dat = c()
    for(ix in 1: (length(report.times)-1)){
      # ix = 1
      
      recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
      lb = report.times[ix]; ub = report.times[ix+1]
      
      imputed = propose_recov_filter2(lb, ub, 
                                      recovers = recovs, 
                                      events.infec = sick_events[-1, ][lb <= sick_events[-1, ]$time & sick_events[-1, ]$time <= ub, c("time","per1"), ], 
                                      nei_infec_temp, 
                                      gam = gam.cur)
      
      # times = c(times, imputed)
      recover.dat = rbind.data.frame(recover.dat, cbind.data.frame(imputed, recovs))
    }
    colnames(recover.dat) = c("time", "per1")
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    
    # Prepare the edges and their time intervals we would like to sample
    edges_to_sample_tmp2 = edges_to_sample
    
    nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
                                       sick_nei = unlist(nei_infec_temp))
    
    interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
    # SI_links2 = cbind.data.frame(rowSort(as.matrix(nei_infec_temp2)), t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
    temp_mat = t(apply(as.matrix(nei_infec_temp2), 1, sort)); rownames(temp_mat) = NULL
    SI_links2= cbind.data.frame(temp_mat, t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
    colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
    temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
    SI_links2 = SI_links2[temp_idx, ]
    interval_idx2 = interval_idx2[temp_idx]
    edges_to_sample_tmp2 = edges_to_sample_tmp2[-(which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2)), ]
    
    extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
    extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                         net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
    extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                         t1 = obs_time[interval_idx2+1],
                                                         net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
    extra_edges2 = as.data.frame(extra_edges2)
    colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
    edges_to_sample_tmp2 = rbind.data.frame(edges_to_sample_tmp2, extra_edges2)
    
    
    edges_to_sample_tmp2$i1 = infect_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$i2 = infect_times[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    
    edges_to_sample_tmp2$i1[edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$i2[edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$r1[edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$r2[edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1] = NA
    homo_time_edges = which(rowSums(is.na(edges_to_sample_tmp2[, 7:10])) == 4)
    edges_to_sample_tmp2_A = edges_to_sample_tmp2[homo_time_edges, ]
    edges_to_sample_tmp2_B = edges_to_sample_tmp2[-homo_time_edges, ]
    
    # We sample the edges where both persons don't change their health status (time-homogeneous sampling)
    # This is easier to vectorize so we do it explicitly
    imputed_net_events_A = c()
    # We first sample some events to make sure the start and end of each edge is the same
    health_status = (infect_times[edges_to_sample_tmp2_A$per1] < edges_to_sample_tmp2_A$t0 & 
                       edges_to_sample_tmp2_A$t0 < recov_times[edges_to_sample_tmp2_A$per1]) + 
      (infect_times[edges_to_sample_tmp2_A$per2] < edges_to_sample_tmp2_A$t0 & 
         edges_to_sample_tmp2_A$t0 < recov_times[edges_to_sample_tmp2_A$per2])
    # Find the pair of agents whose network status at the start and the end are not the same
    tmp_idx = which(edges_to_sample_tmp2_A$net0 != edges_to_sample_tmp2_A$net1) 
    first_type = ifelse(edges_to_sample_tmp2_A$net0[tmp_idx] == 0, 3, 6) + health_status[tmp_idx]
    imputed_net_events_A = rbind.data.frame(imputed_net_events_A,
                                            data.frame(time = sam_texp(length(tmp_idx), 
                                                                       rate = params.cur[first_type], 
                                                                       a = edges_to_sample_tmp2_A$t0[tmp_idx],
                                                                       b = edges_to_sample_tmp2_A$t1[tmp_idx]),
                                                       event = first_type,
                                                       per1 = edges_to_sample_tmp2_A$per1[tmp_idx],
                                                       per2 = edges_to_sample_tmp2_A$per2[tmp_idx])
    )
    
    edges_to_sample_tmp2_A$t0[tmp_idx] = imputed_net_events_A$time
    edges_to_sample_tmp2_A$net0[tmp_idx] = edges_to_sample_tmp2_A$net1[tmp_idx]
    while (nrow(edges_to_sample_tmp2_A) > 0) {
      next_event = ifelse(as.logical(edges_to_sample_tmp2_A$net0), 6, 3) + health_status
      proposed_time = rexp(nrow(edges_to_sample_tmp2_A), rate = params.cur[next_event]) + edges_to_sample_tmp2_A$t0
      accept_idx = which(proposed_time <= edges_to_sample_tmp2_A$t1)
      
      imputed_net_events_A = rbind.data.frame(imputed_net_events_A,
                                              data.frame(time = proposed_time[accept_idx],
                                                         event = next_event[accept_idx],
                                                         per1 = edges_to_sample_tmp2_A$per1[accept_idx],
                                                         per2 = edges_to_sample_tmp2_A$per2[accept_idx]))
      next_event = ifelse(next_event[accept_idx] %in% 3: 5, 3, -3) + next_event[accept_idx]
      
      proposed_time2 = sam_texp(length(accept_idx), 
                                rate = params.cur[next_event], 
                                a = proposed_time[accept_idx],
                                b = edges_to_sample_tmp2_A$t1[accept_idx])
      imputed_net_events_A = rbind.data.frame(imputed_net_events_A,
                                              data.frame(time = sam_texp(length(accept_idx), 
                                                                         rate = params.cur[next_event], 
                                                                         a = proposed_time2,
                                                                         b = edges_to_sample_tmp2_A$t1[accept_idx]),
                                                         event = next_event,
                                                         per1 = edges_to_sample_tmp2_A$per1[accept_idx],
                                                         per2 = edges_to_sample_tmp2_A$per2[accept_idx]))
      edges_to_sample_tmp2_A$t0[accept_idx] = proposed_time2
      edges_to_sample_tmp2_A = edges_to_sample_tmp2_A[accept_idx, ]
      health_status = health_status[accept_idx]
    }
    # Now we sample the edges and time intervals where the health status of both persons change
    imputed_net_events_B = do.call("rbind", mclapply(1: nrow(edges_to_sample_tmp2_B), function(i){
      temp_edge = edges_to_sample_tmp2_B[i, ] %>% as.numeric()
      interval_imputer(sick_link = temp_edge[1: 2],
                       sampling_intvl = temp_edge[3: 4],
                       epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
                       net_status = temp_edge[5: 6],
                       params.cur)
    }, mc.cores = 1))
    
    
    events_temp = rbind.data.frame(sick_events, 
                                   imputed_net_events_A,
                                   imputed_net_events_B)
    events_temp = events_temp[order(events_temp$time), ]
    
    
    PA = parse_augment_cpp(G0, I0, events_temp[-1, ], recover.dat)
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    imp_event_log[it, ] = event.counts
    big.sums = PA$big.sums
    
    if(sum(big.sums<0) > 0) break
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  return(params)
}
infer_miss_recov13_cpp = function(dat, priors, obs_time, force_one = TRUE,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42){
  
  # dat = survey_dat; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T;
  # samples = nsample; burn = nburn; thin = 1; seed = 1; init.params = NULL; timing = T; force_one = TRUE
  
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; 
  # events = dat$events
  reports = dat$health_report; report.times = as.integer(rownames(reports))
  window_length = obs_time[2] - obs_time[1]
  sick_events = dat$sick_events
  net_snapshots = dat$net_snapshots
  interact_report = dat$interaction_report
  
  infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
  infect_times[sick_events$per1] = sick_events$time
  
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
  for (i in 1: nrow(sick_events)) {
    # i = 2
    per = sick_events$per1[i]
    
    ub = report.times[which(reports[, per] == -1)[1]]
    lb = report.times[which(reports[, per] == -1)[1]-1]
    intervals[i, ] = c(lb, ub)
  }
  rownames(intervals) = sick_events$per1
  colnames(intervals) = c("lb", "ub")
  
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = list()
  for (i in 2: nrow(sick_events)) {
    # i = 2
    sick_per = sick_events$per1[i]
    sick_time =  sick_events$time[i]
    time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
    
    # Infectious people at time of sick_per's infection
    other_patients = sick_events$per1[sick_events$time < sick_time]
    other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
    
    # list of persons reporting their interaction with sick_per
    interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
    
    # Identify the set of people who might be the virus transmitter
    # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
    # and those who ever connected with sick_per and infectious
    nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
                             interactd_per[interactd_per %in% other_patients]) 
  }
  names(nei_infec) = sick_events$per1[-1]
  
  
  edges_to_sample = do.call("rbind", lapply(1: (length(obs_time)-1), function(i){
    cbind.data.frame(interact_report[[i]], t0 = i, t1 = i+1)
  }))
  edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
                                                 net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
  edges_to_sample$t0 = obs_time[edges_to_sample$t0]
  edges_to_sample$t1 = obs_time[edges_to_sample$t1]
  
  edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
  
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # Identify network events and intervals to sample
  # run iterations
  S = samples * thin + burn
  imp_event_log = array(NA, dim = c(2, 8, nrow = S))
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # Select the infectious neighbours
    nei_infec_temp = lapply(nei_infec, function(x){
      # x = nei_infec[[40]]
      if (length(x) > 1){
        # sample(x, sample(1: length(x), 1), replace = FALSE)
        sample(x, length(x), replace = FALSE)[1]
      } else{
        x
      }
    })
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    recover.dat = c()
    for(ix in 1: (length(report.times)-1)){
      # ix = 1
      
      recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
      lb = report.times[ix]; ub = report.times[ix+1]
      
      imputed = propose_recov_filter2(lb, ub, 
                                      recovers = recovs, 
                                      events.infec = sick_events[-1, ][lb <= sick_events[-1, ]$time & sick_events[-1, ]$time <= ub, c("time","per1"), ], 
                                      nei_infec_temp, 
                                      gam = gam.cur)
      
      # times = c(times, imputed)
      recover.dat = rbind.data.frame(recover.dat, cbind.data.frame(imputed, recovs))
    }
    colnames(recover.dat) = c("time", "per1")
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    
    # Prepare the edges and their time intervals we would like to sample
    # Pick patient-healthy per pair and note an connection when the healthy person gets sick
    {
      edges_to_sample_tmp2 = edges_to_sample
      
      nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
                                         sick_nei = unlist(nei_infec_temp))
      interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
      temp_mat = t(apply(as.matrix(nei_infec_temp2), 1, sort)); rownames(temp_mat) = NULL
      SI_links2= cbind.data.frame(temp_mat, t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
      colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
      temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
      SI_links2 = SI_links2[temp_idx, ]
      interval_idx2 = interval_idx2[temp_idx]
      edges_to_sample_tmp2 = edges_to_sample_tmp2[-(which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2)), ]
      # cut the interval at the time the healthy person gets sick 
      extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
      extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                           net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
      extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                           t1 = obs_time[interval_idx2+1],
                                                           net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
      extra_edges2 = as.data.frame(extra_edges2)
      colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
      edges_to_sample_tmp2 = rbind.data.frame(edges_to_sample_tmp2, extra_edges2)
    }
    
    # Summarize the agent pair, interval, end point network status, 
    edges_to_sample_tmp2$i1 = infect_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$i2 = infect_times[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$i1[edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$i2[edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$r1[edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$r2[edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1] = NA
    
    
    homo_time_edges = which(rowSums(is.na(edges_to_sample_tmp2[, 7:10])) == 4)
    edges_to_sample_tmp2_A = edges_to_sample_tmp2[homo_time_edges, ]
    edges_to_sample_tmp2_B = edges_to_sample_tmp2[-homo_time_edges, ]
    
    if (force_one){ # We force a pair events happening at the edge disconnected at both endpoints
      temp_idx = which(edges_to_sample_tmp2_A$net0 == edges_to_sample_tmp2_A$net1 & edges_to_sample_tmp2_A$net0 == 0)
      edges_to_sample_tmp2_B = rbind.data.frame(edges_to_sample_tmp2_B, edges_to_sample_tmp2_A[temp_idx, ])
      edges_to_sample_tmp2_A = edges_to_sample_tmp2_A[-temp_idx, ]
    }
    # We sample the edges where both persons don't change their health status (time-homogeneous sampling)
    # This is easier to vectorize so we do it explicitly
    imputed_net_events_A = c()
    # We first sample some events to make sure the start and end of each edge is the same
    health_status = (infect_times[edges_to_sample_tmp2_A$per1] < edges_to_sample_tmp2_A$t0 & 
                       edges_to_sample_tmp2_A$t0 < recov_times[edges_to_sample_tmp2_A$per1]) + 
      (infect_times[edges_to_sample_tmp2_A$per2] < edges_to_sample_tmp2_A$t0 & 
         edges_to_sample_tmp2_A$t0 < recov_times[edges_to_sample_tmp2_A$per2])
    
    # Find the pair of agents whose network status at the start and the end are not the same
    tmp_idx = which(edges_to_sample_tmp2_A$net0 != edges_to_sample_tmp2_A$net1) 
    first_type = ifelse(edges_to_sample_tmp2_A$net0[tmp_idx] == 0, 3, 6) + health_status[tmp_idx]
    imputed_net_events_A = rbind.data.frame(imputed_net_events_A,
                                            data.frame(time = sam_texp(length(tmp_idx), 
                                                                       rate = params.cur[first_type], 
                                                                       a = edges_to_sample_tmp2_A$t0[tmp_idx],
                                                                       b = edges_to_sample_tmp2_A$t1[tmp_idx]),
                                                       event = first_type,
                                                       per1 = edges_to_sample_tmp2_A$per1[tmp_idx],
                                                       per2 = edges_to_sample_tmp2_A$per2[tmp_idx])
    )
    
    edges_to_sample_tmp2_A$t0[tmp_idx] = imputed_net_events_A$time
    edges_to_sample_tmp2_A$net0[tmp_idx] = edges_to_sample_tmp2_A$net1[tmp_idx]
    
    while (nrow(edges_to_sample_tmp2_A) > 0) {
      next_event = ifelse(as.logical(edges_to_sample_tmp2_A$net0), 6, 3) + health_status
      proposed_time = rexp(nrow(edges_to_sample_tmp2_A), rate = params.cur[next_event]) + edges_to_sample_tmp2_A$t0
      accept_idx = which(proposed_time <= edges_to_sample_tmp2_A$t1)
      
      imputed_net_events_A = rbind.data.frame(imputed_net_events_A,
                                              data.frame(time = proposed_time[accept_idx],
                                                         event = next_event[accept_idx],
                                                         per1 = edges_to_sample_tmp2_A$per1[accept_idx],
                                                         per2 = edges_to_sample_tmp2_A$per2[accept_idx]))
      next_event = ifelse(next_event[accept_idx] %in% 3: 5, 3, -3) + next_event[accept_idx]
      
      proposed_time2 = sam_texp(length(accept_idx), 
                                rate = params.cur[next_event], 
                                a = proposed_time[accept_idx],
                                b = edges_to_sample_tmp2_A$t1[accept_idx])
      imputed_net_events_A = rbind.data.frame(imputed_net_events_A,
                                              data.frame(time = sam_texp(length(accept_idx), 
                                                                         rate = params.cur[next_event], 
                                                                         a = proposed_time2,
                                                                         b = edges_to_sample_tmp2_A$t1[accept_idx]),
                                                         event = next_event,
                                                         per1 = edges_to_sample_tmp2_A$per1[accept_idx],
                                                         per2 = edges_to_sample_tmp2_A$per2[accept_idx]))
      edges_to_sample_tmp2_A$t0[accept_idx] = proposed_time2
      edges_to_sample_tmp2_A = edges_to_sample_tmp2_A[accept_idx, ]
      health_status = health_status[accept_idx]
    }
    # Now we sample the edges and time intervals where the health status of both persons change
    imputed_net_events_B = do.call("rbind", mclapply(1: nrow(edges_to_sample_tmp2_B), function(i){
      # i = nrow(edges_to_sample_tmp2_B)
      temp_edge = edges_to_sample_tmp2_B[i, ] %>% as.numeric()
      interval_imputer(sick_link = temp_edge[1: 2],
                       sampling_intvl = temp_edge[3: 4],
                       epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
                       net_status = temp_edge[5: 6],
                       params.cur, force_one = F)
    }, mc.cores = 1))
    
    
    events_temp = rbind.data.frame(sick_events, 
                                   imputed_net_events_A,
                                   imputed_net_events_B)
    events_temp = events_temp[order(events_temp$time), ]
    
    
    PA = parse_augment_cpp(G0, I0, events_temp[-1, ], recover.dat)
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums<0) > 0) break
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  return(list(params, sample_hyperparm = imp_event_log))
}
# infer_miss_recov14_cpp = function(dat, priors, obs_time, force_one = TRUE,
#                                   init.params = NULL,
#                                   verbose = T, plot = T, output.sams = 100, 
#                                   samples = 1000, burn = 100, thin = 1,
#                                   timing = T, seed=42){
#   
#   # dat = survey_dat2; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T;
#   # samples = nsample; burn = nburn; thin = 1; seed = 1; init.params = NULL; timing = T; force_one = TRUE
#   
#   if(timing){ time.st = Sys.time()}
#   
#   set.seed(seed)
#   
#   # preparations
#   G0 = dat$G0; I0 = dat$I0; 
#   # events = dat$events
#   reports = dat$health_report; report.times = as.integer(rownames(reports))
#   window_length = obs_time[2] - obs_time[1]
#   sick_events = dat$sick_events
#   net_snapshots = dat$net_snapshots
#   interact_report = dat$interaction_report
#   
#   infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
#   infect_times[sick_events$per1] = sick_events$time
#   
#   if("truth" %in% names(dat)){
#     true.params = dat$truth
#   }else{
#     true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
#   }
#   
#   
#   # get time intervals to operate on, and those who need exact recovery times imputed
#   intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
#   for (i in 1: nrow(sick_events)) {
#     # i = 8
#     per = sick_events$per1[i]
#     
#     ub = report.times[which(reports[, per] == -1)[1]]
#     lb = report.times[which(reports[, per] == -1)[1]-1]
#     intervals[i, ] = c(lb, ub)
#   }
#   rownames(intervals) = sick_events$per1
#   colnames(intervals) = c("lb", "ub")
#   # Add this to handle those who haven't recovered
#   intervals[which(is.na(intervals[,2])), ] = Inf
#   
#   # get neighborhood info for all infection cases at their infection times
#   nei_infec = list()
#   for (i in 2: nrow(sick_events)) {
#     # i = nrow(sick_events)
#     sick_per = sick_events$per1[i]
#     sick_time =  sick_events$time[i]
#     time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
#     
#     # Infectious people at time of sick_per's infection
#     other_patients = sick_events$per1[sick_events$time < sick_time]
#     other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
#     
#     # list of persons reporting their interaction with sick_per
#     interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
#     
#     # Identify the set of people who might be the virus transmitter
#     # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
#     # and those who ever connected with sick_per and infectious
#     # intersect(interactd_per[interactd_per != sick_per], other_patients)
#     nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
#                              interactd_per[interactd_per %in% other_patients]) 
#     # print(identical(sort(intersect(interactd_per[interactd_per != sick_per], other_patients)), sort(nei_infec[[i-1]])))
#   }
#   names(nei_infec) = sick_events$per1[-1]
#   
#   
#   edges_to_sample = do.call("rbind", lapply(1: (length(obs_time)-1), function(i){
#     cbind.data.frame(interact_report[[i]], t0 = i, t1 = i+1)
#   }))
#   edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
#                                                  net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
#   edges_to_sample$t0 = obs_time[edges_to_sample$t0]
#   edges_to_sample$t1 = obs_time[edges_to_sample$t1]
#   
#   edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
#   
#   
#   # process the priors and initialize param values
#   if(nrow(priors)==4){
#     a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
#     avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
#     b.pr = a.pr/avg.pr
#   }else{
#     a.pr = priors$count
#     b.pr = a.pr/priors$avg
#   }
#   if(is.null(init.params)){
#     # draw from priors
#     params.cur = rgamma(8, shape = a.pr, rate = b.pr)
#   }else{
#     params.cur = init.params
#   }
#   
#   # parameter values storage
#   params = matrix(ncol=8, nrow=samples)
#   vars = c("beta","gamma",
#            "alpha.SS","alpha.SI","alpha.II",
#            "omega.SS","omega.SI","omega.II")
#   colnames(params) = vars
#   
#   # Identify network events and intervals to sample
#   # run iterations
#   S = samples * thin + burn
#   imp_event_log = array(NA, dim = c(2, 8, nrow = S))
#   for(it in 1:S){
#     # it = 1
#     if(verbose){ cat("\nIteration",it,"..\n") }
#     
#     # Select the infectious neighbours
#     nei_infec_temp = lapply(nei_infec, function(x){
#       # x = nei_infec[[40]]
#       if (length(x) > 1){
#         # sample(x, sample(1: length(x), 1), replace = FALSE)
#         sample(x, length(x), replace = FALSE)[1]
#       } else{
#         x
#       }
#     })
#     
#     # (1) propose recovery times
#     gam.cur = params.cur[2]
#     recover.dat = c()
#     for(ix in 1: (length(report.times)-1)){
#       # ix = 1
#       
#       recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
#       lb = report.times[ix]; ub = report.times[ix+1]
#       
#       imputed = propose_recov_filter2(lb, ub, 
#                                       recovers = recovs, 
#                                       events.infec = sick_events[-1, ][lb <= sick_events[-1, ]$time & sick_events[-1, ]$time <= ub, c("time","per1"), ], 
#                                       nei_infec_temp, 
#                                       gam = gam.cur)
#       
#       # times = c(times, imputed)
#       recover.dat = rbind.data.frame(recover.dat, cbind.data.frame(imputed, recovs))
#     }
#     colnames(recover.dat) = c("time", "per1")
#     if(verbose){ cat("Recovery times and network events imputation done.\n") }
#     recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
#     recov_times[recover.dat$per1] = recover.dat$time
#     
#     
#     # Prepare the edges and their time intervals we would like to sample
#     # Pick patient-healthy per pair and note an connection when the healthy person gets sick
#     {
#       edges_to_sample_tmp2 = edges_to_sample
#       
#       nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
#                                          sick_nei = unlist(nei_infec_temp))
#       interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
#       temp_mat = t(apply(as.matrix(nei_infec_temp2), 1, sort)); rownames(temp_mat) = NULL
#       SI_links2= cbind.data.frame(temp_mat, t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
#       colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
#       temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
#       if (length(temp_idx) > 0){
#         SI_links2 = SI_links2[temp_idx, ]
#         interval_idx2 = interval_idx2[temp_idx]
#         edges_to_sample_tmp2 = edges_to_sample_tmp2[-(which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2)), ]
#         # cut the interval at the time the healthy person gets sick 
#         extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
#         extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
#                                                              net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
#         extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
#                                                              t1 = obs_time[interval_idx2+1],
#                                                              net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
#         extra_edges2 = as.data.frame(extra_edges2)
#         colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
#         edges_to_sample_tmp2 = rbind.data.frame(edges_to_sample_tmp2, extra_edges2)
#       }else{
#         edges_to_sample_tmp2 = edges_to_sample
#       }
#     }
#     
#     # Summarize the agent pair, interval, end point network status, 
#     edges_to_sample_tmp2$i1 = infect_times[edges_to_sample_tmp2$per1]
#     edges_to_sample_tmp2$i2 = infect_times[edges_to_sample_tmp2$per2]
#     edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
#     edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
#     edges_to_sample_tmp2$i1[edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
#                               edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1] = NA
#     edges_to_sample_tmp2$i2[edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
#                               edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1] = NA
#     edges_to_sample_tmp2$r1[edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
#                               edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1] = NA
#     edges_to_sample_tmp2$r2[edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
#                               edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1] = NA
#     
#     
#     
#     imputed_net_events_B = do.call("rbind", mclapply(1: nrow(edges_to_sample_tmp2), function(i){
#       # i = nrow(edges_to_sample_tmp2_B)
#       
#       temp_edge = edges_to_sample_tmp2[i, ] %>% as.numeric()
#       interval_imputer2(sick_link = temp_edge[1: 2],
#                         sampling_intvl = temp_edge[3: 4],
#                         epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
#                         net_status = temp_edge[5: 6],
#                         params.cur, force_one = force_one)
#       
#       
#     }, mc.cores = 1))
#     
#     
#     events_temp = rbind.data.frame(sick_events, 
#                                    imputed_net_events_B)
#     events_temp = events_temp[order(events_temp$time), ]
#     
#     
#     PA = parse_augment_cpp(G0, I0, events_temp[-1, ], recover.dat)
#     # PA
#     # 
#     # events = events_temp[-1, ]
#     # a = unique(events[events$event %in% 3: 8, c("per1", "per2")])
#     # 
#     # interact_report[[1]][!duplicated(rbind.data.frame(a, interact_report[[1]]))[-c(1: nrow(a))],]
#     # events_temp2 = rbind.data.frame(sick_events,
#     #                                imputed_net_events_A,
#     #                                imputed_net_events_B,
#     #                                cbind.data.frame(time = recover.dat$time,
#     #                                                 event = 2,
#     #                                                 per1 = recover.dat$per1,
#     #                                                 per2 = NA))
#     # 
#     # events_temp2 = events_temp2[order(events_temp2$time), ]
#     # # parse_augment4(G0, I0, events_temp2[-1, ][1:479, ])
#     # 
#     # temp = parse_augment4(G0, I0, events_temp2[-1, ])
#     # bsum_tab = big.sumer(G0, I0, events_temp2[-1, ])
#     # 
#     # colSums(bsum_tab[,1: 8] * bsum_tab$time_diff)
#     # 
#     # bsum_tab = big.sumer(G0, I0, events_temp2[-1, ][1:480, ])
#     # colSums(bsum_tab[,1: 8] * bsum_tab$time_diff)
#     # 
#     # bsum_tab = big.sumer(G0, I0, events_temp2[-1, ])
#     # colSums(bsum_tab[,1: 8] * bsum_tab$time_diff)
#     
#     event.counts = rep(0, 8)
#     names(event.counts) = 1: 8
#     event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
#     
#     big.sums = PA$big.sums
#     imp_event_log[, ,it] = rbind(event.counts, big.sums)
#     if(sum(big.sums < 0) > 0) stop(paste0("Issue at", it, "-th iter"))
#     # (3) sample params
#     a.post = a.pr + event.counts
#     b.post = b.pr + big.sums
#     params.cur = rgamma(8, shape = a.post, rate = b.post)
#     
#     if(verbose){ 
#       cat("Parameter values sampled:",params.cur,"\n")
#     }
#     
#     ## record this sample after burn-in and thinning
#     if(it > burn & (it - burn) %% thin == 0){
#       s = (it - burn)/thin
#       params[s,] = params.cur
#     }
#     
#   }
#   
#   
#   # timing it
#   if(timing){ 
#     time.en = Sys.time(); 
#     cat("\n\n")
#     print(time.en - time.st)
#     #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
#   }
#   return(list(params, sample_hyperparm = imp_event_log))
# }

infer_miss_recov14_sem = function(dat, priors, obs_time, force_one = TRUE,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42){
  
  # dat = survey_dat; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T;
  # samples = nsample; burn = nburn; thin = 1; seed = 1; init.params = NULL; timing = T; force_one = TRUE
  
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; 
  # events = dat$events
  reports = dat$health_report; report.times = as.integer(rownames(reports))
  window_length = obs_time[2] - obs_time[1]
  sick_events = dat$sick_events
  net_snapshots = dat$net_snapshots
  interact_report = dat$interaction_report
  
  infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
  infect_times[sick_events$per1] = sick_events$time
  
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
  for (i in 1: nrow(sick_events)) {
    # i = 2
    per = sick_events$per1[i]
    
    ub = report.times[which(reports[, per] == -1)[1]]
    lb = report.times[which(reports[, per] == -1)[1]-1]
    intervals[i, ] = c(lb, ub)
  }
  rownames(intervals) = sick_events$per1
  colnames(intervals) = c("lb", "ub")
  
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = list()
  for (i in 2: nrow(sick_events)) {
    # i = 10
    sick_per = sick_events$per1[i]
    sick_time =  sick_events$time[i]
    time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
    
    # Infectious people at time of sick_per's infection
    other_patients = sick_events$per1[sick_events$time < sick_time]
    other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
    
    # list of persons reporting their interaction with sick_per
    interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
    
    # Identify the set of people who might be the virus transmitter
    # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
    # and those who ever connected with sick_per and infectious
    # intersect(interactd_per[interactd_per != sick_per], other_patients)
    nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
                             interactd_per[interactd_per %in% other_patients]) 
    # print(identical(sort(intersect(interactd_per[interactd_per != sick_per], other_patients)), sort(nei_infec[[i-1]])))
  }
  names(nei_infec) = sick_events$per1[-1]
  
  
  edges_to_sample = do.call("rbind", lapply(1: (length(obs_time)-1), function(i){
    cbind.data.frame(interact_report[[i]], t0 = i, t1 = i+1)
  }))
  edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
                                                 net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
  edges_to_sample$t0 = obs_time[edges_to_sample$t0]
  edges_to_sample$t1 = obs_time[edges_to_sample$t1]
  
  edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
  
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # Identify network events and intervals to sample
  # run iterations
  S = samples * thin + burn
  imp_event_log = array(NA, dim = c(2, 8, nrow = S))
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # Select the infectious neighbours
    nei_infec_temp = lapply(nei_infec, function(x){
      # x = nei_infec[[40]]
      if (length(x) > 1){
        # sample(x, sample(1: length(x), 1), replace = FALSE)
        sample(x, length(x), replace = FALSE)[1]
      } else{
        x
      }
    })
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    recover.dat = c()
    for(ix in 1: (length(report.times)-1)){
      # ix = 1
      
      recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
      lb = report.times[ix]; ub = report.times[ix+1]
      
      imputed = propose_recov_filter2(lb, ub, 
                                      recovers = recovs, 
                                      events.infec = sick_events[-1, ][lb <= sick_events[-1, ]$time & sick_events[-1, ]$time <= ub, c("time","per1"), ], 
                                      nei_infec_temp, 
                                      gam = gam.cur)
      
      # times = c(times, imputed)
      recover.dat = rbind.data.frame(recover.dat, cbind.data.frame(imputed, recovs))
    }
    colnames(recover.dat) = c("time", "per1")
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    
    # Prepare the edges and their time intervals we would like to sample
    # Pick patient-healthy per pair and note an connection when the healthy person gets sick
    {
      edges_to_sample_tmp2 = edges_to_sample
      
      nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
                                         sick_nei = unlist(nei_infec_temp))
      interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
      temp_mat = t(apply(as.matrix(nei_infec_temp2), 1, sort)); rownames(temp_mat) = NULL
      SI_links2= cbind.data.frame(temp_mat, t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
      colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
      temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
      SI_links2 = SI_links2[temp_idx, ]
      interval_idx2 = interval_idx2[temp_idx]
      edges_to_sample_tmp2 = edges_to_sample_tmp2[-(which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2)), ]
      # cut the interval at the time the healthy person gets sick 
      extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
      extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                           net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
      extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                           t1 = obs_time[interval_idx2+1],
                                                           net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
      extra_edges2 = as.data.frame(extra_edges2)
      colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
      edges_to_sample_tmp2 = rbind.data.frame(edges_to_sample_tmp2, extra_edges2)
    }
    
    # Summarize the agent pair, interval, end point network status, 
    edges_to_sample_tmp2$i1 = infect_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$i2 = infect_times[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$i1[edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$i2[edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$r1[edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$r2[edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1] = NA
    
    
    
    imputed_net_events_B = do.call("rbind", mclapply(1: nrow(edges_to_sample_tmp2), function(i){
      # i = nrow(edges_to_sample_tmp2_B)
      
      temp_edge = edges_to_sample_tmp2[i, ] %>% as.numeric()
      interval_imputer2(sick_link = temp_edge[1: 2],
                        sampling_intvl = temp_edge[3: 4],
                        epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
                        net_status = temp_edge[5: 6],
                        params.cur, force_one = force_one)
      
      
    }, mc.cores = 1))
    
    
    events_temp = rbind.data.frame(sick_events, 
                                   imputed_net_events_B)
    events_temp = events_temp[order(events_temp$time), ]
    
    
    PA = parse_augment_cpp(G0, I0, events_temp[-1, ], recover.dat)
    # PA
    # 
    # events = events_temp[-1, ]
    # a = unique(events[events$event %in% 3: 8, c("per1", "per2")])
    # 
    # interact_report[[1]][!duplicated(rbind.data.frame(a, interact_report[[1]]))[-c(1: nrow(a))],]
    # events_temp2 = rbind.data.frame(sick_events,
    #                                imputed_net_events_A,
    #                                imputed_net_events_B,
    #                                cbind.data.frame(time = recover.dat$time,
    #                                                 event = 2,
    #                                                 per1 = recover.dat$per1,
    #                                                 per2 = NA))
    # 
    # events_temp2 = events_temp2[order(events_temp2$time), ]
    # # parse_augment4(G0, I0, events_temp2[-1, ][1:479, ])
    # 
    # temp = parse_augment4(G0, I0, events_temp2[-1, ])
    # bsum_tab = big.sumer(G0, I0, events_temp2[-1, ])
    # 
    # colSums(bsum_tab[,1: 8] * bsum_tab$time_diff)
    # 
    # bsum_tab = big.sumer(G0, I0, events_temp2[-1, ][1:480, ])
    # colSums(bsum_tab[,1: 8] * bsum_tab$time_diff)
    # 
    # bsum_tab = big.sumer(G0, I0, events_temp2[-1, ])
    # colSums(bsum_tab[,1: 8] * bsum_tab$time_diff)
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at", it, "-th iter"))
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    # params.cur = rgamma(8, shape = a.post, rate = b.post)
    params.cur = a.post / b.post
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  return(list(params, sample_hyperparm = imp_event_log))
}

infer_miss_recov14_cpp = function(dat, priors, obs_time, force_one = TRUE,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42, remove_no_neighbour = TRUE){
  
  # dat = survey_dat2; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T;
  # samples = nsample; burn = nburn; thin = 1; seed = 1; init.params = NULL; timing = T; force_one = TRUE
  
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dat$G0; I0 = dat$I0; 
  # events = dat$events
  reports = dat$health_report; report.times = as.integer(rownames(reports))
  window_length = obs_time[2] - obs_time[1]
  sick_events = dat$sick_events
  net_snapshots = dat$net_snapshots
  interact_report = dat$interaction_report
  
  infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
  infect_times[sick_events$per1] = sick_events$time
  
  if("truth" %in% names(dat)){
    true.params = dat$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
  for (i in 1: nrow(sick_events)) {
    # i = 8
    per = sick_events$per1[i]
    
    ub = report.times[which(reports[, per] == -1)[1]]
    lb = report.times[which(reports[, per] == -1)[1]-1]
    intervals[i, ] = c(lb, ub)
  }
  rownames(intervals) = sick_events$per1
  colnames(intervals) = c("lb", "ub")
  # Add this to handle those who haven't recovered
  intervals[which(is.na(intervals[,2])), ] = Inf
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = list()
  for (i in 2: nrow(sick_events)) {
    # i = nrow(sick_events)
    sick_per = sick_events$per1[i]
    sick_time =  sick_events$time[i]
    time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
    
    # Infectious people at time of sick_per's infection
    other_patients = sick_events$per1[sick_events$time < sick_time]
    other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
    
    # list of persons reporting their interaction with sick_per
    interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
    
    # Identify the set of people who might be the virus transmitter
    # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
    # and those who ever connected with sick_per and infectious
    # intersect(interactd_per[interactd_per != sick_per], other_patients)
    nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
                             interactd_per[interactd_per %in% other_patients]) 
    # print(identical(sort(intersect(interactd_per[interactd_per != sick_per], other_patients)), sort(nei_infec[[i-1]])))
  }
  names(nei_infec) = sick_events$per1[-1]
  
  if (remove_no_neighbour){
    to_remove_idx = names(nei_infec)[sapply(nei_infec, length) == 0]
    nei_infec = nei_infec[!names(nei_infec) %in% to_remove_idx]
    intervals=intervals[!rownames(intervals) %in% to_remove_idx, ]
    sick_events = sick_events[!sick_events$per1 %in% as.integer(to_remove_idx), ]
  }
  
  
  edges_to_sample = do.call("rbind", lapply(1: (length(obs_time)-1), function(i){
    cbind.data.frame(interact_report[[i]], t0 = i, t1 = i+1)
  }))
  edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
                                                 net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
  edges_to_sample$t0 = obs_time[edges_to_sample$t0]
  edges_to_sample$t1 = obs_time[edges_to_sample$t1]
  
  edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
  
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # Identify network events and intervals to sample
  # run iterations
  S = samples * thin + burn
  imp_event_log = array(NA, dim = c(2, 8, nrow = S))
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # Select the infectious neighbours
    nei_infec_temp = lapply(nei_infec, function(x){
      # x = nei_infec[[40]]
      if (length(x) > 1){
        # sample(x, sample(1: length(x), 1), replace = FALSE)
        sample(x, length(x), replace = FALSE)[1]
      } else{
        x
      }
    })
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    recover.dat = c()
    for(ix in 1: (length(report.times)-1)){
      # ix = 2
      
      recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
      lb = report.times[ix]; ub = report.times[ix+1]
      
      imputed = propose_recov_filter2(lb, ub, 
                                      recovers = recovs, 
                                      events.infec = sick_events[-1, ][lb <= sick_events[-1, ]$time & sick_events[-1, ]$time <= ub, c("time","per1"), ], 
                                      nei_infec_temp, 
                                      gam = gam.cur)
      
      # times = c(times, imputed)
      recover.dat = rbind.data.frame(recover.dat, cbind.data.frame(imputed, recovs))
    }
    colnames(recover.dat) = c("time", "per1")
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    
    # Prepare the edges and their time intervals we would like to sample
    # Pick patient-healthy per pair and note an connection when the healthy person gets sick
    {
      edges_to_sample_tmp2 = edges_to_sample
      
      nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
                                         sick_nei = unlist(nei_infec_temp))
      interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
      temp_mat = t(apply(as.matrix(nei_infec_temp2), 1, sort)); rownames(temp_mat) = NULL
      SI_links2= cbind.data.frame(temp_mat, t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
      colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
      temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
      if (length(temp_idx) > 0){
        SI_links2 = SI_links2[temp_idx, ]
        interval_idx2 = interval_idx2[temp_idx]
        edges_to_sample_tmp2 = edges_to_sample_tmp2[-(which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2)), ]
        # cut the interval at the time the healthy person gets sick 
        extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
        extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                             net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
        extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                             t1 = obs_time[interval_idx2+1],
                                                             net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
        extra_edges2 = as.data.frame(extra_edges2)
        colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
        edges_to_sample_tmp2 = rbind.data.frame(edges_to_sample_tmp2, extra_edges2)
      }else{
        edges_to_sample_tmp2 = edges_to_sample
      }
    }
    
    # Summarize the agent pair, interval, end point network status, 
    edges_to_sample_tmp2$i1 = infect_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$i2 = infect_times[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$i1[edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$i2[edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$r1[edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1] = NA
    edges_to_sample_tmp2$r2[edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                              edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1] = NA
    
    
    
    imputed_net_events_B = do.call("rbind", mclapply(1: nrow(edges_to_sample_tmp2), function(i){
      # i = nrow(edges_to_sample_tmp2_B)
      
      temp_edge = edges_to_sample_tmp2[i, ] %>% as.numeric()
      interval_imputer2(sick_link = temp_edge[1: 2],
                        sampling_intvl = temp_edge[3: 4],
                        epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
                        net_status = temp_edge[5: 6],
                        params.cur, force_one = force_one)
      
      
    }, mc.cores = 1))
    
    
    events_temp = rbind.data.frame(sick_events, 
                                   imputed_net_events_B)
    events_temp = events_temp[order(events_temp$time), ]
    
    
    PA = parse_augment_cpp(G0, I0, events_temp[-1, ], recover.dat)
    # PA
    # 
    # events = events_temp[-1, ]
    # a = unique(events[events$event %in% 3: 8, c("per1", "per2")])
    # 
    # interact_report[[1]][!duplicated(rbind.data.frame(a, interact_report[[1]]))[-c(1: nrow(a))],]
    # events_temp2 = rbind.data.frame(sick_events,
    #                                imputed_net_events_A,
    #                                imputed_net_events_B,
    #                                cbind.data.frame(time = recover.dat$time,
    #                                                 event = 2,
    #                                                 per1 = recover.dat$per1,
    #                                                 per2 = NA))
    # 
    # events_temp2 = events_temp2[order(events_temp2$time), ]
    # # parse_augment4(G0, I0, events_temp2[-1, ][1:479, ])
    # 
    # temp = parse_augment4(G0, I0, events_temp2[-1, ])
    # bsum_tab = big.sumer(G0, I0, events_temp2[-1, ])
    # 
    # colSums(bsum_tab[,1: 8] * bsum_tab$time_diff)
    # 
    # bsum_tab = big.sumer(G0, I0, events_temp2[-1, ][1:480, ])
    # colSums(bsum_tab[,1: 8] * bsum_tab$time_diff)
    # 
    # bsum_tab = big.sumer(G0, I0, events_temp2[-1, ])
    # colSums(bsum_tab[,1: 8] * bsum_tab$time_diff)
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at", it, "-th iter"))
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
    }
    
  }
  
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  if (remove_no_neighbour){
    return(list(params, sample_hyperparm = imp_event_log, removed = to_remove_idx))
  }else{
    return(list(params, sample_hyperparm = imp_event_log))
  }
}
