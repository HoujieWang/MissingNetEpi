# 05/15/2019
# utility functions for inference on incomplete data

library(Matrix)
#library(igraph)
#library(ggplot2)

source("./Truncated_Exponential.R")

# 05/20/24 Added by Houjie to drop events with lease interaction time

drop_interact = function(obj, obs_time, missing_prop = 0.1){
  # obj = survey_dat; obs_time; missing_prop = 0.2
  interaction_report_missing = obj$interaction_report
  for (t in 1: length(obj$interaction_report)) {
    # t = 1
    report_t = obj$interaction_report[[t]]
    
    net_events_t = net_events[obs_time[t] <= net_events$time & net_events$time <= obs_time[t+1], ]
    idx = which(net_events_t$per1 > net_events_t$per2)
    net_events_t[idx, 3: 4] = net_events_t[idx, 4: 3]
    time_interact = vector(mode = 'numeric', length = nrow(report_t))
    for (i in 1: nrow(report_t)) {
      # i = 10
      pair_i = as.numeric(report_t[i, ])
      
      net_status_i = obj$net_snapshots[pair_i[1], pair_i[2], t:(t+1)]
      
      tmp = net_events_t[net_events_t$per1 == pair_i[1] & net_events_t$per2 == pair_i[2], 1: 2]
      tmp$event = ifelse(tmp$event %in% 3: 5, 1, 0)
      all_status = rbind(c(obs_time[t], net_status_i[1]), # get the list of events happen at this agent pair
                         cbind(net_events_t$time[net_events_t$per1 == pair_i[1] & net_events_t$per2 == pair_i[2]], tmp$event), 
                         c(obs_time[t+1], net_status_i[2]))
      
      time_diff = all_status[-1, 1] - all_status[-nrow(all_status), 1]
      
      time_interact[i] = sum(time_diff[as.logical(all_status[-nrow(all_status), 2])])
    }
    interaction_report_missing[[t]] = interaction_report_missing[[t]][-order(time_interact)[1: ceiling(nrow(report_t) * missing_prop)], ]
  }
  
  return(interaction_report_missing)
}

# 03/05/23  Added by Houjie to adapt to the new version of survey data
propose_recov_filter2 = function(lb, ub, recovers, events.infec, nei_infec, gam=0.2){
  # recovers = recovs; events.infec = sick_events[-1, ][lb <= sick_events[-1, ]$time & sick_events[-1, ]$time <= ub, c("time","per1"), ]
  # nei_infec = nei_infec_temp; gam = 0.2
  
  
  # st = min(which(events$time > lb)); en = max(which(events$time <= ub))
  # 
  # # pull up infection log in this interval 
  # events.infec = events[st:en,]
  # events.infec = events.infec[events.infec$event == 1,c("time","per1")]
  
  # if events.infec non-empty:
  # obtain adjusted feasible sampling lower bounds for each person in `recovers`
  # if about-to-recover people are the only sick neighbors of someone infected at t, 
  # then randomly select one of them to mandately recover after t
  bounds = rep(lb, length(recovers))
  if(nrow(events.infec) > 0){
    for(r in 1:nrow(events.infec)){
      # r = 1
      p1 = events.infec$per1[r]
      if(p1 %in% recovers){
        t = events.infec$time[r]
        bounds[recovers==p1] = max(bounds[recovers==p1],t)
      }
      nei = nei_infec[[as.character(p1)]]
      poi = recovers[recovers %in% nei]
      if(length(poi)==length(nei)){
        #cat("POIs for individual",p1,":\n",poi,"\n")
        t = events.infec$time[r]
        if(length(poi)==1){
          p = poi
        }else{
          p = sample(poi, 1)
        }
        bounds[recovers==p] = max(bounds[recovers==p],t)
      }
    }
  }
  
  # sample recovery times under the adjusted bounds
  cands = sam_texp(length(recovers), gam, bounds-lb, ub-lb) + lb
  
  # cat("Interval: [",lb,",",ub,"]\n")
  # cat("To recover:",recovers,"\n")
  # cat("Feasible lower bounds:\n",bounds,"\n")
  # cat("Proposed recovery times:\n",cands,"\n")
  
  return(cands)
}
# 03/05/23  Added by Houjie to coarsen the complete data. The main difference to net_coarsener is that an accurate connection time
# before the times of infection is not recorded
net_coarsener4 = function(net_events, G0, obs_time, more_var = TRUE){
  # net_events = miss_dats1$events[miss_dats1$events$event %in% 3: 8, ]; G0 = miss_dats1$G0
  net_events[net_events$per1 >= net_events$per2, c("per1", "per2")] = 
    net_events[net_events$per1 >= net_events$per2, c("per2", "per1")]
  G0 = as.matrix(G0)
  all_net_snapshots = array(0, dim = c(nrow(G0), nrow(G0), length(obs_time)))
  all_net_snapshots[, , 1] = G0
  # network = net_events$time
  interaction_report = list()
  extra_edges = list()
  for (i in 1: (length(obs_time)-1)) {
    # i = 1
    G1 = G0
    t0 = obs_time[i]
    t1 = obs_time[(i+1)]
    temp_net_events = net_events[t0 < net_events$time & net_events$time <= t1, ] # Filter for the network events happening between t0 and t1
    
    interaction_report[[i]] = unique(temp_net_events[, c("per1", "per2")]) 
    
    for (j in 1: nrow(temp_net_events)) { # Evolve the network from t0 to t1
      # j = 123
      temp_events_net_j = temp_net_events[j, ]
      
      # Set the i,j entry in the adjacency matrix to be 1 for connection and 0 for disconnection
      G1[temp_events_net_j$per1, temp_events_net_j$per2] <- 
        G1[temp_events_net_j$per2, temp_events_net_j$per1] <- as.integer(temp_events_net_j$event %in% 3: 5)
      
    }
    if (more_var){ # if the two agents remained connected during the entire interval, then we still impute for possible net events
      edges_connected = which(G1 == G0 & G1 == 1, arr.ind = TRUE)
      edges_connected = edges_connected[edges_connected[, 1] < edges_connected[, 2], ]
      colnames(edges_connected) = c("per1", "per2")
      edges_connected = edges_connected[!duplicated(rbind(edges_connected, interaction_report[[i]]), fromLast = TRUE)[1: nrow(edges_connected)], ]
      interaction_report[[i]] = rbind.data.frame(interaction_report[[i]], edges_connected)
    }
    
    all_net_snapshots[, , (i+1)] = G1
    # Setting the network at t1 as G0, then we move to the second time window
    G0 = G1
  }
  list(net_snapshots = all_net_snapshots,
       interaction_report = interaction_report)
}


rtrunc_exp = function(intvl, param_temp){
  # param_temp = param[possbile_types]*c(1000, 1,1,1,1)
  # probs = (1 - exp(-param_temp*diff(intvl)))*exp(-param_temp*intvl[-length(intvl)])
  # probs = probs/sum(probs)
  
  lprobs = log1mexp(param_temp*diff(intvl)) - param_temp*intvl[-length(intvl)]
  lprobs = lprobs - logSumExp(lprobs)
  
  interval_idx = sample.int(length(intvl)-1, 1, prob = exp(lprobs))
  c(sam_texp(1, param_temp[interval_idx], intvl[interval_idx], intvl[interval_idx+1]), 
    interval_idx)
}
interval_imputer = function(sick_link, sampling_intvl, epi_times, net_status, params.cur, force_one = TRUE){
  #
  # sick_link = temp_edge[1: 2];
  # sampling_intvl = temp_edge[3: 4];
  # epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]);
  # net_status = temp_edge[5: 6];
  # params.cur = params.cur
  
  # when the two agents report to have net events while remained unconnected at both end points,
  # we force there to be a pair of net events
  # when they remain connected at both endpoints, we do not force the network event
  sick_per = sick_link[1]
  sick_nei = sick_link[2]
  names(epi_times) = c("i1", "r1", "i2", "r2")
  names(sampling_intvl) = c("t0", "t1")
  epi_status = c(1, -1, 1, -1)
  included_epi_idx = which(sampling_intvl[1] <= epi_times & epi_times <= sampling_intvl[2])
  sampling_info = rbind(c(sampling_intvl[1],sort(epi_times[included_epi_idx]), sampling_intvl[2]),
                        c(c(0, epi_status[included_epi_idx][order(epi_times[included_epi_idx])]), NA))
  
  intvl = sampling_info[1, ]
  intvl_idx = 1 # This variable tells which interval we are at
  added_events = c()
  if(net_status[1] != net_status[2]){
    CD_type = ifelse((net_status[2] - net_status[1]) == 1, 3, 6) +
      (epi_times["i1"] < sampling_intvl["t0"] & sampling_intvl["t0"] < epi_times["r1"]) + 
      (epi_times["i2"] < sampling_intvl["t0"] & sampling_intvl["t0"] < epi_times["r2"])
    
    if (intvl_idx == (ncol(sampling_info)-1)){
      possbile_types = CD_type
    }else{
      possbile_types = c(CD_type, cumsum(sampling_info[2, -ncol(sampling_info)][(intvl_idx+1): (ncol(sampling_info)-1)]) + CD_type) # Find the types the inverse event can take
    }
    param_temp = params.cur[possbile_types]
    sampled_event = rtrunc_exp(intvl, param_temp) %>% as.numeric()
    
    added_events = rbind.data.frame(added_events,
                                    data.frame(time = sampled_event[1],
                                               event = possbile_types[sampled_event[2]] %>% as.integer(),
                                               per1 = sick_per,
                                               per2 = sick_nei))
    intvl_idx = intvl_idx + sampled_event[2] - 1
    t_low = sampled_event[1]
    t_up = sampling_info[1, sampling_info[1, ] > t_low][1]
  } else{
    t_low = intvl[1]; t_up = intvl[2]  
  }
  if (force_one &
      net_status[1] == net_status[2]  
      & net_status[2]  == 0
      ){ # If we need at least one pair (dis)conn, we use truncated exp
    CD_type = ifelse(net_status[2] == 0, 3, 6)  + # the future event type in the current interval
      (epi_times["i1"] <= t_low & t_low < epi_times["r1"]) + 
      (epi_times["i2"] <= t_low & t_low < epi_times["r2"])
    
    # We force to have a pair of events
    if (intvl_idx == (ncol(sampling_info)-1)){ # If we are at the last subinterval, it's 
      possbile_types = CD_type
    }else{ # the future event type in the intervals after the current interval
      possbile_types = c(CD_type, cumsum(sampling_info[2, -ncol(sampling_info)][(intvl_idx+1): (ncol(sampling_info)-1)]) + CD_type) # Find the types the inverse event can take
    }
    param_temp = params.cur[possbile_types]
    (sampled_event = rtrunc_exp(intvl, param_temp)) # Here we do the truncated exponential sampling
    t_temp = sampled_event[1]
    next_event_type = possbile_types[sampled_event[2]]
    added_events = rbind.data.frame(added_events, # We add the first sampled event
                                    data.frame(time = t_temp,
                                               event = next_event_type,
                                               per1 = sick_per,
                                               per2 = sick_nei))
    intvl_idx = min(which(intvl > t_temp))-1
    
    CD_type = ifelse(CD_type %in% 6: 8, 3, 6) + 
      (epi_times["i1"] <= t_temp & t_temp < epi_times["r1"]) + 
      (epi_times["i2"] <= t_temp & t_temp < epi_times["r2"])
    # Determine if the inverse event is connection/disconnection
    intvl = c(t_temp, sampling_info[1, which(sampling_info[1, ] > t_temp)]) # Find the intervals the inverse event can happen
    if (intvl_idx == (ncol(sampling_info)-1)){
      possbile_types = CD_type
    }else{
      possbile_types = c(CD_type, cumsum(sampling_info[2, -ncol(sampling_info)][(intvl_idx+1): (ncol(sampling_info)-1)]) + CD_type) # Find the types the inverse event can take
    }
    param_temp = params.cur[possbile_types]
    (sampled_event = rtrunc_exp(intvl, param_temp)) # sample the time and type of the inverse event
    added_events = rbind.data.frame(added_events,
                                    data.frame(time = sampled_event[1],
                                               event = possbile_types[sampled_event[2]],
                                               per1 = sick_per,
                                               per2 = sick_nei))
    
    t_low = sampled_event[1]
    t_up = intvl[intvl > t_low][1]
    intvl_idx = intvl_idx + sampled_event[2] - 1
  }
  
  while (intvl_idx <= (ncol(sampling_info)-1)) {
    CD_type = ifelse(net_status[2] == 0, 3, 6)  + 
      (epi_times["i1"] <= t_low & t_low < epi_times["r1"]) + 
      (epi_times["i2"] <= t_low & t_low < epi_times["r2"])
    next_event_type = CD_type
    param_temp = params.cur[next_event_type]
    t_temp = min(rexp(1, param_temp) + t_low, t_up)
    
    # Determine if the conn/disconn happen before t_up. If happen after, 
    # it means that no event happen before t_up then we move to the second interval
    if (t_temp < t_up){
      # Add the first event
      added_events = rbind.data.frame(added_events,
                                      data.frame(time = t_temp,
                                                 event = next_event_type,
                                                 per1 = sick_per,
                                                 per2 = sick_nei))
      # Add the inverse event
      CD_type = ifelse(CD_type %in% 6: 8, 3, 6) + 
        (epi_times["i1"] <= t_temp & t_temp < epi_times["r1"]) + 
        (epi_times["i2"] <= t_temp & t_temp < epi_times["r2"])
      # Determine if the inverse event is connection/disconnection
      intvl = c(t_temp, sampling_info[1, which(sampling_info[1, ] > t_temp)]) # Find the intervals the inverse event can happen
      if (intvl_idx == (ncol(sampling_info)-1)){
        possbile_types = CD_type
      }else{
        possbile_types = c(CD_type, cumsum(sampling_info[2, -ncol(sampling_info)][(intvl_idx+1): (ncol(sampling_info)-1)]) + CD_type) # Find the types the inverse event can take
      }
      param_temp = params.cur[possbile_types]
      (sampled_event = rtrunc_exp(intvl, param_temp)) # sample the time and type of the inverse event
      added_events = rbind.data.frame(added_events,
                                      data.frame(time = sampled_event[1],
                                                 event = possbile_types[sampled_event[2]],
                                                 per1 = sick_per,
                                                 per2 = sick_nei))
      
      t_low = sampled_event[1]
      t_up = intvl[intvl > t_low][1]
      intvl_idx = intvl_idx + sampled_event[2] - 1
    }else{ # No event happen in the current interval, we go to the next one
      t_low = t_temp
      t_up = intvl[intvl > t_low][1]
      intvl_idx = intvl_idx + 1
    }
    
  }
  added_events  
}

inverse_event = function(net, intvl, param, epi_times, sick_link){
  # net = net_status
  # intvl = sampling_info[1, ]
  # param = params.cur
  
  intvl = c(intvl[1], intvl[intvl > intvl[1]]) # make sure no t0=epi_times
  
  possbile_types = ifelse(net[1] == 0, 3, 6) + 
    (epi_times["i1"] <= intvl[-length(intvl)] & intvl[-length(intvl)] < epi_times["r1"]) + 
    (epi_times["i2"] <= intvl[-length(intvl)] & intvl[-length(intvl)] < epi_times["r2"])
  
  sampled_event = rtrunc_exp(intvl, param[possbile_types]) %>% as.numeric()
  
  data.frame(time = sampled_event[1],
             event = possbile_types[sampled_event[2]] %>% as.integer(),
             per1 = sick_link[1],
             per2 = sick_link[2])
}

interval_imputer2 = function(sick_link, sampling_intvl, epi_times, net_status, params.cur, force_one = TRUE){
  # set.seed(1001)
  # sick_link = temp_edge[1: 2];
  # sampling_intvl = temp_edge[3: 4];
  # epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]);
  # net_status = temp_edge[5: 6];
  # params.cur = params.cur
  # 
  # when the two agents report to have net events while remained unconnected at both end points,
  # we force there to be a pair of net events
  # when they remain connected at both endpoints, we do not force the network event
  
  sick_per = sick_link[1]
  sick_nei = sick_link[2]
  names(epi_times) = c("i1", "r1", "i2", "r2")
  names(sampling_intvl) = c("t0", "t1")
  epi_status = c(1, -1, 1, -1)
  included_epi_idx = which(sampling_intvl[1] <= epi_times & epi_times <= sampling_intvl[2])
  sampling_info = rbind(c(sampling_intvl[1],sort(epi_times[included_epi_idx]), sampling_intvl[2]),
                        c(c(0, epi_status[included_epi_idx][order(epi_times[included_epi_idx])]), NA))
  
  intvl = sampling_info[1, ]
  intvl_idx = 1 # This variable tells which interval we are at
  
  
  added_events = c()
  net.cur = net_status[1]
  
  if(net_status[1] != net_status[2]){
    new_event = inverse_event(net.cur, intvl, params.cur, epi_times, sick_link)
    added_events = rbind.data.frame(added_events, new_event)
    
    intvl = c(new_event$time, intvl[intvl > new_event$time])
    net.cur = ifelse(new_event$event %in% 3: 5, 1, 0)
    
    t.cur =new_event$time
  } else if (force_one & net_status[1] == net_status[2] 
             # & net_status[2]  == 0
              ){ 
    new_event = inverse_event(net.cur, intvl, params.cur, epi_times, sick_link)
    added_events = rbind.data.frame(added_events, new_event)
    intvl = c(new_event$time, intvl[intvl > new_event$time])
    net.cur = ifelse(new_event$event %in% 3: 5, 1, 0)
    
    
    new_event = inverse_event(net.cur, intvl, params.cur, epi_times, sick_link)
    added_events = rbind.data.frame(added_events, new_event)
    intvl = c(new_event$time, intvl[intvl > new_event$time])
    net.cur = ifelse(new_event$event %in% 3: 5, 1, 0)
    
    t.cur =new_event$time
  } else{
    t.cur = intvl[1]
  }
  
  t_up = intvl[intvl > t.cur][1]
  while (t.cur < max(intvl)) {
    
    # Determine the next network event type
    next_type = ifelse(net.cur == 0, 3, 6) + 
      (epi_times["i1"] <= intvl[1] & intvl[1] < epi_times["r1"]) + 
      (epi_times["i2"] <= intvl[1] & intvl[1] < epi_times["r2"])
    
    t.cur = min(rexp(1, params.cur[next_type]) + t.cur, t_up)
    if (t.cur < t_up){ 
      # If there is an network event happen at current subinterval, we record this event
      new_event = data.frame(time = t.cur,
                             event = next_type,
                             per1 = sick_per,
                             per2 = sick_nei)
      added_events = rbind.data.frame(added_events, new_event)
      intvl = c(new_event$time, intvl[intvl > new_event$time])
      net.cur = ifelse(new_event$event %in% 3: 5, 1, 0)
      
      # Then, we add the inverse event
      new_event = inverse_event(net.cur, intvl, params.cur, epi_times, sick_link)
      added_events = rbind.data.frame(added_events, new_event)
      intvl = c(new_event$time, intvl[intvl > new_event$time])
      net.cur = ifelse(new_event$event %in% 3: 5, 1, 0)
      
      t.cur = new_event$time
      t_up = intvl[intvl > t.cur][1]
    }else{ # if no event happens at current subinterval, we move to next subinterval
      t.cur = t_up
      t_up = intvl[intvl > t.cur][1]
    }
  }
  added_events
}

# 01/05/23
# Added by Houjie to return th network at each observation time
net_clipper = function(net_events, G0, obs_time){
  # missing_dats = miss_dats1; window_length = 7
  G0 = as.matrix(G0)
  all_net_snapshots = array(0, dim = c(nrow(G0), nrow(G0), length(obs_time)))
  all_net_snapshots[, , 1] = G0
  network = net_events$time
  for (i in 1: (length(obs_time)-1)) {
    # i = 1
    G1 = G0
    t0 = obs_time[i]
    t1 = obs_time[(i+1)]
    temp_net_events = net_events[t0 < net_events$time & net_events$time <= t1, ] # Filter for the network events happening between t0 and t1
    
    
    for (j in 1: nrow(temp_net_events)) { # Evolve the network from t0 to t1
      # j = 123
      temp_events_net_j = temp_net_events[j, ]
      
      # Set the i,j entry in the adjacency matrix to be 1 for connection and 0 for disconnection
      G1[temp_events_net_j$per1, temp_events_net_j$per2] <- 
        G1[temp_events_net_j$per2, temp_events_net_j$per1] <- as.integer(temp_events_net_j$event %in% 3: 5)
      
    }
    
    all_net_snapshots[, , (i+1)] = G1
    # Setting the network at t1 as G0, then we move to the second time window
    G0 = G1
  }
  all_net_snapshots
}
# 12/28/22
# Added by Houjie to impute some new network events for a given interval with network status at endpoints known.
net_event_imputer = function(input_event, ub, params){
  # input_event: a matrix with columns of "time", "event", "per1" and "per2". 
  # ub: the time upper bound of all the network events to be sampled. Must be as long as the number of rows of input_event
  # params: a vector of parameters of length 8
  idx_sub = 1: nrow(input_event)
  all_new_events = c()
  while (length(idx_sub) > 0) {
    new_events = input_event[idx_sub, ]
    new_events$event = new_events$event + ifelse(new_events$event %in% 3: 5, 3, -3) # Connect a non-link or disconnect a link
    new_events$time = sam_texp(nrow(new_events), params[new_events$event], new_events$time, ub[idx_sub]) # Sample the time of the new events
    
    
    new_inverse_events = new_events
    new_inverse_events$event = new_inverse_events$event + ifelse(new_inverse_events$event %in% 3: 5, 3, -3)
    new_inverse_events$time = sam_texp(nrow(new_inverse_events), params[new_inverse_events$event], new_inverse_events$time, ub[idx_sub])
    
    idx_sub = idx_sub[which(runif(nrow(new_inverse_events)) <= dpois(1, params[new_inverse_events$event] * (ub[idx_sub] - new_inverse_events$time)))]
    # (1 - exp(-params[new_inverse_events$event] * (ub[idx_sub] - new_inverse_events$time))))]
    
    all_new_events = rbind(all_new_events, new_events, new_inverse_events)
  }
  all_new_events
}
# 11/22/22 
# Added by Houjie to remove reconnections and redisconnections for data without recovery times
net_coarsener = function(missing_dats, window_length){
  # missing_dats = miss_dats1; window_length = 7
  G0 = as.matrix(missing_dats$G0); I0 = missing_dats$I0; events = missing_dats$events
  report = missing_dats$report; report.times = missing_dats$report.times
  report_window = report.times[2] - report.times[1]
  
  network_time = events$time
  obs_time = c(0, unique(ceiling(network_time/window_length)*window_length)) #Discrete time points we observe the network
  events_new<- c()
  for (i in 1: (length(obs_time)-1)) {
    # i = 19
    G1 = G0
    t0 = obs_time[i]
    t1 = obs_time[(i+1)]
    temp_events = events[events$time >= t0 & events$time <= t1, ] # Pick the events happening between t0 and t1
    temp_events_net = temp_events[which(!temp_events$event %in% c(1, 2)), ] # Filter for the network events happening between t0 and t1
    
    sick_per1_time = temp_events[temp_events$event == 1, c(1, 3)] # Pick the infection events between t0 and t1
    
    # This vector will record the connection happening before a patient is infected
    keep_real_time = vector("logical", length = nrow(temp_events_net)) 
    for (j in 1: nrow(temp_events_net)) { # Evolve the network from t0 to t1
      temp_events_net_j = temp_events_net[j, ]
      
      # Set the i,j entry in the adjacency matrix to be 1 for connection and 0 for disconnection
      G1[temp_events_net_j$per1, temp_events_net_j$per2] <- 
        G1[temp_events_net_j$per2, temp_events_net_j$per1] <- as.integer(temp_events_net_j$event %in% 3: 5)
      
      # Identify the sick time of future patients involved in the current event
      temp_time = sick_per1_time$time[sick_per1_time$per1 %in% c(temp_events_net_j$per1, temp_events_net_j$per2)]
      
      # If this event is related to future patients, if this event is a connection, 
      # and if this event happens before either of the involved patients get sick, this net event is kept
      if (length(temp_time) > 0 ){
        keep_real_time[j] = temp_events_net_j$time <= max(temp_time) & temp_events_net_j$event %in% 3: 5
      }
    }
    
    # Find the health status of all persons at last report time
    epi_status = report[(floor(t0/report_window)) + 1, ]
    
    # Update health status some persons healthy at last report time but getting sick sometime before t1
    epi_status[events$per1[events$time > report.times[(floor(t0/report_window)) + 1] & 
                             events$time <= t1 & events$event == 1]] = 1
    epi_status[epi_status == -1] = 0
    epi_status[temp_events$per1[temp_events$event == 1]] = 1
    
    # Find the links to disappear at t1
    disconnect = which(G1 - G0 == -1, arr.ind = T)
    
    # Identify the link type
    event_type_discon = rowSums(cbind(epi_status[disconnect[, 1]], epi_status[disconnect[, 2]])) + 6
    
    # Find the links to appear at t1
    connect = which(G1 - G0 == 1, arr.ind = T)
    event_type_connect = rowSums(cbind(epi_status[connect[, 1]], epi_status[connect[, 2]])) + 3
    
    # Summarize the net events obtained by comparing the network at t1 and t0
    events_new_temp = cbind.data.frame(time = t1, 
                                       event = c(event_type_connect, event_type_discon), 
                                       per1 = c(connect[, 1], disconnect[, 1]),
                                       per2 = c(connect[, 2], disconnect[, 2]))
    events_new_temp[events_new_temp$per1 >= events_new_temp$per2, c("per1", "per2")] = 
      events_new_temp[events_new_temp$per1 >= events_new_temp$per2, c("per2", "per1")]
    events_new_temp = unique(events_new_temp) 
    
    # Summarize the connections happening before person getting sick
    events_to_keep = rbind(temp_events[which(temp_events$event == 1), 1: 4], 
                           temp_events_net[which(keep_real_time), 1: 4])
    events_to_keep = events_to_keep[order(events_to_keep$time), ] # reorder the connections are sorted based on time
    
    # Combine the network events
    events_new = rbind.data.frame(events_new, 
                                  rbind.data.frame(events_to_keep, events_new_temp))
    
    # Setting the network at t1 as G0, then we move to the second time window
    G0 = G1
  }
  
  # line 77~83 make sure there's no reconnections before related patients get sick 
  # (i.e. we only ask when did you meet that person for the first time)
  kept_net = events_new[events_new$event == 3 & events_new$time %% 1 != 0, ]
  kept_net[kept_net$per1 >= kept_net$per2, c("per1", "per2")] = 
    kept_net[kept_net$per1 >= kept_net$per2, c("per2", "per1")]
  duplicated_idx = which(rownames(events_new) == rownames(kept_net[duplicated(kept_net[, c("per1", "per2")]), ]))
  if(length(duplicated_idx) > 0){
    events_new = events_new[-duplicated_idx, ] 
  }
  
  missing_dats_out = missing_dats
  missing_dats_out$events = events_new
  missing_dats_out
}

# 12/05/2022
# Added by Houjie to remove reconnections and redisconnections for complete data
net_coarsener2 = function(dats, window_length){
  # dats = full_dats; window_length = 1;
  G0 = as.matrix(dats$G0); I0 = dats$I0; events = dats$events
  
  network_time = events$time
  obs_time = c(0, unique(ceiling(network_time/window_length)*window_length)) #Discrete time points we observe the network
  events_new<- c()
  epi_status = rep(0, nrow(G0))
  epi_status[I0] = 1
  for (i in 1: (length(obs_time)-1)) {
    # i = 52
    G1 = G0
    t0 = obs_time[i]
    t1 = obs_time[(i+1)]
    temp_events = events[events$time >= t0 & events$time <= t1, ] # Pick the events happening between t0 and t1
    temp_events_net = temp_events[which(!temp_events$event %in% c(1, 2)), ] # Filter for the network events happening between t0 and t1
    
    sick_per1_time = temp_events[temp_events$event == 1, c(1, 3)] # Pick the infection events between t0 and t1
    
    # This vector will record the connection happening before a patient is infected
    keep_real_time = vector("logical", length = nrow(temp_events_net)) 
    for (j in 1: nrow(temp_events_net)) { # Evolve the network from t0 to t1
      temp_events_net_j = temp_events_net[j, ]
      
      # Set the i,j entry in the adjacency matrix to be 1 for connection and 0 for disconnection
      G1[temp_events_net_j$per1, temp_events_net_j$per2] <- 
        G1[temp_events_net_j$per2, temp_events_net_j$per1] <- as.integer(temp_events_net_j$event %in% 3: 5)
      
      # Identify the sick time of future patients involved in the current event
      temp_time = sick_per1_time$time[sick_per1_time$per1 %in% c(temp_events_net_j$per1, temp_events_net_j$per2)]
      
      # If this event is related to future patients, if this event is a connection, 
      # and if this event happens before either of the involved patients get sick, this net event is kept
      if (length(temp_time) > 0 ){
        keep_real_time[j] = temp_events_net_j$time <= max(temp_time) & temp_events_net_j$event %in% 3: 5
      }
    }
    
    # Find the health status of all persons at t1 to determine the link type
    temp_events_epi = temp_events[temp_events$event %in% 1: 2, ]
    temp_events_epi$event[temp_events_epi$event == 2] = -1
    for (j in 1: nrow(temp_events_epi)) {
      # j = 1
      temp_events_epi_j = temp_events_epi[j, ]
      epi_status[temp_events_epi_j$per1] = epi_status[temp_events_epi_j$per1] + temp_events_epi_j$event
    }
    
    # Find the links to disappear at t1
    disconnect = which(G1 - G0 == -1, arr.ind = T)
    
    # Identify the link type
    event_type_discon = rowSums(cbind(epi_status[disconnect[, 1]], epi_status[disconnect[, 2]])) + 6
    
    # Find the links to appear at t1
    connect = which(G1 - G0 == 1, arr.ind = T)
    event_type_connect = rowSums(cbind(epi_status[connect[, 1]], epi_status[connect[, 2]])) + 3
    
    # Summarize the net events obtained by comparing the network at t1 and t0
    if (length(c(event_type_connect, event_type_discon)) != 0){
      events_new_temp = cbind.data.frame(time = t1, 
                                         event = c(event_type_connect, event_type_discon), 
                                         per1 = c(connect[, 1], disconnect[, 1]),
                                         per2 = c(connect[, 2], disconnect[, 2]))
      events_new_temp[events_new_temp$per1 >= events_new_temp$per2, c("per1", "per2")] = 
        events_new_temp[events_new_temp$per1 >= events_new_temp$per2, c("per2", "per1")]
      events_new_temp = unique(events_new_temp) 
    } else{
      events_new_temp = c()
    }
    
    # Summarize the connections happening before person getting sick and epi events
    events_to_keep = rbind(temp_events[which(!temp_events$event %in% 3: 8), 1: 4], 
                           temp_events_net[which(keep_real_time), 1: 4])
    events_to_keep = events_to_keep[order(events_to_keep$time), ] # reorder the connections are sorted based on time
    
    # Combine the network events
    events_new = rbind.data.frame(events_new, 
                                  rbind.data.frame(events_to_keep, events_new_temp))
    
    # Setting the network at t1 as G0, then we move to the second time window
    G0 = G1
  }
  
  # line 77~83 make sure there's no reconnections before related patients get sick 
  # (i.e. we only ask when did you meet that person for the first time)
  kept_net = events_new[events_new$event == 3 & events_new$time %% 1 != 0, ]
  kept_net[kept_net$per1 >= kept_net$per2, c("per1", "per2")] = 
    kept_net[kept_net$per1 >= kept_net$per2, c("per2", "per1")]
  duplicated_idx = which(rownames(events_new) == rownames(kept_net[duplicated(kept_net[, c("per1", "per2")]), ]))
  if(length(duplicated_idx) > 0){
    events_new = events_new[-duplicated_idx, ] 
  }
  
  dats_out = dats
  dats_out$events = events_new
  dats_out$infer = NULL
  dats_out$infer2 = NULL
  dats_out
}

# 12/29/2022
# Bug fixed for net_coarsener. 
net_coarsener3 = function(missing_dats = miss_dats1, window_length = 7){
  # missing_dats = miss_dats1; window_length = 7
  G0 = as.matrix(missing_dats$G0); I0 = missing_dats$I0; events = missing_dats$events
  report = missing_dats$report; report.times = missing_dats$report.times
  report_window = report.times[2] - report.times[1]
  
  network_time = events$time
  obs_time = c(0, unique(ceiling(network_time/window_length)*window_length)) #Discrete time points we observe the network
  events_new<- c()
  for (i in 1: (length(obs_time)-1)) {
    # i = 1
    G1 = G0
    t0 = obs_time[i]
    t1 = obs_time[(i+1)]
    temp_events = events[events$time >= t0 & events$time <= t1, ] # Pick the events happening between t0 and t1
    temp_events_net = temp_events[which(!temp_events$event %in% c(1, 2)), ] # Filter for the network events happening between t0 and t1
    
    sick_per1_time = temp_events[temp_events$event == 1, c(1, 3)] # Pick the infection events between t0 and t1
    
    # This vector will record the connection happening before a patient is infected
    keep_real_time = vector("logical", length = nrow(temp_events_net)) 
    for (j in 1: nrow(temp_events_net)) { # Evolve the network from t0 to t1
      # j = 123
      temp_events_net_j = temp_events_net[j, ]
      
      # Set the i,j entry in the adjacency matrix to be 1 for connection and 0 for disconnection
      G1[temp_events_net_j$per1, temp_events_net_j$per2] <- 
        G1[temp_events_net_j$per2, temp_events_net_j$per1] <- as.integer(temp_events_net_j$event %in% 3: 5)
      
      # Identify the sick time of future patients involved in the current event
      temp_time = sick_per1_time$time[sick_per1_time$per1 %in% c(temp_events_net_j$per1, temp_events_net_j$per2)]
      
      # If this event is related to future patients, if this event is a connection, 
      # and if this event happens before either of the involved patients get sick, this net event is kept
      if (length(temp_time) > 0 ){
        keep_real_time[j] = temp_events_net_j$time <= max(temp_time) & temp_events_net_j$event %in% 3: 5
        
        if(keep_real_time[j]){
          # Evolve G0 as well because we will find the coarsened events by G1-G0
          G0[temp_events_net_j$per1, temp_events_net_j$per2] <- 
            G0[temp_events_net_j$per2, temp_events_net_j$per1] <- as.integer(temp_events_net_j$event %in% 3: 5)
        }
      }
    }
    
    # Find the health status of all persons at the report time smaller and nearest to t1 from below
    epi_status = report[max(which(report.times <= t1)), ]
    # epi_status = report[(floor(t0/report_window)) + 1, ]
    
    # Add the changes of health status happening between the small nearest report time and t1
    epi_status[temp_events$per1[temp_events$event == 1 & temp_events$time >= report.times[max(which(report.times <= t1))]]] = 1
    epi_status[temp_events$per1[temp_events$event == 2 & temp_events$time >= report.times[max(which(report.times <= t1))]]] = -1
    
    # let 0 denote healthy and 1 denote infection 
    epi_status[epi_status == -1] = 0
    
    # Find the links to disappear at t1
    disconnect = which(G1 - G0 == -1, arr.ind = T)
    
    # Identify the link type
    event_type_discon = rowSums(cbind(epi_status[disconnect[, 1]], epi_status[disconnect[, 2]])) + 6
    
    # Find the links to appear at t1
    connect = which(G1 - G0 == 1, arr.ind = T)
    event_type_connect = rowSums(cbind(epi_status[connect[, 1]], epi_status[connect[, 2]])) + 3
    
    # Summarize the net events obtained by comparing the network at t1 and t0
    events_new_temp = cbind.data.frame(time = t1, 
                                       event = c(event_type_connect, event_type_discon), 
                                       per1 = c(connect[, 1], disconnect[, 1]),
                                       per2 = c(connect[, 2], disconnect[, 2]))
    events_new_temp[events_new_temp$per1 >= events_new_temp$per2, c("per1", "per2")] = 
      events_new_temp[events_new_temp$per1 >= events_new_temp$per2, c("per2", "per1")]
    events_new_temp = unique(events_new_temp) 
    
    # Make sure there's no reconnections before related patients get sick 
    # (i.e. we only ask when did you meet that person for the first time)
    kept_net = temp_events_net[which(keep_real_time), 1: 4]
    kept_net[kept_net$per1 >= kept_net$per2, c("per1", "per2")] = 
      kept_net[kept_net$per1 >= kept_net$per2, c("per2", "per1")]
    # kept_net = rbind(kept_net, c(6.01, 3, 47, 56))
    kept_net = kept_net[!duplicated(kept_net[, c("per1", "per2")]), ]
    
    # Summarize the connections happening before person getting sick
    events_to_keep = rbind(temp_events[which(temp_events$event == 1), 1: 4], # infection events
                           kept_net) # connections to keep
    events_to_keep = events_to_keep[order(events_to_keep$time), ] # reorder the connections are sorted based on time
    
    # Combine the network events
    events_new = rbind.data.frame(events_new, 
                                  rbind.data.frame(events_to_keep, events_new_temp))
    
    # Setting the network at t1 as G0, then we move to the second time window
    G0 = G1
  }
  
  missing_dats_out = missing_dats
  missing_dats_out$events = events_new
  missing_dats_out
}


# 1. utility functions for evaluating likelihood
# 1.1 trace system status to a given time point
# default model: SIR
obtain_system <- function(G0, I0, events, t.obs, 
                          model="SIR", quarantine=T){
  tmax =  max(events$time)
  if(t.obs >  tmax){
    cat("Required time out of range, use tmax =",tmax, "instead.\n")
    t.obs = tmax
  }
  rowmax = min(which(events$time >= t.obs))
  
  N = nrow(G0); epid = rep(0,N)
  adjmat =  G0; infected =  I0; epid[infected] = 1
  
  for(r in 1:rowmax){
    z =  events[r,]$event
    if (z==1){
      # infection
      p1 = events[r,]$per1
      epid[p1] = 1
    }else if (z==2){
      # recovery
      p1 = events[r,]$per1
      if(model=="SIR"){
        epid[p1] = -1
      }else{
        epid[p1] = 0
      }
    }else{
      # some edge stuff
      p1 = events[r,]$per1
      p2 = events[r,]$per2
      if(quarantine){
        if(z %in% c(3:5)){
          # reconnection
          adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        }else{
          # disconnection
          adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        }
      }else{
        if(z==3){
          adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        }else{
          adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        }
      }
    }
    # report progress
    cat("Processing...", r/rowmax,"\r")
  }
  
  # also return useful summary statistics
  N.I = sum(epid==1); N.S = sum(epid==0)
  # cat("At time", t.obs, ":\n",
  #     "epid vector =", epid, "\n")
  if(quarantine){
    susceptible = which(epid==0)
    infected = which(epid==1)
    Mmax.SS = N.S*(N.S-1)/2
    Mmax.SI = N.S * N.I
    Mmax.II = N.I*(N.I-1)/2
    if(Mmax.SS==0){
      M.SS = 0
    }else{
      M.SS = nnzero(adjmat[susceptible,susceptible])/2
    }
    if(Mmax.SI==0){
      M.SI = 0
    }else{
      M.SI = nnzero(adjmat[susceptible,infected])
    }
    if(Mmax.II==0){
      M.II = 0
    }else{
      M.II = nnzero(adjmat[infected,infected])/2
    }
    
    # stats = c(N.S, N.I, 
    #           Mmax.SS, Mmax.SI, Mmax.SS,
    #           M.SS, M.SI, M.II)
    # names(stats) =  c("N.S","N.I","Mmax.SS","Mmax.SI","Mmax.SS",
    #                   "M.SS","M.SI","M.II")
    
    stats = data.frame(N.S = N.S, N.I = N.I, 
                       Mmax.SS=Mmax.SS, Mmax.SI=Mmax.SI, Mmax.II=Mmax.II,
                       M.SS=M.SS, M.SI=M.SI, M.II=M.II)
    
  }else{
    Mmax = N*(N-1)/2
    M =  nnzero(adjmat)
    
    if(N.S > 0 & N.I > 0){
      N.SI = nnzero(adjmat[susceptible,infected])
    }else{
      N.SI = 0
    }
    
    # stats  = c(N.S, N.I, Mmax, M)
    # names(stats) = c("N.S","N.I","Mmax", "M")
    
    stats  = data.frame(N.S = N.S, N.I = N.I, N.SI = N.SI, 
                        Mmax = Mmax, M = M)
  }
  
  return(list(adjmat = adjmat, epid = epid,
              stats = stats))
}


# 1.2 evaluate the log likelihood of events in a given interval [st,en)
## use log-likelihood for better computation accuracy
## only deal with the coupled case for now
## 05/09/2019
## augmentation: calculate important stats (for inference) as well
eval_interval_loglik <- function(dat, st, en, model="SIR", quarantine=T,
                                 bet = 0.03, gam = 0.15, 
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05),
                                 cal_stats = T){
  
  params = c(bet, gam, alpha.r, alpha.d)
  if(quarantine & length(params) < 8){
    stop("Type-dependent edge rates should be specified! What are you using instead??\n")
  }
  if(!quarantine & length(params) > 4){
    stop("Type-independent edge rates should be specified! What are you using instead??\n")
  }
  
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  tmax =  max(events$time)
  if(st >= en){
    stop("Invalid time interval.\n")
  }
  if(en >  tmax){
    cat("Required ending time out of range, use tmax =",tmax, "instead.\n")
    en = tmax
  }
  
  rowmin = min(which(events$time > st))
  rowmax = max(which(events$time <= en))
  
  if(rowmin > 1){
    sys = obtain_system(G0, I0, events, st, model, quarantine)
    adjmat = sys$adjmat; epid = sys$epid; stats = sys$stats
    
    t.cur = events$time[rowmin-1]
  }else{
    adjmat = G0; epid = rep(0,nrow(G0)); epid[I0] = 1
    susceptible = which(epid==0); infected = which(epid==1)
    N.S = sum(epid==0); N.I = sum(epid==1)
    stats = data.frame(N.S = N.S, N.I = N.I,
                       Mmax.SS = N.S*(N.S-1)/2, 
                       Mmax.SI = N.S * N.I, 
                       Mmax.II = N.I*(N.I-1)/2,
                       M.SS = nnzero(adjmat[susceptible,susceptible])/2, 
                       M.SI = nnzero(adjmat[susceptible,infected]), 
                       M.II = nnzero(adjmat[infected,infected])/2)
    
    t.cur = 0
  }
  
  N.S = stats$N.S; N.I = stats$N.I
  if(quarantine){
    # basic edge stats
    N.SI = stats$M.SI
    M = c(stats$M.SS, stats$M.SI, stats$M.II)
    Mmax= c(stats$Mmax.SS, stats$Mmax.SI, stats$Mmax.II)
    #M.d = Mmax - M
    cat("At the start...\n",
        "N.S:",N.S, "N.I:", N.I,"\n",
        "M:",M,"\n",
        "Mmax:",Mmax,"\n",
        "Md:",Mmax-M,"\n")
  }else{
    stop("This function doesn't deal with the decoupled case for now.\n")
  }
  
  n.E = sum(events$event[rowmin:rowmax] == 1)
  n.R = sum(events$event[rowmin:rowmax] == 2)
  C = rep(0,3); D = rep(0,3)
  
  logsum = 0; expo = 0
  if(cal_stats){
    # data storage for "big.sums"
    big.sums = numeric(length(params))
  }else{
    big.sums = 0
  }
  
  for(r in rowmin:rowmax){
    t.next = events$time[r]
    del.t = t.next - t.cur
    
    
    # update "big.sums"
    if(cal_stats){
      big.sums = big.sums + params * c(N.SI, N.I, Mmax-M, M) * del.t
    }
    # update the exponent part (excluding "-")
    expo = expo + 
      sum(params * c(N.SI, N.I, Mmax-M, M)) * del.t
    
    z =  events$event[r]
    if (z==1){
      # infection
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # update logsum term
      logsum = logsum + log(I.p1)
      # bookkeeping
      epid[p1] = 1
      N.S = N.S - 1
      N.I = N.I + 1
      Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
      M[1] = M[1] - S.p1
      M[2] = M[2] + S.p1 - I.p1
      M[3] = M[3] + I.p1
      N.SI = M[2]
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      N.I = N.I - 1
      if(model=="SIR"){
        epid[p1] = -1
        M[1] = M[1]
        M[2] = M[2] - S.p1
        M[3] = M[3] - I.p1
      }else{
        epid[p1] = 0
        N.S = N.S + 1
        M[1] = M[1] + S.p1
        M[2] = M[2] - S.p1 + I.p1
        M[3] = M[3] - I.p1
      }
      Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
      N.SI = M[2]
    }else{
      # some edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]

      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        # bookkeeping
        if (epid[p1]==0 & epid[p2]==0){
          # S-S type
          logsum = logsum + log(Mmax[1]-M[1])
          C[1] = C[1] + 1
          M[1] = M[1] + 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          logsum = logsum + log(Mmax[3]-M[3])
          C[3] = C[3] + 1
          M[3] = M[3] + 1
        }else{
          # S-I type
          logsum = logsum + log(Mmax[2]-M[2])
          C[2] = C[2] + 1
          M[2] = M[2] + 1
          N.SI = M[2]
        }
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        # bookkeeping
        if (epid[p1]==0 & epid[p2]==0){
          # S-S type
          logsum = logsum + log(M[1])
          D[1] = D[1] + 1
          M[1] = M[1] - 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          logsum = logsum + log(M[3])
          D[3] = D[3] + 1
          M[3] = M[3] - 1
        }else{
          # S-I type
          logsum = logsum + log(M[2])
          D[2] = D[2] + 1
          M[2] = M[2] - 1
          N.SI = M[2]
        }
      }
      
    }
    # report progress
    cat("Processing...", 
        (r-rowmin)/(rowmax-rowmin),"\r")
    
    t.cur = t.next
    
  }
  cat("\nDone.\n")
  cat("Event counts:",c(n.E, n.R, C, D),"\n",
      "logsum:",logsum,"\n",
      "expo:",expo,"\n")
  
  # calculate log-lik
  ll = logsum + 
    sum(log(params)*c(n.E, n.R, C, D))-
    expo
  
  res = list(ll=ll, event.counts = c(n.E, n.R, C, D), big.sums = big.sums)
  
  return(res)
}


# 2.3 a function that augments the observed data w/ imputed recovery times
# and outputs event counts & all the summations
# assume a coupled process
# **deprecated version**
# it disregards all the edge changes related to R people
parse_augment <- function(G0, I0, events, recovers, model="SIR"){
  # dats: a list consisting of everything observed
  # recovers: a dataset w/ variables `time` & `per1`, all imputed recovery times

  N = nrow(G0); epid = rep(0,N)
  adjmat =  G0; infected =  I0; epid[infected] = 1
  
  susceptible = which(epid==0); infected = which(epid==1)
  N.S = sum(epid==0); N.I = sum(epid==1)
  Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
  M = c(nnzero(adjmat[susceptible,susceptible])/2, 
        nnzero(adjmat[susceptible,infected]), 
        nnzero(adjmat[infected,infected])/2)
  N.SI = M[2]
  
  # data storage
  n.E = sum(events$event == 1)
  n.R = sum(events$event == 2) + nrow(recovers)
  C = rep(0,3); D = rep(0,3)
  
  # FOR NOW
  # assume the case of "quarantine" - coupled case
  big.sums = numeric(8) 
  
  # combine the two datasets to make a "complete" dataset
  # and order by ascending time
  events = events[,-c(5:6)]
  recovers$per2 = NA
  recovers$event = 2
  recovers = recovers[,names(events)]
  events = rbind(events, recovers)
  events = events[order(events$time),]
  
  # parse through the augmented data
  t.cur = 0
  for(r in 1:nrow(events)){
    t.next = events$time[r]
    del.t = t.next - t.cur
    
    # update "big.sums"
    big.sums = big.sums + c(N.SI, N.I, Mmax-M, M) * del.t

    z =  events$event[r]
    if (z==1){
      # infection
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      epid[p1] = 1
      N.S = N.S - 1
      N.I = N.I + 1
      Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
      M[1] = M[1] - S.p1
      M[2] = M[2] + S.p1 - I.p1
      M[3] = M[3] + I.p1
      N.SI = M[2]
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      N.I = N.I - 1
      if(model=="SIR"){
        epid[p1] = -1
        M[1] = M[1]
        M[2] = M[2] - S.p1
        M[3] = M[3] - I.p1
      }else{
        epid[p1] = 0
        N.S = N.S + 1
        M[1] = M[1] + S.p1
        M[2] = M[2] - S.p1 + I.p1
        M[3] = M[3] - I.p1
      }
      Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
      N.SI = M[2]
    }else{
      # some edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]
     
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        # bookkeeping
        if (epid[p1]==0 & epid[p2]==0){
          # S-S type
          C[1] = C[1] + 1
          M[1] = M[1] + 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          C[3] = C[3] + 1
          M[3] = M[3] + 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type: neither is recovered
          C[2] = C[2] + 1
          M[2] = M[2] + 1
          N.SI = M[2]
        }
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        # bookkeeping
        if (epid[p1]==0 & epid[p2]==0){
          # S-S type
          D[1] = D[1] + 1
          M[1] = M[1] - 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          D[3] = D[3] + 1
          M[3] = M[3] - 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
          N.SI = M[2]
        }
      }
      
    }
    # report progress
    # cat("Processing augmented dataset...", 
    #     r/nrow(events),"\r")
    
    t.cur = t.next

  }
  #cat("\nDone.\n")
  
  res = list(event.counts = c(n.E, n.R, C, D), big.sums = big.sums)
  
  return(res)
}


# 05/23/2019
# the actually used version
# modify the `parse_augment` function to include R-related edge changes
# 06/02/2019
# debugged
parse_augment2 <- function(G0, I0, events, recovers, model="SIR"){
  # dats: a list consisting of everything observed
  # recovers: a dataset w/ variables `time` & `per1`, all imputed recovery times
  
  N = nrow(G0); epid = rep(0,N)
  adjmat =  G0; epid[I0] = 1
  
  susceptible = which(epid==0); infected = which(epid==1)
  N.S = sum(epid==0); N.I = sum(epid==1)
  if(model == "SIR"){ N.R = 0 }
  Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
  M = c(nnzero(adjmat[susceptible,susceptible])/2, 
        nnzero(adjmat[susceptible,infected]), 
        nnzero(adjmat[infected,infected])/2)
  N.SI = M[2]
  
  # data storage
  n.E = sum(events$event == 1)
  n.R = sum(events$event == 2) + nrow(recovers)
  C = rep(0,3); D = rep(0,3)
  
  # FOR NOW
  # assume the coupled case 
  big.sums = numeric(8) 
  
  # combine the two datasets to make a "complete" dataset
  # and order by ascending time
  events = events[,-c(5:6)]
  recovers$per2 = NA
  recovers$event = 2
  recovers = recovers[,names(events)]
  events = rbind(events, recovers)
  events = events[order(events$time),]
  
  # parse through the augmented data
  t.cur = 0
  for(r in 1:nrow(events)){
    t.next = events$time[r]
    del.t = t.next - t.cur
    
    # update "big.sums"
    big.sums = big.sums + c(N.SI, N.I, Mmax-M, M) * del.t
    
    z = events$event[r]
    if (z==1){
      # infection
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      epid[p1] = 1
      N.S = N.S - 1
      N.I = N.I + 1
      # update counts
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        #Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        Mmax = c((N-N.I)*(N-N.I-1)/2, (N-N.I) * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1 - R.p1
        M[2] = M[2] + S.p1 + R.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = N.SI + S.p1 - I.p1
      }else{
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1
        M[2] = M[2] + S.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = M[2]
      }
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      N.I = N.I - 1
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        epid[p1] = -1
        N.R = N.R + 1
        M[1] = M[1] + S.p1 + R.p1
        M[2] = M[2] - S.p1 - R.p1 + I.p1
        M[3] = M[3] - I.p1
        #Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        Mmax = c((N-N.I)*(N-N.I-1)/2, (N-N.I) * N.I, N.I*(N.I-1)/2)
        N.SI = N.SI - S.p1
      }else{
        epid[p1] = 0
        N.S = N.S + 1
        M[1] = M[1] + S.p1
        M[2] = M[2] - S.p1 + I.p1
        M[3] = M[3] - I.p1
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        N.SI = M[2]
      }
    }else{
      # some edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]
  
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          C[1] = C[1] + 1
          M[1] = M[1] + 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          C[3] = C[3] + 1
          M[3] = M[3] + 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type: neither is recovered
          C[2] = C[2] + 1
          M[2] = M[2] + 1
          N.SI = N.SI + 1
        }else{
          # R-I type
          C[2] = C[2] + 1
          M[2] = M[2] + 1
        }
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          D[1] = D[1] + 1
          M[1] = M[1] - 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          D[3] = D[3] + 1
          M[3] = M[3] - 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
          N.SI = N.SI - 1
        }else{
          # R-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
        }
      }
      
    }
    # report progress
    # cat("Processing augmented dataset...", 
    #     r/nrow(events),"\r")
    
    t.cur = t.next
    
  }
  #cat("\nDone.\n")
  
  res = list(event.counts = c(n.E, n.R, C, D), big.sums = big.sums)
  
  return(res)
}

# A helper function to output the counts of link types at each network event
parse_augment2_test = function(G0_temp = G0, I0, temp_events = full_events_w_recov[full_events_w_recov$event != 2, ], 
                               temp_recovers = full_events_w_recov[full_events_w_recov$event == 2, c(1, 3)], model="SIR"){
  # dats: a list consisting of everything observed
  # recovers: a dataset w/ variables `time` & `per1`, all imputed recovery times
  
  N = nrow(G0_temp); epid = rep(0,N)
  adjmat =  G0_temp; epid[I0] = 1
  
  susceptible = which(epid==0); infected = which(epid==1)
  N.S = sum(epid==0); N.I = sum(epid==1)
  if(model == "SIR"){ N.R = 0 }
  Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
  M = c(nnzero(adjmat[susceptible,susceptible])/2, 
        nnzero(adjmat[susceptible,infected]), 
        nnzero(adjmat[infected,infected])/2)
  N.SI = M[2]
  
  # data storage
  n.E = sum(temp_events$event == 1)
  n.R = sum(temp_events$event == 2) + nrow(temp_recovers)
  C = rep(0,3); D = rep(0,3)
  
  # FOR NOW
  # assume the coupled case 
  big.sums = numeric(8) 
  
  # combine the two datasets to make a "complete" dataset
  # and order by ascending time
  temp_events = temp_events[,-c(5:6)]
  temp_recovers$per2 = NA
  temp_recovers$event = 2
  temp_recovers = temp_recovers[,names(temp_events)]
  temp_events = rbind(temp_events, temp_recovers)
  temp_events = temp_events[order(temp_events$time),]
  
  # parse through the augmented data
  t.cur = 0
  all_big_sums = c()
  for(r in 1:nrow(temp_events)){
    # r = 158
    t.next = temp_events$time[r]
    del.t = t.next - t.cur
    
    # update "big.sums"
    big.sums = big.sums + c(N.SI, N.I, Mmax-M, M) * del.t
    all_big_sums = rbind(all_big_sums, c(N.SI, N.I, Mmax-M, M, del.t))
    z = temp_events$event[r]
    if (z==1){
      # infection
      p1 = temp_events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      epid[p1] = 1
      N.S = N.S - 1
      N.I = N.I + 1
      # update counts
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        #Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        Mmax = c((N-N.I)*(N-N.I-1)/2, (N-N.I) * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1 - R.p1
        M[2] = M[2] + S.p1 + R.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = N.SI + S.p1 - I.p1
      }else{
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1
        M[2] = M[2] + S.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = M[2]
      }
    }else if (z==2){
      # recovery
      p1 = temp_events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      N.I = N.I - 1
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        epid[p1] = -1
        N.R = N.R + 1
        M[1] = M[1] + S.p1 + R.p1
        M[2] = M[2] - S.p1 - R.p1 + I.p1
        M[3] = M[3] - I.p1
        #Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        Mmax = c((N-N.I)*(N-N.I-1)/2, (N-N.I) * N.I, N.I*(N.I-1)/2)
        N.SI = N.SI - S.p1
      }else{
        epid[p1] = 0
        N.S = N.S + 1
        M[1] = M[1] + S.p1
        M[2] = M[2] - S.p1 + I.p1
        M[3] = M[3] - I.p1
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        N.SI = M[2]
      }
    }else{
      # some edge stuff
      p1 = temp_events$per1[r]
      p2 = temp_events$per2[r]
      
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          C[1] = C[1] + 1
          M[1] = M[1] + 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          C[3] = C[3] + 1
          M[3] = M[3] + 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type: neither is recovered
          C[2] = C[2] + 1
          M[2] = M[2] + 1
          N.SI = N.SI + 1
        }else{
          # R-I type
          C[2] = C[2] + 1
          M[2] = M[2] + 1
        }
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          D[1] = D[1] + 1
          M[1] = M[1] - 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          D[3] = D[3] + 1
          M[3] = M[3] - 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
          N.SI = N.SI - 1
        }else{
          # R-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
        }
      }
      
    }
    # report progress
    # cat("Processing augmented dataset...", 
    #     r/nrow(events),"\r")
    
    t.cur = t.next
    
  }
  #cat("\nDone.\n")
  
  res = list(event.counts = c(n.E, n.R, C, D), big.sums = big.sums, big.sum.tab = all_big_sums)
  
  return(res)
}
sourceCpp("big_sumer_cpp.cpp")
parse_augment_cpp = function(G0_temp = G0, I0, temp_events = events_temp[-1, ], 
                             temp_recovers = recover.dat, model="SIR"){
  # dats: a list consisting of everything observed
  # recovers: a dataset w/ variables `time` & `per1`, all imputed recovery times
  # temp_events = events_temp[-1, ]; temp_recovers = recover.dat
  
  # combine the two datasets to make a "complete" dataset
  # and order by ascending time
  temp_events = temp_events[,-c(5:6)]
  temp_recovers$per2 = NA
  temp_recovers$event = 2
  temp_recovers = temp_recovers[,names(temp_events)]
  temp_events = rbind(temp_events, temp_recovers)
  temp_events = temp_events[order(temp_events$time),]
  temp_events$per2[is.na(temp_events$per2)] = 0
  
  
  
  big.sums = big_sumer_cpp(temp_events[, 1],
                           temp_events[, 2],
                           temp_events[, 3],
                           temp_events[, 4], 
                           as.matrix(G0_temp),
                           I0)
  
  res = list(event.counts = table(temp_events$event), big.sums = c(big.sums))
  
  return(res)
}

# 11/29/2022
# Added by Houjie to compute the event counts and big sums for imputed data with coarsened net. 
parse_augment3 <- function(G0, I0, events, recover.dat, model="SIR"){
  # dats: a list consisting of everything observed
  # recover.dat: a dataset w/ variables `time` & `per1`, all imputed recovery times
  
  N = nrow(G0); epid = rep(0,N)
  adjmat =  G0; epid[I0] = 1
  obs_time =  c(0, unique(events$time[events$time %% 1 == 0]))
  susceptible = which(epid==0); infected = which(epid==1)
  N.S = sum(epid==0); N.I = sum(epid==1)
  if(model == "SIR"){ N.R = 0 }
  Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
  M = c(nnzero(adjmat[susceptible,susceptible])/2, 
        nnzero(adjmat[susceptible,infected]), 
        nnzero(adjmat[infected,infected])/2)
  N.SI = M[2]
  
  # data storage
  n.E = sum(events$event == 1)
  n.R = sum(events$event == 2) + nrow(recover.dat)
  C = rep(0,3); D = rep(0,3)
  
  # FOR NOW
  # assume the coupled case 
  big.sums = numeric(8) 
  
  # combine the two datasets to make a "complete" dataset
  # and order by ascending time
  events = events[,-c(5:6)]
  recover.dat$per2 = NA
  recover.dat$event = 2
  recover.dat = recover.dat[,names(events)]
  
  events_imputed = rbind(events, recover.dat)
  events_imputed = events_imputed[order(events_imputed$time),]
  events_imputed = unique(events_imputed)
  
  epi_events = events_imputed[events_imputed$event %in% 1: 2, ]
  epi_events = rbind(c(0, 1, I0, NA), epi_events)
  sick_status = matrix(0, nrow = (length(obs_time)-1), ncol = N)
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
  
  events = events_imputed
  # parse through the augmented data
  t.cur = 0
  for(r in 1:nrow(events)){
    t.next = events$time[r]
    del.t = t.next - t.cur
    
    # update "big.sums"
    big.sums = big.sums + c(N.SI, N.I, Mmax-M, M) * del.t
    
    z = events$event[r]
    if (z==1){
      # infection
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      epid[p1] = 1
      N.S = N.S - 1
      N.I = N.I + 1
      # update counts
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        #Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        Mmax = c((N-N.I)*(N-N.I-1)/2, (N-N.I) * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1 - R.p1
        M[2] = M[2] + S.p1 + R.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = N.SI + S.p1 - I.p1
      }else{
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1
        M[2] = M[2] + S.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = M[2]
      }
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      N.I = N.I - 1
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        epid[p1] = -1
        N.R = N.R + 1
        M[1] = M[1] + S.p1 + R.p1
        M[2] = M[2] - S.p1 - R.p1 + I.p1
        M[3] = M[3] - I.p1
        #Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        Mmax = c((N-N.I)*(N-N.I-1)/2, (N-N.I) * N.I, N.I*(N.I-1)/2)
        N.SI = N.SI - S.p1
      }else{
        epid[p1] = 0
        N.S = N.S + 1
        M[1] = M[1] + S.p1
        M[2] = M[2] - S.p1 + I.p1
        M[3] = M[3] - I.p1
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        N.SI = M[2]
      }
    }else{
      # some edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]
      
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          C[1] = C[1] + 1
          M[1] = M[1] + 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          C[3] = C[3] + 1
          M[3] = M[3] + 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type: neither is recovered
          C[2] = C[2] + 1
          M[2] = M[2] + 1
          N.SI = N.SI + 1
        }else{
          # R-I type
          C[2] = C[2] + 1
          M[2] = M[2] + 1
        }
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          D[1] = D[1] + 1
          M[1] = M[1] - 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          D[3] = D[3] + 1
          M[3] = M[3] - 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
          N.SI = N.SI - 1
        }else{
          # R-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
        }
      }
      
    }
    # report progress
    # cat("Processing augmented dataset...", 
    #     r/nrow(events),"\r")
    
    t.cur = t.next
    
  }
  #cat("\nDone.\n")
  
  res = list(event.counts = c(n.E, n.R, C, D), big.sums = big.sums, full_events = events_imputed)
  
  return(res)
}

# 12/04/2022
# Added by Houjie to compute the event counts and big sums for an complete list of events
parse_augment4 <- function(G0, I0, events, model="SIR"){
  # G0 = G0; I0 =I0; events = events_temp2[-1, ]
  
  N = nrow(G0); epid = rep(0,N)
  adjmat =  G0; epid[I0] = 1
  obs_time =  c(0, unique(events$time[events$time %% 1 == 0]))
  susceptible = which(epid==0); infected = which(epid==1)
  N.S = sum(epid==0); N.I = sum(epid==1)
  if(model == "SIR"){ N.R = 0 }
  Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
  M = c(nnzero(adjmat[susceptible,susceptible])/2, 
        nnzero(adjmat[susceptible,infected]), 
        nnzero(adjmat[infected,infected])/2)
  N.SI = M[2]
  
  # data storage
  n.E = sum(events$event == 1)
  n.R = sum(events$event == 2)
  C = rep(0,3); D = rep(0,3)
  
  # FOR NOW
  # assume the coupled case 
  big.sums = numeric(8) 
  
  events = events[,-c(5:6)]
  # parse through the augmented data
  t.cur = 0
  all_M = matrix(NA, nrow = nrow(events), ncol = 2*length(M)+1)
  for(r in 1:nrow(events)){
    # r = 557
    t.next = events$time[r]
    del.t = t.next - t.cur
    
    # update "big.sums"
    big.sums = big.sums + c(N.SI, N.I, Mmax-M, M) * del.t
    all_M[r,] = c(M, NA, c(sum(adjmat[which(epid!=1), which(epid!=1)])/2, 
                           sum(adjmat[which(epid==0 |epid==-1), which(epid==1)]), 
                           sum(adjmat[which(epid==1), which(epid==1)])/2))
    z = events$event[r]
    if (z==1){
      # infection
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      epid[p1] = 1
      N.S = N.S - 1
      N.I = N.I + 1
      # update counts
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        #Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        Mmax = c((N-N.I)*(N-N.I-1)/2, (N-N.I) * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1 - R.p1
        M[2] = M[2] + S.p1 + R.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = N.SI + S.p1 - I.p1
      }else{
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1
        M[2] = M[2] + S.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = M[2]
      }
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      N.I = N.I - 1
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        epid[p1] = -1
        N.R = N.R + 1
        M[1] = M[1] + S.p1 + R.p1
        M[2] = M[2] - S.p1 - R.p1 + I.p1
        M[3] = M[3] - I.p1
        #Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        Mmax = c((N-N.I)*(N-N.I-1)/2, (N-N.I) * N.I, N.I*(N.I-1)/2)
        N.SI = N.SI - S.p1
      }else{
        epid[p1] = 0
        N.S = N.S + 1
        M[1] = M[1] + S.p1
        M[2] = M[2] - S.p1 + I.p1
        M[3] = M[3] - I.p1
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        N.SI = M[2]
      }
    }else{
      # some edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]
      
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          C[1] = C[1] + 1
          M[1] = M[1] + 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          C[3] = C[3] + 1
          M[3] = M[3] + 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type: neither is recovered
          C[2] = C[2] + 1
          M[2] = M[2] + 1
          N.SI = N.SI + 1
        }else{
          # R-I type
          C[2] = C[2] + 1
          M[2] = M[2] + 1
        }
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          D[1] = D[1] + 1
          M[1] = M[1] - 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          D[3] = D[3] + 1
          M[3] = M[3] - 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
          N.SI = N.SI - 1
        }else{
          # R-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
        }
      }
      
    }
    # report progress
    # cat("Processing augmented dataset...", 
    #     r/nrow(events),"\r")
    
    t.cur = t.next
    
  }
  #cat("\nDone.\n")
  
  res = list(event.counts = c(n.E, n.R, C, D), big.sums = big.sums, all_M)
  
  return(res)
}

# 12/28/2022 A function to compute the big sums
big.sumer = function(G0, I0, events){
  N = nrow(G0)
  # summary_table = matrix(0, nrow = (nrow(events)+1), ncol = length(unique(events$event)))
  big_event_table = matrix(0, nrow = (nrow(events)+1), ncol = 6)
  # summary_table[1] = 1
  t_pre = 0
  
  epi_status = rep(0, N)
  epi_status[I0] = 1
  SI <- I <- H <- M_HH <- M_HI <- M_II <- rep(0, (nrow(events)+1))
  SI[1] = sum(G0[epi_status == 0, epi_status == 1])
  I[1] = sum(epi_status == 1)
  H[1] = sum(epi_status != 1)
  M_HH[1] = sum(G0[epi_status != 1, epi_status != 1])/2
  M_HI[1] = sum(G0[epi_status != 1, epi_status == 1])
  M_II[1] = sum(G0[epi_status == 1, epi_status == 1])/2
  for (ii in 1: nrow(events)) {
    # ii = 9
    
    # t_cur = events$time[ii]
    # event_pre = summary_table[ii, ]
    # event_pre[events$event[ii]] = event_pre[events$event[ii]] + 1
    # summary_table[(ii+1), ] = event_pre
    
    if (events$event[ii] == 1){
      epi_status[events$per1[[ii]]] = 1
    }
    if (events$event[ii] == 2){
      epi_status[events$per1[[ii]]] = -1
    }
    if (events$event[ii] %in% 3: 5){
      G0[events$per1[[ii]], events$per2[[ii]]] = 1
      G0[events$per2[[ii]], events$per1[[ii]]] = 1
    }
    if (events$event[ii] %in% 6: 8){
      G0[events$per1[[ii]], events$per2[[ii]]] = 0
      G0[events$per2[[ii]], events$per1[[ii]]] = 0
    }
    
    SI[(ii+1)] = sum(G0[epi_status == 0, epi_status == 1])
    I[(ii+1)] = sum(epi_status == 1)
    H[(ii+1)] = sum(epi_status != 1)
    M_HH[(ii+1)] = sum(G0[epi_status != 1, epi_status != 1])/2
    M_HI[(ii+1)] = sum(G0[epi_status != 1, epi_status == 1])
    M_II[(ii+1)] = sum(G0[epi_status == 1, epi_status == 1])/2
    
  }
  # summary_table = as.data.frame(summary_table)
  # summary_table$time = c(0, events$time)
  # summary_table$time_diff = c(0, summary_table$time[-1] - summary_table$time[-nrow(summary_table)])
  
  big_sum_table = cbind.data.frame(beta = SI,
                                   gamma = I,
                                   alpha.SS = H*(H-1)/2 - M_HH,
                                   alpha.SI = H*I - M_HI,
                                   alpha.II = I*(I-1)/2-M_II,
                                   omega.SS =  M_HH,
                                   omega.SI = M_HI,
                                   omega.II = M_II,
                                   time = c(0, events$time))
  
  # big_sum_table$time_diff = c(0, big_sum_table$time[-1] - 
  #                               big_sum_table$time[-nrow(big_sum_table)])
  big_sum_table$time_diff = c(big_sum_table$time[-1] - 
                                big_sum_table$time[-nrow(big_sum_table)], 0)
  big_sum_table
}
# 2. Impute missing recovery times on intervals with missingness
# 2.1 a function to obtain time intervals to do imputation on, and involved individuals
## return a table of interval boundaries, (st, en)'s, and a list of vectors of individual labels
get_miss_recov <- function(report, times, events){
  # data storage
  lbs = NULL; ubs = NULL
  miss_recov = list()
  # go through all intervals and check
  nt = length(times)
  for(ix in 2:nt){
    lb = times[ix-1]; ub = times[ix]
    epid.change = report[ix,] - report[ix-1,]
    recovered = which(epid.change < 0)
    
    cat("Interval:",lb,ub,"\n")
    cat("Recovered:",recovered,"\n")
    
    if(length(recovered) > 0){
      # check if they are included in the event log already
      st = min(which(events$time > lb)); en = max(which(events$time <= ub))
      events.sel = events[st:en,]
      events.sel = events[events$event == 2,]
      if(nrow(events.sel) > 0){
        exact.recovered = events.sel$per1
        recovered = recovered[! recovered %in% exact.recovered]
      }
      # record if 'recovered' non-empty
      if(length(recovered) > 0){
        lbs = c(lbs,lb); ubs = c(ubs, ub)
        miss_recov = append(miss_recov, list(recovered))
      }
    }
  }
  return(list(intervals = data.frame(lb=lbs, ub=ubs), recover = miss_recov))
}


# 2.2 a function to obtain the INFECTED neighbors AT TIME OF INFECTION for each individual
# and return a list for permanent storage
# "Infected": up to the latest status report + new infections since then
# option: only return results for a subset of individuals
get_nei_infection <- function(G0, events, report, times, subset=NULL){
  adjmat = G0
  events = events[events$event != 2,]
  nei_infec = list()
  new_infecs = NULL
  lb = 0
 
  # updat adjmat on the fly
  # and record whenever a POI gets infected
  for(r in 1:nrow(events)){
    z =  events$event[r]
    if(z==1){
      # infection: record neighborhood
      p1 = events$per1[r]
      t = events$time[r]
      # keep track of new infections since latest status report
      lb.t = max(which(times <= t))
      if(lb.t == lb){
        new_infecs = c(new_infecs, p1)
      }else{
        new_infecs = p1
        lb = lb.t
      }
      if(is.null(subset) | (!is.null(subset) & p1 %in% subset)){
        nei = which(adjmat[p1,]==1)
        ix = max(which(times <= t))
        infected.pre = which(report[ix,]==1)
        infected = union(infected.pre,new_infecs)
        nei_infec[[as.character(p1)]] = intersect(nei, infected)
      }
      # if(!is.null(subset) & p1 %in% subset){
      #   nei_infec[[as.character(p1)]] = which(adjmat[p1,]==1)
      # }
      # if(is.null(subset)){
      #   nei_infec[[as.character(p1)]] = which(adjmat[p1,]==1)
      # }
    }else{
      # edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]
      # assume type-dependent edge evolution
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
      }
    }
  }
  return(nei_infec)
}



# 2.3 propose recovery times for a given time interval
# 2.3.a "REJECT": propose times from TE and keep rejecting if not consistent with infection trajectory
propose_recov_rej <- function(lb, ub, recovers, events, nei_infec, gam=0.2){
  
  st = min(which(events$time > lb)); en = max(which(events$time <= ub))
  
  # pull up infection log in this interval 
  events.infec = events[st:en,]
  events.infec = events.infec[events.infec$event == 1,c("time","per1")]
  
  # propose initial candicate times
  cands = sam_texp(length(recovers), gam, 0, ub-lb) + lb
  
  #cat("Interval: [",lb,",",ub,"]\n")
  #cat("To recover:",recovers,"\n")
  #cat("infection log:\n")
  #print(events.infec)
  #cat("Proposed recovery times:\n",cands,"\n")
  
  # a sub-function for consistency-checking
  check_consist <- function(cands){
    for(r in 1:nrow(events.infec)){
      p1 = events.infec$per1[r]
      # check self: can't recover before infection
      if(p1 %in% recovers){
        t = events.infec$time[r]
        if(cands[recovers == p1] <= t){
          return(FALSE)
        }
      }
      # then check neighborhood
      nei = nei_infec[[as.character(p1)]]
      poi = recovers[recovers %in% nei]
      if(length(poi)==length(nei)){
        t = events.infec$time[r]
        poi.t = cands[recovers %in% poi]
        if(all(poi.t <= t)){
          return(FALSE)
        }
      }
    }
    return(TRUE)
  }
  
  # if no infection cases to accommodate at all, directly return the proposal
  if(nrow(events.infec)==0){
    return(cands)
  }
  
  # otherwise, rejection sampling
  while (!check_consist(cands)) {
    #cat("Rejected!\n")
    cands = sam_texp(length(recovers), gam, 0, ub-lb) + lb
    #cat("Proposed recovery times:\n",cands,"\n")
  }
  #cat("Accepted!\n")
  return(cands)
}


# 2.3.b "MH": propose times from TE and keep the previous sample if not consistent with infection trajectory
# return something if it is viable, otherwise return a vector of NA's
propose_recov_MH <- function(lb, ub, recovers, events, nei_infec, gam=0.2){
  
  st = min(which(events$time > lb)); en = max(which(events$time <= ub))
  
  # pull up infection log in this interval 
  events.infec = events[st:en,]
  events.infec = events.infec[events.infec$event == 1,c("time","per1")]
  
  # propose candicate times
  cands = sam_texp(length(recovers), gam, 0, ub-lb) + lb
  
  #cat("Interval: [",lb,",",ub,"]\n")
  #cat("To recover:",recovers,"\n")
  #cat("infection log:\n")
  #print(events.infec)
  #cat("Proposed recovery times:\n",cands,"\n")
  
  # a sub-function for consistency-checking
  check_consist <- function(cands){
    for(r in 1:nrow(events.infec)){
      p1 = events.infec$per1[r]
      # check self: can't recover before infection
      if(p1 %in% recovers){
        t = events.infec$time[r]
        if(cands[recovers == p1] <= t){
          return(FALSE)
        }
      }
      # then check neighborhood
      nei = nei_infec[[as.character(p1)]]
      poi = recovers[recovers %in% nei]
      if(length(poi)==length(nei)){
        t = events.infec$time[r]
        poi.t = cands[recovers %in% poi]
        if(all(poi.t <= t)){
          return(FALSE)
        }
      }
    }
    return(TRUE)
  }
  
  # if no infection cases to accommodate at all, directly return the proposal
  if(nrow(events.infec)==0){
    return(cands)
  }
  
  # otherwise, return cands only if they are consistent with observed data
  if(check_consist(cands)){
    return(cands)
  }else{
    return(rep(NA, length(recovers)))
  }
}

# 2.3.c "CHEWBACCA": filter through each infected person's neighborhood to ensure consistency
propose_recov_filter <- function(lb, ub, recovers, events, nei_infec, gam=0.2){
  
  st = min(which(events$time > lb)); en = max(which(events$time <= ub))
  
  # pull up infection log in this interval 
  events.infec = events[st:en,]
  events.infec = events.infec[events.infec$event == 1,c("time","per1")]
  
  # if events.infec non-empty:
  # obtain adjusted feasible sampling lower bounds for each person in `recovers`
  # if about-to-recover people are the only sick neighbors of someone infected at t, 
  # then randomly select one of them to mandately recover after t
  bounds = rep(lb, length(recovers))
  if(nrow(events.infec) > 0){
    for(r in 1:nrow(events.infec)){
      p1 = events.infec$per1[r]
      if(p1 %in% recovers){
        t = events.infec$time[r]
        bounds[recovers==p1] = max(bounds[recovers==p1],t)
      }
      nei = nei_infec[[as.character(p1)]]
      poi = recovers[recovers %in% nei]
      if(length(poi)==length(nei)){
        #cat("POIs for individual",p1,":\n",poi,"\n")
        t = events.infec$time[r]
        if(length(poi)==1){
          p = poi
        }else{
          p = sample(poi, 1)
        }
        bounds[recovers==p] = max(bounds[recovers==p],t)
      }
    }
  }
  
  # sample recovery times under the adjusted bounds
  cands = sam_texp(length(recovers), gam, bounds-lb, ub-lb) + lb
  
  # cat("Interval: [",lb,",",ub,"]\n")
  # cat("To recover:",recovers,"\n")
  # cat("Feasible lower bounds:\n",bounds,"\n")
  # cat("Proposed recovery times:\n",cands,"\n")
  
  return(cands)
}


count_SItype = function(G0, I0, events){
  epid = rep(0, nrow(G0))
  epid[I0] = 1
  
  SI_summary = cbind.data.frame(sick = I0,
                                sus = which(G0[I0, ] == 1 & epid == 0),
                                st = 0,
                                ed = Inf,
                                type = 1)
  
  
  
  for(i in 1:nrow(events)){
    # i = 91
    event_i = events[i, ]
    if (event_i$event == 1){
      temp_idx = which(SI_summary$sus == event_i$per1 & SI_summary$ed == Inf)
      SI_summary$ed[temp_idx] = sapply(SI_summary$ed[temp_idx], function(x){min(x, event_i$time)})
      
      
      epid[event_i$per1] = 1
      if (length(which(G0[event_i$per1, ] == 1 & epid == 0)) > 0){
        temp_summary = cbind.data.frame(sick = event_i$per1,
                                        sus = which(G0[event_i$per1, ] == 1 & epid == 0),
                                        st = event_i$time,
                                        ed = Inf,
                                        type = 1)
        SI_summary = rbind.data.frame(SI_summary,temp_summary)
      }
      
      
      
      
      
    }
    if (event_i$event == 2){
      temp_idx = which(SI_summary$sick == event_i$per1 & SI_summary$ed == Inf)
      
      SI_summary$ed[temp_idx] = sapply(SI_summary$ed[temp_idx], function(x){min(x, event_i$time)})
      epid[event_i$per1] = -1
    }
    if (event_i$event == 4){
      if (all(sort(epid[c(event_i$per1, event_i$per2)]) == c(0, 1))){
        temp_summary = cbind.data.frame(sick = c(event_i$per1, event_i$per2)[epid[c(event_i$per1, event_i$per2)] == 1],
                                        sus =  c(event_i$per1, event_i$per2)[epid[c(event_i$per1, event_i$per2)] == 0],
                                        st = event_i$time,
                                        ed = Inf,
                                        type = 2)
        SI_summary = rbind.data.frame(SI_summary,temp_summary)
      }
      
      
    }
    if (event_i$event == 7){
      sick = c(event_i$per1, event_i$per2)[epid[c(event_i$per1, event_i$per2)] == 1]
      sus =  c(event_i$per1, event_i$per2)[epid[c(event_i$per1, event_i$per2)] == 0]
      
      temp_idx = which(SI_summary$sus == sus & SI_summary$sick == sick & SI_summary$ed == Inf)
      
      SI_summary$ed[temp_idx] = sapply(SI_summary$ed[temp_idx], function(x){min(x, event_i$time)})
    }
    
    
    if (event_i$event %in% 3: 5){
      G0[event_i$per1, event_i$per2] = 1
      G0[event_i$per2, event_i$per1] = 1
    }
    if (event_i$event %in% 6: 8){
      G0[event_i$per1, event_i$per2] = 0
      G0[event_i$per2, event_i$per1] = 0
    }
    
  }
  SI_duration = unlist(SI_summary$ed) - SI_summary$st
  p_start = table(SI_summary$type) / nrow(SI_summary)
  a = sum(SI_duration[SI_summary$type == 1])
  b = sum(SI_duration[SI_summary$type == 2])
  p_duration = c(a, b) / (a+b)
  
  return(rbind(p_start, p_duration))
}


# ## bench mark the two functions (REJECT v.s. CHEWBACCA)
# bp_res =
# bench::press(
#   ix = c(1:length(recovers)),
#   {
#     bench::mark(
#       length(propose_recov_rej(lb=intervals$lb[ix],ub = intervals$ub[ix], recovers = recovers[[ix]],
#                                   events = miss1$events, nei_infec = nei_infec_miss)),
#       length(propose_recov_filter(lb=intervals$lb[ix],ub = intervals$ub[ix], recovers = recovers[[ix]],
#                            events = miss1$events, nei_infec = nei_infec_miss))
#     )
#   }
# )
# 
# bp_tab = bp_res %>% dplyr::select(ix, min, median)
# bp_tab$num_recov = unlist(bp_res$result)
# bp_tab$method = rep(c("rejection","domain"),5)
# bp_tab$min = as.character(bp_tab$min)
# bp_tab$median = as.character(bp_tab$median)
# library(xtable)
# xtable(bp_tab)
# ## not much difference;
# ## but when #(imputation) is big or contraints are complex,
# ## the 'filter' function is more efficient
