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
rowSort_abs <- function(mat) {
  # mat = health_change_time
  # Get the sorting order based on absolute values
  order_indices <- rowOrder(abs(mat), descending = FALSE, stable = FALSE)
  
  # Use matrix indexing to reorder elements row-wise
  sorted_mat <- mat[cbind(rep(1:nrow(mat), each = ncol(mat)), as.vector(t(order_indices)))]
  
  # Reshape the sorted vector back into a matrix
  sorted_mat <- matrix(sorted_mat, nrow = nrow(mat), byrow = TRUE)
  
  return(sorted_mat)
}
impute_trunc = function(edge_summary, params = params.cur, 
                        infect_times = infect_times, recov_times = recov_times){
  # edge_summary = edges_to_sample_tmp2[edges_to_sample_tmp2$net0 != edges_to_sample_tmp2$net1, ]
  health_t0 = (infect_times[edge_summary$per1] <= edge_summary$t0 & 
                 edge_summary$t0 < recov_times[edge_summary$per1]) + 
    (infect_times[edge_summary$per2] <= edge_summary$t0 & 
       edge_summary$t0 < recov_times[edge_summary$per2])
  
  edge_summary$event_idx = ifelse(as.logical(edge_summary$net0), 6, 3) + health_t0
  
  
  edge_summary_diff = edge_summary[edge_summary$net0 != edge_summary$net1, ]
  
  num_breaks = 5-rowsums(is.na(edge_summary_diff[,c(3, 7:10, 4)]))
  
  row_idx_by_break_num = lapply(1: 5, function(i){which(num_breaks == i)})
  
  new_events = c()
  edge_summary_diff_updated = list()
  for (i in 1: 5){
    # i = 2
    if (length(row_idx_by_break_num[[i]]) > 0){
      edge_summary_diff_i = edge_summary_diff[row_idx_by_break_num[[i]], ]
      edge_summary_diff_i$r1 = -edge_summary_diff_i$r1
      edge_summary_diff_i$r2 = -edge_summary_diff_i$r2
      
      if (i > 1){
        health_change_time = edge_summary_diff_i[,7:10]
        health_change_time = matrix(t(health_change_time)[!is.na(t(health_change_time))], ncol = i-1, byrow = T)
        health_change_time = rowSort_abs(health_change_time)
        
        # break points in each interval
        intvls = cbind(edge_summary_diff_i$t0, abs(health_change_time), edge_summary_diff_i$t1)
        # length of each sub-interval by the break points
        intvls_len = intvls[,-1] - intvls[,-ncol(intvls)]
        
        # identify the parameters
        params_idx = t(apply(cbind(edge_summary_diff_i[,ncol(edge_summary_diff_i)], sign(health_change_time)), 1, cumsum))
        params_grid = matrix(params[as.vector(params_idx)], ncol = ncol(params_idx))
        
        # probabilities for each interval
        lprobs_intvl = log1mexp(params_grid*intvls_len) - params_grid * intvls[,-ncol(intvls)]
        probs_intvl = exp(lprobs_intvl - rowMaxs(lprobs_intvl, value = T))
        probs_intvl = probs_intvl / rowsums(probs_intvl)
        
        probs_intvl = probs_intvl %*% upper.tri(diag(ncol(probs_intvl)), diag = TRUE)
        
        # sample the sub-interval where the event will happen
        intvl_idx_sampled = rowMaxs(runif(nrow(probs_intvl)) < probs_intvl)
        
        tmp_idx1 = cbind(1: length(intvl_idx_sampled), intvl_idx_sampled)
        tmp_idx2 = cbind(1: length(intvl_idx_sampled), intvl_idx_sampled+1)
        
        event_times_sampled = sam_texp(n = nrow(intvls), 
                                       rate = params_grid[tmp_idx1],
                                       a = intvls[tmp_idx1],
                                       b = intvls[tmp_idx2])
        new_events = rbind.data.frame(new_events, 
                                      data.frame(time = event_times_sampled,
                                                 event = params_idx[tmp_idx1],
                                                 per1 = edge_summary_diff_i$per1,
                                                 per2 = edge_summary_diff_i$per2,
                                                 edge_idx = rownames(edge_summary_diff_i)))
        
        edge_summary_diff_i$t0 = event_times_sampled
        edge_summary_diff_i$net0 = as.integer(params_idx[tmp_idx1] %in% 3: 5)
        
        edge_summary_diff_i[, 7:10][edge_summary_diff_i[, 7:10] < event_times_sampled] = NA
        edge_summary_diff_updated[[i]] = edge_summary_diff_i[, -ncol(edge_summary_diff_i)]
        
      }else{
        intvls = cbind(edge_summary_diff_i$t0, edge_summary_diff_i$t1)
        intvls_len = intvls[,-1] - intvls[,-ncol(intvls)]
        
        params_idx = edge_summary_diff_i[,ncol(edge_summary_diff_i)]
        params_grid = matrix(params[params_idx], ncol = 1)
        
        event_times_sampled = sam_texp(n = nrow(intvls), 
                                       rate = c(params_grid),
                                       a = intvls[,1],
                                       b = intvls[,2])
        new_events = rbind.data.frame(new_events, 
                                      data.frame(time = event_times_sampled,
                                                 event = params_idx,
                                                 per1 = edge_summary_diff_i$per1,
                                                 per2 = edge_summary_diff_i$per2,
                                                 edge_idx = rownames(edge_summary_diff_i)))
        
        
        edge_summary_diff_i$t0 = event_times_sampled
        edge_summary_diff_i$net0 = as.integer(params_idx %in% 3: 5)
        edge_summary_diff_updated[[i]] = edge_summary_diff_i[,-ncol(edge_summary_diff_i)]
      }
    }
    
    
  }
  
  return(list(new_events = new_events, 
              edge_summary = rbind.data.frame(do.call('rbind', edge_summary_diff_updated))
  ))
}
impute_free = function(edge_summary, params = params.cur, 
                       infect_times = infect_times, recov_times = recov_times){
  
  # edge_summary = edges_to_sample_tmp2
  health_t0 = (infect_times[edge_summary$per1] <= edge_summary$t0 & 
                 edge_summary$t0 < recov_times[edge_summary$per1]) + 
    (infect_times[edge_summary$per2] <= edge_summary$t0 & 
       edge_summary$t0 < recov_times[edge_summary$per2])
  
  edge_summary$event_idx = ifelse(as.logical(edge_summary$net0), 6, 3) + health_t0
  
  num_breaks = 5-rowsums(is.na(edge_summary[,c(3, 7:10, 4)]))
  
  row_idx_by_break_num = lapply(1: 5, function(i){which(num_breaks == i)})
  
  new_events = c()
  edge_summary_updated = list()
  for (i in 1: 5){
    # i = 2
    if (length(row_idx_by_break_num[[i]]) > 0){
      edge_summary_i = edge_summary[row_idx_by_break_num[[i]], ]
      edge_summary_i$r1 = -edge_summary_i$r1
      edge_summary_i$r2 = -edge_summary_i$r2
      
      if (i > 1){
        health_change_time = edge_summary_i[,7:10]
        health_change_time = matrix(t(health_change_time)[!is.na(t(health_change_time))], ncol = i-1, byrow = T)
        health_change_time = rowSort_abs(health_change_time)
        
        # break points in each interval
        intvls = cbind(edge_summary_i$t0, abs(health_change_time), Inf)
        # length of each sub-interval by the break points
        intvls_len = intvls[,-1] - intvls[,-ncol(intvls)]
        
        # identify the parameters
        params_idx = t(apply(cbind(edge_summary_i[,ncol(edge_summary_i)], sign(health_change_time)), 1, cumsum))
        params_grid = matrix(params[as.vector(params_idx)], ncol = ncol(params_idx))
        
        # probabilities for each interval
        lprobs_intvl = log1mexp(params_grid*intvls_len) - params_grid * intvls[,-ncol(intvls)]
        probs_intvl = exp(lprobs_intvl - rowMaxs(lprobs_intvl, value = T))
        probs_intvl = probs_intvl / rowsums(probs_intvl)
        
        probs_intvl = probs_intvl %*% upper.tri(diag(ncol(probs_intvl)), diag = TRUE)
        
        # sample the sub-interval where the event will happen
        intvl_idx_sampled = rowMaxs(runif(nrow(probs_intvl)) < probs_intvl)
        
        tmp_idx1 = cbind(1: length(intvl_idx_sampled), intvl_idx_sampled)
        tmp_idx2 = cbind(1: length(intvl_idx_sampled), intvl_idx_sampled+1)
        
        event_times_sampled = sam_texp(n = nrow(intvls), 
                                       rate = params_grid[tmp_idx1],
                                       a = intvls[tmp_idx1],
                                       b = intvls[tmp_idx2])
        
        
        happening_idx = which(event_times_sampled < edge_summary_i$t1)
        if (length(happening_idx) > 0){
          
          
          new_events = rbind.data.frame(new_events, 
                                        data.frame(time = event_times_sampled[happening_idx],
                                                   event = params_idx[tmp_idx1][happening_idx],
                                                   per1 = edge_summary_i$per1[happening_idx],
                                                   per2 = edge_summary_i$per2[happening_idx],
                                                   edge_idx = rownames(edge_summary_i[happening_idx,])))
          
          edge_summary_i = edge_summary_i[happening_idx, ,drop = FALSE]
          
          edge_summary_i$t0 = event_times_sampled[happening_idx]
          edge_summary_i$net0 = as.integer(params_idx[tmp_idx1][happening_idx] %in% 3: 5)
          edge_summary_i[, 7:10][edge_summary_i[, 7:10] < event_times_sampled[happening_idx]] = NA
          edge_summary_updated[[i]] = edge_summary_i[,-ncol(edge_summary_i)]
          
          
          
        }else{
          edge_summary_updated[[i]] = c()
        }
        
        
      }else{
        intvls = cbind(edge_summary_i$t0, edge_summary_i$t1)
        intvls_len = intvls[,-1] - intvls[,-ncol(intvls)]
        
        params_idx = edge_summary_i[,ncol(edge_summary_i)]
        params_grid = matrix(params[params_idx], ncol = 1)
        
        event_times_sampled = sam_texp(n = nrow(intvls), 
                                       rate = c(params_grid),
                                       a = intvls[,1],
                                       b = Inf)
        happening_idx = which(event_times_sampled < intvls[,2])
        new_events = rbind.data.frame(new_events, 
                                      data.frame(time = event_times_sampled[happening_idx],
                                                 event = params_idx[happening_idx],
                                                 per1 = edge_summary_i$per1[happening_idx],
                                                 per2 = edge_summary_i$per2[happening_idx],
                                                 edge_idx = rownames(edge_summary_i[happening_idx,])))
        
        edge_summary_i = edge_summary_i[happening_idx, ,drop = FALSE]
        
        edge_summary_i$t0 = event_times_sampled[happening_idx]
        edge_summary_i$net0 = as.integer( params_idx[happening_idx] %in% 3: 5)
        edge_summary_updated[[i]] = edge_summary_i[,-ncol(edge_summary_i)]
      }
    }
    
    
  }
  
  return(list(new_events = new_events, 
              edge_summary = rbind.data.frame(do.call('rbind', edge_summary_updated))
  ))
}

impute_all = function(edge_summary, params = params.cur, 
                      infect_times = infect_times, recov_times = recov_times){
  # edge_summary = edges_to_sample_tmp2
  all_new_events = list()
  k = 1
  while (TRUE) {
    # print(k)
    intvl_end_not_match = which(edge_summary$net0 != edge_summary$net1)
    
    if( length(intvl_end_not_match) > 0){
      tmp_dat = impute_trunc(edge_summary[edge_summary$net0 != edge_summary$net1, ], 
                             params, infect_times, recov_times)
      edge_summary = rbind.data.frame(edge_summary[edge_summary$net0 == edge_summary$net1,], tmp_dat$edge_summary)
      all_new_events = append(all_new_events, list(tmp_dat$new_events))
    }
    
    # BUG: the input in impute_all requires net0 = net1
    tmp_dat = impute_free(edge_summary, params, infect_times, recov_times)
    edge_summary = tmp_dat$edge_summary
    all_new_events = append(all_new_events, list(tmp_dat$new_events))
    
    
    if (nrow(edge_summary) == 0){
      break
    }
    k = k + 1
  }
  return(all_new_events)
}
infer_miss_recov17_cpp = function(dat, priors, obs_time, force_one = TRUE, fix_prop = 0.3,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42, remove_no_neighbour = FALSE){
  
  # dat = survey_dat2; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T;
  # samples = 500; burn = 100; thin = 1; seed = 1; init.params = NULL; timing = T; force_one = TRUE
  # remove_no_neighbour = FALSE; fix_prop = 0.9
  
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
  rownames(edges_to_sample) = 1: nrow(edges_to_sample)
  edges_could_fix = rownames(edges_to_sample)
  
  sick_pers = rownames(intervals)
  sick_edges = list()
  for (i in 1: nrow(intervals)){
    # i = 1
    row_i = intervals[i, ]
    sick_edges_i = rownames(edges_to_sample)[((edges_to_sample$per1 == sick_pers[i]) | 
                                                (edges_to_sample$per2 == sick_pers[i])) & 
                                               ((row_i[1] <= edges_to_sample$t0 & edges_to_sample$t0 <= row_i[2]) | 
                                                  (row_i[1] <= edges_to_sample$t1 & edges_to_sample$t1  <= row_i[2]) )]
    if (length(sick_edges_i) > 0){
      sick_edges[[sick_pers[i]]] = sick_edges_i
    }
  }
  sick_edges_vec = unlist(sick_edges)
  names(sick_edges_vec) = unlist(sapply(names(sick_edges), function(x){rep(x, length(sick_edges[[x]]))}))
  
  
  
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
  
  recover.datA <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(recover.datA) = c('time', 'per1')
  
  imputed_net_events_A = list()
  for(it in 1:S){
    # it = 2
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
      
      intervals_i = intervals[!rownames(intervals) %in% as.character(recover.datA$per1),]
      recovs =  rownames(intervals_i)[intervals_i[, 1] == report.times[ix]] %>% as.integer()
      
      # recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
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
    recover.dat = rbind.data.frame(recover.dat, recover.datA)
    recover.dat = recover.dat[order(recover.dat$time),]
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    # Prepare the edges and their time intervals we would like to sample
    # Pick patient-healthy per pair and note an connection when the healthy person gets sick
    
    edges_to_sample_tmp2 = edges_to_sample
    
    {
      nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
                                         sick_nei = unlist(nei_infec_temp))
      interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
      temp_mat = t(apply(as.matrix(nei_infec_temp2), 1, sort)); rownames(temp_mat) = NULL
      SI_links2= cbind.data.frame(temp_mat, t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
      colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
      temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
      
      if (length(temp_idx) > 0){
        SI_links2 = SI_links2[temp_idx, ]
        SI_links2['temp_idx'] = rownames(SI_links2)
        
        interval_idx2 = interval_idx2[temp_idx]
        
        edge_idx_to_split = (which(rbind(SI_links2[,1:4], edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2))
        
        
        SI_links2 = SI_links2 %>% arrange(t0, per1, per2)
        rownames(SI_links2) = SI_links2$temp_idx
        SI_links2 = SI_links2[,1: 4]
        
        tmp = edges_to_sample_tmp2[edge_idx_to_split, ] 
        tmp['edge_id'] = rownames(tmp)
        tmp = tmp %>% arrange(t0, per1, per2)
        edge_to_split = tmp$edge_id
        edge_to_fix = rownames(edges_to_sample_tmp2)[!rownames(edges_to_sample_tmp2) %in% edge_to_split]
        
        
        edges_to_sample_tmp2 = edges_to_sample_tmp2[-edge_idx_to_split, ]
        
        
        
        # cut the interval at the time the healthy person gets sick
        extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
        extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                             net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
        extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                             t1 = obs_time[interval_idx2+1],
                                                             net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
        extra_edges2 = as.data.frame(extra_edges2)
        colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
        
        rownames(extra_edges2)[seq(1, 2*nrow(SI_links2), 2)] = paste0(edge_to_split, "_1")
        rownames(extra_edges2)[seq(2, 2*nrow(SI_links2), 2)] = paste0(edge_to_split, "_2")
        edges_to_sample_tmp2 = rbind.data.frame(edges_to_sample_tmp2, extra_edges2)
      }else{
        edges_to_sample_tmp2 = edges_to_sample
      }
    }
    
    # Summarize the agent pair, interval, end point network status, 
    infect_times2 = infect_times
    infect_times2[abs(infect_times2 - 0) < 1e-4] = NA
    edges_to_sample_tmp2$i1 = infect_times2[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$i2 = infect_times2[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    
    edges_to_sample_tmp2$i1[(edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$i2[(edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$r1[(edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$r2[(edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1)] = NA
    
    idx_to_impute = (!rownames(edges_to_sample_tmp2) %in% imputed_net_events_A$edge_idx)
    edges_to_sample_tmp2= edges_to_sample_tmp2[idx_to_impute,]
    
    
    # initial health status at the beginning of each interval
    st = Sys.time()
    imputed_net_events_B = do.call("rbind", impute_all(edges_to_sample_tmp2, params.cur, infect_times, recov_times))
    ed = Sys.time()
    ed - st
    
    # print(nrow(edges_to_sample_tmp2))
    # st = Sys.time()
    # imputed_net_events_B = mclapply(1: nrow(edges_to_sample_tmp2), function(i){
    #   # i = nrow(edges_to_sample_tmp2)
    #   temp_edge = edges_to_sample_tmp2[i, ] %>% as.numeric()
    #   interval_imputer2(sick_link = temp_edge[1: 2],
    #                     sampling_intvl = temp_edge[3: 4],
    #                     epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
    #                     net_status = temp_edge[5: 6],
    #                     params.cur, force_one = force_one)
    # }, mc.cores = 1)
    # ed = Sys.time()
    # print(ed - st)
    # names(imputed_net_events_B) = rownames(edges_to_sample_tmp2)
    
    if (it == 1){
      imputed_net_events = imputed_net_events_B
    }else{
      imputed_net_events_A = imputed_net_events_A[!imputed_net_events_A$edge_idx %in% edge_to_split, ] 
      imputed_net_events = rbind.data.frame(imputed_net_events_A, imputed_net_events_B)
    }
    
    events_temp = rbind.data.frame(sick_events, imputed_net_events[, -5])
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[-1, ]
    
    
    # net_event_rows = which(events_temp$event > 2)
    # net_events = events_temp[net_event_rows, ]
    # net_events$event = ifelse(net_events$event < 6, 3, 6) +
    #   (infect_times[net_events$per1] <= net_events$time & net_events$time <= recov_times[net_events$per1]) +
    #   (infect_times[net_events$per2] <= net_events$time & net_events$time <= recov_times[net_events$per2])
    # 
    # events_temp[net_event_rows, ] = net_events
    
    PA = parse_augment_cpp(G0, I0, events_temp, recover.dat)
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at ", it, "-th iter"))
    
    
    edge_ids_to_fix = sample(edge_to_fix, length(edge_to_fix)*fix_prop, replace = FALSE)
    
    patient_id_to_fix = sample(sick_pers, length(sick_pers)*fix_prop, replace = FALSE)
    recover.datA = recover.dat[recover.dat$per1 %in% patient_id_to_fix, ]
    # we only fix edges unrelated to selected patients' recovery
    edge_ids_to_fix = edge_ids_to_fix[!edge_ids_to_fix %in% sick_edges_vec[!names(sick_edges_vec) %in% patient_id_to_fix]] 
    imputed_net_events_A = imputed_net_events[imputed_net_events$edge_idx %in% edge_ids_to_fix, ]
    
    
    
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
infer_miss_recov19_cpp = function(dat, priors, obs_time, force_one = TRUE, fix_prop = 0.3,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42, remove_no_neighbour = FALSE){
  
  # dat = survey_dat2; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T;
  # samples = 500; burn = 100; thin = 1; seed = 1; init.params = NULL; timing = T; force_one = TRUE
  # remove_no_neighbour = FALSE; fix_prop = 0
  P = function(t, net1, net2, alpha, omega){
    # t = infect_times[i] - net_status_time; net1 = 
    mat1 = matrix(c(alpha, omega, 
                    alpha, omega), nrow = 2, byrow = T)
    
    mat2 = matrix(c(omega, -omega, 
                    -alpha, alpha), nrow = 2, byrow = T)
    
    ((mat1[net1,net2] + mat2[net1,net2] * exp(-(alpha+omega)*t)) / (omega + alpha))
    
  }
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
  rownames(edges_to_sample) = 1: nrow(edges_to_sample)
  edges_could_fix = rownames(edges_to_sample)
  
  sick_pers = rownames(intervals)
  sick_edges = list()
  for (i in 1: nrow(intervals)){
    # i = 1
    row_i = intervals[i, ]
    sick_edges_i = rownames(edges_to_sample)[((edges_to_sample$per1 == sick_pers[i]) | 
                                                (edges_to_sample$per2 == sick_pers[i])) & 
                                               ((row_i[1] <= edges_to_sample$t0 & edges_to_sample$t0 <= row_i[2]) | 
                                                  (row_i[1] <= edges_to_sample$t1 & edges_to_sample$t1  <= row_i[2]) )]
    if (length(sick_edges_i) > 0){
      sick_edges[[sick_pers[i]]] = sick_edges_i
    }
  }
  sick_edges_vec = unlist(sick_edges)
  names(sick_edges_vec) = unlist(sapply(names(sick_edges), function(x){rep(x, length(sick_edges[[x]]))}))
  
  
  
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
  
  recover.datA <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(recover.datA) = c('time', 'per1')
  
  imputed_net_events_A = list()
  for(it in 1:S){
    # it = 2
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    if (it == 1){
      nei_infec_temp = lapply(nei_infec, function(x){
        # x = nei_infec[[40]]
        if (length(x) > 1){
          # sample(x, sample(1: length(x), 1), replace = FALSE)
          sample(x, length(x), replace = FALSE)[1]
        } else{
          x
        }
      })
      
      
    }else{
      nei_infec_temp = lapply(as.integer(names(nei_infec)), function(i){
        # for (i in as.integer(names(nei_infec))){
        # i = 6
        x = nei_infec[[as.character(i)]]
        if (length(x) > 1){
          
          x = x[recov_times[x] > infect_times[i]]
          
          
          closest_net = infect_times[i] %/% window_length
          net_status_time = closest_net * window_length
          net_status_i = ifelse(net_snapshots[i, x, closest_net+1] == 0, 2, 1)
          
          all_p = sapply(1: length(net_status_i), function(ii){
            # ii = 1
            if (infect_times[x[ii]] <= net_status_time){
              
              P(infect_times[i] - net_status_time, net_status_i[ii], 1, params.cur[4], params.cur[7])
            }else{
              
              p_0 = P(infect_times[x[ii]] - net_status_time, net_status_i[ii], 2, params.cur[3], params.cur[6])
              p_1 = P(infect_times[x[ii]] - net_status_time, net_status_i[ii], 1, params.cur[3], params.cur[6])
              
              
              p_01 = P(infect_times[i] - infect_times[x[ii]], 2, 1, params.cur[4], params.cur[7])
              p_11 = P(infect_times[i] - infect_times[x[ii]], 1, 1, params.cur[4], params.cur[7])
              
              
              p_0*p_01 + p_1*p_11
            }
          })
          nei_infec_i = c()
          while (length(nei_infec_i) == 0) {
            nei_infec_i = x[runif(length(x)) <= all_p]
          }
          if (length(nei_infec_i) > 1){
            sample(nei_infec_i, size=1)
          }else{
            nei_infec_i
          }
        } else{
          x
        }
      })
      names(nei_infec_temp) = names(nei_infec)
    }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    recover.dat = c()
    for(ix in 1: (length(report.times)-1)){
      # ix = 2
      
      intervals_i = intervals[!rownames(intervals) %in% as.character(recover.datA$per1),]
      recovs =  rownames(intervals_i)[intervals_i[, 1] == report.times[ix]] %>% as.integer()
      
      # recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
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
    recover.dat = rbind.data.frame(recover.dat, recover.datA)
    recover.dat = recover.dat[order(recover.dat$time),]
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    # Prepare the edges and their time intervals we would like to sample
    # Pick patient-healthy per pair and note an connection when the healthy person gets sick
    
    edges_to_sample_tmp2 = edges_to_sample
    
    {
      nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
                                         sick_nei = unlist(nei_infec_temp))
      interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
      temp_mat = t(apply(as.matrix(nei_infec_temp2), 1, sort)); rownames(temp_mat) = NULL
      SI_links2= cbind.data.frame(temp_mat, t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
      colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
      temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
      
      if (length(temp_idx) > 0){
        SI_links2 = SI_links2[temp_idx, ]
        SI_links2['temp_idx'] = rownames(SI_links2)
        
        interval_idx2 = interval_idx2[temp_idx]
        
        edge_idx_to_split = (which(rbind(SI_links2[,1:4], edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2))
        
        
        SI_links2 = SI_links2 %>% arrange(t0, per1, per2)
        rownames(SI_links2) = SI_links2$temp_idx
        SI_links2 = SI_links2[,1: 4]
        
        tmp = edges_to_sample_tmp2[edge_idx_to_split, ] 
        tmp['edge_id'] = rownames(tmp)
        tmp = tmp %>% arrange(t0, per1, per2)
        edge_to_split = tmp$edge_id
        edge_to_fix = rownames(edges_to_sample_tmp2)[!rownames(edges_to_sample_tmp2) %in% edge_to_split]
        
        
        edges_to_sample_tmp2 = edges_to_sample_tmp2[-edge_idx_to_split, ]
        
        # cut the interval at the time the healthy person gets sick
        extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
        extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                             net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
        extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                             t1 = obs_time[interval_idx2+1],
                                                             net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
        extra_edges2 = as.data.frame(extra_edges2)
        colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
        
        rownames(extra_edges2)[seq(1, 2*nrow(SI_links2), 2)] = paste0(edge_to_split, "_1")
        rownames(extra_edges2)[seq(2, 2*nrow(SI_links2), 2)] = paste0(edge_to_split, "_2")
        edges_to_sample_tmp2 = rbind.data.frame(edges_to_sample_tmp2, extra_edges2)
      }else{
        edges_to_sample_tmp2 = edges_to_sample
      }
    }
    
    
    # Summarize the agent pair, interval, end point network status, 
    infect_times2 = infect_times
    infect_times2[abs(infect_times2 - 0) < 1e-4] = NA
    edges_to_sample_tmp2$i1 = infect_times2[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$i2 = infect_times2[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    
    edges_to_sample_tmp2$i1[(edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$i2[(edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$r1[(edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$r2[(edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1)] = NA
    
    idx_to_impute = (!rownames(edges_to_sample_tmp2) %in% imputed_net_events_A$edge_idx)
    edges_to_sample_tmp2= edges_to_sample_tmp2[idx_to_impute,]
    
    
    # initial health status at the beginning of each interval
    st = Sys.time()
    imputed_net_events_B = do.call("rbind", impute_all(edges_to_sample_tmp2, params.cur, infect_times, recov_times))
    ed = Sys.time()
    ed - st
    # print(nrow(edges_to_sample_tmp2))
    # st = Sys.time()
    # imputed_net_events_B = mclapply(1: nrow(edges_to_sample_tmp2), function(i){
    #   # i = nrow(edges_to_sample_tmp2)
    #   temp_edge = edges_to_sample_tmp2[i, ] %>% as.numeric()
    #   interval_imputer2(sick_link = temp_edge[1: 2],
    #                     sampling_intvl = temp_edge[3: 4],
    #                     epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
    #                     net_status = temp_edge[5: 6],
    #                     params.cur, force_one = force_one)
    # }, mc.cores = 1)
    # ed = Sys.time()
    # print(ed - st)
    # names(imputed_net_events_B) = rownames(edges_to_sample_tmp2)
    
    if (it == 1){
      imputed_net_events = imputed_net_events_B
    }else{
      imputed_net_events_A = imputed_net_events_A[!imputed_net_events_A$edge_idx %in% edge_to_split, ] 
      imputed_net_events = rbind.data.frame(imputed_net_events_A, imputed_net_events_B)
    }
    
    events_temp = rbind.data.frame(sick_events, imputed_net_events[, -5])
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[-1, ]
    
    
    # net_event_rows = which(events_temp$event > 2)
    # net_events = events_temp[net_event_rows, ]
    # net_events$event = ifelse(net_events$event < 6, 3, 6) +
    #   (infect_times[net_events$per1] <= net_events$time & net_events$time <= recov_times[net_events$per1]) +
    #   (infect_times[net_events$per2] <= net_events$time & net_events$time <= recov_times[net_events$per2])
    # 
    # events_temp[net_event_rows, ] = net_events
    
    PA = parse_augment_cpp(G0, I0, events_temp, recover.dat)
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at ", it, "-th iter"))
    
    
    edge_ids_to_fix = sample(edge_to_fix, length(edge_to_fix)*fix_prop, replace = FALSE)
    
    patient_id_to_fix = sample(sick_pers, length(sick_pers)*fix_prop, replace = FALSE)
    recover.datA = recover.dat[recover.dat$per1 %in% patient_id_to_fix, ]
    # we only fix edges unrelated to selected patients' recovery
    edge_ids_to_fix = edge_ids_to_fix[!edge_ids_to_fix %in% sick_edges_vec[!names(sick_edges_vec) %in% patient_id_to_fix]] 
    imputed_net_events_A = imputed_net_events[imputed_net_events$edge_idx %in% edge_ids_to_fix, ]
    
    
    
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
infer_miss_recov21_cpp = function(dat, priors, obs_time, force_one = TRUE, fix_prop = 0.3,
                                  init.params = NULL, recov_times = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42, remove_no_neighbour = FALSE){
  # 
  # dat = survey_dat2; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T; recov_times = true_recov_times
  # samples = 500; burn = 100; thin = 1; seed = 1; init.params = survey_dat2$truth; timing = T; force_one = TRUE
  # remove_no_neighbour = FALSE; fix_prop = 0
  P = function(t, net1, net2, alpha, omega){
    # t = infect_times[i] - net_status_time; net1 = 
    mat1 = matrix(c(alpha, omega, 
                    alpha, omega), nrow = 2, byrow = T)
    
    mat2 = matrix(c(omega, -omega, 
                    -alpha, alpha), nrow = 2, byrow = T)
    
    ((mat1[net1,net2] + mat2[net1,net2] * exp(-(alpha+omega)*t)) / (omega + alpha))
    
  }
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
  rownames(edges_to_sample) = 1: nrow(edges_to_sample)
  edges_could_fix = rownames(edges_to_sample)
  
  sick_pers = rownames(intervals)
  sick_edges = list()
  for (i in 1: nrow(intervals)){
    # i = 1
    row_i = intervals[i, ]
    sick_edges_i = rownames(edges_to_sample)[((edges_to_sample$per1 == sick_pers[i]) | 
                                                (edges_to_sample$per2 == sick_pers[i])) & 
                                               ((row_i[1] <= edges_to_sample$t0 & edges_to_sample$t0 <= row_i[2]) | 
                                                  (row_i[1] <= edges_to_sample$t1 & edges_to_sample$t1  <= row_i[2]) )]
    if (length(sick_edges_i) > 0){
      sick_edges[[sick_pers[i]]] = sick_edges_i
    }
  }
  sick_edges_vec = unlist(sick_edges)
  names(sick_edges_vec) = unlist(sapply(names(sick_edges), function(x){rep(x, length(sick_edges[[x]]))}))
  
  
  
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
  
  recover.datA <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(recover.datA) = c('time', 'per1')
  
  imputed_net_events_A = list()
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    if (it == 1){
      nei_infec_temp = lapply(nei_infec, function(x){
        # x = nei_infec[[40]]
        if (length(x) > 1){
          # sample(x, sample(1: length(x), 1), replace = FALSE)
          sample(x, length(x), replace = FALSE)[1]
        } else{
          x
        }
      })
      
      
    }else{
      nei_infec_temp = lapply(as.integer(names(nei_infec)), function(i){
        # for (i in as.integer(names(nei_infec))){
        # i = 6
        x = nei_infec[[as.character(i)]]
        if (length(x) > 1){
          
          x = x[recov_times[x] > infect_times[i]]
          
          
          closest_net = infect_times[i] %/% window_length
          net_status_time = closest_net * window_length
          net_status_i = ifelse(net_snapshots[i, x, closest_net+1] == 0, 2, 1)
          
          all_p = sapply(1: length(net_status_i), function(ii){
            # ii = 1
            if (infect_times[x[ii]] <= net_status_time){
              
              P(infect_times[i] - net_status_time, net_status_i[ii], 1, params.cur[4], params.cur[7])
            }else{
              
              p_0 = P(infect_times[x[ii]] - net_status_time, net_status_i[ii], 2, params.cur[3], params.cur[6])
              p_1 = P(infect_times[x[ii]] - net_status_time, net_status_i[ii], 1, params.cur[3], params.cur[6])
              
              
              p_01 = P(infect_times[i] - infect_times[x[ii]], 2, 1, params.cur[4], params.cur[7])
              p_11 = P(infect_times[i] - infect_times[x[ii]], 1, 1, params.cur[4], params.cur[7])
              
              
              p_0*p_01 + p_1*p_11
            }
          })
          nei_infec_i = c()
          while (length(nei_infec_i) == 0) {
            nei_infec_i = x[runif(length(x)) <= all_p]
          }
          if (length(nei_infec_i) > 1){
            sample(nei_infec_i, size=1)
          }else{
            nei_infec_i
          }
        } else{
          x
        }
      })
      names(nei_infec_temp) = names(nei_infec)
    }
    
    # (1) propose recovery times
    if (is.null(recov_times)){
      gam.cur = params.cur[2]
      recover.dat = c()
      for(ix in 1: (length(report.times)-1)){
        # ix = 2
        
        intervals_i = intervals[!rownames(intervals) %in% as.character(recover.datA$per1),]
        recovs =  rownames(intervals_i)[intervals_i[, 1] == report.times[ix]] %>% as.integer()
        
        # recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
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
      recover.dat = rbind.data.frame(recover.dat, recover.datA)
      recover.dat = recover.dat[order(recover.dat$time),]
      if(verbose){ cat("Recovery times and network events imputation done.\n") }
      recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
      recov_times[recover.dat$per1] = recover.dat$time
    }
    
    # Prepare the edges and their time intervals we would like to sample
    # Pick patient-healthy per pair and note an connection when the healthy person gets sick
    
    edges_to_sample_tmp2 = edges_to_sample
    
    {
      nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
                                         sick_nei = unlist(nei_infec_temp))
      interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
      temp_mat = t(apply(as.matrix(nei_infec_temp2), 1, sort)); rownames(temp_mat) = NULL
      SI_links2= cbind.data.frame(temp_mat, t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
      colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
      temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
      
      if (length(temp_idx) > 0){
        SI_links2 = SI_links2[temp_idx, ]
        SI_links2['temp_idx'] = rownames(SI_links2)
        
        interval_idx2 = interval_idx2[temp_idx]
        
        edge_idx_to_split = (which(rbind(SI_links2[,1:4], edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2))
        
        
        SI_links2 = SI_links2 %>% arrange(t0, per1, per2)
        rownames(SI_links2) = SI_links2$temp_idx
        SI_links2 = SI_links2[,1: 4]
        
        tmp = edges_to_sample_tmp2[edge_idx_to_split, ] 
        tmp['edge_id'] = rownames(tmp)
        tmp = tmp %>% arrange(t0, per1, per2)
        edge_to_split = tmp$edge_id
        edge_to_fix = rownames(edges_to_sample_tmp2)[!rownames(edges_to_sample_tmp2) %in% edge_to_split]
        
        
        edges_to_sample_tmp2 = edges_to_sample_tmp2[-edge_idx_to_split, ]
        
        # cut the interval at the time the healthy person gets sick
        extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
        extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                             net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
        extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                             t1 = obs_time[interval_idx2+1],
                                                             net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
        extra_edges2 = as.data.frame(extra_edges2)
        colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
        
        rownames(extra_edges2)[seq(1, 2*nrow(SI_links2), 2)] = paste0(edge_to_split, "_1")
        rownames(extra_edges2)[seq(2, 2*nrow(SI_links2), 2)] = paste0(edge_to_split, "_2")
        
        extra_edges2$zero_allowed = FALSE
        extra_edges2$zero_allowed[seq(2, 2*nrow(SI_links2), 2)] = TRUE 
        
        edges_to_sample_tmp2$zero_allowed = FALSE
        edges_to_sample_tmp2 = rbind.data.frame(edges_to_sample_tmp2, extra_edges2)
      }else{
        edges_to_sample_tmp2 = edges_to_sample
        edges_to_sample_tmp2$zero_allowed = FALSE
      }
    }
    
    
    # Summarize the agent pair, interval, end point network status, 
    infect_times2 = infect_times
    infect_times2[abs(infect_times2 - 0) < 1e-4] = NA
    edges_to_sample_tmp2$i1 = infect_times2[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$i2 = infect_times2[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    
    edges_to_sample_tmp2$i1[(edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$i2[(edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$r1[(edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$r2[(edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1)] = NA
    
    idx_to_impute = (!rownames(edges_to_sample_tmp2) %in% imputed_net_events_A$edge_idx)
    edges_to_sample_tmp2= edges_to_sample_tmp2[idx_to_impute,]
    
    
    edges_to_sample_tmp2$health_at_t0 = (infect_times[edges_to_sample_tmp2$per1] <= edges_to_sample_tmp2$t0 & 
                                           edges_to_sample_tmp2$t0 < recov_times[edges_to_sample_tmp2$per1]) + 
      (infect_times[edges_to_sample_tmp2$per2] <= edges_to_sample_tmp2$t0 & 
         edges_to_sample_tmp2$t0 < recov_times[edges_to_sample_tmp2$per2])
    
    # initial health status at the beginning of each interval
    
    # st = Sys.time()
    if (!is.null(init.params)) params.cur = c(params.cur[1:2], init.params[-c(1,2)])
    imputed_net_events_B = do.call("rbind", impute_all(edges_to_sample_tmp2, params.cur, infect_times, recov_times))
    
    # ed = Sys.time()
    # ed - st
    # print(nrow(edges_to_sample_tmp2))
    # st = Sys.time()
    # imputed_net_events_B = mclapply(1: nrow(edges_to_sample_tmp2), function(i){
    #   
    # # for (i in 1: nrow(edges_to_sample_tmp2)){
    #   # i = 8
    #   temp_edge = edges_to_sample_tmp2[i, ]
    #   
    #   breakpoints = rbind(temp_edge[7:10], c(1, 1, -1, -1))
    #   breakpoints = breakpoints[, !is.na(breakpoints[1, ]), drop=F]
    #   
    #   breakpoints = breakpoints[, order(as.numeric(breakpoints[1,])), drop=F]
    #   
    #   intvl = c(temp_edge$t0, as.numeric(breakpoints[1, ]), temp_edge$t1)
    #   health_status = cumsum(c(temp_edge$health_at_t0, as.integer(breakpoints[2,])))
    #   
    #   
    #   
    #   alpha = params.cur[3:5][health_status+1]
    #   omega = params.cur[6:8][health_status+1]
    #   
    #   
    #   if (ncol(breakpoints) == 0){
    #     breakpoints = c()
    #   }else{
    #     breakpoints = as.numeric(breakpoints[1,])
    #   }
    #   while (TRUE) {
    #     sims = simulate_ctmc(alpha, omega, breakpoints - temp_edge$t0, temp_edge$t1-temp_edge$t0, temp_edge$net0, temp_edge$net1)
    #     if (temp_edge$net0 != temp_edge$net1 | 
    #         (temp_edge$net0 == temp_edge$net1 & nrow(sims) > 1) | 
    #         (temp_edge$net0 == temp_edge$net1 & grepl("_", rownames(temp_edge)))){
    #       break
    #     }
    #   }
    #   sims$time = sims$time + temp_edge$t0
    #   
    #   sims = cbind.data.frame(sims, per1 = temp_edge$per1, per2 = temp_edge$per2)
    #   sims = sims[-1,]
    #   
    #   
    #   sims$event = health_status[as.integer(cut(sims$time, intvl, 1: (length(intvl)-1)))] + 
    #     ifelse(sims$event == 1, 3, 6)
    #   
    #   sims
    #   
    # }, mc.cores = 1)
    # imputed_net_events_B = do.call("rbind", imputed_net_events_B)
    # # sum(imputed_net_events_B$event == 8)
    # ed = Sys.time()
    # print(ed - st)
    # names(imputed_net_events_B) = rownames(edges_to_sample_tmp2)
    
    if (it == 1){
      imputed_net_events = imputed_net_events_B
    }else{
      imputed_net_events_A = imputed_net_events_A[!imputed_net_events_A$edge_idx %in% edge_to_split, ] 
      imputed_net_events = rbind.data.frame(imputed_net_events_A, imputed_net_events_B)
    }
    
    events_temp = rbind.data.frame(sick_events, imputed_net_events[, -5])
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[-1, ]
    
    
    # net_event_rows = which(events_temp$event > 2)
    # net_events = events_temp[net_event_rows, ]
    # net_events$event = ifelse(net_events$event < 6, 3, 6) +
    #   (infect_times[net_events$per1] <= net_events$time & net_events$time <= recov_times[net_events$per1]) +
    #   (infect_times[net_events$per2] <= net_events$time & net_events$time <= recov_times[net_events$per2])
    # 
    # events_temp[net_event_rows, ] = net_events
    recover.dat = cbind.data.frame(time = recov_times[recov_times < max(recov_times)],
                                   per1 = which(recov_times < max(recov_times)))
    recover.dat = recover.dat[order(recover.dat$time),]
    
    PA = parse_augment_cpp(G0, I0, events_temp, recover.dat)
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at ", it, "-th iter"))
    
    
    edge_ids_to_fix = sample(edge_to_fix, length(edge_to_fix)*fix_prop, replace = FALSE)
    
    patient_id_to_fix = sample(sick_pers, length(sick_pers)*fix_prop, replace = FALSE)
    recover.datA = recover.dat[recover.dat$per1 %in% patient_id_to_fix, ]
    # we only fix edges unrelated to selected patients' recovery
    edge_ids_to_fix = edge_ids_to_fix[!edge_ids_to_fix %in% sick_edges_vec[!names(sick_edges_vec) %in% patient_id_to_fix]] 
    imputed_net_events_A = imputed_net_events[imputed_net_events$edge_idx %in% edge_ids_to_fix, ]
    
    
    
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
infer_miss_recov23_cpp = function(dat, net_status_infect, priors, obs_time, force_one = TRUE, fix_prop = 0.3,
                                  init.params = NULL, recov_times = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42, remove_no_neighbour = FALSE, relax = TRUE){
  # recov_times = true_recov_times
  # dat = survey_dat2; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T
  # samples = 500; burn = 100; thin = 1; seed = 1; init.params = survey_dat2$truth; timing = T; force_one = TRUE
  # remove_no_neighbour = FALSE; fix_prop = 0; timing = TRUE
  # net_status_infect = net_status_infect
  {
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
    interact_report_ind = array(0, dim = c(nrow(G0), nrow(G0), length(interact_report)))
    
    for (t in 1: length(interact_report)) {
      # t = 1
      interact_report_ind[cbind(as.matrix(interact_report[[t]]), t)] = 1
      interact_report_ind[cbind(as.matrix(interact_report[[t]])[,2:1], t)] = 1
    }
    
    infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
    infect_times[sick_events$per1] = sick_events$time
    
    if("truth" %in% names(dat)){
      true.params = dat$truth
    }else{
      true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
    }
    
    # Get recovery intervals for all patients
    {
      intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
      for (i in 1: nrow(sick_events)) {
        per = sick_events$per1[i]
        
        ub = report.times[which(reports[, per] == -1)[1]]
        lb = report.times[which(reports[, per] == -1)[1]-1]
        intervals[i, ] = c(lb, ub)
      }
      rownames(intervals) = sick_events$per1
      colnames(intervals) = c("lb", "ub")
      # Add this to handle those who haven't recovered
      intervals[which(is.na(intervals[,2])), ] = Inf
      intervals[,1] = rowMaxs(cbind(intervals[,1], infect_times[sick_events$per1]),value = TRUE)
      intervals = intervals[order(as.integer(rownames(intervals))), ]
    }
    # Get neighbourhood for all patients
    {
      nei_infec = list()
      for (i in 2: nrow(sick_events)) {
        # i = nrow(sick_events)
        sick_per = sick_events$per1[i]
        sick_time =  sick_events$time[i]
        time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
        
        # Infectious people at time of sick_per's infection
        other_patients = sick_events$per1[sick_events$time < sick_time]
        other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
        
        # The rule to determine if patient j belongs to the infectious neighbourhood of newly infected patient k
        # is to check from the intervals when j is sick to the one when i is sick, i and j have ever interacted
        all_infect_nei = c()
        for (j in 1: length(other_patients)) {
          # j = 34
          infect_nei = other_patients[j]
          
          
          st_idx = (infect_times[infect_nei] %/% window_length + 1)
          ed_idx = (infect_times[sick_per] %/% window_length + 1)
          intvl_idx = st_idx: ed_idx
          
          if (any(interact_report_ind[infect_nei, sick_per, intvl_idx] == 1) |
              any(net_snapshots[infect_nei, sick_per, c(intvl_idx, ed_idx+1)] == 1)){
            all_infect_nei = c(all_infect_nei, infect_nei)
          }
          
        }
        nei_infec[[i-1]] = all_infect_nei
        # list of persons reporting their interaction with sick_per
        # interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
        
        # Identify the set of people who might be the virus transmitter
        # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
        # and those who ever connected with sick_per and infectious
        # intersect(interactd_per[interactd_per != sick_per], other_patients)
        # nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
        #                          interactd_per[interactd_per %in% other_patients]) 
        # print(identical(sort(intersect(interactd_per[interactd_per != sick_per], other_patients)), sort(nei_infec[[i-1]])))
      }
      names(nei_infec) = sick_events$per1[-1]
      
      if (remove_no_neighbour){
        to_remove_idx = names(nei_infec)[sapply(nei_infec, length) == 0]
        nei_infec = nei_infec[!names(nei_infec) %in% to_remove_idx]
        intervals=intervals[!rownames(intervals) %in% to_remove_idx, ]
        sick_events = sick_events[!sick_events$per1 %in% as.integer(to_remove_idx), ]
      }
    }
    
    tmp = cbind(1: nrow(G0), infect_times)
    tmp = tmp[tmp[,1] %in% sick_events$per1, ]
    sick_order = tmp[order(tmp[,2]), 1]
    nei_infec = list()
    for (ii in 2: length(sick_order)) {
      nei_infec[[as.character(sick_order[ii])]] = sick_order[1: (ii-1)]
    }
    ik_intvl_idx = cbind(infect_times %/% window_length + 1,
                         infect_times %/% window_length + 2)
    
    non_SI_interact_report_ind = interact_report_ind
    
    for (i in as.integer(names(nei_infec))) {
      # i = 111
      nei_infec_i = nei_infec[[as.character(i)]]
      for (j in nei_infec_i) {
        non_SI_interact_report_ind[rbind(cbind(i, j, ik_intvl_idx[j]: ik_intvl_idx[i]),
                                         cbind(j, i, ik_intvl_idx[j]: ik_intvl_idx[i]))] = 0
      }
    }
    
    
    edges_to_sample = which(non_SI_interact_report_ind == 1, arr.ind = TRUE)
    edges_to_sample = cbind.data.frame(edges_to_sample, edges_to_sample[,3]+1)
    edges_to_sample = edges_to_sample[edges_to_sample[,1] < edges_to_sample[,2],]
    colnames(edges_to_sample) = c('per1', 'per2', 't0', 't1')
    edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
                                                   net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
    edges_to_sample$t0 = obs_time[edges_to_sample$t0]
    edges_to_sample$t1 = obs_time[edges_to_sample$t1]
    
    edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
    rownames(edges_to_sample) = 1: nrow(edges_to_sample)
    # edges_could_fix = rownames(edges_to_sample)
    
    
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
    
    recover.datA <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(recover.datA) = c('time', 'per1')
  }
  
  # imputed_net_events_A = list()
  for(it in 1:S){
    # it = 4
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    ncount = 0
    while (TRUE){
      recov_times_tmp = sam_texp(nrow(intervals), params.cur[2], 0, intervals[,2]-intervals[,1]) + intervals[,1]
      recov_times_tmp[is.nan(recov_times_tmp)] = Inf
      if (all(sapply(names(net_status_infect)[-1], function(sick_per){
        # sick_per = '106'
        conn_nei = net_status_infect[[sick_per]]
        
        any(infect_times[conn_nei] < infect_times[as.integer(sick_per)] & 
            recov_times_tmp[conn_nei] > infect_times[as.integer(sick_per)])
      }))){
        recov_times = rep(Inf, nrow(G0))
        recov_times[as.integer(rownames(intervals))] = as.vector(recov_times_tmp)
        print(paste0('sampling recovery times failed ', ncount, ' times'))
        break
      }else{
        ncount = ncount + 1
      }
    }
    # print('recovery times for: ')
    # print(recov_times[c(20, 106)])
    SI_links = SI_events_imputer(sick_order,
                                 nei_infec,
                                 recov_times,
                                 infect_times,
                                 net_status_infect,
                                 window_length,
                                 obs_time,
                                 net_snapshots,
                                 interact_report_ind,
                                 params.cur, relax=relax)
    
    edges_to_sample_tmp2 = edges_to_sample
    # Summarize the agent pair, interval, end point network status,
    infect_times2 = infect_times
    infect_times2[abs(infect_times2 - 0) < 1e-4] = NA
    edges_to_sample_tmp2$i1 = infect_times2[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$i2 = infect_times2[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    
    edges_to_sample_tmp2$i1[(edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$i2[(edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$r1[(edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$r2[(edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1)] = NA
    
    # idx_to_impute = (!rownames(edges_to_sample_tmp2) %in% imputed_net_events_A$edge_idx)
    # edges_to_sample_tmp2= edges_to_sample_tmp2[idx_to_impute,]
    
    
    edges_to_sample_tmp2$health_at_t0 = (infect_times[edges_to_sample_tmp2$per1] <= edges_to_sample_tmp2$t0 &
                                           edges_to_sample_tmp2$t0 < recov_times[edges_to_sample_tmp2$per1]) +
      (infect_times[edges_to_sample_tmp2$per2] <= edges_to_sample_tmp2$t0 &
         edges_to_sample_tmp2$t0 < recov_times[edges_to_sample_tmp2$per2])
    
    # # initial health status at the beginning of each interval
    # 
    # # st = Sys.time()
    # if (!is.null(init.params)) params.cur = c(params.cur[1:2], init.params[-c(1,2)])
    imputed_net_events_B = do.call("rbind", impute_all(edges_to_sample_tmp2, params.cur, infect_times, recov_times))
    imputed_net_events = rbind.data.frame(SI_links, imputed_net_events_B[,1:4])
    
    
    
    events_temp = rbind.data.frame(sick_events, imputed_net_events)
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[-1, ]
    
    # 
    # net_event_rows = which(events_temp$event > 2)
    # net_events = events_temp[net_event_rows, ]
    # net_events$event = ifelse(net_events$event < 6, 3, 6) +
    #   (infect_times[net_events$per1] <= net_events$time & net_events$time <= recov_times[net_events$per1]) +
    #   (infect_times[net_events$per2] <= net_events$time & net_events$time <= recov_times[net_events$per2])
    # events_temp[net_event_rows, ] = net_events
    recover.dat = cbind.data.frame(time = recov_times[recov_times < max(recov_times)],
                                   per1 = which(recov_times < max(recov_times)))
    recover.dat = recover.dat[order(recover.dat$time),]
    # # 
    # table(events_temp$event)
    PA = parse_augment_cpp(G0, I0, events_temp, recover.dat)
    print(PA)
    # {
    #   net_events_temp = events_temp[events_temp$event %in% 3: 8,]
    #   G0 = dats$G0
    #   for (tt in 1: nrow(net_events_temp)) {
    #     # tt = 4585
    #     if (tt %% 100 == 0) print(tt)
    #     event_tt = net_events_temp[tt, ]
    # 
    #     G0[event_tt$per1, event_tt$per2] <- G0[event_tt$per2, event_tt$per1] <- G0[event_tt$per2, event_tt$per1] + ifelse(event_tt$event%in% 3:5, 1, -1)
    #     if (!G0[event_tt$per2, event_tt$per1] %in% 0: 1){
    #       print('bad')
    #       break
    #     }
    #   }
    # }
    
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at ", it, "-th iter"))
    
    
    # edge_ids_to_fix = sample(edge_to_fix, length(edge_to_fix)*fix_prop, replace = FALSE)
    
    # patient_id_to_fix = sample(sick_pers, length(sick_pers)*fix_prop, replace = FALSE)
    # recover.datA = recover.dat[recover.dat$per1 %in% patient_id_to_fix, ]
    # we only fix edges unrelated to selected patients' recovery
    # edge_ids_to_fix = edge_ids_to_fix[!edge_ids_to_fix %in% sick_edges_vec[!names(sick_edges_vec) %in% patient_id_to_fix]] 
    # imputed_net_events_A = imputed_net_events[imputed_net_events$edge_idx %in% edge_ids_to_fix, ]
    
    
    
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
# A faster version of infer_miss_recov23_cpp
infer_miss_recov24_cpp = function(dat, net_status_infect, priors, obs_time, force_one = TRUE, fix_prop = 0.3,
                                  init.params = NULL, recov_times = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42, remove_no_neighbour = FALSE){
  
  # recov_times = true_recov_times
  # set.seed(1)
  # dat = survey_dat2; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T
  # samples = 500; burn = 100; thin = 1; seed = 1; init.params = survey_dat2$truth; timing = T; force_one = TRUE
  # remove_no_neighbour = FALSE; fix_prop = 0; timing = TRUE; relax= TRUE
  # net_status_infect = net_status_infect
  {
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
    interact_report_ind = array(0, dim = c(nrow(G0), nrow(G0), length(interact_report)))
    
    for (t in 1: length(interact_report)) {
      # t = 1
      interact_report_ind[cbind(as.matrix(interact_report[[t]]), t)] = 1
      interact_report_ind[cbind(as.matrix(interact_report[[t]])[,2:1], t)] = 1
    }
    
    infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
    infect_times[sick_events$per1] = sick_events$time
    
    if("truth" %in% names(dat)){
      true.params = dat$truth
    }else{
      true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
    }
    
    # Get recovery intervals for all patients
    {
      intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
      for (i in 1: nrow(sick_events)) {
        per = sick_events$per1[i]
        
        ub = report.times[which(reports[, per] == -1)[1]]
        lb = report.times[which(reports[, per] == -1)[1]-1]
        intervals[i, ] = c(lb, ub)
      }
      rownames(intervals) = sick_events$per1
      colnames(intervals) = c("lb", "ub")
      # Add this to handle those who haven't recovered
      intervals[which(is.na(intervals[,2])), ] = Inf
      intervals[,1] = rowMaxs(cbind(intervals[,1], infect_times[sick_events$per1]),value = TRUE)
      intervals = intervals[order(as.integer(rownames(intervals))), ]
    }
    # Get neighbourhood for all patients
    {
      nei_infec = list()
      for (i in 2: nrow(sick_events)) {
        # i = nrow(sick_events)
        sick_per = sick_events$per1[i]
        sick_time =  sick_events$time[i]
        time_idx = cut(sick_time, obs_time, labels = 1: (length(obs_time)-1))
        
        # Infectious people at time of sick_per's infection
        other_patients = sick_events$per1[sick_events$time < sick_time]
        other_patients = other_patients[intervals[as.character(other_patients), "ub"] > sick_time]
        
        # The rule to determine if patient j belongs to the infectious neighbourhood of newly infected patient k
        # is to check from the intervals when j is sick to the one when i is sick, i and j have ever interacted
        all_infect_nei = c()
        for (j in 1: length(other_patients)) {
          # j = 34
          infect_nei = other_patients[j]
          
          
          st_idx = (infect_times[infect_nei] %/% window_length + 1)
          ed_idx = (infect_times[sick_per] %/% window_length + 1)
          intvl_idx = st_idx: ed_idx
          
          if (any(interact_report_ind[infect_nei, sick_per, intvl_idx] == 1) |
              any(net_snapshots[infect_nei, sick_per, c(intvl_idx, ed_idx+1)] == 1)){
            all_infect_nei = c(all_infect_nei, infect_nei)
          }
          
        }
        nei_infec[[i-1]] = all_infect_nei
        # list of persons reporting their interaction with sick_per
        # interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
        
        # Identify the set of people who might be the virus transmitter
        # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0, 
        # and those who ever connected with sick_per and infectious
        # intersect(interactd_per[interactd_per != sick_per], other_patients)
        # nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
        #                          interactd_per[interactd_per %in% other_patients]) 
        # print(identical(sort(intersect(interactd_per[interactd_per != sick_per], other_patients)), sort(nei_infec[[i-1]])))
      }
      names(nei_infec) = sick_events$per1[-1]
      
      if (remove_no_neighbour){
        to_remove_idx = names(nei_infec)[sapply(nei_infec, length) == 0]
        nei_infec = nei_infec[!names(nei_infec) %in% to_remove_idx]
        intervals=intervals[!rownames(intervals) %in% to_remove_idx, ]
        sick_events = sick_events[!sick_events$per1 %in% as.integer(to_remove_idx), ]
      }
    }
    
    tmp = cbind(1: nrow(G0), infect_times)
    tmp = tmp[tmp[,1] %in% sick_events$per1, ]
    sick_order = tmp[order(tmp[,2]), 1]
    nei_infec = list()
    for (ii in 2: length(sick_order)) {
      nei_infec[[as.character(sick_order[ii])]] = sick_order[1: (ii-1)]
    }
    ik_intvl_idx = cbind(infect_times %/% window_length + 1,
                         infect_times %/% window_length + 2)
    
    non_SI_interact_report_ind = interact_report_ind
    
    for (i in as.integer(names(nei_infec))) {
      # i = 111
      nei_infec_i = nei_infec[[as.character(i)]]
      for (j in nei_infec_i) {
        non_SI_interact_report_ind[rbind(cbind(i, j, ik_intvl_idx[j]: ik_intvl_idx[i]),
                                         cbind(j, i, ik_intvl_idx[j]: ik_intvl_idx[i]))] = 0
      }
    }
    
    
    edges_to_sample = which(non_SI_interact_report_ind == 1, arr.ind = TRUE)
    edges_to_sample = cbind.data.frame(edges_to_sample, edges_to_sample[,3]+1)
    edges_to_sample = edges_to_sample[edges_to_sample[,1] < edges_to_sample[,2],]
    colnames(edges_to_sample) = c('per1', 'per2', 't0', 't1')
    edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])], 
                                                   net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
    edges_to_sample$t0 = obs_time[edges_to_sample$t0]
    edges_to_sample$t1 = obs_time[edges_to_sample$t1]
    
    edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
    if (nrow(edges_to_sample) > 0) rownames(edges_to_sample) = 1: nrow(edges_to_sample)
    
    # edges_could_fix = rownames(edges_to_sample)
    
    
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
    
    recover.datA <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(recover.datA) = c('time', 'per1')
  }
  
  # imputed_net_events_A = list()
  for(it in 1:S){
    # it = 1
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    ncount = 0
    while (TRUE){
      recov_times_tmp = rep(Inf, nrow(G0))
      recov_times_tmp[as.integer(rownames(intervals))] = sam_texp(nrow(intervals), params.cur[2], 0, intervals[,2]-intervals[,1]) + intervals[,1]
      # recov_times_tmp = sam_texp(nrow(intervals), params.cur[2], 0, intervals[,2]-intervals[,1]) + intervals[,1]
      recov_times_tmp[is.nan(recov_times_tmp)] = Inf
      if (all(sapply(names(net_status_infect)[-1], function(sick_per){
        # sick_per = '28'
        conn_nei = net_status_infect[[sick_per]]
        
        any(infect_times[conn_nei] < infect_times[as.integer(sick_per)] & 
            recov_times_tmp[conn_nei] > infect_times[as.integer(sick_per)])
      }))){
        # recov_times = rep(Inf, nrow(G0))
        # recov_times[as.integer(rownames(intervals))] = as.vector(recov_times_tmp)
        recov_times = as.vector(recov_times_tmp)
        print(paste0('sampling recovery times failed ', ncount, ' times'))
        break
      }else{
        ncount = ncount + 1
      }
    }
    # print('recovery times for: ')
    # print(recov_times[c(20, 106)])
    # set.seed(1)
    # st = Sys.time()
    # SI_links = SI_events_imputer(sick_order,
    #                              nei_infec,
    #                              recov_times,
    #                              infect_times,
    #                              net_status_infect,
    #                              window_length,
    #                              obs_time,
    #                              net_snapshots,
    #                              interact_report_ind,
    #                              params.cur, relax=relax)
    # 
    # ed = Sys.time()
    # ed - st
    
    st = Sys.time()
    SI_links = SI_events_imputer2(sick_order,
                                  nei_infec,
                                  recov_times,
                                  infect_times,
                                  net_status_infect,
                                  window_length,
                                  obs_time,
                                  net_snapshots,
                                  interact_report_ind,
                                  params.cur)
    ed = Sys.time()
    ed - st
    
    edges_to_sample_tmp2 = edges_to_sample
    # Summarize the agent pair, interval, end point network status,
    infect_times2 = infect_times
    infect_times2[abs(infect_times2 - 0) < 1e-4] = NA
    edges_to_sample_tmp2$i1 = infect_times2[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$i2 = infect_times2[edges_to_sample_tmp2$per2]
    edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    
    edges_to_sample_tmp2$i1[(edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$i2[(edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$r1[(edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1)] = NA
    
    edges_to_sample_tmp2$r2[(edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                               edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1)] = NA
    
    # idx_to_impute = (!rownames(edges_to_sample_tmp2) %in% imputed_net_events_A$edge_idx)
    # edges_to_sample_tmp2= edges_to_sample_tmp2[idx_to_impute,]
    
    
    edges_to_sample_tmp2$health_at_t0 = (infect_times[edges_to_sample_tmp2$per1] <= edges_to_sample_tmp2$t0 &
                                           edges_to_sample_tmp2$t0 < recov_times[edges_to_sample_tmp2$per1]) +
      (infect_times[edges_to_sample_tmp2$per2] <= edges_to_sample_tmp2$t0 &
         edges_to_sample_tmp2$t0 < recov_times[edges_to_sample_tmp2$per2])
    
    # # initial health status at the beginning of each interval
    # 
    # # st = Sys.time()
    # if (!is.null(init.params)) params.cur = c(params.cur[1:2], init.params[-c(1,2)])
    imputed_net_events_B = do.call("rbind", impute_all(edges_to_sample_tmp2, params.cur, infect_times, recov_times))
    imputed_net_events = rbind.data.frame(SI_links, imputed_net_events_B[,1:4])
    
    
    
    events_temp = rbind.data.frame(sick_events, imputed_net_events)
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[-1, ]
    
    # 
    # net_event_rows = which(events_temp$event > 2)
    # net_events = events_temp[net_event_rows, ]
    # net_events$event = ifelse(net_events$event < 6, 3, 6) +
    #   (infect_times[net_events$per1] <= net_events$time & net_events$time <= recov_times[net_events$per1]) +
    #   (infect_times[net_events$per2] <= net_events$time & net_events$time <= recov_times[net_events$per2])
    # events_temp[net_event_rows, ] = net_events
    recover.dat = cbind.data.frame(time = recov_times[recov_times < max(recov_times)],
                                   per1 = which(recov_times < max(recov_times)))
    recover.dat = recover.dat[order(recover.dat$time),]
    # # 
    # table(events_temp$event)
    PA = parse_augment_cpp(G0, I0, events_temp, recover.dat)
    print(PA)
    # {
    #   net_events_temp = events_temp[events_temp$event %in% 3: 8,]
    #   G0 = dats$G0
    #   for (tt in 1: nrow(net_events_temp)) {
    #     # tt = 4585
    #     if (tt %% 100 == 0) print(tt)
    #     event_tt = net_events_temp[tt, ]
    # 
    #     G0[event_tt$per1, event_tt$per2] <- G0[event_tt$per2, event_tt$per1] <- G0[event_tt$per2, event_tt$per1] + ifelse(event_tt$event%in% 3:5, 1, -1)
    #     if (!G0[event_tt$per2, event_tt$per1] %in% 0: 1){
    #       print('bad')
    #       break
    #     }
    #   }
    # }
    
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at ", it, "-th iter"))
    
    
    # edge_ids_to_fix = sample(edge_to_fix, length(edge_to_fix)*fix_prop, replace = FALSE)
    
    # patient_id_to_fix = sample(sick_pers, length(sick_pers)*fix_prop, replace = FALSE)
    # recover.datA = recover.dat[recover.dat$per1 %in% patient_id_to_fix, ]
    # we only fix edges unrelated to selected patients' recovery
    # edge_ids_to_fix = edge_ids_to_fix[!edge_ids_to_fix %in% sick_edges_vec[!names(sick_edges_vec) %in% patient_id_to_fix]] 
    # imputed_net_events_A = imputed_net_events[imputed_net_events$edge_idx %in% edge_ids_to_fix, ]
    
    
    
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
infer_miss_recov25_cpp = function(dat, net_status_infect, priors, obs_time, force_one = TRUE, fix_prop = 0.3,
                                  init.params = NULL, recov_times = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42, remove_no_neighbour = FALSE){
  process_network_edges <- function(sampled_edges, infection_times, recovery_times) {
    # sampled_edges = edges_to_sample_HI[38,];
    # infection_times = infect_times
    # recovery_times = recov_times
    # Create a temporary copy
    temp_edges <- sampled_edges
    
    # Summarize the agent pair, interval, and endpoint network status
    adjusted_infect_times <- infection_times
    adjusted_infect_times[abs(adjusted_infect_times - 0) < 1e-4] <- NA
    
    temp_edges$i1 <- adjusted_infect_times[temp_edges$per1]
    temp_edges$i2 <- adjusted_infect_times[temp_edges$per2]
    temp_edges$r1 <- recovery_times[temp_edges$per1]
    temp_edges$r2 <- recovery_times[temp_edges$per2]
    
    # Adjust values based on time constraints
    temp_edges$i1[(temp_edges$i1 < temp_edges$t0 |
                     temp_edges$i1 > temp_edges$t1)] <- NA
    
    temp_edges$i2[(temp_edges$i2 < temp_edges$t0 |
                     temp_edges$i2 > temp_edges$t1)] <- NA
    
    temp_edges$r1[(temp_edges$r1 < temp_edges$t0 |
                     temp_edges$r1 > temp_edges$t1)] <- NA
    
    temp_edges$r2[(temp_edges$r2 < temp_edges$t0 |
                     temp_edges$r2 > temp_edges$t1)] <- NA
    
    # Compute health status at t0
    temp_edges$health_at_t0 <- (infection_times[temp_edges$per1] <= temp_edges$t0 &
                                  temp_edges$t0 < recovery_times[temp_edges$per1]) + 
      (infection_times[temp_edges$per2] <= temp_edges$t0 &
         temp_edges$t0 < recovery_times[temp_edges$per2])
    
    return(temp_edges)
  }
  
  
  extract_filtered_edges <- function(interaction_matrix, network_data, time_vector) {
    
    
    # Find edges where interaction report indicator is 1
    
    filtered_edges <- which(interaction_matrix == 1, arr.ind = TRUE)
    
    # Add a new column for time shift
    filtered_edges <- cbind.data.frame(filtered_edges, filtered_edges[,3] + 1)
    
    # Keep only edges where node1 < node2 to avoid duplicates
    filtered_edges <- filtered_edges[filtered_edges[,1] < filtered_edges[,2],]
    
    # Assign column names
    colnames(filtered_edges) <- c('per1', 'per2', 't0', 't1')
    
    # Extract network snapshots for t0 and t1
    filtered_edges <- cbind(filtered_edges, 
                            cbind(net0 = network_data[as.matrix(filtered_edges[, -4])], 
                                  net1 = network_data[as.matrix(filtered_edges[, -3])]))
    
    # Map time indices to actual observation times
    filtered_edges$t0 <- time_vector[filtered_edges$t0]
    filtered_edges$t1 <- time_vector[filtered_edges$t1]
    
    # Remove duplicate rows
    filtered_edges <- filtered_edges[!duplicated(filtered_edges), ]
    
    # Reset row names if there are rows
    if (nrow(filtered_edges) > 0) rownames(filtered_edges) <- 1:nrow(filtered_edges)
    
    return(filtered_edges)
  }
  # 
  # recov_times = true_recov_times
  # set.seed(1)
  # dat = survey_dat2; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T
  # samples = 500; burn = 100; thin = 1; seed = 1; init.params = survey_dat2$truth; timing = T; force_one = TRUE
  # remove_no_neighbour = FALSE; fix_prop = 0; timing = TRUE; relax= TRUE; verbose = TRUE
  # net_status_infect = net_status_infect
  {
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
    interact_report_ind = array(0, dim = c(nrow(G0), nrow(G0), length(interact_report)))
    healthy_agents = c(1: nrow(G0))[-sick_events$per1]
    for (t in 1: length(interact_report)) {
      # t = 1
      interact_report_ind[cbind(as.matrix(interact_report[[t]]), t)] = 1
      interact_report_ind[cbind(as.matrix(interact_report[[t]])[,2:1], t)] = 1
    }
    
    infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
    infect_times[sick_events$per1] = sick_events$time
    
    if("truth" %in% names(dat)){
      true.params = dat$truth
    }else{
      true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
    }
    
    # Get recovery intervals for all patients
    {
      intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
      for (i in 1: nrow(sick_events)) {
        # i = 1
        per = sick_events$per1[i]
        
        ub = report.times[which(reports[, per] == -1)[1]]
        lb = report.times[which(reports[, per] == -1)[1]-1]
        intervals[i, ] = c(lb, ub)
      }
      rownames(intervals) = sick_events$per1
      colnames(intervals) = c("lb", "ub")
      # Add this to handle those who haven't recovered
      intervals = intervals[rowSums(!is.na(intervals)) == 2, ]
      intervals[,1] = rowMaxs(cbind(intervals[,1], infect_times[as.integer(rownames(intervals))]),value = TRUE)
      intervals = intervals[order(as.integer(rownames(intervals))), ]
    }
    intervals[intervals[,2] > obs_time[length(obs_time)], 2] = obs_time[length(obs_time)]
    
    # Get neighbourhood for all patients
    # {
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
    #     # The rule to determine if patient j belongs to the infectious neighbourhood of newly infected patient k
    #     # is to check from the intervals when j is sick to the one when i is sick, i and j have ever interacted
    #     all_infect_nei = c()
    #     for (j in 1: length(other_patients)) {
    #       # j = 34
    #       infect_nei = other_patients[j]
    # 
    # 
    #       st_idx = (infect_times[infect_nei] %/% window_length + 1)
    #       ed_idx = (infect_times[sick_per] %/% window_length + 1)
    #       intvl_idx = st_idx: ed_idx
    # 
    #       if (any(interact_report_ind[infect_nei, sick_per, intvl_idx] == 1) |
    #           any(net_snapshots[infect_nei, sick_per, c(intvl_idx, ed_idx+1)] == 1)){
    #         all_infect_nei = c(all_infect_nei, infect_nei)
    #       }
    # 
    #     }
    #     nei_infec[[i-1]] = all_infect_nei
    #     # list of persons reporting their interaction with sick_per
    #     # interactd_per = interact_report[[time_idx]][rowSums(interact_report[[time_idx]] == sick_per) > 0, ] %>% unlist()
    # 
    #     # Identify the set of people who might be the virus transmitter
    #     # They're the infectious people at time of sick_per's infection who are connected with sick_per at t0,
    #     # and those who ever connected with sick_per and infectious
    #     # intersect(interactd_per[interactd_per != sick_per], other_patients)
    #     # nei_infec[[i-1]] = union(other_patients[as.logical(net_snapshots[, , time_idx][cbind(sick_per, other_patients)])],
    #     #                          interactd_per[interactd_per %in% other_patients])
    #     # print(identical(sort(intersect(interactd_per[interactd_per != sick_per], other_patients)), sort(nei_infec[[i-1]])))
    #   }
    #   names(nei_infec) = sick_events$per1[-1]
    # 
    #   if (remove_no_neighbour){
    #     to_remove_idx = names(nei_infec)[sapply(nei_infec, length) == 0]
    #     nei_infec = nei_infec[!names(nei_infec) %in% to_remove_idx]
    #     intervals=intervals[!rownames(intervals) %in% to_remove_idx, ]
    #     sick_events = sick_events[!sick_events$per1 %in% as.integer(to_remove_idx), ]
    #   }
    # }
    
    tmp = cbind(1: nrow(G0), infect_times)
    tmp = tmp[tmp[,1] %in% sick_events$per1, ]
    
    sick_order = tmp[order(tmp[,2]), 1]
    nei_infec = list()
    for (ii in 2: length(sick_order)) {
      nei_infec[[as.character(sick_order[ii])]] = sick_order[1: (ii-1)]
    }
    ik_intvl_idx = cbind(infect_times %/% window_length + 1,
                         infect_times %/% window_length + 2)
    
    
    non_SI_interact_report_ind = interact_report_ind
    
    for (i in as.integer(names(nei_infec))) {
      # i = 111
      nei_infec_i = nei_infec[[as.character(i)]]
      for (j in nei_infec_i) {
        non_SI_interact_report_ind[rbind(cbind(i, j, ik_intvl_idx[j]: ik_intvl_idx[i]),
                                         cbind(j, i, ik_intvl_idx[j]: ik_intvl_idx[i]))] = 0
      }
    }
    
    edges_to_sample_HI_HH = extract_filtered_edges(non_SI_interact_report_ind, net_snapshots, obs_time)
    
    report.times[report.times >= obs_time[length(obs_time)]] = obs_time[length(obs_time)]-1e-5
    sick_range = cbind(infect_times[sick_order] %/% window_length + 1,
                       apply(reports[, sick_order], 2, function(x){
                         if (-1 %in% x){
                           report.times[which(x == -1)[1]]
                         }else{
                           report.times[length(x)]
                         }
                       })
    )
    rownames(sick_range) = sick_order
    sick_range[,2] = sick_range[,2] %/% window_length + 1
    
    tmp_idx_set = do.call('rbind', lapply(1: nrow(sick_range), function(ii){
      cbind(sick_order[ii], sick_range[ii,1]: sick_range[ii,2])
    }))
    zero_out_idx = do.call('rbind', lapply(healthy_agents, function(ii){
      cbind(ii, tmp_idx_set)
    }))
    zero_out_idx = rbind(zero_out_idx, zero_out_idx[,c(2,1,3)])
    non_SI_interact_report_ind[zero_out_idx] = 0
    
    edges_to_sample_HH = extract_filtered_edges(non_SI_interact_report_ind, net_snapshots, obs_time)
    
    edges_to_sample_HI = edges_to_sample_HI_HH[!duplicated(rbind.data.frame(edges_to_sample_HH, edges_to_sample_HI_HH))[-c(1: nrow(edges_to_sample_HH))],]
    
    #
    # edges_to_sample = which(non_SI_interact_report_ind == 1, arr.ind = TRUE)
    # edges_to_sample = cbind.data.frame(edges_to_sample, edges_to_sample[,3]+1)
    # edges_to_sample = edges_to_sample[edges_to_sample[,1] < edges_to_sample[,2],]
    # colnames(edges_to_sample) = c('per1', 'per2', 't0', 't1')
    # edges_to_sample = cbind(edges_to_sample, cbind(net0 = net_snapshots[as.matrix(edges_to_sample[, -4])],
    #                                                net1 = net_snapshots[as.matrix(edges_to_sample[, -3])]))
    # edges_to_sample$t0 = obs_time[edges_to_sample$t0]
    # edges_to_sample$t1 = obs_time[edges_to_sample$t1]
    #
    # edges_to_sample = edges_to_sample[!duplicated(edges_to_sample), ]
    # if (nrow(edges_to_sample) > 0) rownames(edges_to_sample) = 1: nrow(edges_to_sample)
    #
    # edges_to_sample_health_sick = edges_to_sample_no2[!duplicated(rbind.data.frame(edges_to_sample, edges_to_sample_no2))[-c(1: nrow(edges_to_sample))],]
    # # edges_could_fix = rownames(edges_to_sample)
    
    
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
    
    recover.datA <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(recover.datA) = c('time', 'per1')
  }
  
  # imputed_net_events_A = list()
  for(it in 1:S){
    # it = 1
    if (it %% 100 == 0) print(it)
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    ncount = 0
    while (TRUE){
      recov_times_tmp = rep(Inf, nrow(G0))
      recov_times_tmp[as.integer(rownames(intervals))] = sam_texp(nrow(intervals), params.cur[2], 0, intervals[,2]-intervals[,1]) + intervals[,1]
      # recov_times_tmp = sam_texp(nrow(intervals), params.cur[2], 0, intervals[,2]-intervals[,1]) + intervals[,1]
      recov_times_tmp[is.nan(recov_times_tmp)] = Inf
      if (all(sapply(names(net_status_infect)[-1], function(sick_per){
        # sick_per = '28'
        conn_nei = net_status_infect[[sick_per]]
        
        any(infect_times[conn_nei] < infect_times[as.integer(sick_per)] &
            recov_times_tmp[conn_nei] > infect_times[as.integer(sick_per)])
      }))){
        # recov_times = rep(Inf, nrow(G0))
        # recov_times[as.integer(rownames(intervals))] = as.vector(recov_times_tmp)
        recov_times = as.vector(recov_times_tmp)
        if (verbose) print(paste0('sampling recovery times failed ', ncount, ' times'))
        break
      }else{
        ncount = ncount + 1
      }
    }
    # 
    # 
    st = Sys.time()
    SI_links = SI_events_imputer2(sick_order,
                                  nei_infec,
                                  recov_times,
                                  infect_times,
                                  net_status_infect,
                                  window_length,
                                  obs_time,
                                  net_snapshots,
                                  interact_report_ind,
                                  params.cur)
    ed = Sys.time()
    ed - st
    
    # st = Sys.time()
    # SI_links2 = SI_events_imputer2_no_sick(sick_order,
    #                                        healthy_agents,
    #                                        sick_range,
    #                                        net_status_infect,
    #                                        recov_times,
    #                                        infect_times,
    #                                        window_length,
    #                                        obs_time,
    #                                        net_snapshots,
    #                                        interact_report_ind,
    #                                        params.cur)
    # ed = Sys.time()
    # ed - st
    
    # edges_to_sample_tmp2 = edges_to_sample_health_sick
    # # Summarize the agent pair, interval, end point network status,
    # infect_times2 = infect_times
    # infect_times2[abs(infect_times2 - 0) < 1e-4] = NA
    # edges_to_sample_tmp2$i1 = infect_times2[edges_to_sample_tmp2$per1]
    # edges_to_sample_tmp2$i2 = infect_times2[edges_to_sample_tmp2$per2]
    # edges_to_sample_tmp2$r1 = recov_times[edges_to_sample_tmp2$per1]
    # edges_to_sample_tmp2$r2 = recov_times[edges_to_sample_tmp2$per2]
    # 
    # edges_to_sample_tmp2$i1[(edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
    #                            edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1)] = NA
    # 
    # edges_to_sample_tmp2$i2[(edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
    #                            edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1)] = NA
    # 
    # edges_to_sample_tmp2$r1[(edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
    #                            edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1)] = NA
    # 
    # edges_to_sample_tmp2$r2[(edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
    #                            edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1)] = NA
    # 
    # # idx_to_impute = (!rownames(edges_to_sample_tmp2) %in% imputed_net_events_A$edge_idx)
    # # edges_to_sample_tmp2= edges_to_sample_tmp2[idx_to_impute,]
    # 
    # edges_to_sample_tmp2$health_at_t0 = (infect_times[edges_to_sample_tmp2$per1] <= edges_to_sample_tmp2$t0 &
    #                                        edges_to_sample_tmp2$t0 < recov_times[edges_to_sample_tmp2$per1]) +
    #   (infect_times[edges_to_sample_tmp2$per2] <= edges_to_sample_tmp2$t0 &
    #      edges_to_sample_tmp2$t0 < recov_times[edges_to_sample_tmp2$per2])
    
    # # initial health status at the beginning of each interval
    # 
    # # st = Sys.time()
    # if (!is.null(init.params)) params.cur = c(params.cur[1:2], init.params[-c(1,2)])
    
    # tmp = impute_all(edges_to_sample_tmp2, params.cur, infect_times, recov_times)
    # 
    # sapply(tmp, function(x){any(x$event == 1)})
    
    edges_to_sample_tmp2_HI = process_network_edges(edges_to_sample_HI, infect_times, recov_times)
    edges_to_sample_tmp2_HH = process_network_edges(edges_to_sample_HH, infect_times, recov_times)
    
    imputed_net_events_HI = do.call("rbind", impute_all(edges_to_sample_tmp2_HI, params.cur + c(rep(0,6), params.cur[1], 0), infect_times, recov_times))
    imputed_net_events_HH = do.call("rbind", impute_all(edges_to_sample_tmp2_HH, params.cur, infect_times, recov_times))
    
    imputed_net_events = rbind.data.frame(SI_links, imputed_net_events_HI[,1:4], imputed_net_events_HH[,1:4])
    # imputed_net_events = rbind.data.frame(SI_links, SI_links2, imputed_net_events_HH[,1:4])
    # 
    # imputed_net_events_B = do.call("rbind", impute_all(edges_to_sample_tmp2, params.cur, infect_times, recov_times))
    # 
    # imputed_net_events_health_sick = do.call("rbind", impute_all(edges_to_sample_tmp2, params.cur + c(rep(0,6), params.cur[1], 0), infect_times, recov_times))
    # 
    # 
    # imputed_net_events = rbind.data.frame(SI_links, SI_links2,imputed_net_events_B[,1:4])
    # imputed_net_events1 = rbind.data.frame(SI_links,imputed_net_events_B[,1:4], imputed_net_events_health_sick[,1:4])
    
    events_temp = rbind.data.frame(sick_events, imputed_net_events)
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[-1, ]
    
    recover.dat = cbind.data.frame(time = recov_times[recov_times < max(recov_times)],
                                   per1 = which(recov_times < max(recov_times)))
    recover.dat = recover.dat[order(recover.dat$time),]
    # # 
    # table(events_temp$event)
    PA = parse_augment_cpp(G0, I0, events_temp, recover.dat)
    if (verbose) print(PA)
    
    
    
    # {
    #   net_events_temp = events_temp[events_temp$event %in% 3: 8,]
    #   G0 = dats$G0
    #   for (tt in 1: nrow(net_events_temp)) {
    #     # tt = 660
    #     if (tt %% 100 == 0) print(tt)
    #     event_tt = net_events_temp[tt, ]
    # 
    #     G0[event_tt$per1, event_tt$per2] <- G0[event_tt$per2, event_tt$per1] <- G0[event_tt$per2, event_tt$per1] + ifelse(event_tt$event%in% 3:5, 1, -1)
    #     if (!G0[event_tt$per2, event_tt$per1] %in% 0: 1){
    #       print('bad')
    #       break
    #     }
    #   }
    # }
    # net_events_temp[net_events_temp$per1 == 52 & net_events_temp$per2 == 3,]
    # tmp = net_snapshots[51, 3,]
    # names(tmp) = obs_time
    # tmp2 = interact_report_ind[51, 3,]
    # names(tmp2) = obs_time[-1]
    # 
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at ", it, "-th iter"))
    
    
    # edge_ids_to_fix = sample(edge_to_fix, length(edge_to_fix)*fix_prop, replace = FALSE)
    
    # patient_id_to_fix = sample(sick_pers, length(sick_pers)*fix_prop, replace = FALSE)
    # recover.datA = recover.dat[recover.dat$per1 %in% patient_id_to_fix, ]
    # we only fix edges unrelated to selected patients' recovery
    # edge_ids_to_fix = edge_ids_to_fix[!edge_ids_to_fix %in% sick_edges_vec[!names(sick_edges_vec) %in% patient_id_to_fix]] 
    # imputed_net_events_A = imputed_net_events[imputed_net_events$edge_idx %in% edge_ids_to_fix, ]
    
    
    
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


infer_miss_recov27_cpp = function(dat, net_status_infect, priors, obs_time, force_one = TRUE,
                                  init.params = NULL, recov_times = NULL,
                                  verbose = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42, remove_no_neighbour = FALSE){
  process_network_edges <- function(sampled_edges, infection_times, recovery_times) {
    # sampled_edges = edges_to_sample_HI[38,];
    # infection_times = infect_times
    # recovery_times = recov_times
    # Create a temporary copy
    temp_edges <- sampled_edges
    
    # Summarize the agent pair, interval, and endpoint network status
    adjusted_infect_times <- infection_times
    adjusted_infect_times[abs(adjusted_infect_times - 0) < 1e-4] <- NA
    
    temp_edges$i1 <- adjusted_infect_times[temp_edges$per1]
    temp_edges$i2 <- adjusted_infect_times[temp_edges$per2]
    temp_edges$r1 <- recovery_times[temp_edges$per1]
    temp_edges$r2 <- recovery_times[temp_edges$per2]
    
    # Adjust values based on time constraints
    temp_edges$i1[(temp_edges$i1 < temp_edges$t0 |
                     temp_edges$i1 > temp_edges$t1)] <- NA
    
    temp_edges$i2[(temp_edges$i2 < temp_edges$t0 |
                     temp_edges$i2 > temp_edges$t1)] <- NA
    
    temp_edges$r1[(temp_edges$r1 < temp_edges$t0 |
                     temp_edges$r1 > temp_edges$t1)] <- NA
    
    temp_edges$r2[(temp_edges$r2 < temp_edges$t0 |
                     temp_edges$r2 > temp_edges$t1)] <- NA
    
    # Compute health status at t0
    temp_edges$health_at_t0 <- (infection_times[temp_edges$per1] <= temp_edges$t0 &
                                  temp_edges$t0 < recovery_times[temp_edges$per1]) + 
      (infection_times[temp_edges$per2] <= temp_edges$t0 &
         temp_edges$t0 < recovery_times[temp_edges$per2])
    
    return(temp_edges)
  }
  
  
  extract_filtered_edges <- function(interaction_matrix, network_data, time_vector) {
    
    
    # Find edges where interaction report indicator is 1
    
    filtered_edges <- which(interaction_matrix == 1, arr.ind = TRUE)
    
    # Add a new column for time shift
    filtered_edges <- cbind.data.frame(filtered_edges, filtered_edges[,3] + 1)
    
    # Keep only edges where node1 < node2 to avoid duplicates
    filtered_edges <- filtered_edges[filtered_edges[,1] < filtered_edges[,2],]
    
    # Assign column names
    colnames(filtered_edges) <- c('per1', 'per2', 't0', 't1')
    
    # Extract network snapshots for t0 and t1
    filtered_edges <- cbind(filtered_edges, 
                            cbind(net0 = network_data[as.matrix(filtered_edges[, -4])], 
                                  net1 = network_data[as.matrix(filtered_edges[, -3])]))
    
    # Map time indices to actual observation times
    filtered_edges$t0 <- time_vector[filtered_edges$t0]
    filtered_edges$t1 <- time_vector[filtered_edges$t1]
    
    # Remove duplicate rows
    filtered_edges <- filtered_edges[!duplicated(filtered_edges), ]
    
    # Reset row names if there are rows
    if (nrow(filtered_edges) > 0) rownames(filtered_edges) <- 1:nrow(filtered_edges)
    
    return(filtered_edges)
  }

  {
    if(timing){ time.st = Sys.time()}
    
    set.seed(seed)
    
    # preparations
    G0 = dat$G0; I0 = dat$I0;
    # events = dat$events
    reports = dat$health_report; report.times = as.integer(rownames(reports))
    window_length = obs_time[2] - obs_time[1]
    sick_events = dat$sick_events
    sick_agents = sick_events$per1
    healthy_agents = c(1: nrow(G0))[-sick_events$per1]
    net_snapshots = dat$net_snapshots
    interact_report = dat$interaction_report
    interact_report_ind = array(0, dim = c(nrow(G0), nrow(G0), length(interact_report)))
    for (t in 1: length(interact_report)) {
      # t = 1
      interact_report_ind[cbind(as.matrix(interact_report[[t]]), t)] = 1
      interact_report_ind[cbind(as.matrix(interact_report[[t]])[,2:1], t)] = 1
    }
    
    interact_report_ind_II = array(0, dim = c(nrow(G0), nrow(G0), length(interact_report)))
    for (t in 1: length(interact_report)) {
      # t = 1
      interact_report_ind_II[cbind(as.matrix(interact_report[[t]][(interact_report[[t]][, 1] %in% sick_agents) & 
                                                                    (interact_report[[t]][, 2] %in% sick_agents),
      ]), t)] = 1
      interact_report_ind_II[cbind(as.matrix(interact_report[[t]][(interact_report[[t]][, 1] %in% sick_agents) & 
                                                                    (interact_report[[t]][, 2] %in% sick_agents),
                                                                  2: 1]), t)] = 1
    }
    
    
    infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
    infect_times[sick_events$per1] = sick_events$time
    
    if("truth" %in% names(dat)){
      true.params = dat$truth
    }else{
      true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
    }
    
    # Get recovery intervals for all patients
    {
      intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
      for (i in 1: nrow(sick_events)) {
        # i = 1
        per = sick_events$per1[i]
        
        ub = report.times[which(reports[, per] == -1)[1]]
        lb = report.times[which(reports[, per] == -1)[1]-1]
        intervals[i, ] = c(lb, ub)
      }
      rownames(intervals) = sick_events$per1
      colnames(intervals) = c("lb", "ub")
      # Add this to handle those who haven't recovered
      intervals = intervals[rowSums(!is.na(intervals)) == 2, ]
      intervals[,1] = rowMaxs(cbind(intervals[,1], infect_times[as.integer(rownames(intervals))]),value = TRUE)
      intervals = intervals[order(as.integer(rownames(intervals))), ]
    }
    intervals[intervals[,2] > obs_time[length(obs_time)], 2] = obs_time[length(obs_time)]
    
    
    tmp = cbind(1: nrow(G0), infect_times)
    tmp = tmp[tmp[,1] %in% sick_events$per1, ]
    
    sick_order = tmp[order(tmp[,2]), 1]
    nei_infec = list()
    for (ii in 2: length(sick_order)) {
      nei_infec[[as.character(sick_order[ii])]] = sick_order[1: (ii-1)]
    }
    ik_intvl_idx = cbind(infect_times %/% window_length + 1,
                         infect_times %/% window_length + 2)
    
    
    edges_to_sample_II = extract_filtered_edges(interact_report_ind_II, net_snapshots, obs_time)
    
    non_SI_interact_report_ind = interact_report_ind
    
    for (i in as.integer(names(nei_infec))) {
      # i = 111
      nei_infec_i = nei_infec[[as.character(i)]]
      for (j in nei_infec_i) {
        non_SI_interact_report_ind[rbind(cbind(i, j, ik_intvl_idx[j]: ik_intvl_idx[i]),
                                         cbind(j, i, ik_intvl_idx[j]: ik_intvl_idx[i]))] = 0
      }
    }
    
    edges_to_sample_HI_HH = extract_filtered_edges(non_SI_interact_report_ind, net_snapshots, obs_time)
    
    report.times[report.times >= obs_time[length(obs_time)]] = obs_time[length(obs_time)]-1e-5
    sick_range = cbind(infect_times[sick_order] %/% window_length + 1,
                       apply(reports[, sick_order], 2, function(x){
                         if (-1 %in% x){
                           report.times[which(x == -1)[1]]
                         }else{
                           report.times[length(x)]
                         }
                       })
    )
    rownames(sick_range) = sick_order
    sick_range[,2] = sick_range[,2] %/% window_length + 1
    
    tmp_idx_set = do.call('rbind', lapply(1: nrow(sick_range), function(ii){
      cbind(sick_order[ii], sick_range[ii,1]: sick_range[ii,2])
    }))
    zero_out_idx = do.call('rbind', lapply(healthy_agents, function(ii){
      cbind(ii, tmp_idx_set)
    }))
    zero_out_idx = rbind(zero_out_idx, zero_out_idx[,c(2,1,3)])
    non_SI_interact_report_ind[zero_out_idx] = 0
    
    edges_to_sample_HH = extract_filtered_edges(non_SI_interact_report_ind, net_snapshots, obs_time)
    
    edges_to_sample_HI = edges_to_sample_HI_HH[!duplicated(rbind.data.frame(edges_to_sample_HH, edges_to_sample_HI_HH))[-c(1: nrow(edges_to_sample_HH))],]
    
    
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
    
    recover.datA <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(recover.datA) = c('time', 'per1')
  }
  
  for(it in 1:S){
    # it = 1
    if (it %% 100 == 0) print(it)
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    #### Sample Recovery Times ####
    ncount = 0
    while (TRUE){
      recov_times_tmp = rep(Inf, nrow(G0))
      recov_times_tmp[as.integer(rownames(intervals))] = sam_texp(nrow(intervals), params.cur[2], 0, intervals[,2]-intervals[,1]) + intervals[,1]
      recov_times_tmp[is.nan(recov_times_tmp)] = Inf
      if (all(sapply(names(net_status_infect)[-1], function(sick_per){
        # sick_per = '28'
        conn_nei = net_status_infect[[sick_per]]
        
        any(infect_times[conn_nei] < infect_times[as.integer(sick_per)] &
            recov_times_tmp[conn_nei] > infect_times[as.integer(sick_per)])
      }))){
        recov_times = as.vector(recov_times_tmp)
        if (verbose) print(paste0('sampling recovery times failed ', ncount, ' times'))
        break
      }else{
        ncount = ncount + 1
      }
    }
    
    #### Sample Missing Interaction Events ####
    SI_links = SI_events_imputer3(sick_order,
                                  nei_infec,
                                  recov_times,
                                  infect_times,
                                  net_status_infect,
                                  window_length,
                                  obs_time,
                                  net_snapshots,
                                  interact_report_ind,
                                  params.cur)

    SI_links2 = SI_events_imputer3_no_sick(sick_order,
                                           healthy_agents,
                                           sick_range,
                                           net_status_infect,
                                           recov_times,
                                           infect_times,
                                           window_length,
                                           obs_time,
                                           net_snapshots,
                                           interact_report_ind,
                                           params.cur)

  
    edges_to_sample_tmp2_HH = process_network_edges(edges_to_sample_HH, infect_times, recov_times)
    
    imputed_net_events_HH = do.call('rbind.data.frame',apply(edges_to_sample_tmp2_HH, 1, function(x){
      x = as.numeric(x)
      while (TRUE) {
        event_tmp = simulate_ctmc3_cpp(params.cur[3], params.cur[6], vector('numeric', 0), 
                                       as.numeric(x[4] - x[3]), 
                                       x[5],
                                       x[6])
        if(nrow(event_tmp) > 0) break
      }
      event_tmp = cbind.data.frame(event_tmp, per1 = x[1], per2 = x[2])
      event_tmp$event = ifelse(event_tmp$event == 1, 3, 6)
      event_tmp$time = event_tmp$time + x[3]
      
      event_tmp
    }))
    
    
    imputed_net_events = rbind.data.frame(SI_links, SI_links2, imputed_net_events_HH[,1:4])
    
    events_temp = rbind.data.frame(sick_events, imputed_net_events)
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[-1, ]
    
    recover.dat = cbind.data.frame(time = recov_times[recov_times < max(recov_times)],
                                   per1 = which(recov_times < max(recov_times)))
    recover.dat = recover.dat[order(recover.dat$time),]

    
    PA = parse_augment_cpp(G0, I0, events_temp, recover.dat)
    if (verbose) print(PA)
    
    #### Sample Model Parameters ####
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at ", it, "-th iter"))
  
    
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

##################################################################################################################
SI_events_imputer4 = function(sick_order2,
                              nei_infec2,
                              recov_times2,
                              infect_times2,
                              net_status_infect2,
                              window_length2,
                              obs_time2,
                              net_snapshots2,
                              interact_report_ind2,
                              params.cur2,
                              force_one){
  
  last_SI_events_and_edge_info = lapply(sick_order2[-1], function(i){
    x = nei_infec2[[as.character(i)]]
    non_nei = x[recov_times2[x] <= infect_times2[i]]
    x = x[recov_times2[x] > infect_times2[i]]
    
    net_status_ik = rep(0, length(x))
    net_status_ik[x %in% net_status_infect2[[as.character(i)]]] = 1
    
    # this is the network status of (j ,k)
    net_status_ij = as.integer(sapply(as.character(x), function(xx){i %in% net_status_infect2[[xx]]}))
    
    intvl_idx = (infect_times2[x] %/% window_length2) + 1
    
    out_t0ij = cbind.data.frame(per1 = x, per2 = i,
                                t0 = obs_time2[intvl_idx],
                                t1 = infect_times2[x],
                                net0 = net_snapshots2[cbind(x, i, intvl_idx)], 
                                net1 = net_status_ij,
                                force1 = interact_report_ind2[cbind(x, i, intvl_idx)])
    out_t0ij$t1[out_t0ij$t1-out_t0ij$t0 == 0] = out_t0ij$t0[out_t0ij$t1-out_t0ij$t0 == 0] + 1e-5
    
    out_ikt1 = cbind.data.frame(per1 = x, per2 = i,
                                t0 = infect_times2[i],
                                t1 = obs_time2[(infect_times2[i] %/% window_length2) + 2],
                                net0 = net_status_ik, 
                                net1 = net_snapshots2[cbind(x, i, (infect_times2[i] %/% window_length2) + 2)],
                                force1 = interact_report_ind2[cbind(x, i, (infect_times2[i] %/% window_length2) + 1)])
    out_ijik = cbind.data.frame(per1 = x, per2 = i, 
                                t0 = infect_times2[x], 
                                t1 = infect_times2[i],
                                net0 = net_status_ij, 
                                net1 = net_status_ik,
                                force1 = interact_report_ind2[cbind(x, i, (infect_times2[i] %/% window_length2) + 1)])
    edge_info_nei = get_SI_link_info(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2)
    
    intvl_idx_k = infect_times2[i] %/% window_length2+1
    edge_info_nonnei = list()
    for (j in non_nei){
      ij = infect_times2[j]; rj = recov_times2[j]
      intvl_idx_j = ij %/% window_length2+1
      
      t0 = obs_time2[intvl_idx_j]
      
      link_type_changes = rbind(c(t0,ij, rj, infect_times2[i], recov_times2[i]), 
                                c(0, 1, -1, 1, -1),
                                c(0, 1, -1, 0, 0),
                                c(intvl_idx_j, rep(NA, 4)),
                                c(net_snapshots2[i,j,intvl_idx_j], 
                                  as.integer(i %in% net_status_infect2[as.character(j)]), 
                                  NA,  as.integer(j %in% net_status_infect2[as.character(i)]), NA))
      
      link_type_changes = cbind(link_type_changes, rbind(obs_time2[(intvl_idx_j+1):(intvl_idx_k+1)], NA, NA, (intvl_idx_j+1):(intvl_idx_k+1), net_snapshots2[i,j, intvl_idx_j: intvl_idx_k]))
      link_type_changes = link_type_changes[, order(link_type_changes[1,])]
      
      link_type_changes = link_type_changes[, link_type_changes[1,] <= obs_time2[(intvl_idx_k+1)]]
      link_type_changes[2, is.na(link_type_changes[2, ])] = 0 
      link_type_changes[2, ] = cumsum(link_type_changes[2,])
      
      link_type_changes[3, is.na(link_type_changes[3, ])] = 0
      link_type_changes[3, ] = cumsum(link_type_changes[3,])
      for (jj in (intvl_idx_j): (intvl_idx_k)) {
        if (interact_report_ind2[j, i, jj] == 1){
          tmp = link_type_changes[, which(link_type_changes[4, ] == jj): which(link_type_changes[4, ] == (jj+1))]
          
          link_info_i = dissect_matrix(tmp)
          
          
          
          edge_info_nonnei = append(edge_info_nonnei, 
                                    (lapply(link_info_i, function(tmp_x){
                                      init_time = tmp_x[1, 1]
                                      tmp_x[1,] = tmp_x[1,] - init_time
                                      n = ncol(tmp_x)
                                      list(params.cur2[tmp_x[2, 1: (n-1)]+3], params.cur2[tmp_x[2, 1: (n-1)]+6] + tmp_x[3, 1: (n-1)]*params.cur2[1],
                                           tmp_x[1, -c(1, n)],
                                           tmp_x[1, n],
                                           tmp_x[5, 1],
                                           tmp_x[5, n],
                                           init_time, c(i, j))
                                    })))
          
        }
      }
    }
    
    
    out = append(edge_info_nei$edge_info, edge_info_nonnei)
    return(list(edge_info_nei$events, out))
  })
  
  last_SI_events = do.call('rbind.data.frame', lapply(last_SI_events_and_edge_info, function(x){x[[1]]}))
  
  SI_edge_info_by_sick_per = lapply(last_SI_events_and_edge_info, function(x){x[[2]]})
  all_edge_info =  do.call('c', lapply(last_SI_events_and_edge_info, function(x){x[[2]]}))
  
  
  # Apply function
  grouped_interactions <- regroup_interactions(all_edge_info, window_length)
  
  tmp = do.call('rbind', lapply(grouped_interactions, function(x){
    num_samples = 1
    while (TRUE) {
      events = do.call('rbind', lapply(x, function(xx){
        tmp_events = simulate_ctmc3_cpp(xx[[1]], xx[[2]], xx[[3]], xx[[4]], xx[[5]], xx[[6]])
        if (nrow(tmp_events) > 0) tmp_events$time = tmp_events$time + xx[[7]]
        tmp_events
      }))
      if (nrow(events) > 0 | num_samples >= force_one) break
      num_samples = num_samples + 1
    }
    if (nrow(events) > 1){
      out = cbind.data.frame(events, per1 = x[[1]][[8]][1], per2 = x[[1]][[8]][2])
    }else{
      out = data.frame(time=c(),event=c(), per1=c(), per2=c())
    }
    out
  }))
  tmp$event = ifelse(tmp$event == 1, 3, 6) + 
    (infect_times2[tmp$per1] <= tmp$time &
       tmp$time < recov_times2[tmp$per1]) +
    (infect_times2[tmp$per2] <= tmp$time &
       tmp$time < recov_times2[tmp$per2])
  
  out_data  = rbind.data.frame(last_SI_events, tmp)
  out_data[order(out_data$time),]
}

SI_events_imputer4_no_sick = function(sick_agents2,
                                      healthy_agents2,
                                      sick_range2,
                                      nei_infec2,
                                      recov_times2,
                                      infect_times2,
                                      window_length2,
                                      obs_time2,
                                      net_snapshots2,
                                      interact_report_ind2,
                                      params.cur2,
                                      force_one){
  
  all_edge_info = do.call('c', lapply(healthy_agents2, function(i){
    edge_info = list()
    for (j in sick_agents2){
      intvl_idx_k = sick_range2[as.character(j),2]
      ij = infect_times2[j]; rj = recov_times2[j]
      intvl_idx_j = ij %/% window_length2+1
      
      t0 = obs_time2[intvl_idx_j]
      
      
      link_type_changes = rbind(c(t0, ij, rj), # time
                                c(0, 1, -1), # health status
                                c(0, 1, -1), # indication of SI link
                                c(intvl_idx_j, rep(NA, 2)), # observational intervel idx 
                                c(net_snapshots2[j,i,intvl_idx_j], as.integer(i %in% nei_infec2[[as.character(j)]]), NA)) # network status
      
      link_type_changes = cbind(link_type_changes, rbind(obs_time2[(intvl_idx_j+1):(intvl_idx_k+1)], NA, NA, (intvl_idx_j+1):(intvl_idx_k+1), net_snapshots2[j,i,(intvl_idx_j+1):(intvl_idx_k+1)]))
      link_type_changes = link_type_changes[, order(link_type_changes[1,])]
      
      link_type_changes = link_type_changes[, link_type_changes[1,] <= obs_time2[(intvl_idx_k+1)]]
      link_type_changes[2, is.na(link_type_changes[2, ])] = 0 
      link_type_changes[2, ] = cumsum(link_type_changes[2,])
      
      link_type_changes[3, is.na(link_type_changes[3, ])] = 0
      link_type_changes[3, ] = cumsum(link_type_changes[3,])
      
      link_type_changes = rbind(link_type_changes, rep(0, ncol(link_type_changes))) 
      nn = ncol(link_type_changes)
      link_type_changes[,-nn][6, !is.na(link_type_changes[5, -nn])] = interact_report_ind2[j, i, link_type_changes[,-nn][1, !is.na(link_type_changes[5, -nn])] %/% window_length2 + 1]
      
      if (link_type_changes[1,1] == link_type_changes[1,2]){
        link_type_changes = link_type_changes[,-1]
      }
      breakpoint_idx = which(!is.na(link_type_changes[5,]))
      for (jj in 1: (ncol(link_type_changes)-1)) {
        if (link_type_changes[6, jj] == 1 & !is.na(link_type_changes[5, jj])){
          
          last_jj = min(breakpoint_idx[breakpoint_idx > jj])
          
          if (jj == (last_jj-1)){
            breakpoints = vector('numeric')
          }else{
            breakpoints = link_type_changes[1,(jj+1) : (last_jj-1)] - link_type_changes[1, jj]
          }
          edge_info = append(edge_info, 
                             list(list(params.cur2[link_type_changes[2,jj: (last_jj-1)]+3], params.cur2[link_type_changes[2,jj: (last_jj-1)]+6] + link_type_changes[3, jj: (last_jj-1)]*params.cur2[1],
                                       breakpoints, # break points
                                       link_type_changes[1, last_jj] - link_type_changes[1, jj], 
                                       link_type_changes[5, jj],
                                       link_type_changes[5, last_jj],
                                       link_type_changes[1, jj], c(j, i))))
        }
      }
    }
    
    edge_info
  }))
  
  # Apply function
  grouped_interactions <- regroup_interactions(all_edge_info, window_length)
  
  tmp = do.call('rbind', lapply(grouped_interactions, function(x){
    num_samples = 1
    while (TRUE) {
      
      events = do.call('rbind', lapply(x, function(xx){
        tmp_events = simulate_ctmc3_cpp(xx[[1]], xx[[2]], xx[[3]], xx[[4]], xx[[5]], xx[[6]])
        if (nrow(tmp_events) > 0) tmp_events$time = tmp_events$time + xx[[7]]
        tmp_events
      }))
      if (nrow(events) > 0 | num_samples >= force_one) break
      num_samples = num_samples + 1
    }
    if (nrow(events) > 1){
      out = cbind.data.frame(events, per1 = x[[1]][[8]][1], per2 = x[[1]][[8]][2])
    }else{
      out = data.frame(time=c(),event=c(), per1=c(), per2=c())
    }
    out
  }))
  tmp$event = ifelse(tmp$event == 1, 3, 6) + 
    (infect_times2[tmp$per1] <= tmp$time &
       tmp$time < recov_times2[tmp$per1]) +
    (infect_times2[tmp$per2] <= tmp$time &
       tmp$time < recov_times2[tmp$per2])
  
  out_data  = tmp[order(tmp$time),]
  out_data
}


infer_miss_recov28_cpp = function(dat, net_status_infect, priors, obs_time, force_one = 100,
                                  init.params = NULL, recov_times = NULL,
                                  verbose = T, 
                                  samples = 1000, burn = 100, thin = 1, seed=42){
  #' Infer Missing Recovery Times and Model Parameters
  #'
  #' This function performs MCMC-based inference to impute missing recovery times and estimate model parameters
  #' for an epidemic process on a dynamic contact network.
  #'
  #' @param dat A list containing data for the model, including the initial network `G0`, initial infection states `I0`,
  #'            health reports, sick events, interaction reports, and network snapshots.
  #' @param net_status_infect A list mapping each sick individual to the set of individuals they could have been infected by (network-based).
  #' @param priors A data frame specifying the prior distribution parameters: columns `count` (shape) and `avg` (mean).
  #'               If it has 4 rows, the function will expand it to 8 parameters assuming repeated values.
  #' @param obs_time A numeric vector of observation times, typically the time range over which the epidemic is observed.
  #' @param force_one int, default 100. This parameter is the least number of times that we need to sample an edge. Some edge can have very low probability
  #' to have events occur within that period, but if we have to sample these edges, we will sample at least `force_one` times. After that, we will stop, 
  #' even we do not have anything for that edge. This is a safeguard to ensure that we do not get stuck.
  #' 
  #' @param init.params Optional. A vector of initial values for the 8 model parameters. If NULL, values will be sampled from priors.
  #' @param recov_times Optional. Predefined recovery times. If NULL, they will be inferred.
  #' @param verbose Logical. If TRUE, intermediate output and progress will be printed during execution.
  #' @param samples Number of posterior samples to return after burn-in and thinning.
  #' @param burn Number of initial MCMC iterations to discard (burn-in).
  #' @param thin Thinning interval for MCMC. Only every `thin`th sample will be kept.
  #' @param seed Random seed for reproducibility.
  #'
  #' @return A list containing:
  #'         - A matrix of sampled model parameters (`params`) with one row per retained iteration.
  #'         - An array `sample_hyperparm` of shape (2, 8, total_iterations), recording event counts and sufficient statistics.
  #'
  #' @examples
  #' result <- infer_miss_recov27_cpp(dat, net_status_infect, priors, obs_time)
  
  process_network_edges <- function(sampled_edges, infection_times, recovery_times) {
    # sampled_edges = edges_to_sample_HI[38,];
    # infection_times = infect_times
    # recovery_times = recov_times
    # Create a temporary copy
    temp_edges <- sampled_edges
    
    # Summarize the agent pair, interval, and endpoint network status
    adjusted_infect_times <- infection_times
    adjusted_infect_times[abs(adjusted_infect_times - 0) < 1e-4] <- NA
    
    temp_edges$i1 <- adjusted_infect_times[temp_edges$per1]
    temp_edges$i2 <- adjusted_infect_times[temp_edges$per2]
    temp_edges$r1 <- recovery_times[temp_edges$per1]
    temp_edges$r2 <- recovery_times[temp_edges$per2]
    
    # Adjust values based on time constraints
    temp_edges$i1[(temp_edges$i1 < temp_edges$t0 |
                     temp_edges$i1 > temp_edges$t1)] <- NA
    
    temp_edges$i2[(temp_edges$i2 < temp_edges$t0 |
                     temp_edges$i2 > temp_edges$t1)] <- NA
    
    temp_edges$r1[(temp_edges$r1 < temp_edges$t0 |
                     temp_edges$r1 > temp_edges$t1)] <- NA
    
    temp_edges$r2[(temp_edges$r2 < temp_edges$t0 |
                     temp_edges$r2 > temp_edges$t1)] <- NA
    
    # Compute health status at t0
    temp_edges$health_at_t0 <- (infection_times[temp_edges$per1] <= temp_edges$t0 &
                                  temp_edges$t0 < recovery_times[temp_edges$per1]) + 
      (infection_times[temp_edges$per2] <= temp_edges$t0 &
         temp_edges$t0 < recovery_times[temp_edges$per2])
    
    return(temp_edges)
  }
  
  
  extract_filtered_edges <- function(interaction_matrix, network_data, time_vector) {
    
    
    # Find edges where interaction report indicator is 1
    
    filtered_edges <- which(interaction_matrix == 1, arr.ind = TRUE)
    
    # Add a new column for time shift
    filtered_edges <- cbind.data.frame(filtered_edges, filtered_edges[,3] + 1)
    
    # Keep only edges where node1 < node2 to avoid duplicates
    filtered_edges <- filtered_edges[filtered_edges[,1] < filtered_edges[,2],]
    
    # Assign column names
    colnames(filtered_edges) <- c('per1', 'per2', 't0', 't1')
    
    # Extract network snapshots for t0 and t1
    filtered_edges <- cbind(filtered_edges, 
                            cbind(net0 = network_data[as.matrix(filtered_edges[, -4])], 
                                  net1 = network_data[as.matrix(filtered_edges[, -3])]))
    
    # Map time indices to actual observation times
    filtered_edges$t0 <- time_vector[filtered_edges$t0]
    filtered_edges$t1 <- time_vector[filtered_edges$t1]
    
    # Remove duplicate rows
    filtered_edges <- filtered_edges[!duplicated(filtered_edges), ]
    
    # Reset row names if there are rows
    if (nrow(filtered_edges) > 0) rownames(filtered_edges) <- 1:nrow(filtered_edges)
    
    return(filtered_edges)
  }
  
  {
    set.seed(seed)
    
    # preparations
    G0 = dat$G0; I0 = dat$I0;
    # events = dat$events
    reports = dat$health_report; report.times = as.integer(rownames(reports))
    window_length = obs_time[2] - obs_time[1]
    sick_events = dat$sick_events
    sick_agents = sick_events$per1
    healthy_agents = c(1: nrow(G0))[-sick_events$per1]
    net_snapshots = dat$net_snapshots
    interact_report = dat$interaction_report
    interact_report_ind = array(0, dim = c(nrow(G0), nrow(G0), length(interact_report)))
    for (t in 1: length(interact_report)) {
      # t = 1
      interact_report_ind[cbind(as.matrix(interact_report[[t]]), t)] = 1
      interact_report_ind[cbind(as.matrix(interact_report[[t]])[,2:1], t)] = 1
    }
    
    interact_report_ind_II = array(0, dim = c(nrow(G0), nrow(G0), length(interact_report)))
    for (t in 1: length(interact_report)) {
      # t = 1
      interact_report_ind_II[cbind(as.matrix(interact_report[[t]][(interact_report[[t]][, 1] %in% sick_agents) & 
                                                                    (interact_report[[t]][, 2] %in% sick_agents),
      ]), t)] = 1
      interact_report_ind_II[cbind(as.matrix(interact_report[[t]][(interact_report[[t]][, 1] %in% sick_agents) & 
                                                                    (interact_report[[t]][, 2] %in% sick_agents),
                                                                  2: 1]), t)] = 1
    }
    
    
    infect_times = rep(obs_time[length(obs_time)]+1, length = nrow(G0))
    infect_times[sick_events$per1] = sick_events$time
    
    if("truth" %in% names(dat)){
      true.params = dat$truth
    }else{
      true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
    }
    
    # Get recovery intervals for all patients
    {
      intervals = matrix(NA, nrow = nrow(sick_events), ncol = 2)
      for (i in 1: nrow(sick_events)) {
        # i = 1
        per = sick_events$per1[i]
        
        ub = report.times[which(reports[, per] == -1)[1]]
        lb = report.times[which(reports[, per] == -1)[1]-1]
        intervals[i, ] = c(lb, ub)
      }
      rownames(intervals) = sick_events$per1
      colnames(intervals) = c("lb", "ub")
      # Add this to handle those who haven't recovered
      intervals = intervals[rowSums(!is.na(intervals)) == 2, ]
      intervals[,1] = rowMaxs(cbind(intervals[,1], infect_times[as.integer(rownames(intervals))]),value = TRUE)
      intervals = intervals[order(as.integer(rownames(intervals))), ]
    }
    intervals[intervals[,2] > obs_time[length(obs_time)], 2] = obs_time[length(obs_time)]
    
    
    tmp = cbind(1: nrow(G0), infect_times)
    tmp = tmp[tmp[,1] %in% sick_events$per1, ]
    
    sick_order = tmp[order(tmp[,2]), 1]
    nei_infec = list()
    for (ii in 2: length(sick_order)) {
      nei_infec[[as.character(sick_order[ii])]] = sick_order[1: (ii-1)]
    }
    ik_intvl_idx = cbind(infect_times %/% window_length + 1,
                         infect_times %/% window_length + 2)
    
    
    edges_to_sample_II = extract_filtered_edges(interact_report_ind_II, net_snapshots, obs_time)
    
    non_SI_interact_report_ind = interact_report_ind
    
    for (i in as.integer(names(nei_infec))) {
      # i = 111
      nei_infec_i = nei_infec[[as.character(i)]]
      for (j in nei_infec_i) {
        non_SI_interact_report_ind[rbind(cbind(i, j, ik_intvl_idx[j]: ik_intvl_idx[i]),
                                         cbind(j, i, ik_intvl_idx[j]: ik_intvl_idx[i]))] = 0
      }
    }
    
    edges_to_sample_HI_HH = extract_filtered_edges(non_SI_interact_report_ind, net_snapshots, obs_time)
    
    report.times[report.times >= obs_time[length(obs_time)]] = obs_time[length(obs_time)]-1e-5
    sick_range = cbind(infect_times[sick_order] %/% window_length + 1,
                       apply(reports[, sick_order], 2, function(x){
                         if (-1 %in% x){
                           report.times[which(x == -1)[1]]
                         }else{
                           report.times[length(x)]
                         }
                       })
    )
    rownames(sick_range) = sick_order
    sick_range[,2] = sick_range[,2] %/% window_length + 1
    
    tmp_idx_set = do.call('rbind', lapply(1: nrow(sick_range), function(ii){
      cbind(sick_order[ii], sick_range[ii,1]: sick_range[ii,2])
    }))
    zero_out_idx = do.call('rbind', lapply(healthy_agents, function(ii){
      cbind(ii, tmp_idx_set)
    }))
    zero_out_idx = rbind(zero_out_idx, zero_out_idx[,c(2,1,3)])
    non_SI_interact_report_ind[zero_out_idx] = 0
    
    edges_to_sample_HH = extract_filtered_edges(non_SI_interact_report_ind, net_snapshots, obs_time)
    
    edges_to_sample_HI = edges_to_sample_HI_HH[!duplicated(rbind.data.frame(edges_to_sample_HH, edges_to_sample_HI_HH))[-c(1: nrow(edges_to_sample_HH))],]
    
    
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
    
    recover.datA <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(recover.datA) = c('time', 'per1')
  }
  
  for(it in 1:S){
    # it = 1
    if (it %% 100 == 0) print(it)
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    #### Sample Recovery Times ####
    ncount = 0
    while (TRUE){
      recov_times_tmp = rep(Inf, nrow(G0))
      recov_times_tmp[as.integer(rownames(intervals))] = sam_texp(nrow(intervals), params.cur[2], 0, intervals[,2]-intervals[,1]) + intervals[,1]
      recov_times_tmp[is.nan(recov_times_tmp)] = Inf
      if (all(sapply(names(net_status_infect)[-1], function(sick_per){
        # sick_per = '28'
        conn_nei = net_status_infect[[sick_per]]
        
        any(infect_times[conn_nei] < infect_times[as.integer(sick_per)] &
            recov_times_tmp[conn_nei] > infect_times[as.integer(sick_per)])
      }))){
        recov_times = as.vector(recov_times_tmp)
        if (verbose) print(paste0('sampling recovery times failed ', ncount, ' times'))
        break
      }else{
        ncount = ncount + 1
      }
    }
    
    #### Sample Missing Interaction Events ####
    SI_links = SI_events_imputer4(sick_order,
                                  nei_infec,
                                  recov_times,
                                  infect_times,
                                  net_status_infect,
                                  window_length,
                                  obs_time,
                                  net_snapshots,
                                  interact_report_ind,
                                  params.cur,
                                  force_one)
    SI_links2 = SI_events_imputer4_no_sick(sick_order,
                                           healthy_agents,
                                           sick_range,
                                           net_status_infect,
                                           recov_times,
                                           infect_times,
                                           window_length,
                                           obs_time,
                                           net_snapshots,
                                           interact_report_ind,
                                           params.cur,
                                           force_one)
    
    
    edges_to_sample_tmp2_HH = process_network_edges(edges_to_sample_HH, infect_times, recov_times)
    
    imputed_net_events_HH = do.call('rbind.data.frame',apply(edges_to_sample_tmp2_HH, 1, function(x){
      x = as.numeric(x)
      num_samples = 1
      while (TRUE) {
        event_tmp = simulate_ctmc3_cpp(params.cur[3], params.cur[6], vector('numeric', 0), 
                                       as.numeric(x[4] - x[3]), 
                                       x[5],
                                       x[6])
        if(nrow(event_tmp) > 0 | num_samples >= force_one) break
        num_samples = num_samples + 1
      }
      if (nrow(event_tmp) > 0){
        event_tmp = cbind.data.frame(event_tmp, per1 = x[1], per2 = x[2])
        event_tmp$event = ifelse(event_tmp$event == 1, 3, 6)
        event_tmp$time = event_tmp$time + x[3]
      }else{
        event_tmp = data.frame(time=c(),event=c(), per1=c(), per2=c())
      }
      
      event_tmp
    }))
    
    
    imputed_net_events = rbind.data.frame(SI_links, SI_links2, imputed_net_events_HH[,1:4])
    
    events_temp = rbind.data.frame(sick_events, imputed_net_events)
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[-1, ]
    
    recover.dat = cbind.data.frame(time = recov_times[recov_times < max(recov_times)],
                                   per1 = which(recov_times < max(recov_times)))
    recover.dat = recover.dat[order(recover.dat$time),]
    
    
    PA = parse_augment_cpp(G0, I0, events_temp, recover.dat)
    if (verbose) print(PA)
    
    #### Sample Model Parameters ####
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at ", it, "-th iter"))
    
    
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
  
  return(list(params, sample_hyperparm = imp_event_log))
}

##################################################################################################################

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


infer_miss_recov16_cpp = function(dat, priors, obs_time, force_one = TRUE, fix_prop = 0.3,
                                  init.params = NULL,
                                  verbose = T, plot = T, output.sams = 100, 
                                  samples = 1000, burn = 100, thin = 1,
                                  timing = T, seed=42, remove_no_neighbour = FALSE){
  
  # dat = survey_dat2; priors = pr; obs_time = obs_time; output.sams = 100; plot = F; verbose = T;
  # samples = 50; burn = 0; thin = 1; seed = 1; init.params = NULL; timing = T; force_one = TRUE
  # remove_no_neighbour = FALSE; fix_prop = 0.9
  
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
  rownames(edges_to_sample) = 1: nrow(edges_to_sample)
  edges_could_fix = rownames(edges_to_sample)
  
  sick_pers = rownames(intervals)
  sick_edges = list()
  for (i in 1: nrow(intervals)){
    # i = 1
    row_i = intervals[i, ]
    sick_edges_i = rownames(edges_to_sample)[((edges_to_sample$per1 == sick_pers[i]) | 
                                                (edges_to_sample$per2 == sick_pers[i])) & 
                                               ((row_i[1] <= edges_to_sample$t0 & edges_to_sample$t0 <= row_i[2]) | 
                                                  (row_i[1] <= edges_to_sample$t1 & edges_to_sample$t1  <= row_i[2]) )]
    if (length(sick_edges_i) > 0){
      sick_edges[[sick_pers[i]]] = sick_edges_i
    }
  }
  sick_edges_vec = unlist(sick_edges)
  names(sick_edges_vec) = unlist(sapply(names(sick_edges), function(x){rep(x, length(sick_edges[[x]]))}))
  
  
  
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
  
  recover.datA <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(recover.datA) = c('time', 'per1')
  
  imputed_net_events_A = list()
  for(it in 1:S){
    # it = 2
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
      
      intervals_i = intervals[!rownames(intervals) %in% as.character(recover.datA$per1),]
      recovs =  rownames(intervals_i)[intervals_i[, 1] == report.times[ix]] %>% as.integer()
      
      # recovs =  rownames(intervals)[intervals[, 1] == report.times[ix]] %>% as.integer()
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
    recover.dat = rbind.data.frame(recover.dat, recover.datA)
    recover.dat = recover.dat[order(recover.dat$time),]
    if(verbose){ cat("Recovery times and network events imputation done.\n") }
    recov_times = rep(max(c(obs_time[length(obs_time)], recover.dat$time))+1, length = nrow(G0))
    recov_times[recover.dat$per1] = recover.dat$time
    
    # Prepare the edges and their time intervals we would like to sample
    # Pick patient-healthy per pair and note an connection when the healthy person gets sick
    
    edges_to_sample_tmp2 = edges_to_sample
    
    {
      nei_infec_temp2 = cbind.data.frame(sick_per = as.integer(names(nei_infec_temp)), 
                                         sick_nei = unlist(nei_infec_temp))
      interval_idx2 = cut(infect_times[nei_infec_temp2$sick_per], obs_time, 1: (length(obs_time)-1)) %>% as.integer()
      temp_mat = t(apply(as.matrix(nei_infec_temp2), 1, sort)); rownames(temp_mat) = NULL
      SI_links2= cbind.data.frame(temp_mat, t0 = obs_time[interval_idx2], t1 = obs_time[interval_idx2+1])
      colnames(SI_links2) = colnames(edges_to_sample_tmp2)[1: 4]
      temp_idx = which(rbind(SI_links2, edges_to_sample_tmp2[, 1: 4]) %>% duplicated(fromLast = TRUE))
      
      if (length(temp_idx) > 0){
        SI_links2 = SI_links2[temp_idx, ]
        SI_links2['temp_idx'] = rownames(SI_links2)
        
        interval_idx2 = interval_idx2[temp_idx]
        
        edge_idx_to_split = (which(rbind(SI_links2[,1:4], edges_to_sample_tmp2[, 1: 4]) %>% duplicated()) - nrow(SI_links2))
        
        
        SI_links2 = SI_links2 %>% arrange(t0, per1, per2)
        rownames(SI_links2) = SI_links2$temp_idx
        SI_links2 = SI_links2[,1: 4]
        
        tmp = edges_to_sample_tmp2[edge_idx_to_split, ] 
        tmp['edge_id'] = rownames(tmp)
        tmp = tmp %>% arrange(t0, per1, per2)
        edge_to_split = tmp$edge_id
        edge_to_fix = rownames(edges_to_sample_tmp2)[!rownames(edges_to_sample_tmp2) %in% edge_to_split]
        
        
        edges_to_sample_tmp2 = edges_to_sample_tmp2[-edge_idx_to_split, ]
        
        
        
        # cut the interval at the time the healthy person gets sick
        extra_edges2 = matrix(NA, nrow = 2*nrow(SI_links2), ncol = ncol(edges_to_sample_tmp2))
        extra_edges2[seq(1, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 3], t1 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                             net0 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2))], net1 = 1) %>% as.matrix()
        extra_edges2[seq(2, 2*nrow(SI_links2), 2), ] = cbind(SI_links2[, 1: 2], t0 = infect_times[nei_infec_temp2$sick_per[as.integer(rownames(SI_links2))]]+1e-5, 
                                                             t1 = obs_time[interval_idx2+1],
                                                             net0 = 1, net1 = net_snapshots[as.matrix(cbind(SI_links2[, 1: 2], interval_idx2+1))])  %>% as.matrix()
        extra_edges2 = as.data.frame(extra_edges2)
        colnames(extra_edges2) = colnames(edges_to_sample_tmp2)
        
        rownames(extra_edges2)[seq(1, 2*nrow(SI_links2), 2)] = paste0(edge_to_split, "_1")
        rownames(extra_edges2)[seq(2, 2*nrow(SI_links2), 2)] = paste0(edge_to_split, "_2")
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
    
    edges_to_sample_tmp2= edges_to_sample_tmp2[!rownames(edges_to_sample_tmp2) %in% names(imputed_net_events_A),]
    print(nrow(edges_to_sample_tmp2))
    st = Sys.time()
    imputed_net_events_B = mclapply(1: nrow(edges_to_sample_tmp2), function(i){
      # i = nrow(edges_to_sample_tmp2_B)
      temp_edge = edges_to_sample_tmp2[i, ] %>% as.numeric()
      interval_imputer2(sick_link = temp_edge[1: 2],
                        sampling_intvl = temp_edge[3: 4],
                        epi_times = c(infect_times[temp_edge[1]], recov_times[temp_edge[1]], infect_times[temp_edge[2]], recov_times[temp_edge[2]]),
                        net_status = temp_edge[5: 6],
                        params.cur, force_one = force_one)
    }, mc.cores = 1)
    ed = Sys.time()
    print(ed - st)
    names(imputed_net_events_B) = rownames(edges_to_sample_tmp2)
    
    imputed_net_events_A = imputed_net_events_A[!names(imputed_net_events_A) %in% edge_to_split]
    imputed_net_events = append(imputed_net_events_A, imputed_net_events_B)
    
    events_temp = rbind.data.frame(sick_events, 
                                   do.call("rbind", imputed_net_events_A),
                                   do.call("rbind", imputed_net_events_B)
    )
    events_temp = events_temp[order(events_temp$time), ]
    events_temp = events_temp[-1, ]
    
    
    net_event_rows = which(events_temp$event > 2)
    net_events = events_temp[net_event_rows, ]
    net_events$event = ifelse(net_events$event < 6, 3, 6) +
      (infect_times[net_events$per1] <= net_events$time & net_events$time <= recov_times[net_events$per1]) +
      (infect_times[net_events$per2] <= net_events$time & net_events$time <= recov_times[net_events$per2])
    
    events_temp[net_event_rows, ] = net_events
    
    PA = parse_augment_cpp(G0, I0, events_temp, recover.dat)
    
    event.counts = rep(0, 8)
    names(event.counts) = 1: 8
    event.counts[as.integer(names(PA$event.counts))] = PA$event.counts
    
    big.sums = PA$big.sums
    imp_event_log[, ,it] = rbind(event.counts, big.sums)
    if(sum(big.sums < 0) > 0) stop(paste0("Issue at ", it, "-th iter"))
    
    
    edge_ids_to_fix = sample(edge_to_fix, length(edge_to_fix)*fix_prop, replace = FALSE)
    
    patient_id_to_fix = sample(sick_pers, length(sick_pers)*fix_prop, replace = FALSE)
    recover.datA = recover.dat[recover.dat$per1 %in% patient_id_to_fix, ]
    # we only fix edges unrelated to selected patients' recovery
    edge_ids_to_fix = edge_ids_to_fix[!edge_ids_to_fix %in% sick_edges_vec[!names(sick_edges_vec) %in% patient_id_to_fix]] 
    imputed_net_events_A = imputed_net_events[edge_ids_to_fix]
    
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

