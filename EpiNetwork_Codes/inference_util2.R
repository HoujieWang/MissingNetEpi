
dissect_matrix = function(mat){
  link_info_i = list()
  l = 1; r = 1
  while (r < ncol(mat)) {
    r = r + 1
    if (!is.na(mat[5, r])){
      link_info_i = append(link_info_i, list(mat[, l: r]))
      l = r
    }
  }
  return(link_info_i)
}
regroup_interactions <- function(interaction_list, window_length) {
  # interaction_list = all_edge_info
  grouped <- list()
  
  for (interaction in interaction_list) {
    # interaction = interaction_list[[1]]
    agent_pair <- interaction[[8]]  # [[8]]-th element (agent pair)
    time_value <- interaction[[7]]  # [[7]]-th element (numeric value)
    
    # Extract agent indices
    agentA <- agent_pair[1]
    agentB <- agent_pair[2]
    
    # Ensure consistent ordering of pairs (min, max) for key
    pair_key <- paste(sort(c(agentA, agentB)), collapse = "-")
    
    # Compute time-based grouping key
    time_key <- floor(time_value / window_length)
    
    # Final grouping key
    final_key <- paste(pair_key, time_key, sep = "_")
    
    if (!final_key %in% names(grouped)) {
      grouped[[final_key]] <- list()
    }
    
    # Append interaction
    grouped[[final_key]] <- append(grouped[[final_key]], list(interaction))
  }
  
  return(grouped)
}
SI_events_imputer3 = function(sick_order2,
                              nei_infec2,
                              recov_times2,
                              infect_times2,
                              net_status_infect2,
                              window_length2,
                              obs_time2,
                              net_snapshots2,
                              interact_report_ind2,
                              params.cur2){
  # sick_order2 = sick_order
  # nei_infec2 = nei_infec
  # recov_times2 = recov_times
  # infect_times2 = infect_times
  # net_status_infect2 = net_status_infect
  # window_length2 = window_length
  # obs_time2 = obs_time
  # net_snapshots2 = net_snapshots
  # interact_report_ind2 = interact_report_ind
  # params.cur2 = params.cur
  
  last_SI_events_and_edge_info = lapply(sick_order2[-1], function(i){
    # for (i in as.integer(names(nei_infec))) {
    # print(i)
    # i = 10
    # set.seed(10001)
    x = nei_infec2[[as.character(i)]]
    non_nei = x[recov_times2[x] <= infect_times2[i]]
    x = x[recov_times2[x] > infect_times2[i]]
    
    net_status_ik = rep(0, length(x))
    net_status_ik[x %in% net_status_infect2[[as.character(i)]]] = 1
    
    # this is the network status of (j ,k)
    net_status_ij = as.integer(sapply(as.character(x), function(xx){i %in% net_status_infect2[[xx]]}))
    
    intvl_idx = (infect_times2[x] %/% window_length2) + 1
    
    # print(paste0('x is ', x))
    out_t0ij = cbind.data.frame(per1 = x, per2 = i,
                                t0 = obs_time2[intvl_idx],
                                t1 = infect_times2[x],
                                net0 = net_snapshots2[cbind(x, i, intvl_idx)], 
                                net1 = net_status_ij,
                                force1 = interact_report_ind2[cbind(x, i, intvl_idx)])
    # print('out_t0ij')
    out_t0ij$t1[out_t0ij$t1-out_t0ij$t0 == 0] = out_t0ij$t0[out_t0ij$t1-out_t0ij$t0 == 0] + 1e-5
    
    out_ikt1 = cbind.data.frame(per1 = x, per2 = i,
                                t0 = infect_times2[i],
                                t1 = obs_time2[(infect_times2[i] %/% window_length2) + 2],
                                net0 = net_status_ik, 
                                net1 = net_snapshots2[cbind(x, i, (infect_times2[i] %/% window_length2) + 2)],
                                force1 = interact_report_ind2[cbind(x, i, (infect_times2[i] %/% window_length2) + 1)])
    # print('out_ikt1')
    out_ijik = cbind.data.frame(per1 = x, per2 = i, 
                                t0 = infect_times2[x], 
                                t1 = infect_times2[i],
                                net0 = net_status_ij, 
                                net1 = net_status_ik,
                                force1 = interact_report_ind2[cbind(x, i, (infect_times2[i] %/% window_length2) + 1)])
    # print('out_ijik')
    edge_info_nei = get_SI_link_info(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2)
    # if (nrow(out_t0ij) > 0){
    #   # events_SI_sick_nei= impute_SI_links2(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2, relax = relax)
    #   # 
    #   # rbenchmark::benchmark(impute_SI_links2(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2, relax = relax),
    #   #                       get_SI_link_info(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2))
    #   
    #   edge_info_nei = get_SI_link_info(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2)
    # }
    # 
    # else{
    #   events_SI_sick_nei = data.frame()
    # }
    
    
    intvl_idx_k = infect_times2[i] %/% window_length2+1
    edge_info_nonnei = list()
    for (j in non_nei){
      # j = 2
      # print(c(j))
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
        # jj = 1
        if (interact_report_ind2[j, i, jj] == 1){
          tmp = link_type_changes[, which(link_type_changes[4, ] == jj): which(link_type_changes[4, ] == (jj+1))]
          
          link_info_i = dissect_matrix(tmp)
          
          
          
          # set.seed(1)
          edge_info_nonnei = append(edge_info_nonnei, 
                                    (lapply(link_info_i, function(tmp_x){
                                      # tmp_x = link_info_i[[1]]
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
    while (TRUE) {
      
      events = do.call('rbind', lapply(x, function(xx){
        tmp_events = simulate_ctmc3_cpp(xx[[1]], xx[[2]], xx[[3]], xx[[4]], xx[[5]], xx[[6]])
        if (nrow(tmp_events) > 0) tmp_events$time = tmp_events$time + xx[[7]]
        tmp_events
      }))
      if (nrow(events) > 0) break
    }
    out = cbind.data.frame(events, per1 = x[[1]][[8]][1], per2 = x[[1]][[8]][2])
    # out$event = ifelse(out$event == 1, 3, 6) +
    #   (infect_times2[out$per1] <= out$time &
    #      out$time < recov_times2[out$per1]) +
    #   (infect_times2[out$per2] <= out$time &
    #      out$time < recov_times2[out$per2])
    out
  }))
  
  # tmp = run_simulations_cpp(all_edge_info)
  tmp$event = ifelse(tmp$event == 1, 3, 6) + 
    (infect_times2[tmp$per1] <= tmp$time &
       tmp$time < recov_times2[tmp$per1]) +
    (infect_times2[tmp$per2] <= tmp$time &
       tmp$time < recov_times2[tmp$per2])
  
  out_data  = rbind.data.frame(last_SI_events, tmp)
  out_data[order(out_data$time),]
}

SI_events_imputer3_no_sick = function(sick_agents2,
                                      healthy_agents2,
                                      sick_range2,
                                      nei_infec2,
                                      recov_times2,
                                      infect_times2,
                                      window_length2,
                                      obs_time2,
                                      net_snapshots2,
                                      interact_report_ind2,
                                      params.cur2){
  
  # sick_agents2 = sick_order
  # healthy_agents2 = healthy_agents
  # sick_range2 = sick_range
  # nei_infec2 = nei_infec
  # recov_times2 = recov_times
  # infect_times2 = infect_times
  # window_length2 = window_length
  # obs_time2 = obs_time
  # net_snapshots2 = net_snapshots
  # interact_report_ind2 = interact_report_ind
  # params.cur2 = params.cur
  
  all_edge_info = do.call('c', lapply(healthy_agents2, function(i){
    # for (i in as.integer(names(nei_infec))) {
    # i = 3
    # set.seed(10001)
    edge_info = list()
    for (j in sick_agents2){
      # j = 51
      # intvl_idx_k = recov_times2[j] %/% window_length2+1
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
    while (TRUE) {
      
      events = do.call('rbind', lapply(x, function(xx){
        tmp_events = simulate_ctmc3_cpp(xx[[1]], xx[[2]], xx[[3]], xx[[4]], xx[[5]], xx[[6]])
        if (nrow(tmp_events) > 0) tmp_events$time = tmp_events$time + xx[[7]]
        tmp_events
      }))
      if (nrow(events) > 0) break
    }
    out = cbind.data.frame(events, per1 = x[[1]][[8]][1], per2 = x[[1]][[8]][2])
    # out$event = ifelse(out$event == 1, 3, 6) +
    #   (infect_times2[out$per1] <= out$time &
    #      out$time < recov_times2[out$per1]) +
    #   (infect_times2[out$per2] <= out$time &
    #      out$time < recov_times2[out$per2])
    out
  }))
  # tmp = run_simulations_cpp(all_edge_info)
  tmp$event = ifelse(tmp$event == 1, 3, 6) + 
    (infect_times2[tmp$per1] <= tmp$time &
       tmp$time < recov_times2[tmp$per1]) +
    (infect_times2[tmp$per2] <= tmp$time &
       tmp$time < recov_times2[tmp$per2])
  
  out_data  = tmp[order(tmp$time),]
  out_data
}

SI_events_imputer2_no_sick = function(sick_agents2,
                                      healthy_agents2,
                                      sick_range2,
                                      nei_infec2,
                                      recov_times2,
                                      infect_times2,
                                      window_length2,
                                      obs_time2,
                                      net_snapshots2,
                                      interact_report_ind2,
                                      params.cur2){
  
  # sick_agents2 = sick_order
  # healthy_agents2 = healthy_agents
  # sick_range2 = sick_range
  # nei_infec2 = nei_infec
  # recov_times2 = recov_times
  # infect_times2 = infect_times
  # window_length2 = window_length
  # obs_time2 = obs_time
  # net_snapshots2 = net_snapshots
  # interact_report_ind2 = interact_report_ind
  # params.cur2 = params.cur
  
  all_edge_info = do.call('c', lapply(healthy_agents2, function(i){
    # for (i in as.integer(names(nei_infec))) {
    # i = 3
    # set.seed(10001)
    edge_info = list()
    for (j in sick_agents2){
      # j = 51
      # intvl_idx_k = recov_times2[j] %/% window_length2+1
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
  
  tmp = run_simulations_cpp(all_edge_info)
  tmp$event = ifelse(tmp$event == 1, 3, 6) + 
    (infect_times2[tmp$per1] <= tmp$time &
       tmp$time < recov_times2[tmp$per1]) +
    (infect_times2[tmp$per2] <= tmp$time &
       tmp$time < recov_times2[tmp$per2])
  
  out_data  = tmp[order(tmp$time),]
  out_data
}




postSI_events_info <- function(recovery_times, out_ikt1_dat, ii, params) {
  # recovery_times = recov_times2; out_ikt1_dat = out_ikt1_dat; ii = 1
  # Extract recovery times for the individuals involved
  recov_times_ii <- recovery_times[c(out_ikt1_dat$per1[ii], out_ikt1_dat$per2[ii])]
  
  # Determine if recovery times are within the specified interval
  recov_in_intvl <- c(
    (out_ikt1_dat$t0[ii] <= recov_times_ii[1] & recov_times_ii[1] <= out_ikt1_dat$t1[ii]),
    (out_ikt1_dat$t0[ii] <= recov_times_ii[2] & recov_times_ii[2] <= out_ikt1_dat$t1[ii])
  )
  
  # Create a data frame of recovery times and filter those within the interval
  recov_times_ii <- cbind.data.frame(recov_times_ii, recov_in_intvl)
  recov_times_ii <- recov_times_ii$recov_times_ii[as.logical(recov_times_ii$recov_in_intvl)]
  
  # Simulate CTMC events based on the number of recovery times in the interval
  if (length(recov_times_ii) == 0) {
    
    list(params[5], params[8], vector('numeric'), out_ikt1_dat$t1[ii] - out_ikt1_dat$t0[ii], 
         out_ikt1_dat$net0[ii], out_ikt1_dat$net1[ii], 
         out_ikt1_dat$t0[ii], c(out_ikt1_dat$per1[ii], out_ikt1_dat$per2[ii]))
    # postSI_events <- simulate_ctmc2(
    #   params[5], params[8], vector('numeric'), out_ikt1_dat$t1[ii] - out_ikt1_dat$t0[ii], 
    #   out_ikt1_dat$net0[ii], out_ikt1_dat$net1[ii]
    # )
  } else if (length(recov_times_ii) == 1) {
    # postSI_events <- simulate_ctmc2(
    #   params[5:4], params[8:7], c(recov_times_ii) - out_ikt1_dat$t0[ii], 
    #   out_ikt1_dat$t1[ii] - out_ikt1_dat$t0[ii], 
    #   out_ikt1_dat$net0[ii], out_ikt1_dat$net1[ii]
    # )
    list(params[5:4], params[8:7], c(recov_times_ii) - out_ikt1_dat$t0[ii], 
         out_ikt1_dat$t1[ii] - out_ikt1_dat$t0[ii], 
         out_ikt1_dat$net0[ii], out_ikt1_dat$net1[ii], 
         out_ikt1_dat$t0[ii], c(out_ikt1_dat$per1[ii], out_ikt1_dat$per2[ii]))
  } else if (length(recov_times_ii) == 2) {
    # postSI_events <- simulate_ctmc2(
    #   params[5:3], params[8:6], c(recov_times_ii) - out_ikt1_dat$t0[ii], 
    #   out_ikt1_dat$t1[ii] - out_ikt1_dat$t0[ii], 
    #   out_ikt1_dat$net0[ii], out_ikt1_dat$net1[ii]
    # )
    list(
      params[5:3], params[8:6], sort(recov_times_ii) - out_ikt1_dat$t0[ii], 
      out_ikt1_dat$t1[ii] - out_ikt1_dat$t0[ii], 
      out_ikt1_dat$net0[ii], out_ikt1_dat$net1[ii], 
      out_ikt1_dat$t0[ii], c(out_ikt1_dat$per1[ii], out_ikt1_dat$per2[ii]))
    
  }
  
  # # Adjust the time of events and add additional information if events exist
  # if (nrow(postSI_events) > 0) {
  #   postSI_events$time <- postSI_events$time + out_ikt1_dat$t0[ii]
  #   postSI_events$event = ifelse(postSI_events$event == 1, 3, 6) + 
  #     3 - as.integer(cut(postSI_events$time, c(out_ikt1_dat$t0[ii], recov_times_ii, Inf)))
  # }
  # 
  # return(postSI_events)
}

get_SI_link_info <- function(out_t0ij_dat, out_ikt1_dat, out_ijik_dat, recov_times2, params = params.cur, net_snapshots, interact_report_ind) { #relax = TRUE if we do not force an event
  #
  # out_t0ij_dat = out_t0ij; out_ikt1_dat = out_ikt1; out_ijik_dat = out_ijik
  # params = params.cur;recov_times2 = recov_times2
  # Ensure input statuses are vectors
  n_processes <- nrow(out_ijik_dat)
  alpha = params[4]
  omega = params[7]
  beta =  params[1]

  # Initialize states at t_b
  # t_b = out_ijik_dat$t1[1]
  # status_tb = out_ijik_dat$net1
  # status <- out_ijik_dat$net1 # This is the changing slice of network status
  # time <- out_ijik_dat$t1[1]
  # events <- data.frame(time = numeric(0),
  #                      event = character(0),
  #                      per1 = integer(0),
  #                      per2 = integer(0)) # Store transitions
  # time_in_1 <- rep(0, n_processes) # Track time spent in 1 for each process

  # Sample the second-to-last event time backward
  # disconnect_rate <- sum(ifelse(status_tb == 1 & out_ijik_dat$force1 == 1, (omega + beta), 0))
  # connect_rate <- sum(ifelse(status_tb == 0 & out_ijik_dat$force1 == 1, alpha, 0))
  # total_rate <- disconnect_rate + connect_rate

  # snapshot_idx = rowMaxs(cbind(obs_time[out_ijik_dat$t1 %/% window_length+1],
  #                              out_ijik_dat$t0))
  # This is the closest observable network status, which is either A(ij) or A(t'),
  # if t' is the largest observation time point smaller than ik

  status_0 = net_snapshots[cbind(out_ijik_dat$per1, out_ijik_dat$per2, out_ijik_dat$t1[1] %/% window_length+1)]
  status_0[out_ijik_dat$t0 > obs_time[out_ijik_dat$t1 %/% window_length+1]] =
    out_ijik_dat$net0[out_ijik_dat$t0 > obs_time[out_ijik_dat$t1 %/% window_length+1]]

  out_ijik_dat$closest_net0 = status_0
  out_ijik_dat$closest_t0 = out_ijik_dat$t0
  st_of_SI_link_before_t1 = which(out_ijik_dat$t0 < obs_time[out_ijik_dat$t1 %/% window_length+1])
  out_ijik_dat$closest_t0[st_of_SI_link_before_t1] = obs_time[out_ijik_dat$t1 %/% window_length+1][st_of_SI_link_before_t1]


  slice_idx = (out_t0ij_dat$t1 %/% window_length+1)
  event_ind = rep(0, nrow(out_t0ij_dat))
  remaining_idx = 1: nrow(out_t0ij_dat)
  while (length(remaining_idx) > 0) {

    event_ind[remaining_idx] = event_ind[remaining_idx] +
      interact_report_ind[cbind(as.matrix(out_t0ij_dat[remaining_idx,1:2]), slice_idx[remaining_idx])]
    slice_idx = slice_idx + 1
    remaining_idx = remaining_idx[slice_idx <= (out_ikt1_dat$t0[1] %/% window_length+1)]
  }


  edge_info = list()
  for (ii in which(event_ind > 0)) {
    # ii=3
    ij_intvl = obs_time[out_ijik_dat$t0[ii] %/% window_length + 1:2]
    cur_t0 = out_ijik_dat$closest_t0[ii]

    if (out_ijik_dat$t1[ii] < ij_intvl[2]){# This means that both ij and ik are in the same observation interval
      # This while loop determines if we need to resample due to the mismatch of number of events
      if (out_ijik_dat$force1[ii] == 1){


        edge_info = append(edge_info,
                           list(list(params[4], params[7]+params[1], vector('numeric'), out_ijik_dat$t1[ii]-out_ijik_dat$closest_t0[ii],
                                     out_ijik_dat$closest_net0[ii], out_ijik_dat$net1[ii],
                                     out_ijik_dat$closest_t0[ii], c(out_ijik_dat$per1[ii], out_ijik_dat$per2[ii]))))

        edge_info = append(edge_info, list(postSI_events_info(recov_times2, out_ikt1_dat, ii, params)))

        edge_info = append(edge_info,
                           list(list(params[3], params[5], vector('numeric'), out_t0ij_dat$t1[ii]-out_t0ij_dat$t0[ii],
                                     out_t0ij_dat$net0[ii], out_t0ij_dat$net1[ii],
                                     out_t0ij_dat$t0[ii], c(out_ijik_dat$per1[ii], out_ijik_dat$per2[ii]))))

      }

    }else{
      # Check if we need to impute in the interval of ik
      if (out_ijik_dat$force1[ii] == 1){

        edge_info = append(edge_info,
                           list(list(params[4], params[7]+params[1], vector('numeric'), out_ijik_dat$t1[ii]-out_ijik_dat$closest_t0[ii],
                                     out_ijik_dat$closest_net0[ii], out_ijik_dat$net1[ii],
                                     out_ijik_dat$closest_t0[ii], c(out_ijik_dat$per1[ii], out_ijik_dat$per2[ii]))))
        edge_info = append(edge_info, list(postSI_events_info(recov_times2, out_ikt1_dat, ii, params)))

      }

      # Check if we need to impute in the interval of ij
      if (out_t0ij_dat$force1[ii] == 1){

        edge_info = append(edge_info,
                           list(list(params[3], params[5], vector('numeric'), out_t0ij_dat$t1[ii]-out_t0ij_dat$t0[ii],
                                     out_t0ij_dat$net0[ii], out_t0ij_dat$net1[ii],
                                     out_t0ij_dat$t0[ii], c(out_t0ij_dat$per1[ii],out_t0ij_dat$per2[ii]))))

        edge_info = append(edge_info,
                           list(list(params[4], params[7]+params[1], vector('numeric'), obs_time[out_t0ij_dat$t1[ii] %/% window_length + 2]-out_t0ij_dat$t1[ii],
                                     out_t0ij_dat$net1[ii], net_snapshots[out_t0ij_dat$per1[ii], out_t0ij_dat$per2[ii], out_t0ij_dat$t1[ii] %/% window_length + 2],
                                     out_t0ij_dat$t1[ii],c(out_t0ij_dat$per1[ii],out_t0ij_dat$per2[ii]))))


      }

      # the SI links ranges over more than two observation times
      first_SI_end_idx = out_t0ij_dat$t1[ii] %/% window_length + 2
      first_SI_end_time = obs_time[first_SI_end_idx]
      last_SI_st_time = obs_time[out_ijik_dat$t1[ii] %/% window_length + 1]
      # middle_events = data.frame(time=numeric(0), event=integer(0),per1 = integer(0),per2 = integer(0))
      if(first_SI_end_time < last_SI_st_time){
        num_intvl_run_over = (last_SI_st_time - first_SI_end_time) / window_length
        for (iii in 0: (num_intvl_run_over-1)) {
          # iii = 0
          # if there is event
          if (interact_report_ind[out_t0ij_dat$per1[ii], out_t0ij_dat$per2[ii], first_SI_end_idx+iii]){
            first_status =  net_snapshots[out_t0ij_dat$per1[ii], out_t0ij_dat$per2[ii], first_SI_end_idx+iii]
            end_status =  net_snapshots[out_t0ij_dat$per1[ii], out_t0ij_dat$per2[ii], first_SI_end_idx+iii+1]

            edge_info = append(edge_info,
                               list(list(params[4], params[1]+params[7], vector('numeric'), window_length, first_status, end_status,
                                         (first_SI_end_idx+iii-1)*window_length, c(out_t0ij_dat$per1[ii], out_t0ij_dat$per2[ii]))))

          }
        }
      }

      # events = rbind.data.frame(events, preSI_SI_events, middle_events, SI_postSI_events)
    }
  }
  list(events = data.frame(time = numeric(0),
                           event = character(0),
                           per1 = integer(0),
                           per2 = integer(0)),
       edge_info = edge_info)

}
SI_events_imputer2 = function(sick_order2,
                              nei_infec2,
                              recov_times2,
                              infect_times2,
                              net_status_infect2,
                              window_length2,
                              obs_time2,
                              net_snapshots2,
                              interact_report_ind2,
                              params.cur2){
  # sick_order2 = sick_order
  # nei_infec2 = nei_infec
  # recov_times2 = recov_times
  # infect_times2 = infect_times
  # net_status_infect2 = net_status_infect
  # window_length2 = window_length
  # obs_time2 = obs_time
  # net_snapshots2 = net_snapshots
  # interact_report_ind2 = interact_report_ind
  # params.cur2 = params.cur
  
  last_SI_events_and_edge_info = lapply(sick_order2[-1], function(i){
    # for (i in as.integer(names(nei_infec))) {
    # print(i)
    # i = 2
    # set.seed(10001)
    x = nei_infec2[[as.character(i)]]
    non_nei = x[recov_times2[x] <= infect_times2[i]]
    x = x[recov_times2[x] > infect_times2[i]]
    
    net_status_ik = rep(0, length(x))
    net_status_ik[x %in% net_status_infect2[[as.character(i)]]] = 1
    
    # this is the network status of (j ,k)
    net_status_ij = as.integer(sapply(as.character(x), function(xx){i %in% net_status_infect2[[xx]]}))
    
    intvl_idx = (infect_times2[x] %/% window_length2) + 1
    
    # print(paste0('x is ', x))
    out_t0ij = cbind.data.frame(per1 = x, per2 = i,
                                t0 = obs_time2[intvl_idx],
                                t1 = infect_times2[x],
                                net0 = net_snapshots2[cbind(x, i, intvl_idx)], 
                                net1 = net_status_ij,
                                force1 = interact_report_ind2[cbind(x, i, intvl_idx)])
    # print('out_t0ij')
    out_t0ij$t1[out_t0ij$t1-out_t0ij$t0 == 0] = out_t0ij$t0[out_t0ij$t1-out_t0ij$t0 == 0] + 1e-5
    
    out_ikt1 = cbind.data.frame(per1 = x, per2 = i,
                                t0 = infect_times2[i],
                                t1 = obs_time2[(infect_times2[i] %/% window_length2) + 2],
                                net0 = net_status_ik, 
                                net1 = net_snapshots2[cbind(x, i, (infect_times2[i] %/% window_length2) + 2)],
                                force1 = interact_report_ind2[cbind(x, i, (infect_times2[i] %/% window_length2) + 1)])
    # print('out_ikt1')
    out_ijik = cbind.data.frame(per1 = x, per2 = i, 
                                t0 = infect_times2[x], 
                                t1 = infect_times2[i],
                                net0 = net_status_ij, 
                                net1 = net_status_ik,
                                force1 = interact_report_ind2[cbind(x, i, (infect_times2[i] %/% window_length2) + 1)])
    # print('out_ijik')
    edge_info_nei = get_SI_link_info(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2)
    # if (nrow(out_t0ij) > 0){
    #   # events_SI_sick_nei= impute_SI_links2(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2, relax = relax)
    #   # 
    #   # rbenchmark::benchmark(impute_SI_links2(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2, relax = relax),
    #   #                       get_SI_link_info(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2))
    #   
    #   edge_info_nei = get_SI_link_info(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2)
    # }
    # 
    # else{
    #   events_SI_sick_nei = data.frame()
    # }
    intvl_idx_k = infect_times2[i] %/% window_length2+1
    edge_info_nonnei = list()
    for (j in non_nei){
      # j = 124
      # print(c(j))
      ij = infect_times2[j]; rj = recov_times2[j]
      intvl_idx_j = ij %/% window_length2+1
      
      t0 = obs_time2[intvl_idx_j]
      
      link_type_changes = rbind(c(t0,ij, rj, infect_times2[i], recov_times2[i]), 
                                c(0, 1, -1, 1, -1),
                                c(0, 1, -1, 0, 0),
                                c(intvl_idx_j, rep(NA, 4)))
      
      link_type_changes = cbind(link_type_changes, rbind(obs_time2[(intvl_idx_j+1):(intvl_idx_k+1)], NA, NA, (intvl_idx_j+1):(intvl_idx_k+1)))
      link_type_changes = link_type_changes[, order(link_type_changes[1,])]
      
      link_type_changes = link_type_changes[, link_type_changes[1,] <= obs_time2[(intvl_idx_k+1)]]
      link_type_changes[2, is.na(link_type_changes[2, ])] = 0 
      link_type_changes[2, ] = cumsum(link_type_changes[2,])
      
      link_type_changes[3, is.na(link_type_changes[3, ])] = 0
      link_type_changes[3, ] = cumsum(link_type_changes[3,])
      for (jj in (intvl_idx_j): (intvl_idx_k)) {
        # jj = 3
        if (interact_report_ind2[j, i, jj] == 1){
          tmp = link_type_changes[, which(link_type_changes[4, ] == jj): which(link_type_changes[4, ] == (jj+1))]
          # set.seed(1)
          edge_info_nonnei = append(edge_info_nonnei, 
                                    list(list(params.cur2[tmp[2,-ncol(tmp)]+3], params.cur2[tmp[2,-ncol(tmp)]+6] + tmp[3,-ncol(tmp)]*params.cur2[1],
                                              tmp[1,-c(1,ncol(tmp))] - tmp[1, 1], window_length2, 
                                              net_snapshots2[j,i,jj],
                                              net_snapshots2[j,i,jj+1], 
                                              tmp[1, 1], c(j, i))))
          
        }
      }
    }
    
    
    out = append(edge_info_nei$edge_info, edge_info_nonnei)
    return(list(edge_info_nei$events, out))
  })
  
  last_SI_events = do.call('rbind.data.frame', lapply(last_SI_events_and_edge_info, function(x){x[[1]]}))
  all_edge_info =  do.call('c', lapply(last_SI_events_and_edge_info, function(x){x[[2]]}))
  
  
  tmp = run_simulations_cpp(all_edge_info)
  tmp$event = ifelse(tmp$event == 1, 3, 6) + 
    (infect_times2[tmp$per1] <= tmp$time &
       tmp$time < recov_times2[tmp$per1]) +
    (infect_times2[tmp$per2] <= tmp$time &
       tmp$time < recov_times2[tmp$per2])
  
  out_data  = rbind.data.frame(last_SI_events, tmp)
  out_data[order(out_data$time),]
}
SI_events_imputer = function(sick_order2,
                             nei_infec2,
                             recov_times2,
                             infect_times2,
                             net_status_infect2,
                             window_length2,
                             obs_time2,
                             net_snapshots2,
                             interact_report_ind2,
                             params.cur2, 
                             relax = TRUE){
  # sick_order2 = sick_order
  # nei_infec2 = nei_infec
  # recov_times2 = recov_times
  # infect_times2 = infect_times
  # net_status_infect2 = net_status_infect
  # window_length2 = window_length
  # obs_time2 = obs_time
  # net_snapshots2 = net_snapshots
  # interact_report_ind2 = interact_report_ind
  # params.cur2 = params.cur
  
  do.call('rbind', lapply(sick_order2[-1], function(i){
    # for (i in as.integer(names(nei_infec))) {
    # print(i)
    # i = 346
    # set.seed(10001)
    x = nei_infec2[[as.character(i)]]
    non_nei = x[recov_times2[x] <= infect_times2[i]]
    x = x[recov_times2[x] > infect_times2[i]]
    
    net_status_ik = rep(0, length(x))
    net_status_ik[x %in% net_status_infect2[[as.character(i)]]] = 1
    
    # this is the network status of (j ,k)
    net_status_ij = as.integer(sapply(as.character(x), function(xx){i %in% net_status_infect2[[xx]]}))
    
    intvl_idx = (infect_times2[x] %/% window_length2) + 1
    
    # print(paste0('x is ', x))
    out_t0ij = cbind.data.frame(per1 = x, per2 = i,
                                t0 = obs_time2[intvl_idx],
                                t1 = infect_times2[x],
                                net0 = net_snapshots2[cbind(x, i, intvl_idx)], 
                                net1 = net_status_ij,
                                force1 = interact_report_ind2[cbind(x, i, intvl_idx)])
    # print('out_t0ij')
    out_t0ij$t1[out_t0ij$t1-out_t0ij$t0 == 0] = out_t0ij$t0[out_t0ij$t1-out_t0ij$t0 == 0] + 1e-5
    
    out_ikt1 = cbind.data.frame(per1 = x, per2 = i,
                                t0 = infect_times2[i],
                                t1 = obs_time2[(infect_times2[i] %/% window_length2) + 2],
                                net0 = net_status_ik, 
                                net1 = net_snapshots2[cbind(x, i, (infect_times2[i] %/% window_length2) + 2)],
                                force1 = interact_report_ind2[cbind(x, i, (infect_times2[i] %/% window_length2) + 1)])
    # print('out_ikt1')
    out_ijik = cbind.data.frame(per1 = x, per2 = i, 
                                t0 = infect_times2[x], 
                                t1 = infect_times2[i],
                                net0 = net_status_ij, 
                                net1 = net_status_ik,
                                force1 = interact_report_ind2[cbind(x, i, (infect_times2[i] %/% window_length2) + 1)])
    # print('out_ijik')
    
    if (nrow(out_t0ij) > 0){
      # set.seed(100)
      # impute_SI_links(out_t0ij, out_ikt1, out_ijik, params.cur, relax = TRUE)
      # 
      events_SI_sick_nei= impute_SI_links2(out_t0ij, out_ikt1, out_ijik, recov_times2, params.cur2, net_snapshots2, interact_report_ind2, relax = relax)
    }else{
      events_SI_sick_nei = data.frame()
    }
    # print('SI')
    intvl_idx_k = infect_times2[i] %/% window_length2+1
    if (length(non_nei) > 0){
      events_SI_nonsick_nei = do.call('rbind', lapply(non_nei, function(j){
        # j = 337
        # print(c(j))
        ij = infect_times2[j]; rj = recov_times2[j]
        intvl_idx_j = ij %/% window_length2+1
        
        t0 = obs_time2[intvl_idx_j]
        # t1 = obs_time[intvl_idx_j+1]
        
        
        # link_type_changes = rbind(c(t0, ij, rj, infect_times2[i], recov_times2[i]), 
        #                           c(0, 1, -1, 1, -1),
        #                           c(0, 1, 0, 0, 0),
        #                           c(intvl_idx_j, rep(NA, 4)))
        # 
        # link_type_changes = cbind(link_type_changes, rbind(obs_time2[(intvl_idx_j+1):(intvl_idx_k+1)], NA, NA, (intvl_idx_j+1):(intvl_idx_k+1)))
        # link_type_changes = link_type_changes[, order(link_type_changes[1,])]
        # 
        # link_type_changes = link_type_changes[, link_type_changes[1,] <= obs_time2[(intvl_idx_k+1)]]
        # link_type_changes[2, is.na(link_type_changes[2, ])] = 0 
        # link_type_changes[2, ] = cumsum(link_type_changes[2,])
        # link_type_changes[3, which(is.na(link_type_changes[3,]))] = link_type_changes[3, which(is.na(link_type_changes[3,]))-1]
        link_type_changes = rbind(c(t0,ij, rj, infect_times2[i], recov_times2[i]), 
                                  c(0, 1, -1, 1, -1),
                                  c(0, 1, -1, 0, 0),
                                  c(intvl_idx_j, rep(NA, 4)))
        
        link_type_changes = cbind(link_type_changes, rbind(obs_time2[(intvl_idx_j+1):(intvl_idx_k+1)], NA, NA, (intvl_idx_j+1):(intvl_idx_k+1)))
        link_type_changes = link_type_changes[, order(link_type_changes[1,])]
        
        link_type_changes = link_type_changes[, link_type_changes[1,] <= obs_time2[(intvl_idx_k+1)]]
        link_type_changes[2, is.na(link_type_changes[2, ])] = 0 
        link_type_changes[2, ] = cumsum(link_type_changes[2,])
        
        link_type_changes[3, is.na(link_type_changes[3, ])] = 0
        link_type_changes[3, ] = cumsum(link_type_changes[3,])
        
        # link_type_changes[3, which(is.na(link_type_changes[3,]))] = link_type_changes[3, which(is.na(link_type_changes[3,]))-1]
        
        events_j = data.frame()
        for (jj in (intvl_idx_j): (intvl_idx_k)) {
          # jj = 3
          if (interact_report_ind2[j, i, jj] == 1){
            tmp = link_type_changes[, which(link_type_changes[4, ] == jj): which(link_type_changes[4, ] == (jj+1))]
            # set.seed(1)
            tmp_events = simulate_ctmc2(params.cur2[tmp[2,-ncol(tmp)]+3], params.cur2[tmp[2,-ncol(tmp)]+6] + tmp[3,-ncol(tmp)]*params.cur2[1],
                                        tmp[1,-c(1,ncol(tmp))] - tmp[1, 1], window_length2, 
                                        net_snapshots2[j,i,jj],
                                        net_snapshots2[j,i,jj+1])
            tmp_events$time = tmp_events$time + tmp[1, 1]
            events_j = rbind.data.frame(events_j,tmp_events)
            
          }
          
        }
        if (nrow(events_j) > 0){
          events_j$event = ifelse(events_j$event == 1, 3, 6) + 
            (ij <= events_j$time & events_j$time <= rj) + 
            (infect_times2[i] <= events_j$time & events_j$time <= recov_times2[i])
          events_j$per1 = j
          events_j$per2 = i
          events_j 
        }else{
          NULL
        }
      }))
    }else{
      events_SI_nonsick_nei = data.frame()
    }
    return(rbind.data.frame(events_SI_sick_nei,events_SI_nonsick_nei))
  }))
  
}
impute_SI_links2 <- function(out_t0ij_dat, out_ikt1_dat, out_ijik_dat, recov_times2, params = params.cur, net_snapshots, interact_report_ind, relax = TRUE) { #relax = TRUE if we do not force an event
  # 
  # out_t0ij_dat = out_t0ij; out_ikt1_dat = out_ikt1; out_ijik_dat = out_ijik
  # params = params.cur;recov_times2 = recov_times2
  # Ensure input statuses are vectors
  n_processes <- nrow(out_ijik_dat)
  alpha = params[4]
  omega = params[7]
  beta =  params[1]
  
  
  
  # Initialize states at t_b
  t_b = out_ijik_dat$t1[1]
  status_tb = out_ijik_dat$net1
  status <- out_ijik_dat$net1 # This is the changing slice of network status
  time <- out_ijik_dat$t1[1]
  events <- data.frame(time = numeric(0), 
                       event = character(0), 
                       per1 = integer(0), 
                       per2 = integer(0)) # Store transitions
  # time_in_1 <- rep(0, n_processes) # Track time spent in 1 for each process
  
  # Sample the second-to-last event time backward
  disconnect_rate <- sum(ifelse(status_tb == 1 & out_ijik_dat$force1 == 1, (omega + beta), 0))
  connect_rate <- sum(ifelse(status_tb == 0 & out_ijik_dat$force1 == 1, alpha, 0))
  total_rate <- disconnect_rate + connect_rate
  
  snapshot_idx = rowMaxs(cbind(obs_time[out_ijik_dat$t1 %/% window_length+1], 
                               out_ijik_dat$t0))
  # This is the closest observable network status, which is either A(ij) or A(t'),
  # if t' is the largest observation time point smaller than ik
  
  status_0 = net_snapshots[cbind(out_ijik_dat$per1, out_ijik_dat$per2, out_ijik_dat$t1[1] %/% window_length+1)]
  status_0[out_ijik_dat$t0 > obs_time[out_ijik_dat$t1 %/% window_length+1]] = 
    out_ijik_dat$net0[out_ijik_dat$t0 > obs_time[out_ijik_dat$t1 %/% window_length+1]]
  
  out_ijik_dat$closest_net0 = status_0
  out_ijik_dat$closest_t0 = out_ijik_dat$t0
  st_of_SI_link_before_t1 = which(out_ijik_dat$t0 < obs_time[out_ijik_dat$t1 %/% window_length+1])
  out_ijik_dat$closest_t0[st_of_SI_link_before_t1] = obs_time[out_ijik_dat$t1 %/% window_length+1][st_of_SI_link_before_t1]
  
  
  # set.seed(10001)
  if (any(c(out_t0ij_dat$force1, out_ijik_dat$force1, out_ikt1_dat$force1) == 1)){ # If at least one edge needs imputation
    
    if (sum(out_ijik_dat$force1) == 0){
      out_ijik_dat$t1 = out_ijik_dat$closest_t0
      second_last_id = -1
    }else{
      lb = max(out_ijik_dat$closest_t0[out_ijik_dat$force1 == 1])
      t_second_last <- t_b - rexp(1, total_rate)
      if (t_second_last < lb & all(status[out_ijik_dat$force1 == 1] == status_0[out_ijik_dat$force1 == 1])){
        out_ijik_dat$t1 = lb
        second_last_id = -1
      }else{
        if (t_second_last < lb & !all(status[out_ijik_dat$force1 == 1] == status_0[out_ijik_dat$force1 == 1])){
          t_second_last <- t_b - sam_texp(1, total_rate, 0, t_b - lb)
        }
        event_type <- sample(0: 1, 1,prob = c(disconnect_rate, connect_rate) / total_rate)
        
        target_indices = which(status == 1-event_type & out_ijik_dat$force1 == 1)
        if (length(target_indices) == 1){
          process_idx = target_indices
        }else{
          process_idx <- sample(target_indices, 1)
        }
        status[process_idx] <- event_type
        
        events <- rbind(events, data.frame(time = t_second_last, event = 4 + 3*event_type, 
                                           per1 = out_ijik_dat$per1[process_idx], 
                                           per2 = out_ijik_dat$per2[process_idx]))
        out_ijik_dat$t1 = t_second_last
        second_last_id = process_idx
        # out_ijik_dat$net1 = status[process_idx]
      }
    }
    
    
    
    for (ii in 1:n_processes) {
      # ii = 2
      # set.seed(1)
      # print(paste0('it stops at ', ii,' sub link'))
      ij_intvl = obs_time[out_ijik_dat$t0[ii] %/% window_length + 1:2]
      cur_t0 = out_ijik_dat$closest_t0[ii]
      
      if (out_ijik_dat$t1[ii] < ij_intvl[2]){# This means that both ij and ik are in the same observation interval
        # This while loop determines if we need to resample due to the mismatch of number of events
        if (out_ijik_dat$force1[ii] == 1){
          while (TRUE) {
            if (out_ijik_dat$t1[ii] -out_ijik_dat$t0[ii] > 0){
              process_events = backward_sampler(state_init = status[ii], 
                                                out_ijik_dat = out_ijik_dat,
                                                ii,
                                                process_time_init = t_second_last,
                                                omega = params[7],
                                                alpha = params[4],
                                                beta = params[1])
            }else{
              process_events = data.frame()
            }
            
            postSI_events = postSI_events_sampler(recov_times2, out_ikt1_dat, ii, params)
            
            
            preSI_events = simulate_ctmc2(params[3], params[5], vector('numeric'), out_t0ij_dat$t1[ii]-out_t0ij_dat$t0[ii], 
                                          out_t0ij_dat$net0[ii], out_t0ij_dat$net1[ii])
            if(nrow(preSI_events) > 0){
              preSI_events$time = preSI_events$time + out_t0ij_dat$t0[ii]
              preSI_events$event = ifelse(preSI_events$event == 1, 3, 6)
            }
            
            preSI_SI_postSI_events = rbind.data.frame(preSI_events, process_events, postSI_events)
            if (nrow(preSI_SI_postSI_events) > 0){
              preSI_SI_postSI_events = cbind.data.frame(preSI_SI_postSI_events, 
                                                        per1 = out_ikt1_dat$per1[ii],
                                                        per2 = out_ikt1_dat$per2[ii])
              # print(preSI_SI_postSI_events)
              events = rbind.data.frame(events, preSI_SI_postSI_events)
              break
            }
            # If the second last event happens at this edge
            if (second_last_id == ii) break
            
            # Here we allow to break with nrow(preSI_SI_postSI_events) == 0
            if (relax) break
            
          }
        }
        # else{
        #   preSI_SI_postSI_events = data.frame(time=numeric(0), event=integer(0),per1 = integer(0),per2 = integer(0))
        # }
        # events = rbind.data.frame(events, preSI_SI_postSI_events)
      }else{
        # Check if we need to impute in the interval of ik
        if (out_ijik_dat$force1[ii] == 1){
          while (TRUE) {
            process_events = backward_sampler(status[ii], 
                                              out_ijik_dat = out_ijik_dat,
                                              ii,
                                              process_time = t_second_last,
                                              omega = params[7],
                                              alpha = params[4],
                                              beta = params[1])
            postSI_events = postSI_events_sampler(recov_times2, out_ikt1_dat, ii, params)
            SI_postSI_events = rbind.data.frame(process_events, postSI_events)
            
            
            if (nrow(SI_postSI_events) > 0){
              SI_postSI_events = cbind.data.frame(SI_postSI_events, 
                                                  per1 = out_ikt1_dat$per1[ii],
                                                  per2 = out_ikt1_dat$per2[ii])
              break
            }
            if (second_last_id == ii) break
            if (relax){
              SI_postSI_events = data.frame(time=numeric(0), event=integer(0),per1 = integer(0),per2 = integer(0))
              break
            }
          }
        }else{
          SI_postSI_events = data.frame(time=numeric(0), event=integer(0),per1 = integer(0),per2 = integer(0))
        }
        
        # Check if we need to impute in the interval of ij
        if (out_t0ij_dat$force1[ii] == 1){
          while (TRUE) {
            preSI_events = simulate_ctmc2(params[3], params[5], vector('numeric'), out_t0ij_dat$t1[ii]-out_t0ij_dat$t0[ii], 
                                          out_t0ij_dat$net0[ii], out_t0ij_dat$net1[ii])
            if(nrow(preSI_events) > 0){
              preSI_events$time = preSI_events$time + out_t0ij_dat$t0[ii]
              preSI_events$event = ifelse(preSI_events$event == 1, 3, 6)
            }else{
              preSI_events = data.frame(time=numeric(0), event=integer(0),per1 = integer(0),per2 = integer(0))
            }
            
            preSI_events2 = simulate_ctmc2(params[4], params[7]+params[1], vector('numeric'), obs_time[out_t0ij_dat$t1[ii] %/% window_length + 2]-out_t0ij_dat$t1[ii], 
                                           out_t0ij_dat$net1[ii], net_snapshots[out_t0ij_dat$per1[ii], out_t0ij_dat$per2[ii], out_t0ij_dat$t1[ii] %/% window_length + 2])
            
            if (nrow(preSI_events2) > 0){
              preSI_events2$time = preSI_events2$time + out_t0ij_dat$t1[ii]
              preSI_events2$event = ifelse(preSI_events2$event == 1, 4, 7)
            }else{
              preSI_events2 = data.frame(time=numeric(0), event=integer(0),per1 = integer(0),per2 = integer(0))
            }
            
            preSI_SI_events = rbind.data.frame(preSI_events, preSI_events2)
            if (nrow(preSI_SI_events) > 0){
              preSI_SI_events = cbind.data.frame(preSI_SI_events,
                                                 per1 = out_ikt1_dat$per1[ii],
                                                 per2 = out_ikt1_dat$per2[ii])
              break
            }
            
            if (relax) break
          }
        }else{
          preSI_SI_events = data.frame(time=numeric(0), event=integer(0),per1 = integer(0),per2 = integer(0))
        }
        
        # the SI links ranges over more than two observation intervals
        first_SI_end_idx = out_t0ij_dat$t1[ii] %/% window_length + 2
        first_SI_end_time = obs_time[first_SI_end_idx]
        last_SI_st_time = obs_time[out_ijik_dat$t1[ii] %/% window_length + 1]
        middle_events = data.frame(time=numeric(0), event=integer(0),per1 = integer(0),per2 = integer(0))
        if(first_SI_end_time < last_SI_st_time){
          num_intvl_run_over = (last_SI_st_time - first_SI_end_time) / window_length
          for (iii in 0: (num_intvl_run_over-1)) {
            # iii = 0
            # if there is event
            if (interact_report_ind[out_t0ij_dat$per1[ii], out_t0ij_dat$per2[ii], first_SI_end_idx+iii]){
              first_status =  net_snapshots[out_t0ij_dat$per1[ii], out_t0ij_dat$per2[ii], first_SI_end_idx+iii]
              end_status =  net_snapshots[out_t0ij_dat$per1[ii], out_t0ij_dat$per2[ii], first_SI_end_idx+iii+1]
              while(TRUE){
                tmp_events = simulate_ctmc2(params[4], params[1]+params[7], vector('numeric'), window_length, first_status, end_status)
                tmp_events$event = ifelse(tmp_events$event == 1, 4, 7)
                if (nrow(tmp_events) > 0){
                  tmp_events = cbind.data.frame(tmp_events, 
                                                per1 = out_t0ij_dat$per1[ii],
                                                per2 = out_t0ij_dat$per2[ii])
                  middle_events = rbind.data.frame(middle_events, tmp_events)
                  break
                }
                if (relax){
                  break
                }
              }
            }
          }
        }
        
        events = rbind.data.frame(events, preSI_SI_events, middle_events, SI_postSI_events)
      }
    }
    # If all processes meet the start condition, return the results
    # total_time_in_1 <- sum(time_in_1)
    # Flip the events
    events <- events[order(events$time), ]
    events$time <- events$time
    return(events)
    
  }else{
    return(data.frame())
  }
}

backward_sampler <- function(state_init, out_ijik_dat, ii, process_time_init, omega, alpha, beta) {
  # state_init = status[ii]; out_ijik_dat; process_time_init = t_second_last
  while (TRUE) {
    # Initialize state and time for individual process
    process_events <- data.frame(time = numeric(0), 
                                 event = character(0),
                                 per1 = integer(0),
                                 per2 = integer(0))
    # time_in_1_process <- 0
    state = state_init
    process_time = process_time_init
    # This while loop keeps going backward until the end of the edge
    texp_sample = FALSE
    while (TRUE) {
      rate <- if (state == 1) omega + beta else alpha
      if (!texp_sample & state != out_ijik_dat$closest_net0[ii]){
        dt <- sam_texp(1, rate, 0, process_time - out_ijik_dat$closest_t0[ii])
        texp_sample = TRUE
      }else{
        dt <- rexp(1, rate)
      }
      process_time <- process_time - dt
      
      if (process_time <= out_ijik_dat$closest_t0[ii]) {
        # Add remaining time if the process ends at 0
        # time_in_1_process <- time_in_1_process + (t_second_last - process_time) * (state == 1)
        break
      }
      
      # Record event
      # time_in_1_process <- time_in_1_process + dt * (state == 1)
      state <- 1 - state
      
      process_events <- rbind(process_events, data.frame(time = process_time, 
                                                         event = ifelse(state == 1, 7, 4))
                              # ,
                              # per1 = out_ijik_dat$per1[ii],
                              # per2 = out_ijik_dat$per2[ii])
      )
    }
    
    # Check if the trajectory matches the initial state
    if (state == out_ijik_dat$closest_net0[ii]) {
      break
    }
  }
  return(process_events)
  
}



postSI_events_sampler <- function(recovery_times, out_ikt1_dat, ii, params) {
  # recovery_times = recov_times2; out_ikt1_dat = out_ikt1_dat; ii = 2
  # Extract recovery times for the individuals involved
  recov_times_ii <- recovery_times[c(out_ikt1_dat$per1[ii], out_ikt1_dat$per2[ii])]
  
  # Determine if recovery times are within the specified interval
  recov_in_intvl <- c(
    (out_ikt1_dat$t0[ii] <= recov_times_ii[1] & recov_times_ii[1] <= out_ikt1_dat$t1[ii]),
    (out_ikt1_dat$t0[ii] <= recov_times_ii[2] & recov_times_ii[2] <= out_ikt1_dat$t1[ii])
  )
  
  # Create a data frame of recovery times and filter those within the interval
  recov_times_ii <- cbind.data.frame(recov_times_ii, recov_in_intvl)
  recov_times_ii <- recov_times_ii$recov_times_ii[as.logical(recov_times_ii$recov_in_intvl)]
  
  # Simulate CTMC events based on the number of recovery times in the interval
  if (length(recov_times_ii) == 0) {
    postSI_events <- simulate_ctmc2(
      params[5], params[8], vector('numeric'), out_ikt1_dat$t1[ii] - out_ikt1_dat$t0[ii], 
      out_ikt1_dat$net0[ii], out_ikt1_dat$net1[ii]
    )
  } else if (length(recov_times_ii) == 1) {
    postSI_events <- simulate_ctmc2(
      params[5:4], params[8:7], c(recov_times_ii) - out_ikt1_dat$t0[ii], 
      out_ikt1_dat$t1[ii] - out_ikt1_dat$t0[ii], 
      out_ikt1_dat$net0[ii], out_ikt1_dat$net1[ii]
    )
  } else if (length(recov_times_ii) == 2) {
    postSI_events <- simulate_ctmc2(
      params[5:3], params[8:6], c(recov_times_ii) - out_ikt1_dat$t0[ii], 
      out_ikt1_dat$t1[ii] - out_ikt1_dat$t0[ii], 
      out_ikt1_dat$net0[ii], out_ikt1_dat$net1[ii]
    )
  }
  
  # Adjust the time of events and add additional information if events exist
  if (nrow(postSI_events) > 0) {
    postSI_events$time <- postSI_events$time + out_ikt1_dat$t0[ii]
    postSI_events$event = ifelse(postSI_events$event == 1, 3, 6) + 
      3 - as.integer(cut(postSI_events$time, c(out_ikt1_dat$t0[ii], recov_times_ii, Inf)))
  }
  
  return(postSI_events)
}
trunc_piecewise_obj = function(x, u, intervals, lambda, T_max){(trunc_piecewise_cdf(x, intervals, lambda, T_max) - u)^2}
trunc_piecewise_cdf = function(x, intervals, lambda, ub){
  (1 - survival_piecewise_const_exp_cpp(x,intervals, lambda)) / (1 - survival_piecewise_const_exp_cpp(ub, intervals, lambda))
}

simulate_ctmc2 <- function(conn_rate, disconn_rate, breakpoints, T_max, first_state, last_state) {
  # alpha, omega: vectors of rates for each sub-interval
  # breakpoints: vector of times defining sub-intervals (excluding 0 and T_max)
  # T_max: total time
  # first_state: initial state
  # last_state: required final state
  
  # conn_rate = params[5:4]; disconn_rate = params[8:7];
  # breakpoints = c(recov_times_ii) - out_ikt1_dat$t0[ii]
  # T_max = out_ikt1_dat$t1[ii] - out_ikt1_dat$t0[ii]; first_state = out_ikt1_dat$net0[ii]; last_state = out_ikt1_dat$net1[ii]
  
  # Append 0 and T_max to the breakpoints for complete intervals
  
  # breakpoints = c(1); T_max = 3; first_state = 0; last_state = 1
  # conn_rate = c(0.5, 1); disconn_rate = c(0.5, 1)*2;
  intervals <- c(0, breakpoints, Inf)
  
  if (first_state == last_state){
    # benchmark(simulate_ctmc(conn_rate, disconn_rate, breakpoints, T_max, first_state, last_state),
    #   simulate_ctmc_cpp(conn_rate, disconn_rate, breakpoints, T_max, first_state, last_state), replications = 1000)
    # 
    # 
    # simulate_ctmc(conn_rate, disconn_rate, breakpoints, T_max, first_state, last_state)
    simulate_ctmc_cpp(conn_rate, disconn_rate, breakpoints, T_max, first_state, last_state)
  }else{
    if (first_state == 0){
      lambda = conn_rate
    }else{
      lambda = disconn_rate
    }
    
    st_time = optimize(trunc_piecewise_obj, u =runif(1), intervals = c(0, breakpoints, Inf), lambda = lambda, T_max = T_max,
                       lower = 0, upper = T_max)$minimum
    
    current_time <- st_time
    current_state <- 1 - first_state
    breakpoints = breakpoints[breakpoints > st_time]
    if (length(breakpoints) == 0){
      breakpoints = vector('numeric')
    } else{
      breakpoints = breakpoints - st_time
    }
    
    nintvl = length(conn_rate)
    
    events = simulate_ctmc_cpp(conn_rate[(nintvl - length(breakpoints)): nintvl], 
                               disconn_rate[(nintvl - length(breakpoints)): nintvl], 
                               breakpoints, T_max-st_time, current_state, last_state)
    if(nrow(events) > 0){
      events$time = events$time + st_time
    }
    rbind.data.frame(data.frame(time = st_time, event = 1 - first_state), events)
  }
}



impute_network_events <- function(
    edges_to_sample2,
    infect_times2,
    recov_times2,
    params.cur2
) {
  # edges_to_sample2 = edges_to_sample
  # infect_times2 = infect_times
  # recov_times2 = recov_times
  # params.cur2 = params.cur

  edges_to_sample_tmp2 = edges_to_sample2
  
  # Summarize the agent pair, interval, end point network status, 
  infect_times2[abs(infect_times2 - 0) < 1e-4] = NA
  edges_to_sample_tmp2$i1 = infect_times2[edges_to_sample_tmp2$per1]
  edges_to_sample_tmp2$i2 = infect_times2[edges_to_sample_tmp2$per2]
  edges_to_sample_tmp2$r1 = recov_times2[edges_to_sample_tmp2$per1]
  edges_to_sample_tmp2$r2 = recov_times2[edges_to_sample_tmp2$per2]
  
  edges_to_sample_tmp2$i1[(edges_to_sample_tmp2$i1 < edges_to_sample_tmp2$t0 |
                             edges_to_sample_tmp2$i1 > edges_to_sample_tmp2$t1)] = NA
  
  edges_to_sample_tmp2$i2[(edges_to_sample_tmp2$i2 < edges_to_sample_tmp2$t0 |
                             edges_to_sample_tmp2$i2 > edges_to_sample_tmp2$t1)] = NA
  
  edges_to_sample_tmp2$r1[(edges_to_sample_tmp2$r1 < edges_to_sample_tmp2$t0 |
                             edges_to_sample_tmp2$r1 > edges_to_sample_tmp2$t1)] = NA
  
  edges_to_sample_tmp2$r2[(edges_to_sample_tmp2$r2 < edges_to_sample_tmp2$t0 |
                             edges_to_sample_tmp2$r2 > edges_to_sample_tmp2$t1)] = NA
  
  edges_to_sample_tmp2$health_at_t0 = (infect_times2[edges_to_sample_tmp2$per1] <= edges_to_sample_tmp2$t0 & 
                                         edges_to_sample_tmp2$t0 < recov_times2[edges_to_sample_tmp2$per1]) + 
    (infect_times2[edges_to_sample_tmp2$per2] <= edges_to_sample_tmp2$t0 & 
       edges_to_sample_tmp2$t0 < recov_times2[edges_to_sample_tmp2$per2])
  
  # initial health status at the beginning of each interval
  
  # st = Sys.time()
  if (!is.null(init.params)) params.cur2 = c(params.cur2[1:2], init.params[-c(1,2)])
  imputed_net_events_B = do.call("rbind", impute_all(edges_to_sample_tmp2, params.cur2, infect_times2, recov_times2))
  imputed_net_events_B
}



