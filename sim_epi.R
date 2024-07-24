rm(list = ls())
library(Rcpp)
library(matrixStats)
library(VGAM)
source("./sim_inference.R")
source("./inference_util.R")
source("./pre_process.R")
source("./inference_miss_recov.R")
library(parallel)

nsample = 4e3; nburn = 1e3; window_length = 12
pr = data.frame(count = rep(1,4), avg = c(0.075, 0.1, 0.005, 0.05))
beta = c(0.05, 0.075)
net_params = rbind(c(.005, 0.001, .005, .05, 0.1, .05),
                   c(.005, 0.001, .005, .05, 0.1, .05)*2)
sim_params = cbind(rep(beta, nrow(net_params)), 0.1, matrix(rep(net_params, each = length(beta)), ncol = 6))

idx = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

all_dats = readRDS("all_dats_2023-11-11.rds")

all_res_coar_net = list()
all_res_full_net = list()
for (j in 1: nrow(sim_params)) {
  # j = 4
  print(paste0("Now running the ", j ,"-th param set"))
  dats = all_dats[[idx]][[j]]

  miss_dats1 = miss_recovery(dats, interval = 7, miss_prop = 1)
  # obs_time = c(0, ceiling(max(miss_dats1$events$time)))
  obs_time = c(0, unique(ceiling(miss_dats1$events$time/window_length)*window_length))
  dats_no_recov_dis_net = net_coarsener3(miss_dats1, window_length)
  
  
  net_events = dats$events[dats$events$event %in% 3: 8, 1: 4]
  survey_dat = net_coarsener4(net_events, dats$G0, obs_time, FALSE)
  survey_dat$G0 = dats$G0
  survey_dat$I0 = dats$I0
  survey_dat$health_report = miss_dats1$report
  survey_dat$sick_events = rbind.data.frame(c(0, 1, dats$I0, NA),
                                            dats$events[dats$events$event == 1, 1: 4])
  survey_dat$truth = dats$truth
  
  
  res.fil3 = infer_miss_recov14_cpp(dat = survey_dat, force_one = TRUE, priors = pr, obs_time = obs_time, output.sams = 100, plot = F, verbose = F,
                                    samples = nsample, burn = nburn, thin = 1, seed = 1)
  
  # res.fil = infer_miss_recov12_cpp(dat = survey_dat, priors = pr, obs_time = obs_time, output.sams = 100, plot = F, verbose = T,
  #                                   samples = nsample, burn = nburn, thin = 1, seed = 1)
  # 
  res.fil2 = infer_miss_recov3(dat = dats_no_recov_dis_net, priors = pr, obs_time = obs_time, output.sams = 100, plot = F, verbose = T,
                               samples = nsample, burn = nburn, thin = 1, impute = "filter", seed = 1)
  res_array = array(NA, dim = c(dim(res.fil3[[1]]), 2))
  res_array[, , 1] = res.fil2
  res_array[, , 2] = res.fil3[[1]]
  
  temp_res = list()
  temp_res[[1]] = res_array
  names(temp_res) = paste0(window_length, "d")
  
  # Inference based on sampled recovery times and complete network observations
  res.fil1 = infer_miss_recov(dats = miss_dats1, priors = pr, output.sams = 100, plot = F, verbose = F,
                              samples = nsample, burn = nburn, thin = 1, impute = "filter", seed = 1)
  
  all_res_coar_net[[j]] = temp_res
  all_res_full_net[[j]] = res.fil1
}





out_result = list(
  full_net_res = all_res_full_net,
  coar_net_res = all_res_coar_net)
saveRDS(out_result, file = paste0("all_model_", as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')), ".rds"))