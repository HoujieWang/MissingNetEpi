rm(list = ls())
library(Rcpp)
source("./sim_inference.R")
source("./inference_util.R")
source("./inference_util2.R")
source("./pre_process.R")
source("./inference_miss_recov.R")
sourceCpp("./inverse_cdf_sampler.cpp")
sourceCpp("./naive_rejection_sampler.cpp")
sourceCpp("./modified_rej_sampler.cpp")
library(parallel)
library(VGAM)
library(matrixStats)
library(Rfast)

nsample = 500; nburn = 250;
pr = data.frame(count = rep(1,4), avg = c(0.075, 0.1, 0.005, 0.05))
beta = c(0.05, 0.075)
net_params = rbind(c(.005, 0.001, .005, .05, 0.1, .05),
                   c(.005, 0.001, .005, .05, 0.1, .05)*2)
sim_params = cbind(rep(beta, nrow(net_params)), 0.1, matrix(rep(net_params, each = length(beta)), ncol = 6))

idx = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# idx = 1
all_dats = readRDS("small_all_dats_2025-02-03.rds")


all_res = list()

j = 4
dats = all_dats[[idx]][[j]]

if (abs(dats$truth[3] - 0.005) < 1e-5) {
  window_length = 5
  obs_time = seq(0, 60, window_length)
  dats$events = dats$events[dats$events$time <= max(obs_time), ] 
}
if (abs(dats$truth[3] - 0.01) < 1e-5) {
  window_length = 5
  obs_time = seq(0, 30, window_length)
  dats$events = dats$events[dats$events$time <= max(obs_time), ] 
}

miss_dats1 = miss_recovery(dats, interval = 5, miss_prop = 1)

dats_no_recov_dis_net = net_coarsener3(miss_dats1, window_length)
net_status_infect = infect_net_recorder(dats$events[dats$events$event %in% c(1,3: 8), 1: 4], dats$G0, dats$I0)
net_events = dats$events[dats$events$event %in% 3: 8, 1: 4]
net_events[net_events$per1 > net_events$per2, 3: 4] = net_events[net_events$per1 > net_events$per2, 4:3]
# survey_dat1 = net_coarsener4(net_events, dats$G0, obs_time, TRUE)
survey_dat2 = net_coarsener4(net_events, dats$G0, obs_time, FALSE)


survey_dat2$G0 = dats$G0
survey_dat2$I0 = dats$I0
survey_dat2$health_report = miss_dats1$report
survey_dat2$sick_events = rbind.data.frame(c(0, 1, dats$I0, NA),
                                           dats$events[dats$events$event == 1, 1: 4])
survey_dat2$truth = dats$truth


# 
# true_recov_times = rep(dats$events$time[length(dats$events$time)]+1, nrow(dats$G0))
# true_recov_times[dats$events$per1[dats$events$event == 2]] = dats$events$time[dats$events$event == 2]
# true_recov_times[colSums(survey_dat2$health_report == 1) == 0] = Inf


# 

res.fil1 = infer_miss_recov27_cpp(dat = survey_dat2, net_status_infect, fix_prop = 0, 
                                  # init.params = survey_dat2$truth, 
                                  # recov_times = true_recov_times, 
                                  force_one = TRUE, priors = pr, obs_time = obs_time, output.sams = 100, plot = F, verbose = T,
                                  samples = nsample, burn = nburn, thin = 1, seed = 1)

# tmp = real_events[real_events$event == 7, ]
# 
# which((tmp$time >= true_recov_times[tmp$per1] & infect_times[tmp$per2] <= tmp$time & tmp$time <= true_recov_times[tmp$per2]) |
#   (tmp$time >= true_recov_times[tmp$per2] & infect_times[tmp$per1] <= tmp$time & tmp$time <= true_recov_times[tmp$per1]))

# colMeans(res.fil2[[1]])
# dats$infer2[nrow(dats$infer2), ]
# dats$infer[(nrow(dats$infer)-7): nrow(dats$infer),]
# 
# real_events = dats$events
# real_events[which(real_events$per2 < real_events$per1), c(3, 4)] = real_events[which(real_events$per2 < real_events$per1), c(4, 3)]
saveRDS(res.fil1, file = paste0("run_sim_small_", j, '_', as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')), ".rds"))
