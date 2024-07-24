rm(list = ls())
library(parallel)
library(Rcpp)
library(Rfast)
library(RcppZiggurat)
library(RcppAlgos)
library(tidyverse)
source("./sim_inference.R")
source("./inference_util.R")
source("./pre_process.R")
source("./inference_miss_recov.R")

priors = data.frame(count = rep(1,4), avg = c(0.075, 0.1, 0.005, 0.05))
N_pop = 500; Tmax = 100; event.max = 1.5e4
all_parm = expand.grid(self_prop = c(0.1, 0.5, 0.9),
                       lag = c(1, 3),
                       altru_d = c(1))
all_dats = list()
seeds = sample(1e8, 10)

idx = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
for (i in 1: nrow(all_parm)) {
  # i = 1
  print(paste0("NOW GENERATING the", i, "-th DATA"))
  parm = as.numeric(all_parm[i, ])
  self_per = c(1: N_pop)[as.logical(rbinom(N_pop, 1, parm[1]))]
  temp_dat = list()
  for (j in 1: length(seeds)) {
    set.seed(seeds[j])
    dats = stochastic_coevolve_infer6(N=N_pop, 
                                      tmax=100, event.max = event.max,
                                      priors = data.frame(count = rep(1,4), avg = c(0.075, 0.1, 0.005, 0.05)),
                                      bet=0.03, gam=0.1, 
                                      init.infec = 1,
                                      alpha.r=c(0.01, 0.001, 0.005, 0.01, 0.01, 0.01), 
                                      alpha.d=c(0.05, parm[3], 0.075, 0.05, 0.05, 0.05),
                                      init.net = NULL, selfish = c(1: N_pop)[as.logical(rbinom(N_pop, 1, parm[1]))], 
                                      lag_period = parm[2],
                                      init.p = 0.1)
    temp_dat[[j]] = dats
  }
  all_dats[[i]] = temp_dat
}
saveRDS(all_dats, file = paste0("sim_data.rds"))

for (lag in unique(all_parm$lag)) {
  # lag = 1
  idx = which(all_parm$lag == lag)
  plot_dat = c()
  line_id = 0
  for (i in idx) {
    # i = 1
    plot_dat = rbind.data.frame(plot_dat,
                                cbind.data.frame(do.call("rbind", lapply(1: length(all_dats[[i]]), 
                                                                         function(ii){cbind(all_dats[[i]][[ii]]$events$time, all_dats[[i]][[ii]]$SIR_record, line_id+ii)})),
                                                 prop = all_parm$self_prop[i])
                                )
    line_id = line_id + length(all_dats[[i]])
  }
  colnames(plot_dat) = c("time", "SI", "S", "I", "R", "id", "selfish%")
  plot_dat$id = as.factor(plot_dat$id)
  plot_dat$`selfish%` = as.factor(plot_dat$`selfish%`)
  plot_dat_temp = plot_dat[plot_dat$id %in% c(5, 15, 25), ]
  fig1 = ggplot(data = plot_dat_temp, aes(x = time, y = SI, group = id,color = `selfish%`)) + geom_line()
  fig2 = ggplot(data = plot_dat_temp, aes(x = time, y = S, group = id,color = `selfish%`)) + geom_line()
  fig3 = ggplot(data = plot_dat_temp, aes(x = time, y = I, group = id,color = `selfish%`)) + geom_line() + xlim(0, 11) + ylim(0, 450)
  fig4 = ggplot(data = plot_dat_temp, aes(x = time, y = R, group = id,color = `selfish%`)) + geom_line()
  
}

# table(all_dats[[3]]$events$event)

# curve_type = c("SI", "S", "I", "R")
# # for (i in 1: 4) {
#   i = 3
#   SI_tab = data.frame()
#   for (ii in 4:6) {
#     # ii = 1
#     SI_tab = rbind.data.frame(SI_tab, 
#                               cbind.data.frame(all_dats[[ii]]$events$time,
#                                 all_dats[[ii]]$SIR_record[, i],
#                                 paste0("(", all_parm[ii,1], "; ",  all_parm[ii,2]))
#                               )
#   }
#   colnames(SI_tab) = c("time", "value", "type")
#   ggplot(data = SI_tab,aes(x = time, y = value, color = type)) + 
#     geom_line() + xlim(0, 25) + 
#     ggtitle(paste0(curve_type[i], " vs. time"))
# # }

