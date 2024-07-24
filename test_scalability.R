rm(list = ls())
library(Rcpp)
source("./sim_inference.R")
source("./inference_util.R")
source("./pre_process.R")
source("./inference_miss_recov.R")
library(parallel)
library(VGAM)
library(matrixStats)
library(ggplot2)
library(extrafont)
nsample = 50; nburn = 0
pr = data.frame(count = rep(1,4), avg = c(0.075, 0.1, 0.005, 0.05))
Nagents = 1000
init.p = 0.003
scaling_factor = c(0.1, 0.5, 1, 1.5)


set.seed(10001)
data_gen_time = c()
tot_events <- running_time <- matrix(NA, nrow = length(scaling_factor), ncol = 2)
for (i in 1: length(scaling_factor)){

  # seed = sample.int(.Machine$integer.max, 1)
  # set.seed(seed)
  # Generate the data (initial network, infector, epid/network events, and inference retults)
  nonstop = TRUE
  st = Sys.time()
  while (nonstop) { 
    dats = stochastic_coevolve_infer4(N=Nagents, bet=0.05, gam = 0.1, model = "SIR",tmax = 50,
                                      quarantine=T, init.p = init.p,
                                      alpha.r = c(.002, .001, .002)*scaling_factor[i],
                                      alpha.d = c(.015, 0.1, .015)*scaling_factor[i], 
                                      Bayes = T, MLE = T, priors = pr, verbose = F, plot = F, infer.verbose = F,
                                      return.infer = T, samples = 1000)
    event_count = table(factor(dats$events$event, levels = 1: 8))
    if(sum(event_count > 1) == length(event_count)){nonstop = FALSE}
  }
  ed = Sys.time()
  
  data_gen_time[i] = ed - st
  
  
  
  # Discard recovery times every 7 days
  miss_dats1 = miss_recovery(dats, interval = 7, miss_prop = 1)
  window_length = 7
  obs_time = c(0, unique(ceiling(miss_dats1$events$time/window_length)*window_length))
  

  net_events = dats$events[dats$events$event %in% 3: 8, 1: 4]
  survey_dat2 = net_coarsener4(net_events, dats$G0, obs_time, FALSE)
  survey_dat2$G0 = dats$G0
  survey_dat2$I0 = dats$I0
  survey_dat2$health_report = miss_dats1$report
  survey_dat2$sick_events = rbind.data.frame(c(0, 1, dats$I0, NA),
                                             dats$events[dats$events$event == 1, 1: 4])
  survey_dat2$truth = dats$truth
  tot_events[i, 1] = sum(sapply(survey_dat2$interaction_report, nrow))
  st = Sys.time()
  res.fil9 = infer_miss_recov14_cpp(dat = survey_dat2, force_one = TRUE, priors = pr, obs_time = obs_time, output.sams = 100, plot = F, verbose = T,
                                    samples = 50, burn = 0, thin = 1, seed = 1)
  ed = Sys.time()
  
  running_time[i, 1] = as.numeric(difftime(ed, st, units = "secs"))
  
  
  
  miss_dats1 = miss_recovery(dats, interval = 7, miss_prop = 1)
  window_length = 14
  obs_time = c(0, unique(ceiling(miss_dats1$events$time/window_length)*window_length))
  
  
  net_events = dats$events[dats$events$event %in% 3: 8, 1: 4]
  survey_dat2 = net_coarsener4(net_events, dats$G0, obs_time, FALSE)
  survey_dat2$G0 = dats$G0
  survey_dat2$I0 = dats$I0
  survey_dat2$health_report = miss_dats1$report
  survey_dat2$sick_events = rbind.data.frame(c(0, 1, dats$I0, NA),
                                             dats$events[dats$events$event == 1, 1: 4])
  survey_dat2$truth = dats$truth
  tot_events[i, 2] = sum(sapply(survey_dat2$interaction_report, nrow))
  st = Sys.time()
  res.fil9 = infer_miss_recov14_cpp(dat = survey_dat2, force_one = TRUE, priors = pr, obs_time = obs_time, output.sams = 100, plot = F, verbose = T,
                                    samples = 50, burn = 0, thin = 1, seed = 1)
  ed = Sys.time()
  
  running_time[i,2] = as.numeric(difftime(ed, st, units = "secs"))
}


plot_data = cbind.data.frame(num_events = (c(tot_events)),
                 running_time = c(running_time / 50),
                 obs_intvl = rep(c("7-Day", "14-Day"), each = 4))

saveRDS(plot_data, file = "plot_scalability.rds")

font_import()
loadfonts('pdf')  # For other devices if using Cairo


p1 = ggplot(data = plot_data, aes(x = num_events, y = running_time, color = obs_intvl)) + 
  geom_line() + geom_point() + xlab("Number of Events") + ylab("Imputation Time (s)") +
  # ggtitle("Imputation time  vs. Number of edges to impute") + 
  labs(color = "Obs. Interval") +
  theme(legend.position = c(0.05, 0.95), # Top left, inside
        legend.justification = c("left", "top"), # Anchor point of legend
        legend.box.margin = margin(-5, -5, -5, -5),
        text = element_text(family = "Times New Roman"),
        # legend.title = "Obs. Interval",
        # legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "white"),  # Set background to white
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("scalability.pdf", plot = p1, width = 6, height = 4, units = "in")

