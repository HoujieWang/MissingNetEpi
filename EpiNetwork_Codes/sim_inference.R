# 05/06/2019
# Simulation & Inference
# work with COMPLETE data


library(Matrix)
library(igraph)
library(ggplot2)
library(dplyr)



# 05/20/2019
# 2. Modified function: for SIR model, recovered people adopt S- rates 
# for their re- and dis- connections
# alpha.R* = alpha.S* and omega.R* = omega.S*
# ALSO: return true parameters too
# 05/21/2019
# add more info to "MLE.tab" to enable plotting over repeated experiments

# 02/04/2020
# add network degree plots to address Rev.1 concern
stochastic_coevolve_infer2 <- 
  function(N=50, tmax=100, 
           bet=0.05, gam=0.2, init.infec = 1, model = "SIS",
           alpha.r=0.005, alpha.d=0.05, quarantine=F,
           init.net = NULL, init.p = 0.1,
           verbose = F, plot = F,
           infer.interval = 100, 
           Bayes = F, priors = NULL, samples = 500,
           infer.verbose = F, infer.plot = F, return.infer = F,
           net.deg.plot = F, net.interval = 100,
           demo.sleep = F){
    ## alpha.* = c(SSrate, SIrate, IIrate) ##
    #set.seed(seed)
    
    # check edge rates validity
    if(quarantine & (length(alpha.r) < 3 | length(alpha.d) < 3)){
      stop("Under quarantine, need to specify 3 values for edge connection/disconnection rates respectively!")
    }
    
    # output basic info
    cat("N =", N, " Model:", model,"\n",
        "beta =", bet, " gamma =", gam, "\n",
        "init.infec =", init.infec, " init.p =", init.p, "\n",
        "alpha.r =", alpha.r, "\n", "alpha.d =", alpha.d, "\n",
        "quarantine?", quarantine,"\n",
        "Bayesian inference?", Bayes, "\n",
        "inference interval =", infer.interval, "\n")
    
    # summarization stats to be recorded
    time = NULL
    event.type = NULL
    preval = NULL
    dens = NULL
    # add involved individual labels 
    per1 = NULL
    per2 = NULL
    
    # more stats to record for inference
    if(quarantine){
      event.counts = numeric(8)
      big.sums = numeric(8)
    }else{
      event.counts = numeric(4)
      big.sums = numeric(4)
    }
    tot.event = 0
    
    if(Bayes){
      if(quarantine & nrow(priors)==4){
        a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
        avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
        b.pr = a.pr/avg.pr
      }else{
        a.pr = priors$count
        b.pr = a.pr/priors$avg
      }
      if(quarantine){
        vars = c("beta","gamma",
                 "alpha.SS","alpha.SI","alpha.II",
                 "omega.SS","omega.SI","omega.II")
      }else{
        vars = c("beta","gamma","alpha","omega")
      }
      Bayes.tab = NULL
    }else{
      MLE.tab = NULL
    }
    params = c(bet, gam, alpha.r, alpha.d)
    
    
    # initialize ER(p) graph if the initial netwok is not provided
    if(is.null(init.net)){
      init.net = sample_gnp(n = N, p = init.p, directed = F)
    }
    # convert to adjacency matrix
    adjmat = get.adjacency(init.net, type=c("both"), names=F)
    adjmat = as(adjmat,"dgTMatrix")
    # save the initial network structure (adj.mat)
    init.adj = adjmat
    
    
    ## check whether adjmat is symmetric and nnzero() returns an even number
    if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
      stop(paste("The initial adjmat isn't symmetric, where nonzero =", 
                 nnzero(adjmat),"\n"))
    }
    
    # initialize epidemic
    infected = sample(N, init.infec)
    epid = rep(0,N); epid[infected] = 1
    cat("Initially infected: ",infected,"\n")
    # save initial cases
    init.infected = infected
    
    # initialize data objects to keep track of
    susceptible = which(epid==0)
    N.I = init.infec; N.S = N - N.I
    t.cur = 0
    if(quarantine){
      N.c = numeric(3); N.d = numeric(3)
      adjmat.ii = adjmat[infected,infected]
      N.c[3] = nnzero(adjmat.ii)/2
      N.d[3] = N.I * (N.I - 1)/2 - N.c[3]
      # both need this:
      adjmat.si = adjmat[susceptible,infected]
      if(model=="SIR"){
        recovered = NULL; N.R = 0
        adjmat.rs = adjmat[union(susceptible,recovered), union(susceptible,recovered)]
        N.c[1] = nnzero(adjmat.rs)/2
        N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
        adjmat.ri = adjmat.si
        N.c[2] = nnzero(adjmat.ri)
        N.d[2] = N.S * N.I - N.c[2]
      }else{
        adjmat.ss = adjmat[susceptible,susceptible]
        N.c[1] = nnzero(adjmat.ss)/2
        N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
        
        N.c[2] = nnzero(adjmat.si)
        N.d[2] = N.S * N.I - N.c[2]
      }
      # at the start: N.R = 0, so this works
      N.SI = N.c[2]
    }else{
      N.c = nnzero(adjmat)/2
      N.d = N * (N-1)/2 - N.c
      
      adjmat.si = adjmat[susceptible,infected]
      N.SI = nnzero(adjmat.si)
    }
    
    
    # iterations
    while (t.cur < tmax) {
      # NOTES #
      # a. check to make sure N.S also > 0
      # b. handle edge cases (N.S = 1 and/or N.I = 1)
      
      # # check whether all N's are non-negative
      # Ns = c(N.SI, N.I, N.d, N.c)
      # if(any(Ns < 0) | all(Ns == 0)){
      #   msg = paste("Weird stuff happened! All the Ns =", paste(Ns,collapse = ", "), "\n")
      #   stop(msg)
      # }
      
      # 1. calculate sum of risks and sample the next event time
      Lambda = bet * N.SI + N.I * gam + sum(alpha.r * N.d) + sum(alpha.d * N.c)
      t.next = t.cur + rexp(1, rate = Lambda)
      
      # 2. determine the event type
      ratios = c(bet * N.SI, N.I * gam, alpha.r * N.d, alpha.d * N.c)
      # 2.1 deal with network change
      if(quarantine){
        
        z = sample(8, 1, prob = ratios)
        
        if(z==3){
          # reconnect a random S-S/S-R/R-R pair
          if(model == "SIS"){
            num.dis = (N.S - 1) - rowSums(adjmat.ss)
            from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
            to.cands = which(adjmat.ss[from.ind,]==0)
          }else{
            num.dis = (N.S + N.R - 1) - rowSums(adjmat.rs)
            from.ind = sample(N.S + N.R, 1, prob = num.dis); from = union(susceptible, recovered)[from.ind]
            to.cands = which(adjmat.rs[from.ind,]==0)
          }
          to.cands = to.cands[to.cands!=from.ind]
          if(length(to.cands) == 1){
            to.ind = to.cands
          }else{
            to.ind = sample(to.cands,1)
          }
          if(model == "SIS"){
            to = susceptible[to.ind]
            adjmat.ss[as.integer(from.ind),as.integer(to.ind)] = 1
            adjmat.ss[as.integer(to.ind),as.integer(from.ind)] = 1
            msg = paste0("reconnect an S-S pair, (",from,",",to,")")
          }else{
            to = union(susceptible, recovered)[to.ind]
            adjmat.rs[as.integer(from.ind),as.integer(to.ind)] = 1
            adjmat.rs[as.integer(to.ind),as.integer(from.ind)] = 1
            msg = paste0("reconnect an S-S/S-R/R-R pair, (",from,",",to,")")
          }
          adjmat[as.integer(from),as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          N.c[1] = N.c[1] + 1; N.d[1] = N.d[1] - 1
        }
        if(z==4){
          # reconnect a random S-I/R-I pair
          if(model == "SIS"){
            if (N.S == 1){
              from = susceptible
              to.ind = sample(which(adjmat.si==0),1); to = infected[to.ind]
            }else if (N.I == 1){
              num.dis = N.I - adjmat.si
              from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
              to = infected
            }else{
              num.dis = N.I - rowSums(adjmat.si)
              from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
              to.ind = sample(which(adjmat.si[from.ind,]==0),1); to = infected[to.ind]
            }
            msg = paste0("reconnect an S-I pair, (",from,",",to,").")
          }else{
            #N.SR = N.S + N.R; 
            # 03/27/2020: debug
            sus.recov = union(susceptible, recovered)
            N.SR = length(sus.recov)
            if (N.SR == 1){
              from = sus.recov
              # 03/27/2020: debug
              num.dis = N.SR - as.vector(adjmat.ri)
              to.ind = sample(N.I, 1, prob = num.dis); to = infected[to.ind]
              #to.ind = sample(which(adjmat.ri==0),1); to = infected[to.ind]
            }else if (N.I == 1){
              num.dis = N.I - as.vector(adjmat.ri)
              from.ind = sample(N.SR,1,prob = num.dis); from = sus.recov[from.ind]
              to = infected
            }else{
              num.dis = N.I - rowSums(adjmat.ri)
              # debugging...
              #cat(N.S, N.R, N.I, '\n')
              #cat(length(susceptible), length(recovered), length(infected),'\n')
              #print(dim(adjmat.ri))
              from.ind = sample(N.SR,1,prob = num.dis); from = sus.recov[from.ind]
              
              # 03/28/2020: debug
              vec.dis = 1 - as.vector(adjmat.ri[from.ind,])
              to.ind = sample(N.I,1,prob = vec.dis); to = infected[to.ind]
              #to.ind = sample(which(adjmat.ri[from.ind,]==0),1); to = infected[to.ind]
            }
            msg = paste0("reconnect an S-I/R-I pair, (",from,",",to,").")
          }
          adjmat[as.integer(from),as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          # 03/27/2020
          # try to fix a bug
          if(length(susceptible)>0){
            adjmat.si = adjmat[susceptible,infected]
            N.SI = nnzero(adjmat.si)
          }else{
            N.SI = 0
          }
          
          if(model == "SIR"){ adjmat.ri = adjmat[sus.recov,infected] }
          
          N.c[2] = N.c[2] + 1; N.d[2] = N.d[2] - 1
          #N.SI = N.c[2]
        }
        if(z==5){
          # reconnect a random I-I pair
          num.dis = (N.I - 1) - rowSums(adjmat.ii)
          from.ind = sample(N.I,1,prob = num.dis); from = infected[from.ind]
          to.cands = which(adjmat.ii[from.ind,]==0)
          to.cands = to.cands[to.cands!=from.ind]
          if(length(to.cands) == 1){
            to.ind = to.cands
          }else{
            to.ind = sample(to.cands,1)
          }
          to = infected[to.ind]
          adjmat.ii[as.integer(from.ind),as.integer(to.ind)] = 1
          adjmat.ii[as.integer(to.ind),as.integer(from.ind)] = 1
          adjmat[as.integer(from),as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          N.c[3] = N.c[3] + 1; N.d[3] = N.d[3] - 1
          
          msg = paste0("reconnect an I-I pair, (",from,",",to,").")
        }
        if(z==6){
          # disconnect a random S-S/R-S/R-R pair
          dis.ind = sample(N.c[1], 1) 
          if(model == "SIS"){
            from.ind = adjmat.ss@i[dis.ind] + 1; to.ind = adjmat.ss@j[dis.ind] + 1
            adjmat.ss[as.integer(from.ind),as.integer(to.ind)] = 0
            adjmat.ss[as.integer(to.ind),as.integer(from.ind)] = 0
            
            from = susceptible[from.ind]; to = susceptible[to.ind]
            
            msg = paste0("disconnect an S-S pair, (",from,",",to,").")
          }else{
            sus.recov = union(susceptible, recovered)
            from.ind = adjmat.rs@i[dis.ind] + 1; to.ind = adjmat.rs@j[dis.ind] + 1
            adjmat.rs[as.integer(from.ind),as.integer(to.ind)] = 0
            adjmat.rs[as.integer(to.ind),as.integer(from.ind)] = 0
            
            from = sus.recov[from.ind]; to = sus.recov[to.ind]
            
            msg = paste0("disconnect an S-S/S-R/R-R pair, (",from,",",to,").")
          }
          adjmat[as.integer(from),as.integer(to)] = 0
          adjmat[as.integer(to),as.integer(from)] = 0
          
          N.c[1] = N.c[1] - 1; N.d[1] = N.d[1] + 1
          
          # ## sometimes from/to is just empty... check it out
          # if(length(from)!=1){
          #   stop("from.ind = ", from.ind, 
          #        " but susceptible length is ", length(susceptible),
          #        " N.c[1] =", N.c[1], " dis.ind =", dis.ind, 
          #        " while adjmat.ss@i has length ", length(adjmat.ss@i))
          # }
          # ##
        }
        if(z==7){
          # disconnect a random S-I/R-I pair
          # dis.ind = sample(N.c[2], 1)
          if(model == "SIS"){
            if (N.S == 1){
              from = susceptible
              to.ind = sample(which(adjmat.si==1),1); to = infected[to.ind]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.si = adjmat[susceptible,infected]
            }else if (N.I == 1){
              num.con = as.vector(adjmat.si)
              from.ind = sample(N.S,1,prob = num.con); from = susceptible[from.ind]
              to = infected
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.si = adjmat[susceptible,infected]
            }else{
              dis.ind = sample(N.c[2], 1)
              from = adjmat.si@i[dis.ind] + 1; to = adjmat.si@j[dis.ind] + 1
              adjmat.si[as.integer(from),as.integer(to)] = 0
              
              from = susceptible[from]; to = infected[to]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
            }
            msg = paste0("disconnect an S-I pair, (",from,",",to,").")
          }else{
            N.SR = N.S + N.R; sus.recov = union(susceptible, recovered)
            if (N.SR == 1){
              from = sus.recov
              to.ind = sample(which(adjmat.ri==1),1); to = infected[to.ind]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.ri = adjmat[sus.recov,infected]
            }else if (N.I == 1){
              num.con = adjmat.ri
              from.ind = sample(N.SR,1,prob = num.con); from = sus.recov[from.ind]
              to = infected
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.ri = adjmat[sus.recov,infected]
            }else{
              dis.ind = sample(N.c[2], 1)
              from = adjmat.ri@i[dis.ind] + 1; to = adjmat.ri@j[dis.ind] + 1
              # debugging...
              #cat('dis.ind=',dis.ind,'\n')
              #cat('from=',from,', to=',to, '\n')
              #cat(dim(adjmat.ri), '\n')
              #cat(length(adjmat.ri@i),  length(adjmat.ri@j), N.c[2], '\n')
              
              if(length(adjmat.ri@i) != N.c[2]){
                cat('Discrepancy when disconnecting H-I...\n')
                stop(paste0("At time ", t.next, ", ",msg,"\n"))
              }
              
              #adjmat.ri[as.integer(from),as.integer(to)] = 0
              
              from = sus.recov[from]; to = infected[to]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.ri = adjmat[sus.recov,infected]
            }
            if(length(susceptible)>0){
              adjmat.si = adjmat[susceptible,infected]
            }else{
              adjmat.si = NULL
            }
            msg = paste0("disconnect an S-I/R-I pair, (",from,",",to,").")
          }
          
          N.c[2] = N.c[2] - 1; N.d[2] = N.d[2] + 1
          N.SI = nnzero(adjmat.si)

        }
        if(z==8){
          # disconnect a random I-I pair
          dis.ind = sample(N.c[3], 1) 
          from = adjmat.ii@i[dis.ind] + 1; to = adjmat.ii@j[dis.ind] + 1
          adjmat.ii[as.integer(from),as.integer(to)] = 0
          adjmat.ii[as.integer(to),as.integer(from)] = 0
          
          from = infected[from]; to = infected[to]
          adjmat[as.integer(from),as.integer(to)] = 0
          adjmat[as.integer(to),as.integer(from)] = 0
          
          N.c[3] = N.c[3] - 1; N.d[3] = N.d[3] + 1
          
          msg = paste0("disconnect an I-I pair, (",from,",",to,").")
        }
      }else{
        z = sample(4, 1, prob = ratios)
        if(z==3){
          # reconnect a random pair
          num.dis = (N - 1) - rowSums(adjmat)
          from = sample(N,1,prob = num.dis)
          to.cands = which(adjmat[from,]==0)
          to.cands = to.cands[to.cands!=from]
          if(length(to.cands) == 1){
            to = to.cands
          }else{
            to= sample(to.cands,1)
          }
          adjmat[as.integer(from), as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          
          N.c = N.c + 1; N.d = N.d - 1
          
          # UPDATE: if any of S or I is empty now
          if(length(susceptible) == 0 | length(infected) == 0){
            adjmat.si = NULL; N.SI = 0
          }else{
            adjmat.si = adjmat[susceptible,infected]
            N.SI = nnzero(adjmat.si)
          }
          
          msg = paste0("reconnect a pair, (",from,",",to,").")
        }
        if(z==4){
          # disconnect a random pair
          dis.ind = sample(N.c, 1)
          from = adjmat@i[dis.ind] + 1; to = adjmat@j[dis.ind] + 1
          adjmat[from,to] = 0; adjmat[to,from] = 0
          N.c = N.c - 1; N.d = N.d + 1
          
          # UPDATE: if any of S or I is empty now
          if(length(susceptible)==0 | length(infected)==0){
            adjmat.si = NULL; N.SI = 0
          }else{
            adjmat.si = adjmat[susceptible,infected]
            N.SI = nnzero(adjmat.si)
          }
          
          msg = paste0("disconnect a pair, (",from,",",to,").")
        }
      }
      
      #2.2 deal with epidemic progress
      if(z==2){
        # recover an infected person
        rec.ind = sample(N.I,1)
        recover = infected[rec.ind]
        infected = infected[-rec.ind]
        N.I = N.I - 1 
        if(model == "SIS"){
          susceptible = c(susceptible,recover)
          epid[recover] = 0
          N.S = N.S + 1
        }else{
          epid[recover] = -1
          if(quarantine){
            recovered = c(recovered, recover)
            N.R = N.R + 1
          }
        }
        
        msg = paste0("individual ", recover,
                     " is recovered. Disease prevalence = ",
                     N.I/N, ".")
      }
      if(z==1){
        # select a susceptible person to infect
        if(N.S == 1){
          # only one person to infect
          new.infec = susceptible
          susceptible = NULL
        }else{
          # multiple susceptible inviduals to choose from
          if(N.I == 1){
            # only one infected person
            num.nei.i = adjmat.si
          }else{
            num.nei.i = rowSums(adjmat.si)
          }
          infec.ind = sample(N.S,1,prob = num.nei.i)
          new.infec = susceptible[infec.ind]
          susceptible = susceptible[-infec.ind]
        }
        epid[new.infec] = 1
        N.S = N.S - 1
        infected = c(infected, new.infec); N.I = N.I + 1
        
        msg = paste0("individual ", new.infec,
                     " is infected. Disease prevalence = ",
                     N.I/N, ".")
      }
      
      # A sanity check: if N.I==0, stop the simulation
      if(N.I==0){
        if(verbose){
          cat(paste0("At time ", t.next, ", "),msg,"\n")
        }
        # still update inference statistics
        del.t = t.next - t.cur
        z = as.integer(z)
        event.counts[z] = event.counts[z] + 1
        tot.event = tot.event + 1
        
        big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
        
        # still have to record the event
        time = c(time, t.next)
        event.type = c(event.type, z)
        preval = c(preval, N.I)
        dens = c(dens, nnzero(adjmat))
        # here: a recovery must be the latest event
        per1 = c(per1, recover)
        per2 = c(per2, NA)
        
        cat("Disease extinct at time", t.next, ", simulation stops.\n")
        break
      }
      
      
      # # check whether adjmat is still symmetric--things have gone wrong somewhere
      # if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
      #   stop(paste("The adjmat isn't symmetric any more, where nonzero =", 
      #              nnzero(adjmat), "and that's after z =", z,
      #              "from.ind =", from.ind, "to.ind = ", to.ind, 
      #              "to.cands =", paste(to.cands,collapse = ","),
      #              "from =", from, "to =", to,"\n"))
      # }
      
      
      # 2.2.* make updates about S-I, S-S, I-I connections if z==1 or 2
      # while handling edge cases of N.S <= 1 and/or N.I == 1
      if(z %in% c(1,2)){
        if(N.S == 0){
          adjmat.si = NULL
          N.SI = 0
          if(quarantine && model == "SIR" && N.R == 0){
            adjmat.ri = NULL
            N.c[2] = 0
          }
          # 03/27/2020: debug
          if(quarantine && model == "SIR" && N.R > 0){
            adjmat.ri = adjmat[recovered,infected]
            N.c[2] = nnzero(adjmat.ri)
          }
        }else{
          adjmat.si = adjmat[susceptible,infected]
          N.SI = nnzero(adjmat.si)
          if(quarantine && model == "SIR"){
            adjmat.ri = adjmat[union(susceptible,recovered),infected]
            N.c[2] = nnzero(adjmat.ri)
          }
        }
        if(quarantine){
          if(model == "SIS"){
            if(N.S <= 1){
              adjmat.ss = 0
              N.c[1] = 0
              N.d[1] = 0
            }else{
              adjmat.ss = adjmat[susceptible,susceptible]
              N.c[1] = nnzero(adjmat.ss)/2
              N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
            }
            N.c[2] = N.SI
            N.d[2] = N.S * N.I - N.c[2]
          }else{
            #N.SR = N.S + N.R; 
            sus.recov = union(susceptible,recovered)
            # 03/27/2020: debug
            N.SR = length(sus.recov)
            
            if(N.SR <= 1){
              adjmat.rs = 0
              N.c[1] = 0
              N.d[1] = 0
            }else{
              adjmat.rs = adjmat[sus.recov,sus.recov]
              N.c[1] = nnzero(adjmat.rs)/2
              N.d[1] = N.SR * (N.SR - 1)/2 - N.c[1]
            }
            N.d[2] = N.SR * N.I - N.c[2]
          }
          if(N.I == 1){
            adjmat.ii = 0
            N.c[3] = 0
            N.d[3] = 0
          }else{
            adjmat.ii = adjmat[infected,infected]
            N.c[3] = nnzero(adjmat.ii)/2
            N.d[3] = N.I * (N.I - 1)/2 - N.c[3]
          }
        }
      }
      
      # 3. update inference statistics and do inference occasionally
      # 3.1 update stats
      del.t = t.next - t.cur
      z = as.integer(z)
      event.counts[z] = event.counts[z] + 1
      tot.event = tot.event + 1
      
      big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
      
      # 3.2 inference if ...
      if(tot.event %% infer.interval == 0){
        if(Bayes){
          # Bayesian inference: Gamma sampling
          a.post = a.pr + event.counts
          b.post = b.pr + big.sums
          cat("With", tot.event,"events, posterior means and variances:\n",
              a.post/b.post, "\n",
              a.post/(b.post^2), "\n")
          # record posterior stats
          if(return.infer){
            m.post = a.post/b.post
            lb.post = qgamma(.025, a.post, b.post)
            ub.post = qgamma(.975, a.post, b.post)
            this.chunk = cbind(m.post, lb.post, ub.post, event.counts, tot.event)
            Bayes.tab = rbind(Bayes.tab, this.chunk)
          }
          # plot
          if(infer.plot){
            par(mfrow=c(2,2))
            for(ix in 1:length(params)){
              v = vars[ix]; p = params[ix]
              sams = rgamma(samples, shape=a.post[ix], rate = b.post[ix])
              xl = range(sams, p)
              hist(sams, main=paste("Posterior samples for",vars[ix]),
                   xlab = vars[ix], xlim=xl, col="lightblue")
              abline(v=p,col="red")
              if(demo.sleep & ix%%4 == 0){
                Sys.sleep(0.2)
              }
            }
          }
        }else{
          # MLE
          MLEs = event.counts/big.sums
          if(infer.verbose){
            cat("With",tot.event, "events, MLEs:","\n", MLEs, "\n")
          }
          MLE.tab = rbind(MLE.tab, c(tot.event, MLEs, event.counts))
        }
      }
      
      # updated 02/04/2020
      # 3.3 examine network degree distribution if
      if(net.deg.plot & (tot.event %% net.interval == 0)){
        degs = rowSums(adjmat)
        avg.deg = mean(degs)
        poi.degs = rpois(N, avg.deg)
        
        deg.dat = data.frame(Degrees = degs, DegPoi = poi.degs)
        print(
        ggplot(deg.dat) + geom_density(aes(x = Degrees, color="Actual degree")) +
          geom_density(aes(x = DegPoi, color="Poisson")) +
          labs(color="Distribution", 
               caption = paste("At t =", round(t.cur, digits = 2))) + 
          theme_bw(base_size = 14) %+replace%
          theme(legend.position = c(0.8,0.8))
        )
      }
      
      
      
      # 4. report the event and summarize and set t.cur
      
      t.cur = t.next
      
      # record the event
      time = c(time, t.cur)
      event.type = c(event.type, z)
      preval = c(preval, N.I)
      dens = c(dens, nnzero(adjmat))
      # record labels of individual(s) involved
      if (z >= 3){
        per1 = c(per1, from)
        per2 = c(per2, to)
      }else if (z==2){
        per1 = c(per1, recover)
        per2 = c(per2, NA)
      }else{
        per1 = c(per1, new.infec)
        per2 = c(per2, NA)
      }
      
      
      if(verbose){
        cat(paste0("At time ", t.cur, ", "),msg,"\n")
      }
      
    }
    
    # return results/summary
    cat("Simulation finished. At the end: \n", 
        "Disease prevalence =", N.I/N, "\n",
        "Network density = ", nnzero(adjmat)/(N*(N-1)),"\n")
    event.log = data.frame(time = time, event = event.type,
                           per1 = per1, per2 = per2,
                           preval = preval/N, dens = dens/(N*(N-1)))
    
    # inference at the end
    if(Bayes){
      # Bayesian inference: Gamma sampling
      a.post = a.pr + event.counts
      b.post = b.pr + big.sums
      cat("Posterior means and variances with", tot.event, "events in total:\n",
          a.post/b.post, "\n",
          a.post/(b.post^2), "\n")
      cat("Event type counts:\n", event.counts,"\n")
      # record posterior stats
      if(return.infer){
        m.post = a.post/b.post
        lb.post = qgamma(.025, a.post, b.post)
        ub.post = qgamma(.975, a.post, b.post)
        this.chunk = cbind(m.post, lb.post, ub.post, event.counts, tot.event)
        Bayes.tab = rbind(Bayes.tab, this.chunk)
      }
      # plot
      if(infer.plot){
        par(mfrow=c(2,2))
        for(ix in 1:length(params)){
          v = vars[ix]; p = params[ix]
          sams = rgamma(samples, shape=a.post[ix], rate = b.post[ix])
          xl = range(sams, p)
          hist(sams, main=paste("Posterior samples for",vars[ix]),
               xlab = vars[ix], xlim=xl, col="lightblue")
          abline(v=p,col="red")
          if(demo.sleep & ix%%4 == 0){
            Sys.sleep(0.2)
          }
        }
      }
    }else{
      # MLE
      MLEs = event.counts/big.sums
      cat("MLEs with", tot.event, "events in total:","\n", MLEs, "\n")
      cat("Event type counts:\n", event.counts,"\n")
      
      MLE.tab = rbind(MLE.tab, c(tot.event, MLEs, event.counts))
      MLE.tab = as.data.frame(MLE.tab)
      if(quarantine){
        names(MLE.tab) = c("event.num","beta","gamma",
                           "alpha.SS","alpha.SI","alpha.II",
                           "omega.SS","omega.SI","omega.II",
                           "infection","recovery",
                           "SS.connect","SI.connect","II.connect",
                           "SS.discon","SI.discon","II.discon")
      }else{
        names(MLE.tab) = c("event.num","beta","gamma","alpha","omega",
                           "infection","recovery","connect","discon")
      }
    }
    
    if(plot){
      subt = paste(model," beta =",bet," gamma =", gam, 
                   " alpha.r =", alpha.r, " alpha.d =", alpha.d)
      plot(preval ~ time, data = event.log, 
           xlab="Time", ylab="Disease prevalence", 
           main = "Disease Prevalence vs. Time", sub = subt,
           ylim = c(0,1), type="l",lwd=2)
      plot(dens ~ time, data = event.log,
           xlab="Time", ylab="Edge density", 
           main = "Network Edge Density vs. Time", sub = subt,
           type="l",lwd=2)
    }
    
    if(infer.plot){
      #params = c(bet, gam, alpha.r, alpha.d)
      if(Bayes){
        par(mfrow=c(1,1))
      }else{
        # a function to calculate the frequentist 95% CI for MLEs
        MLE.CI <- function(x){
          l.hat = x[1]; n.obs = x[2]
          if (is.na(l.hat) || l.hat == 0){
            lb = 0; ub = 0
          }else if (n.obs >= 20){
            lb = l.hat * (1 + qnorm(.025)/sqrt(n.obs))
            ub = l.hat * (1 + qnorm(.975)/sqrt(n.obs))
          }else{
            lb = l.hat * qchisq(.025, 2*n.obs) / (2*n.obs)
            ub = l.hat * qchisq(.975, 2*n.obs) / (2*n.obs)
          }
          return(c(lb,ub))
        }
        
        lcol = length(names(MLE.tab)[-1])
        vars = names(MLE.tab)[-1][1:(lcol/2)]
        counts = names(MLE.tab)[-1][(lcol/2 + 1):lcol]
        e.num = MLE.tab$event.num
        if(quarantine){
          par(mfrow=c(2,2))
        }else{
          par(mfrow=c(2,2))
        }
        for(ix in 1:length(vars)){
          v = vars[ix]
          if(all(is.na(MLE.tab[[v]]))){
            plot.new()
          }else{
            tab.v = MLE.tab[,c(v,counts[ix])]
            bounds.v = t(apply(tab.v, 1, MLE.CI))
            bounds.v = as.data.frame(bounds.v)
            names(bounds.v) = c("lb","ub")
            #s.m = mean(MLE.tab[[v]], na.rm=T)
            #s.sd = sd(MLE.tab[[v]], na.rm=T)
            yl = range(min(bounds.v$lb), max(bounds.v$ub), params[ix])
            plot(MLE.tab[[v]]~e.num, xlab="# events", ylab=v,
                 main = paste("MLE for", v, "v.s. # events"), 
                 type="l", ylim=yl)
            lines(bounds.v$lb~e.num, lty=2, col="gray")
            lines(bounds.v$ub~e.num, lty=2, col="gray")
          }
          abline(h=params[ix], col="red")
        }
        par(mfrow=c(1,1))
      }
    }
    
    # add: print out big.sums if Bayesian
    if(Bayes){
      cat("big.sums: \n", big.sums,"\n")
    }
    
    # 02/04/2020: degree distribution at the end
    if(net.deg.plot){
      degs = rowSums(adjmat)
      avg.deg = mean(degs)
      poi.degs = rpois(N, avg.deg)
      
      deg.dat = data.frame(Degrees = degs, DegPoi = poi.degs)
      print(
      ggplot(deg.dat) + geom_density(aes(x = Degrees, color="Actual degree")) +
        geom_density(aes(x = DegPoi, color="Poisson")) +
        labs(color="Distribution", 
             caption = paste("At the end, t =", round(t.cur, digits = 2))) + 
        theme_bw(base_size = 14) %+replace%
        theme(legend.position = c(0.8,0.8))
      )
    }
    
    
    # return event data
    # MODIFICATION: return true parametes too
    if(return.infer){
      if(Bayes){
        Bayes.tab = as.data.frame(Bayes.tab)
        Bayes.tab$var = vars
        infer.tab = Bayes.tab
      }else{
        infer.tab = MLE.tab
      }
      return(list(G0 = init.adj, I0 = init.infected,
                  events = event.log, infer=infer.tab,
                  truth = params))
    }else{
      return(list(G0 = init.adj, I0 = init.infected,
                  events = event.log,truth = params))
    }
    
  }

stochastic_coevolve_infer3 <- 
  function(N=50, tmax=100, 
           bet=0.05, gam=0.2, init.infec = 1, model = "SIS",
           alpha.r=0.005, alpha.d=0.05, quarantine=F,
           init.net = NULL, init.p = 0.1,
           verbose = F, plot = F,
           infer.interval = 100, 
           Bayes = T, MLE = T, priors = NULL, samples = 500,
           infer.verbose = F, infer.plot = F, return.infer = F,
           net.deg.plot = F, net.interval = 100,
           demo.sleep = F){
    ## alpha.* = c(SSrate, SIrate, IIrate) ##
    #set.seed(seed)
    
    # check edge rates validity
    if(quarantine & (length(alpha.r) < 3 | length(alpha.d) < 3)){
      stop("Under quarantine, need to specify 3 values for edge connection/disconnection rates respectively!")
    }
    
    # output basic info
    cat("N =", N, " Model:", model,"\n",
        "beta =", bet, " gamma =", gam, "\n",
        "init.infec =", init.infec, " init.p =", init.p, "\n",
        "alpha.r =", alpha.r, "\n", "alpha.d =", alpha.d, "\n",
        "quarantine?", quarantine,"\n",
        "Bayesian inference?", Bayes, "\n",
        "inference interval =", infer.interval, "\n")
    
    # summarization stats to be recorded
    time = NULL
    event.type = NULL
    preval = NULL
    dens = NULL
    # add involved individual labels 
    per1 = NULL
    per2 = NULL
    
    # more stats to record for inference
    if(quarantine){
      event.counts = numeric(8)
      big.sums = numeric(8)
    }else{
      event.counts = numeric(4)
      big.sums = numeric(4)
    }
    tot.event = 0
    
    if(Bayes){
      if(quarantine & nrow(priors)==4){
        a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
        avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
        b.pr = a.pr/avg.pr
      }else{
        a.pr = priors$count
        b.pr = a.pr/priors$avg
      }
      if(quarantine){
        vars = c("beta","gamma",
                 "alpha.SS","alpha.SI","alpha.II",
                 "omega.SS","omega.SI","omega.II")
      }else{
        vars = c("beta","gamma","alpha","omega")
      }
      Bayes.tab = NULL
    }
    if (MLE){
      MLE.tab = NULL
    }
    params = c(bet, gam, alpha.r, alpha.d)
    
    
    # initialize ER(p) graph if the initial netwok is not provided
    if(is.null(init.net)){
      init.net = sample_gnp(n = N, p = init.p, directed = F)
    }
    # convert to adjacency matrix
    adjmat = get.adjacency(init.net, type=c("both"), names=F)
    adjmat = as(adjmat,"dgTMatrix")
    # save the initial network structure (adj.mat)
    init.adj = adjmat
    
    
    ## check whether adjmat is symmetric and nnzero() returns an even number
    if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
      stop(paste("The initial adjmat isn't symmetric, where nonzero =", 
                 nnzero(adjmat),"\n"))
    }
    
    # initialize epidemic
    infected = sample(N, init.infec)
    epid = rep(0,N); epid[infected] = 1
    cat("Initially infected: ",infected,"\n")
    # save initial cases
    init.infected = infected
    
    # initialize data objects to keep track of
    susceptible = which(epid==0)
    N.I = init.infec; N.S = N - N.I
    t.cur = 0
    if(quarantine){
      N.c = numeric(3); N.d = numeric(3)
      adjmat.ii = adjmat[infected,infected]
      N.c[3] = nnzero(adjmat.ii)/2
      N.d[3] = N.I * (N.I - 1)/2 - N.c[3]
      # both need this:
      adjmat.si = adjmat[susceptible,infected]
      if(model=="SIR"){
        recovered = NULL; N.R = 0
        adjmat.rs = adjmat[union(susceptible,recovered), union(susceptible,recovered)]
        N.c[1] = nnzero(adjmat.rs)/2
        N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
        adjmat.ri = adjmat.si
        N.c[2] = nnzero(adjmat.ri)
        N.d[2] = N.S * N.I - N.c[2]
      }else{
        adjmat.ss = adjmat[susceptible,susceptible]
        N.c[1] = nnzero(adjmat.ss)/2
        N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
        
        N.c[2] = nnzero(adjmat.si)
        N.d[2] = N.S * N.I - N.c[2]
      }
      # at the start: N.R = 0, so this works
      N.SI = N.c[2]
    }else{
      N.c = nnzero(adjmat)/2
      N.d = N * (N-1)/2 - N.c
      
      adjmat.si = adjmat[susceptible,infected]
      N.SI = nnzero(adjmat.si)
    }
    
    
    # iterations
    while (t.cur < tmax) {
      # NOTES #
      # a. check to make sure N.S also > 0
      # b. handle edge cases (N.S = 1 and/or N.I = 1)
      
      # # check whether all N's are non-negative
      # Ns = c(N.SI, N.I, N.d, N.c)
      # if(any(Ns < 0) | all(Ns == 0)){
      #   msg = paste("Weird stuff happened! All the Ns =", paste(Ns,collapse = ", "), "\n")
      #   stop(msg)
      # }
      
      # 1. calculate sum of risks and sample the next event time
      Lambda = bet * N.SI + N.I * gam + sum(alpha.r * N.d) + sum(alpha.d * N.c)
      t.next = t.cur + rexp(1, rate = Lambda)
      
      # 2. determine the event type
      ratios = c(bet * N.SI, N.I * gam, alpha.r * N.d, alpha.d * N.c)
      # 2.1 deal with network change
      if(quarantine){
        
        z = sample(8, 1, prob = ratios)
        
        if(z==3){
          # reconnect a random S-S/S-R/R-R pair
          if(model == "SIS"){
            num.dis = (N.S - 1) - rowSums(adjmat.ss)
            from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
            to.cands = which(adjmat.ss[from.ind,]==0)
          }else{
            num.dis = (N.S + N.R - 1) - rowSums(adjmat.rs)
            from.ind = sample(N.S + N.R, 1, prob = num.dis); from = union(susceptible, recovered)[from.ind]
            to.cands = which(adjmat.rs[from.ind,]==0)
          }
          to.cands = to.cands[to.cands!=from.ind]
          if(length(to.cands) == 1){
            to.ind = to.cands
          }else{
            to.ind = sample(to.cands,1)
          }
          if(model == "SIS"){
            to = susceptible[to.ind]
            adjmat.ss[as.integer(from.ind),as.integer(to.ind)] = 1
            adjmat.ss[as.integer(to.ind),as.integer(from.ind)] = 1
            msg = paste0("reconnect an S-S pair, (",from,",",to,")")
          }else{
            to = union(susceptible, recovered)[to.ind]
            adjmat.rs[as.integer(from.ind),as.integer(to.ind)] = 1
            adjmat.rs[as.integer(to.ind),as.integer(from.ind)] = 1
            msg = paste0("reconnect an S-S/S-R/R-R pair, (",from,",",to,")")
          }
          adjmat[as.integer(from),as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          N.c[1] = N.c[1] + 1; N.d[1] = N.d[1] - 1
        }
        if(z==4){
          # reconnect a random S-I/R-I pair
          if(model == "SIS"){
            if (N.S == 1){
              from = susceptible
              to.ind = sample(which(adjmat.si==0),1); to = infected[to.ind]
            }else if (N.I == 1){
              num.dis = N.I - adjmat.si
              from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
              to = infected
            }else{
              num.dis = N.I - rowSums(adjmat.si)
              from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
              to.ind = sample(which(adjmat.si[from.ind,]==0),1); to = infected[to.ind]
            }
            msg = paste0("reconnect an S-I pair, (",from,",",to,").")
          }else{
            #N.SR = N.S + N.R; 
            # 03/27/2020: debug
            sus.recov = union(susceptible, recovered)
            N.SR = length(sus.recov)
            if (N.SR == 1){
              from = sus.recov
              # 03/27/2020: debug
              num.dis = N.SR - as.vector(adjmat.ri)
              to.ind = sample(N.I, 1, prob = num.dis); to = infected[to.ind]
              #to.ind = sample(which(adjmat.ri==0),1); to = infected[to.ind]
            }else if (N.I == 1){
              num.dis = N.I - as.vector(adjmat.ri)
              from.ind = sample(N.SR,1,prob = num.dis); from = sus.recov[from.ind]
              to = infected
            }else{
              num.dis = N.I - rowSums(adjmat.ri)
              # debugging...
              #cat(N.S, N.R, N.I, '\n')
              #cat(length(susceptible), length(recovered), length(infected),'\n')
              #print(dim(adjmat.ri))
              from.ind = sample(N.SR,1,prob = num.dis); from = sus.recov[from.ind]
              
              # 03/28/2020: debug
              vec.dis = 1 - as.vector(adjmat.ri[from.ind,])
              to.ind = sample(N.I,1,prob = vec.dis); to = infected[to.ind]
              #to.ind = sample(which(adjmat.ri[from.ind,]==0),1); to = infected[to.ind]
            }
            msg = paste0("reconnect an S-I/R-I pair, (",from,",",to,").")
          }
          adjmat[as.integer(from),as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          # 03/27/2020
          # try to fix a bug
          if(length(susceptible)>0){
            adjmat.si = adjmat[susceptible,infected]
            N.SI = nnzero(adjmat.si)
          }else{
            N.SI = 0
          }
          
          if(model == "SIR"){ adjmat.ri = adjmat[sus.recov,infected] }
          
          N.c[2] = N.c[2] + 1; N.d[2] = N.d[2] - 1
          #N.SI = N.c[2]
        }
        if(z==5){
          # reconnect a random I-I pair
          num.dis = (N.I - 1) - rowSums(adjmat.ii)
          from.ind = sample(N.I,1,prob = num.dis); from = infected[from.ind]
          to.cands = which(adjmat.ii[from.ind,]==0)
          to.cands = to.cands[to.cands!=from.ind]
          if(length(to.cands) == 1){
            to.ind = to.cands
          }else{
            to.ind = sample(to.cands,1)
          }
          to = infected[to.ind]
          adjmat.ii[as.integer(from.ind),as.integer(to.ind)] = 1
          adjmat.ii[as.integer(to.ind),as.integer(from.ind)] = 1
          adjmat[as.integer(from),as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          N.c[3] = N.c[3] + 1; N.d[3] = N.d[3] - 1
          
          msg = paste0("reconnect an I-I pair, (",from,",",to,").")
        }
        if(z==6){
          # disconnect a random S-S/R-S/R-R pair
          dis.ind = sample(N.c[1], 1) 
          if(model == "SIS"){
            from.ind = adjmat.ss@i[dis.ind] + 1; to.ind = adjmat.ss@j[dis.ind] + 1
            adjmat.ss[as.integer(from.ind),as.integer(to.ind)] = 0
            adjmat.ss[as.integer(to.ind),as.integer(from.ind)] = 0
            
            from = susceptible[from.ind]; to = susceptible[to.ind]
            
            msg = paste0("disconnect an S-S pair, (",from,",",to,").")
          }else{
            sus.recov = union(susceptible, recovered)
            from.ind = adjmat.rs@i[dis.ind] + 1; to.ind = adjmat.rs@j[dis.ind] + 1
            adjmat.rs[as.integer(from.ind),as.integer(to.ind)] = 0
            adjmat.rs[as.integer(to.ind),as.integer(from.ind)] = 0
            
            from = sus.recov[from.ind]; to = sus.recov[to.ind]
            
            msg = paste0("disconnect an S-S/S-R/R-R pair, (",from,",",to,").")
          }
          adjmat[as.integer(from),as.integer(to)] = 0
          adjmat[as.integer(to),as.integer(from)] = 0
          
          N.c[1] = N.c[1] - 1; N.d[1] = N.d[1] + 1
          
          # ## sometimes from/to is just empty... check it out
          # if(length(from)!=1){
          #   stop("from.ind = ", from.ind, 
          #        " but susceptible length is ", length(susceptible),
          #        " N.c[1] =", N.c[1], " dis.ind =", dis.ind, 
          #        " while adjmat.ss@i has length ", length(adjmat.ss@i))
          # }
          # ##
        }
        if(z==7){
          # disconnect a random S-I/R-I pair
          # dis.ind = sample(N.c[2], 1)
          if(model == "SIS"){
            if (N.S == 1){
              from = susceptible
              to.ind = sample(which(adjmat.si==1),1); to = infected[to.ind]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.si = adjmat[susceptible,infected]
            }else if (N.I == 1){
              num.con = as.vector(adjmat.si)
              from.ind = sample(N.S,1,prob = num.con); from = susceptible[from.ind]
              to = infected
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.si = adjmat[susceptible,infected]
            }else{
              dis.ind = sample(N.c[2], 1)
              from = adjmat.si@i[dis.ind] + 1; to = adjmat.si@j[dis.ind] + 1
              adjmat.si[as.integer(from),as.integer(to)] = 0
              
              from = susceptible[from]; to = infected[to]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
            }
            msg = paste0("disconnect an S-I pair, (",from,",",to,").")
          }else{
            N.SR = N.S + N.R; sus.recov = union(susceptible, recovered)
            if (N.SR == 1){
              from = sus.recov
              to.ind = sample(which(adjmat.ri==1),1); to = infected[to.ind]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.ri = adjmat[sus.recov,infected]
            }else if (N.I == 1){
              num.con = adjmat.ri
              from.ind = sample(N.SR,1,prob = num.con); from = sus.recov[from.ind]
              to = infected
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.ri = adjmat[sus.recov,infected]
            }else{
              dis.ind = sample(N.c[2], 1)
              from = adjmat.ri@i[dis.ind] + 1; to = adjmat.ri@j[dis.ind] + 1
              # debugging...
              #cat('dis.ind=',dis.ind,'\n')
              #cat('from=',from,', to=',to, '\n')
              #cat(dim(adjmat.ri), '\n')
              #cat(length(adjmat.ri@i),  length(adjmat.ri@j), N.c[2], '\n')
              
              if(length(adjmat.ri@i) != N.c[2]){
                cat('Discrepancy when disconnecting H-I...\n')
                stop(paste0("At time ", t.next, ", ",msg,"\n"))
              }
              
              #adjmat.ri[as.integer(from),as.integer(to)] = 0
              
              from = sus.recov[from]; to = infected[to]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.ri = adjmat[sus.recov,infected]
            }
            if(length(susceptible)>0){
              adjmat.si = adjmat[susceptible,infected]
            }else{
              adjmat.si = NULL
            }
            msg = paste0("disconnect an S-I/R-I pair, (",from,",",to,").")
          }
          
          N.c[2] = N.c[2] - 1; N.d[2] = N.d[2] + 1
          N.SI = nnzero(adjmat.si)
          
        }
        if(z==8){
          # disconnect a random I-I pair
          dis.ind = sample(N.c[3], 1) 
          from = adjmat.ii@i[dis.ind] + 1; to = adjmat.ii@j[dis.ind] + 1
          adjmat.ii[as.integer(from),as.integer(to)] = 0
          adjmat.ii[as.integer(to),as.integer(from)] = 0
          
          from = infected[from]; to = infected[to]
          adjmat[as.integer(from),as.integer(to)] = 0
          adjmat[as.integer(to),as.integer(from)] = 0
          
          N.c[3] = N.c[3] - 1; N.d[3] = N.d[3] + 1
          
          msg = paste0("disconnect an I-I pair, (",from,",",to,").")
        }
      }else{
        z = sample(4, 1, prob = ratios)
        if(z==3){
          # reconnect a random pair
          num.dis = (N - 1) - rowSums(adjmat)
          from = sample(N,1,prob = num.dis)
          to.cands = which(adjmat[from,]==0)
          to.cands = to.cands[to.cands!=from]
          if(length(to.cands) == 1){
            to = to.cands
          }else{
            to= sample(to.cands,1)
          }
          adjmat[as.integer(from), as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          
          N.c = N.c + 1; N.d = N.d - 1
          
          # UPDATE: if any of S or I is empty now
          if(length(susceptible) == 0 | length(infected) == 0){
            adjmat.si = NULL; N.SI = 0
          }else{
            adjmat.si = adjmat[susceptible,infected]
            N.SI = nnzero(adjmat.si)
          }
          
          msg = paste0("reconnect a pair, (",from,",",to,").")
        }
        if(z==4){
          # disconnect a random pair
          dis.ind = sample(N.c, 1)
          from = adjmat@i[dis.ind] + 1; to = adjmat@j[dis.ind] + 1
          adjmat[from,to] = 0; adjmat[to,from] = 0
          N.c = N.c - 1; N.d = N.d + 1
          
          # UPDATE: if any of S or I is empty now
          if(length(susceptible)==0 | length(infected)==0){
            adjmat.si = NULL; N.SI = 0
          }else{
            adjmat.si = adjmat[susceptible,infected]
            N.SI = nnzero(adjmat.si)
          }
          
          msg = paste0("disconnect a pair, (",from,",",to,").")
        }
      }
      
      #2.2 deal with epidemic progress
      if(z==2){
        # recover an infected person
        rec.ind = sample(N.I,1)
        recover = infected[rec.ind]
        infected = infected[-rec.ind]
        N.I = N.I - 1 
        if(model == "SIS"){
          susceptible = c(susceptible,recover)
          epid[recover] = 0
          N.S = N.S + 1
        }else{
          epid[recover] = -1
          if(quarantine){
            recovered = c(recovered, recover)
            N.R = N.R + 1
          }
        }
        
        msg = paste0("individual ", recover,
                     " is recovered. Disease prevalence = ",
                     N.I/N, ".")
      }
      if(z==1){
        # select a susceptible person to infect
        if(N.S == 1){
          # only one person to infect
          new.infec = susceptible
          susceptible = NULL
        }else{
          # multiple susceptible inviduals to choose from
          if(N.I == 1){
            # only one infected person
            num.nei.i = adjmat.si
          }else{
            num.nei.i = rowSums(adjmat.si)
          }
          infec.ind = sample(N.S,1,prob = num.nei.i)
          new.infec = susceptible[infec.ind]
          susceptible = susceptible[-infec.ind]
        }
        epid[new.infec] = 1
        N.S = N.S - 1
        infected = c(infected, new.infec); N.I = N.I + 1
        
        msg = paste0("individual ", new.infec,
                     " is infected. Disease prevalence = ",
                     N.I/N, ".")
      }
      
      # A sanity check: if N.I==0, stop the simulation
      if(N.I==0){
        if(verbose){
          cat(paste0("At time ", t.next, ", "),msg,"\n")
        }
        # still update inference statistics
        del.t = t.next - t.cur
        z = as.integer(z)
        event.counts[z] = event.counts[z] + 1
        tot.event = tot.event + 1
        
        big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
        
        # still have to record the event
        time = c(time, t.next)
        event.type = c(event.type, z)
        preval = c(preval, N.I)
        dens = c(dens, nnzero(adjmat))
        # here: a recovery must be the latest event
        per1 = c(per1, recover)
        per2 = c(per2, NA)
        
        cat("Disease extinct at time", t.next, ", simulation stops.\n")
        break
      }
      
      
      # # check whether adjmat is still symmetric--things have gone wrong somewhere
      # if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
      #   stop(paste("The adjmat isn't symmetric any more, where nonzero =", 
      #              nnzero(adjmat), "and that's after z =", z,
      #              "from.ind =", from.ind, "to.ind = ", to.ind, 
      #              "to.cands =", paste(to.cands,collapse = ","),
      #              "from =", from, "to =", to,"\n"))
      # }
      
      
      # 2.2.* make updates about S-I, S-S, I-I connections if z==1 or 2
      # while handling edge cases of N.S <= 1 and/or N.I == 1
      if(z %in% c(1,2)){
        if(N.S == 0){
          adjmat.si = NULL
          N.SI = 0
          if(quarantine && model == "SIR" && N.R == 0){
            adjmat.ri = NULL
            N.c[2] = 0
          }
          # 03/27/2020: debug
          if(quarantine && model == "SIR" && N.R > 0){
            adjmat.ri = adjmat[recovered,infected]
            N.c[2] = nnzero(adjmat.ri)
          }
        }else{
          adjmat.si = adjmat[susceptible,infected]
          N.SI = nnzero(adjmat.si)
          if(quarantine && model == "SIR"){
            adjmat.ri = adjmat[union(susceptible,recovered),infected]
            N.c[2] = nnzero(adjmat.ri)
          }
        }
        if(quarantine){
          if(model == "SIS"){
            if(N.S <= 1){
              adjmat.ss = 0
              N.c[1] = 0
              N.d[1] = 0
            }else{
              adjmat.ss = adjmat[susceptible,susceptible]
              N.c[1] = nnzero(adjmat.ss)/2
              N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
            }
            N.c[2] = N.SI
            N.d[2] = N.S * N.I - N.c[2]
          }else{
            #N.SR = N.S + N.R; 
            sus.recov = union(susceptible,recovered)
            # 03/27/2020: debug
            N.SR = length(sus.recov)
            
            if(N.SR <= 1){
              adjmat.rs = 0
              N.c[1] = 0
              N.d[1] = 0
            }else{
              adjmat.rs = adjmat[sus.recov,sus.recov]
              N.c[1] = nnzero(adjmat.rs)/2
              N.d[1] = N.SR * (N.SR - 1)/2 - N.c[1]
            }
            N.d[2] = N.SR * N.I - N.c[2]
          }
          if(N.I == 1){
            adjmat.ii = 0
            N.c[3] = 0
            N.d[3] = 0
          }else{
            adjmat.ii = adjmat[infected,infected]
            N.c[3] = nnzero(adjmat.ii)/2
            N.d[3] = N.I * (N.I - 1)/2 - N.c[3]
          }
        }
      }
      
      # 3. update inference statistics and do inference occasionally
      # 3.1 update stats
      del.t = t.next - t.cur
      z = as.integer(z)
      event.counts[z] = event.counts[z] + 1
      tot.event = tot.event + 1
      
      big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
      
      # 3.2 inference if ...
      if(tot.event %% infer.interval == 0){
        if(Bayes){
          # Bayesian inference: Gamma sampling
          a.post = a.pr + event.counts
          b.post = b.pr + big.sums
          cat("With", tot.event,"events, posterior means and variances:\n",
              a.post/b.post, "\n",
              a.post/(b.post^2), "\n")
          # record posterior stats
          if(return.infer){
            m.post = a.post/b.post
            lb.post = qgamma(.025, a.post, b.post)
            ub.post = qgamma(.975, a.post, b.post)
            this.chunk = cbind(m.post, lb.post, ub.post, event.counts, tot.event)
            Bayes.tab = rbind(Bayes.tab, this.chunk)
          }
          # plot
          if(infer.plot){
            par(mfrow=c(2,2))
            for(ix in 1:length(params)){
              v = vars[ix]; p = params[ix]
              sams = rgamma(samples, shape=a.post[ix], rate = b.post[ix])
              xl = range(sams, p)
              hist(sams, main=paste("Posterior samples for",vars[ix]),
                   xlab = vars[ix], xlim=xl, col="lightblue")
              abline(v=p,col="red")
              if(demo.sleep & ix%%4 == 0){
                Sys.sleep(0.2)
              }
            }
          }
        }
        if (MLE){
          # MLE
          MLEs = event.counts/big.sums
          if(infer.verbose){
            cat("With",tot.event, "events, MLEs:","\n", MLEs, "\n")
          }
          MLE.tab = rbind(MLE.tab, c(tot.event, MLEs, event.counts))
        }
      }
      
      # updated 02/04/2020
      # 3.3 examine network degree distribution if
      if(net.deg.plot & (tot.event %% net.interval == 0)){
        degs = rowSums(adjmat)
        avg.deg = mean(degs)
        poi.degs = rpois(N, avg.deg)
        
        deg.dat = data.frame(Degrees = degs, DegPoi = poi.degs)
        print(
          ggplot(deg.dat) + geom_density(aes(x = Degrees, color="Actual degree")) +
            geom_density(aes(x = DegPoi, color="Poisson")) +
            labs(color="Distribution", 
                 caption = paste("At t =", round(t.cur, digits = 2))) + 
            theme_bw(base_size = 14) %+replace%
            theme(legend.position = c(0.8,0.8))
        )
      }
      
      
      
      # 4. report the event and summarize and set t.cur
      
      t.cur = t.next
      
      # record the event
      time = c(time, t.cur)
      event.type = c(event.type, z)
      preval = c(preval, N.I)
      dens = c(dens, nnzero(adjmat))
      # record labels of individual(s) involved
      if (z >= 3){
        per1 = c(per1, from)
        per2 = c(per2, to)
      }else if (z==2){
        per1 = c(per1, recover)
        per2 = c(per2, NA)
      }else{
        per1 = c(per1, new.infec)
        per2 = c(per2, NA)
      }
      
      
      if(verbose){
        cat(paste0("At time ", t.cur, ", "),msg,"\n")
      }
      
    }
    
    # return results/summary
    cat("Simulation finished. At the end: \n", 
        "Disease prevalence =", N.I/N, "\n",
        "Network density = ", nnzero(adjmat)/(N*(N-1)),"\n")
    event.log = data.frame(time = time, event = event.type,
                           per1 = per1, per2 = per2,
                           preval = preval/N, dens = dens/(N*(N-1)))
    
    # inference at the end
    if(Bayes){
      # Bayesian inference: Gamma sampling
      a.post = a.pr + event.counts
      b.post = b.pr + big.sums
      cat("Posterior means and variances with", tot.event, "events in total:\n",
          a.post/b.post, "\n",
          a.post/(b.post^2), "\n")
      cat("Event type counts:\n", event.counts,"\n")
      # record posterior stats
      if(return.infer){
        m.post = a.post/b.post
        lb.post = qgamma(.025, a.post, b.post)
        ub.post = qgamma(.975, a.post, b.post)
        this.chunk = cbind(m.post, lb.post, ub.post, event.counts, tot.event)
        Bayes.tab = rbind(Bayes.tab, this.chunk)
      }
      # plot
      if(infer.plot){
        par(mfrow=c(2,2))
        for(ix in 1:length(params)){
          v = vars[ix]; p = params[ix]
          sams = rgamma(samples, shape=a.post[ix], rate = b.post[ix])
          xl = range(sams, p)
          hist(sams, main=paste("Posterior samples for",vars[ix]),
               xlab = vars[ix], xlim=xl, col="lightblue")
          abline(v=p,col="red")
          if(demo.sleep & ix%%4 == 0){
            Sys.sleep(0.2)
          }
        }
      }
    }
    if(MLE){
      # MLE
      MLEs = event.counts/big.sums
      cat("MLEs with", tot.event, "events in total:","\n", MLEs, "\n")
      cat("Event type counts:\n", event.counts,"\n")
      
      MLE.tab = rbind(MLE.tab, c(tot.event, MLEs, event.counts))
      MLE.tab = as.data.frame(MLE.tab)
      if(quarantine){
        names(MLE.tab) = c("event.num","beta","gamma",
                           "alpha.SS","alpha.SI","alpha.II",
                           "omega.SS","omega.SI","omega.II",
                           "infection","recovery",
                           "SS.connect","SI.connect","II.connect",
                           "SS.discon","SI.discon","II.discon")
      }else{
        names(MLE.tab) = c("event.num","beta","gamma","alpha","omega",
                           "infection","recovery","connect","discon")
      }
    }
    
    if(plot){
      subt = paste(model," beta =",bet," gamma =", gam, 
                   " alpha.r =", alpha.r, " alpha.d =", alpha.d)
      plot(preval ~ time, data = event.log, 
           xlab="Time", ylab="Disease prevalence", 
           main = "Disease Prevalence vs. Time", sub = subt,
           ylim = c(0,1), type="l",lwd=2)
      plot(dens ~ time, data = event.log,
           xlab="Time", ylab="Edge density", 
           main = "Network Edge Density vs. Time", sub = subt,
           type="l",lwd=2)
    }
    
    if(infer.plot){
      #params = c(bet, gam, alpha.r, alpha.d)
      if(Bayes){
        par(mfrow=c(1,1))
      }else{
        # a function to calculate the frequentist 95% CI for MLEs
        MLE.CI <- function(x){
          l.hat = x[1]; n.obs = x[2]
          if (is.na(l.hat) || l.hat == 0){
            lb = 0; ub = 0
          }else if (n.obs >= 20){
            lb = l.hat * (1 + qnorm(.025)/sqrt(n.obs))
            ub = l.hat * (1 + qnorm(.975)/sqrt(n.obs))
          }else{
            lb = l.hat * qchisq(.025, 2*n.obs) / (2*n.obs)
            ub = l.hat * qchisq(.975, 2*n.obs) / (2*n.obs)
          }
          return(c(lb,ub))
        }
        
        lcol = length(names(MLE.tab)[-1])
        vars = names(MLE.tab)[-1][1:(lcol/2)]
        counts = names(MLE.tab)[-1][(lcol/2 + 1):lcol]
        e.num = MLE.tab$event.num
        if(quarantine){
          par(mfrow=c(2,2))
        }else{
          par(mfrow=c(2,2))
        }
        for(ix in 1:length(vars)){
          v = vars[ix]
          if(all(is.na(MLE.tab[[v]]))){
            plot.new()
          }else{
            tab.v = MLE.tab[,c(v,counts[ix])]
            bounds.v = t(apply(tab.v, 1, MLE.CI))
            bounds.v = as.data.frame(bounds.v)
            names(bounds.v) = c("lb","ub")
            #s.m = mean(MLE.tab[[v]], na.rm=T)
            #s.sd = sd(MLE.tab[[v]], na.rm=T)
            yl = range(min(bounds.v$lb), max(bounds.v$ub), params[ix])
            plot(MLE.tab[[v]]~e.num, xlab="# events", ylab=v,
                 main = paste("MLE for", v, "v.s. # events"), 
                 type="l", ylim=yl)
            lines(bounds.v$lb~e.num, lty=2, col="gray")
            lines(bounds.v$ub~e.num, lty=2, col="gray")
          }
          abline(h=params[ix], col="red")
        }
        par(mfrow=c(1,1))
      }
    }
    
    # add: print out big.sums if Bayesian
    if(Bayes){
      cat("big.sums: \n", big.sums,"\n")
    }
    
    # 02/04/2020: degree distribution at the end
    if(net.deg.plot){
      degs = rowSums(adjmat)
      avg.deg = mean(degs)
      poi.degs = rpois(N, avg.deg)
      
      deg.dat = data.frame(Degrees = degs, DegPoi = poi.degs)
      print(
        ggplot(deg.dat) + geom_density(aes(x = Degrees, color="Actual degree")) +
          geom_density(aes(x = DegPoi, color="Poisson")) +
          labs(color="Distribution", 
               caption = paste("At the end, t =", round(t.cur, digits = 2))) + 
          theme_bw(base_size = 14) %+replace%
          theme(legend.position = c(0.8,0.8))
      )
    }
    
    
    # return event data
    # MODIFICATION: return true parametes too
    if(return.infer){
      if(Bayes){
        Bayes.tab = as.data.frame(Bayes.tab)
        Bayes.tab$var = vars
        infer.tab = Bayes.tab
      }else{
        infer.tab = MLE.tab
      }
      return(list(G0 = init.adj, I0 = init.infected,
                  events = event.log, infer=infer.tab,
                  infer2 = MLE.tab,
                  truth = params))
    }else{
      return(list(G0 = init.adj, I0 = init.infected,
                  events = event.log,truth = params))
    }
    
  }

# Rolling output suppressed
stochastic_coevolve_infer4 <- 
  function(N=50, tmax=100, 
           bet=0.05, gam=0.2, init.infec = 1, model = "SIS",
           alpha.r=0.005, alpha.d=0.05, quarantine=F,
           init.net = NULL, init.p = 0.1,
           verbose = F, plot = F,
           infer.interval = 100, 
           Bayes = T, MLE = T, priors = NULL, samples = 500,
           infer.verbose = F, infer.plot = F, return.infer = F,
           net.deg.plot = F, net.interval = 100,
           demo.sleep = F){
    ## alpha.* = c(SSrate, SIrate, IIrate) ##
    #set.seed(seed)
    
    # check edge rates validity
    if(quarantine & (length(alpha.r) < 3 | length(alpha.d) < 3)){
      stop("Under quarantine, need to specify 3 values for edge connection/disconnection rates respectively!")
    }
    
    # output basic info
    # cat("N =", N, " Model:", model,"\n",
    #     "beta =", bet, " gamma =", gam, "\n",
    #     "init.infec =", init.infec, " init.p =", init.p, "\n",
    #     "alpha.r =", alpha.r, "\n", "alpha.d =", alpha.d, "\n",
    #     "quarantine?", quarantine,"\n",
    #     "Bayesian inference?", Bayes, "\n",
    #     "inference interval =", infer.interval, "\n")
    
    # summarization stats to be recorded
    time = NULL
    event.type = NULL
    preval = NULL
    dens = NULL
    # add involved individual labels 
    per1 = NULL
    per2 = NULL
    
    # more stats to record for inference
    if(quarantine){
      event.counts = numeric(8)
      big.sums = numeric(8)
    }else{
      event.counts = numeric(4)
      big.sums = numeric(4)
    }
    tot.event = 0
    
    if(Bayes){
      if(quarantine & nrow(priors)==4){
        a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
        avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
        b.pr = a.pr/avg.pr
      }else{
        a.pr = priors$count
        b.pr = a.pr/priors$avg
      }
      if(quarantine){
        vars = c("beta","gamma",
                 "alpha.SS","alpha.SI","alpha.II",
                 "omega.SS","omega.SI","omega.II")
      }else{
        vars = c("beta","gamma","alpha","omega")
      }
      Bayes.tab = NULL
    }
    if (MLE){
      MLE.tab = NULL
    }
    params = c(bet, gam, alpha.r, alpha.d)
    
    
    # initialize ER(p) graph if the initial netwok is not provided
    if(is.null(init.net)){
      init.net = sample_gnp(n = N, p = init.p, directed = F)
    }
    # convert to adjacency matrix
    adjmat = get.adjacency(init.net, type=c("both"), names=F)
    adjmat = as(adjmat,"dgTMatrix")
    # save the initial network structure (adj.mat)
    init.adj = adjmat
    
    
    ## check whether adjmat is symmetric and nnzero() returns an even number
    if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
      stop(paste("The initial adjmat isn't symmetric, where nonzero =", 
                 nnzero(adjmat),"\n"))
    }
    
    # initialize epidemic
    infected = sample(N, init.infec)
    epid = rep(0,N); epid[infected] = 1
    cat("Initially infected: ",infected,"\n")
    # save initial cases
    init.infected = infected
    
    # initialize data objects to keep track of
    susceptible = which(epid==0)
    N.I = init.infec; N.S = N - N.I
    t.cur = 0
    if(quarantine){
      N.c = numeric(3); N.d = numeric(3)
      adjmat.ii = adjmat[infected,infected]
      N.c[3] = nnzero(adjmat.ii)/2
      N.d[3] = N.I * (N.I - 1)/2 - N.c[3]
      # both need this:
      adjmat.si = adjmat[susceptible,infected]
      if(model=="SIR"){
        recovered = NULL; N.R = 0
        adjmat.rs = adjmat[union(susceptible,recovered), union(susceptible,recovered)]
        N.c[1] = nnzero(adjmat.rs)/2
        N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
        adjmat.ri = adjmat.si
        N.c[2] = nnzero(adjmat.ri)
        N.d[2] = N.S * N.I - N.c[2]
      }else{
        adjmat.ss = adjmat[susceptible,susceptible]
        N.c[1] = nnzero(adjmat.ss)/2
        N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
        
        N.c[2] = nnzero(adjmat.si)
        N.d[2] = N.S * N.I - N.c[2]
      }
      # at the start: N.R = 0, so this works
      N.SI = N.c[2]
    }else{
      N.c = nnzero(adjmat)/2
      N.d = N * (N-1)/2 - N.c
      
      adjmat.si = adjmat[susceptible,infected]
      N.SI = nnzero(adjmat.si)
    }
    
    
    # iterations
    counter = 10
    while (t.cur < tmax) {
      
      # NOTES #
      # a. check to make sure N.S also > 0
      # b. handle edge cases (N.S = 1 and/or N.I = 1)
      
      # # check whether all N's are non-negative
      # Ns = c(N.SI, N.I, N.d, N.c)
      # if(any(Ns < 0) | all(Ns == 0)){
      #   msg = paste("Weird stuff happened! All the Ns =", paste(Ns,collapse = ", "), "\n")
      #   stop(msg)
      # }
      
      # 1. calculate sum of risks and sample the next event time
      Lambda = bet * N.SI + N.I * gam + sum(alpha.r * N.d) + sum(alpha.d * N.c)
      t.next = t.cur + rexp(1, rate = Lambda)
      
      if ((t.cur < counter) & (t.next > counter)) {
        print(paste0(counter, " Days have passed"))
        counter = counter + 10
      }
        
      # 2. determine the event type
      ratios = c(bet * N.SI, N.I * gam, alpha.r * N.d, alpha.d * N.c)
      # 2.1 deal with network change
      if(quarantine){
        
        z = sample(8, 1, prob = ratios)
        
        if(z==3){
          # reconnect a random S-S/S-R/R-R pair
          if(model == "SIS"){
            num.dis = (N.S - 1) - rowSums(adjmat.ss)
            from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
            to.cands = which(adjmat.ss[from.ind,]==0)
          }else{
            num.dis = (N.S + N.R - 1) - rowSums(adjmat.rs)
            from.ind = sample(N.S + N.R, 1, prob = num.dis); from = union(susceptible, recovered)[from.ind]
            to.cands = which(adjmat.rs[from.ind,]==0)
          }
          to.cands = to.cands[to.cands!=from.ind]
          if(length(to.cands) == 1){
            to.ind = to.cands
          }else{
            to.ind = sample(to.cands,1)
          }
          if(model == "SIS"){
            to = susceptible[to.ind]
            adjmat.ss[as.integer(from.ind),as.integer(to.ind)] = 1
            adjmat.ss[as.integer(to.ind),as.integer(from.ind)] = 1
            msg = paste0("reconnect an S-S pair, (",from,",",to,")")
          }else{
            to = union(susceptible, recovered)[to.ind]
            adjmat.rs[as.integer(from.ind),as.integer(to.ind)] = 1
            adjmat.rs[as.integer(to.ind),as.integer(from.ind)] = 1
            msg = paste0("reconnect an S-S/S-R/R-R pair, (",from,",",to,")")
          }
          adjmat[as.integer(from),as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          N.c[1] = N.c[1] + 1; N.d[1] = N.d[1] - 1
        }
        if(z==4){
          # reconnect a random S-I/R-I pair
          if(model == "SIS"){
            if (N.S == 1){
              from = susceptible
              to.ind = sample(which(adjmat.si==0),1); to = infected[to.ind]
            }else if (N.I == 1){
              num.dis = N.I - adjmat.si
              from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
              to = infected
            }else{
              num.dis = N.I - rowSums(adjmat.si)
              from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
              to.ind = sample(which(adjmat.si[from.ind,]==0),1); to = infected[to.ind]
            }
            msg = paste0("reconnect an S-I pair, (",from,",",to,").")
          }else{
            #N.SR = N.S + N.R; 
            # 03/27/2020: debug
            sus.recov = union(susceptible, recovered)
            N.SR = length(sus.recov)
            if (N.SR == 1){
              from = sus.recov
              # 03/27/2020: debug
              num.dis = N.SR - as.vector(adjmat.ri)
              to.ind = sample(N.I, 1, prob = num.dis); to = infected[to.ind]
              #to.ind = sample(which(adjmat.ri==0),1); to = infected[to.ind]
            }else if (N.I == 1){
              num.dis = N.I - as.vector(adjmat.ri)
              from.ind = sample(N.SR,1,prob = num.dis); from = sus.recov[from.ind]
              to = infected
            }else{
              num.dis = N.I - rowSums(adjmat.ri)
              # debugging...
              #cat(N.S, N.R, N.I, '\n')
              #cat(length(susceptible), length(recovered), length(infected),'\n')
              #print(dim(adjmat.ri))
              from.ind = sample(N.SR,1,prob = num.dis); from = sus.recov[from.ind]
              
              # 03/28/2020: debug
              vec.dis = 1 - as.vector(adjmat.ri[from.ind,])
              to.ind = sample(N.I,1,prob = vec.dis); to = infected[to.ind]
              #to.ind = sample(which(adjmat.ri[from.ind,]==0),1); to = infected[to.ind]
            }
            msg = paste0("reconnect an S-I/R-I pair, (",from,",",to,").")
          }
          adjmat[as.integer(from),as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          # 03/27/2020
          # try to fix a bug
          if(length(susceptible)>0){
            adjmat.si = adjmat[susceptible,infected]
            N.SI = nnzero(adjmat.si)
          }else{
            N.SI = 0
          }
          
          if(model == "SIR"){ adjmat.ri = adjmat[sus.recov,infected] }
          
          N.c[2] = N.c[2] + 1; N.d[2] = N.d[2] - 1
          #N.SI = N.c[2]
        }
        if(z==5){
          # reconnect a random I-I pair
          num.dis = (N.I - 1) - rowSums(adjmat.ii)
          from.ind = sample(N.I,1,prob = num.dis); from = infected[from.ind]
          to.cands = which(adjmat.ii[from.ind,]==0)
          to.cands = to.cands[to.cands!=from.ind]
          if(length(to.cands) == 1){
            to.ind = to.cands
          }else{
            to.ind = sample(to.cands,1)
          }
          to = infected[to.ind]
          adjmat.ii[as.integer(from.ind),as.integer(to.ind)] = 1
          adjmat.ii[as.integer(to.ind),as.integer(from.ind)] = 1
          adjmat[as.integer(from),as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          N.c[3] = N.c[3] + 1; N.d[3] = N.d[3] - 1
          
          msg = paste0("reconnect an I-I pair, (",from,",",to,").")
        }
        if(z==6){
          # disconnect a random S-S/R-S/R-R pair
          dis.ind = sample(N.c[1], 1) 
          if(model == "SIS"){
            from.ind = adjmat.ss@i[dis.ind] + 1; to.ind = adjmat.ss@j[dis.ind] + 1
            adjmat.ss[as.integer(from.ind),as.integer(to.ind)] = 0
            adjmat.ss[as.integer(to.ind),as.integer(from.ind)] = 0
            
            from = susceptible[from.ind]; to = susceptible[to.ind]
            
            msg = paste0("disconnect an S-S pair, (",from,",",to,").")
          }else{
            sus.recov = union(susceptible, recovered)
            from.ind = adjmat.rs@i[dis.ind] + 1; to.ind = adjmat.rs@j[dis.ind] + 1
            adjmat.rs[as.integer(from.ind),as.integer(to.ind)] = 0
            adjmat.rs[as.integer(to.ind),as.integer(from.ind)] = 0
            
            from = sus.recov[from.ind]; to = sus.recov[to.ind]
            
            msg = paste0("disconnect an S-S/S-R/R-R pair, (",from,",",to,").")
          }
          adjmat[as.integer(from),as.integer(to)] = 0
          adjmat[as.integer(to),as.integer(from)] = 0
          
          N.c[1] = N.c[1] - 1; N.d[1] = N.d[1] + 1
          
          # ## sometimes from/to is just empty... check it out
          # if(length(from)!=1){
          #   stop("from.ind = ", from.ind, 
          #        " but susceptible length is ", length(susceptible),
          #        " N.c[1] =", N.c[1], " dis.ind =", dis.ind, 
          #        " while adjmat.ss@i has length ", length(adjmat.ss@i))
          # }
          # ##
        }
        if(z==7){
          # disconnect a random S-I/R-I pair
          # dis.ind = sample(N.c[2], 1)
          if(model == "SIS"){
            if (N.S == 1){
              from = susceptible
              to.ind = sample(which(adjmat.si==1),1); to = infected[to.ind]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.si = adjmat[susceptible,infected]
            }else if (N.I == 1){
              num.con = as.vector(adjmat.si)
              from.ind = sample(N.S,1,prob = num.con); from = susceptible[from.ind]
              to = infected
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.si = adjmat[susceptible,infected]
            }else{
              dis.ind = sample(N.c[2], 1)
              from = adjmat.si@i[dis.ind] + 1; to = adjmat.si@j[dis.ind] + 1
              adjmat.si[as.integer(from),as.integer(to)] = 0
              
              from = susceptible[from]; to = infected[to]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
            }
            msg = paste0("disconnect an S-I pair, (",from,",",to,").")
          }else{
            N.SR = N.S + N.R; sus.recov = union(susceptible, recovered)
            if (N.SR == 1){
              from = sus.recov
              to.ind = sample(which(adjmat.ri==1),1); to = infected[to.ind]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.ri = adjmat[sus.recov,infected]
            }else if (N.I == 1){
              num.con = adjmat.ri
              from.ind = sample(N.SR,1,prob = num.con); from = sus.recov[from.ind]
              to = infected
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.ri = adjmat[sus.recov,infected]
            }else{
              dis.ind = sample(N.c[2], 1)
              from = adjmat.ri@i[dis.ind] + 1; to = adjmat.ri@j[dis.ind] + 1
              # debugging...
              #cat('dis.ind=',dis.ind,'\n')
              #cat('from=',from,', to=',to, '\n')
              #cat(dim(adjmat.ri), '\n')
              #cat(length(adjmat.ri@i),  length(adjmat.ri@j), N.c[2], '\n')
              
              if(length(adjmat.ri@i) != N.c[2]){
                cat('Discrepancy when disconnecting H-I...\n')
                stop(paste0("At time ", t.next, ", ",msg,"\n"))
              }
              
              #adjmat.ri[as.integer(from),as.integer(to)] = 0
              
              from = sus.recov[from]; to = infected[to]
              adjmat[as.integer(from),as.integer(to)] = 0
              adjmat[as.integer(to),as.integer(from)] = 0
              adjmat.ri = adjmat[sus.recov,infected]
            }
            if(length(susceptible)>0){
              adjmat.si = adjmat[susceptible,infected]
            }else{
              adjmat.si = NULL
            }
            msg = paste0("disconnect an S-I/R-I pair, (",from,",",to,").")
          }
          
          N.c[2] = N.c[2] - 1; N.d[2] = N.d[2] + 1
          N.SI = nnzero(adjmat.si)
          
        }
        if(z==8){
          # disconnect a random I-I pair
          dis.ind = sample(N.c[3], 1) 
          from = adjmat.ii@i[dis.ind] + 1; to = adjmat.ii@j[dis.ind] + 1
          adjmat.ii[as.integer(from),as.integer(to)] = 0
          adjmat.ii[as.integer(to),as.integer(from)] = 0
          
          from = infected[from]; to = infected[to]
          adjmat[as.integer(from),as.integer(to)] = 0
          adjmat[as.integer(to),as.integer(from)] = 0
          
          N.c[3] = N.c[3] - 1; N.d[3] = N.d[3] + 1
          
          msg = paste0("disconnect an I-I pair, (",from,",",to,").")
        }
      }else{
        z = sample(4, 1, prob = ratios)
        if(z==3){
          # reconnect a random pair
          num.dis = (N - 1) - rowSums(adjmat)
          from = sample(N,1,prob = num.dis)
          to.cands = which(adjmat[from,]==0)
          to.cands = to.cands[to.cands!=from]
          if(length(to.cands) == 1){
            to = to.cands
          }else{
            to= sample(to.cands,1)
          }
          adjmat[as.integer(from), as.integer(to)] = 1
          adjmat[as.integer(to),as.integer(from)] = 1
          
          N.c = N.c + 1; N.d = N.d - 1
          
          # UPDATE: if any of S or I is empty now
          if(length(susceptible) == 0 | length(infected) == 0){
            adjmat.si = NULL; N.SI = 0
          }else{
            adjmat.si = adjmat[susceptible,infected]
            N.SI = nnzero(adjmat.si)
          }
          
          msg = paste0("reconnect a pair, (",from,",",to,").")
        }
        if(z==4){
          # disconnect a random pair
          dis.ind = sample(N.c, 1)
          from = adjmat@i[dis.ind] + 1; to = adjmat@j[dis.ind] + 1
          adjmat[from,to] = 0; adjmat[to,from] = 0
          N.c = N.c - 1; N.d = N.d + 1
          
          # UPDATE: if any of S or I is empty now
          if(length(susceptible)==0 | length(infected)==0){
            adjmat.si = NULL; N.SI = 0
          }else{
            adjmat.si = adjmat[susceptible,infected]
            N.SI = nnzero(adjmat.si)
          }
          
          msg = paste0("disconnect a pair, (",from,",",to,").")
        }
      }
      
      #2.2 deal with epidemic progress
      if(z==2){
        # recover an infected person
        rec.ind = sample(N.I,1)
        recover = infected[rec.ind]
        infected = infected[-rec.ind]
        N.I = N.I - 1 
        if(model == "SIS"){
          susceptible = c(susceptible,recover)
          epid[recover] = 0
          N.S = N.S + 1
        }else{
          epid[recover] = -1
          if(quarantine){
            recovered = c(recovered, recover)
            N.R = N.R + 1
          }
        }
        
        msg = paste0("individual ", recover,
                     " is recovered. Disease prevalence = ",
                     N.I/N, ".")
      }
      if(z==1){
        # select a susceptible person to infect
        if(N.S == 1){
          # only one person to infect
          new.infec = susceptible
          susceptible = NULL
        }else{
          # multiple susceptible inviduals to choose from
          if(N.I == 1){
            # only one infected person
            num.nei.i = adjmat.si
          }else{
            num.nei.i = rowSums(adjmat.si)
          }
          infec.ind = sample(N.S,1,prob = num.nei.i)
          new.infec = susceptible[infec.ind]
          susceptible = susceptible[-infec.ind]
        }
        epid[new.infec] = 1
        N.S = N.S - 1
        infected = c(infected, new.infec); N.I = N.I + 1
        
        msg = paste0("individual ", new.infec,
                     " is infected. Disease prevalence = ",
                     N.I/N, ".")
      }
      
      # A sanity check: if N.I==0, stop the simulation
      if(N.I==0){
        if(verbose){
          cat(paste0("At time ", t.next, ", "),msg,"\n")
        }
        # still update inference statistics
        del.t = t.next - t.cur
        z = as.integer(z)
        event.counts[z] = event.counts[z] + 1
        tot.event = tot.event + 1
        
        big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
        
        # still have to record the event
        time = c(time, t.next)
        event.type = c(event.type, z)
        preval = c(preval, N.I)
        dens = c(dens, nnzero(adjmat))
        # here: a recovery must be the latest event
        per1 = c(per1, recover)
        per2 = c(per2, NA)
        
        cat("Disease extinct at time", t.next, ", simulation stops.\n")
        break
      }
      
      
      # # check whether adjmat is still symmetric--things have gone wrong somewhere
      # if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
      #   stop(paste("The adjmat isn't symmetric any more, where nonzero =", 
      #              nnzero(adjmat), "and that's after z =", z,
      #              "from.ind =", from.ind, "to.ind = ", to.ind, 
      #              "to.cands =", paste(to.cands,collapse = ","),
      #              "from =", from, "to =", to,"\n"))
      # }
      
      
      # 2.2.* make updates about S-I, S-S, I-I connections if z==1 or 2
      # while handling edge cases of N.S <= 1 and/or N.I == 1
      if(z %in% c(1,2)){
        if(N.S == 0){
          adjmat.si = NULL
          N.SI = 0
          if(quarantine && model == "SIR" && N.R == 0){
            adjmat.ri = NULL
            N.c[2] = 0
          }
          # 03/27/2020: debug
          if(quarantine && model == "SIR" && N.R > 0){
            adjmat.ri = adjmat[recovered,infected]
            N.c[2] = nnzero(adjmat.ri)
          }
        }else{
          adjmat.si = adjmat[susceptible,infected]
          N.SI = nnzero(adjmat.si)
          if(quarantine && model == "SIR"){
            adjmat.ri = adjmat[union(susceptible,recovered),infected]
            N.c[2] = nnzero(adjmat.ri)
          }
        }
        if(quarantine){
          if(model == "SIS"){
            if(N.S <= 1){
              adjmat.ss = 0
              N.c[1] = 0
              N.d[1] = 0
            }else{
              adjmat.ss = adjmat[susceptible,susceptible]
              N.c[1] = nnzero(adjmat.ss)/2
              N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
            }
            N.c[2] = N.SI
            N.d[2] = N.S * N.I - N.c[2]
          }else{
            #N.SR = N.S + N.R; 
            sus.recov = union(susceptible,recovered)
            # 03/27/2020: debug
            N.SR = length(sus.recov)
            
            if(N.SR <= 1){
              adjmat.rs = 0
              N.c[1] = 0
              N.d[1] = 0
            }else{
              adjmat.rs = adjmat[sus.recov,sus.recov]
              N.c[1] = nnzero(adjmat.rs)/2
              N.d[1] = N.SR * (N.SR - 1)/2 - N.c[1]
            }
            N.d[2] = N.SR * N.I - N.c[2]
          }
          if(N.I == 1){
            adjmat.ii = 0
            N.c[3] = 0
            N.d[3] = 0
          }else{
            adjmat.ii = adjmat[infected,infected]
            N.c[3] = nnzero(adjmat.ii)/2
            N.d[3] = N.I * (N.I - 1)/2 - N.c[3]
          }
        }
      }
      
      # 3. update inference statistics and do inference occasionally
      # 3.1 update stats
      del.t = t.next - t.cur
      z = as.integer(z)
      event.counts[z] = event.counts[z] + 1
      tot.event = tot.event + 1
      
      big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
      
      # 3.2 inference if ...
      if(tot.event %% infer.interval == 0){
        if(Bayes){
          # Bayesian inference: Gamma sampling
          a.post = a.pr + event.counts
          b.post = b.pr + big.sums
          if(infer.verbose){
            cat("With", tot.event,"events, posterior means and variances:\n",
                a.post/b.post, "\n",
                a.post/(b.post^2), "\n")
          }
          # record posterior stats
          if(return.infer){
            m.post = a.post/b.post
            lb.post = qgamma(.025, a.post, b.post)
            ub.post = qgamma(.975, a.post, b.post)
            this.chunk = cbind(m.post, lb.post, ub.post, event.counts, tot.event)
            Bayes.tab = rbind(Bayes.tab, this.chunk)
          }
          # plot
          if(infer.plot){
            par(mfrow=c(2,2))
            for(ix in 1:length(params)){
              v = vars[ix]; p = params[ix]
              sams = rgamma(samples, shape=a.post[ix], rate = b.post[ix])
              xl = range(sams, p)
              hist(sams, main=paste("Posterior samples for",vars[ix]),
                   xlab = vars[ix], xlim=xl, col="lightblue")
              abline(v=p,col="red")
              if(demo.sleep & ix%%4 == 0){
                Sys.sleep(0.2)
              }
            }
          }
        }
        if (MLE){
          # MLE
          MLEs = event.counts/big.sums
          if(infer.verbose){
            cat("With",tot.event, "events, MLEs:","\n", MLEs, "\n")
          }
          MLE.tab = rbind(MLE.tab, c(tot.event, MLEs, event.counts))
        }
      }
      
      # updated 02/04/2020
      # 3.3 examine network degree distribution if
      if(net.deg.plot & (tot.event %% net.interval == 0)){
        degs = rowSums(adjmat)
        avg.deg = mean(degs)
        poi.degs = rpois(N, avg.deg)
        
        deg.dat = data.frame(Degrees = degs, DegPoi = poi.degs)
        print(
          ggplot(deg.dat) + geom_density(aes(x = Degrees, color="Actual degree")) +
            geom_density(aes(x = DegPoi, color="Poisson")) +
            labs(color="Distribution", 
                 caption = paste("At t =", round(t.cur, digits = 2))) + 
            theme_bw(base_size = 14) %+replace%
            theme(legend.position = c(0.8,0.8))
        )
      }
      
      
      
      # 4. report the event and summarize and set t.cur
      
      t.cur = t.next
      
      # record the event
      time = c(time, t.cur)
      event.type = c(event.type, z)
      preval = c(preval, N.I)
      dens = c(dens, nnzero(adjmat))
      # record labels of individual(s) involved
      if (z >= 3){
        per1 = c(per1, from)
        per2 = c(per2, to)
      }else if (z==2){
        per1 = c(per1, recover)
        per2 = c(per2, NA)
      }else{
        per1 = c(per1, new.infec)
        per2 = c(per2, NA)
      }
      
      
      if(verbose){
        cat(paste0("At time ", t.cur, ", "),msg,"\n")
      }
      
    }
    
    # return results/summary
    cat("Simulation finished. At the end: \n", 
        "Disease prevalence =", N.I/N, "\n",
        "Network density = ", nnzero(adjmat)/(N*(N-1)),"\n")
    event.log = data.frame(time = time, event = event.type,
                           per1 = per1, per2 = per2,
                           preval = preval/N, dens = dens/(N*(N-1)))
    
    # inference at the end
    if(Bayes){
      # Bayesian inference: Gamma sampling
      a.post = a.pr + event.counts
      b.post = b.pr + big.sums
      cat("Posterior means and variances with", tot.event, "events in total:\n",
          a.post/b.post, "\n",
          a.post/(b.post^2), "\n")
      cat("Event type counts:\n", event.counts,"\n")
      # record posterior stats
      if(return.infer){
        m.post = a.post/b.post
        lb.post = qgamma(.025, a.post, b.post)
        ub.post = qgamma(.975, a.post, b.post)
        this.chunk = cbind(m.post, lb.post, ub.post, event.counts, tot.event)
        Bayes.tab = rbind(Bayes.tab, this.chunk)
      }
      # plot
      if(infer.plot){
        par(mfrow=c(2,2))
        for(ix in 1:length(params)){
          v = vars[ix]; p = params[ix]
          sams = rgamma(samples, shape=a.post[ix], rate = b.post[ix])
          xl = range(sams, p)
          hist(sams, main=paste("Posterior samples for",vars[ix]),
               xlab = vars[ix], xlim=xl, col="lightblue")
          abline(v=p,col="red")
          if(demo.sleep & ix%%4 == 0){
            Sys.sleep(0.2)
          }
        }
      }
    }
    if(MLE){
      # MLE
      MLEs = event.counts/big.sums
      cat("MLEs with", tot.event, "events in total:","\n", MLEs, "\n")
      cat("Event type counts:\n", event.counts,"\n")
      
      MLE.tab = rbind(MLE.tab, c(tot.event, MLEs, event.counts))
      MLE.tab = as.data.frame(MLE.tab)
      if(quarantine){
        names(MLE.tab) = c("event.num","beta","gamma",
                           "alpha.SS","alpha.SI","alpha.II",
                           "omega.SS","omega.SI","omega.II",
                           "infection","recovery",
                           "SS.connect","SI.connect","II.connect",
                           "SS.discon","SI.discon","II.discon")
      }else{
        names(MLE.tab) = c("event.num","beta","gamma","alpha","omega",
                           "infection","recovery","connect","discon")
      }
    }
    
    if(plot){
      subt = paste(model," beta =",bet," gamma =", gam, 
                   " alpha.r =", alpha.r, " alpha.d =", alpha.d)
      plot(preval ~ time, data = event.log, 
           xlab="Time", ylab="Disease prevalence", 
           main = "Disease Prevalence vs. Time", sub = subt,
           ylim = c(0,1), type="l",lwd=2)
      plot(dens ~ time, data = event.log,
           xlab="Time", ylab="Edge density", 
           main = "Network Edge Density vs. Time", sub = subt,
           type="l",lwd=2)
    }
    
    if(infer.plot){
      #params = c(bet, gam, alpha.r, alpha.d)
      if(Bayes){
        par(mfrow=c(1,1))
      }else{
        # a function to calculate the frequentist 95% CI for MLEs
        MLE.CI <- function(x){
          l.hat = x[1]; n.obs = x[2]
          if (is.na(l.hat) || l.hat == 0){
            lb = 0; ub = 0
          }else if (n.obs >= 20){
            lb = l.hat * (1 + qnorm(.025)/sqrt(n.obs))
            ub = l.hat * (1 + qnorm(.975)/sqrt(n.obs))
          }else{
            lb = l.hat * qchisq(.025, 2*n.obs) / (2*n.obs)
            ub = l.hat * qchisq(.975, 2*n.obs) / (2*n.obs)
          }
          return(c(lb,ub))
        }
        
        lcol = length(names(MLE.tab)[-1])
        vars = names(MLE.tab)[-1][1:(lcol/2)]
        counts = names(MLE.tab)[-1][(lcol/2 + 1):lcol]
        e.num = MLE.tab$event.num
        if(quarantine){
          par(mfrow=c(2,2))
        }else{
          par(mfrow=c(2,2))
        }
        for(ix in 1:length(vars)){
          v = vars[ix]
          if(all(is.na(MLE.tab[[v]]))){
            plot.new()
          }else{
            tab.v = MLE.tab[,c(v,counts[ix])]
            bounds.v = t(apply(tab.v, 1, MLE.CI))
            bounds.v = as.data.frame(bounds.v)
            names(bounds.v) = c("lb","ub")
            #s.m = mean(MLE.tab[[v]], na.rm=T)
            #s.sd = sd(MLE.tab[[v]], na.rm=T)
            yl = range(min(bounds.v$lb), max(bounds.v$ub), params[ix])
            plot(MLE.tab[[v]]~e.num, xlab="# events", ylab=v,
                 main = paste("MLE for", v, "v.s. # events"), 
                 type="l", ylim=yl)
            lines(bounds.v$lb~e.num, lty=2, col="gray")
            lines(bounds.v$ub~e.num, lty=2, col="gray")
          }
          abline(h=params[ix], col="red")
        }
        par(mfrow=c(1,1))
      }
    }
    
    # add: print out big.sums if Bayesian
    if(Bayes){
      cat("big.sums: \n", big.sums,"\n")
    }
    
    # 02/04/2020: degree distribution at the end
    if(net.deg.plot){
      degs = rowSums(adjmat)
      avg.deg = mean(degs)
      poi.degs = rpois(N, avg.deg)
      
      deg.dat = data.frame(Degrees = degs, DegPoi = poi.degs)
      print(
        ggplot(deg.dat) + geom_density(aes(x = Degrees, color="Actual degree")) +
          geom_density(aes(x = DegPoi, color="Poisson")) +
          labs(color="Distribution", 
               caption = paste("At the end, t =", round(t.cur, digits = 2))) + 
          theme_bw(base_size = 14) %+replace%
          theme(legend.position = c(0.8,0.8))
      )
    }
    
    
    # return event data
    # MODIFICATION: return true parametes too
    if(return.infer){
      if(Bayes){
        Bayes.tab = as.data.frame(Bayes.tab)
        Bayes.tab$var = vars
        infer.tab = Bayes.tab
      }else{
        infer.tab = MLE.tab
      }
      return(list(G0 = init.adj, I0 = init.infected,
                  events = event.log, infer=infer.tab,
                  infer2 = MLE.tab,
                  truth = params))
    }else{
      return(list(G0 = init.adj, I0 = init.infected,
                  events = event.log,truth = params))
    }
    
  }

# Added 04/23/2023 by Houjie. 
stochastic_coevolve_infer5 <- 
  function(N=75, 
           tmax=100, 
           priors,
           bet=0.1, gam=0.15, 
           init.infec = 1,
           alpha.r=c(0.015, 0.0095, 0.0125, 0.004, 0.007, 0.0100), 
           alpha.d=c(0.10, 0.15, 0.07, 0.20, 0.12, 0.04),
           init.net = NULL, selfish = c(1: N)[as.logical(rbinom(N, 1, 0.5))], lag_period = 7,
           init.p = 0.1){
    # set.seed(10001)
    # N=75; tmax=100; bet=0.1; gam=0.15; init.infec = 1; init.net = NULL; init.p = 0.1
    # selfish = c(1: N)[as.logical(rbinom(N, 1, 0.5))]; init.net = NULL; lag_period = NULL
    # 
    # priors = data.frame(count = rep(1,4), avg = c(0.075, 0.1, 0.005, 0.05))
    # alpha.r = c(0.015, 0.0095, 0.0125, 0.004, 0.007, 0.0100) # corresponds to c(0.005, 0.002, 0.0075)
    # alpha.d = c(0.10, 0.15, 0.07, 0.20, 0.12, 0.04) # corresponds to c(0.02, 0.1, 0.05)
    
    all_status = expand.grid(lapply(1: 4, function(i){c(0, 1)})) %>% as.matrix()
    colnames(all_status) = c("Ai", "xi", "Aj", "xj")
    Ai = all_status[, 1]
    xi = all_status[, 2]
    Aj = all_status[, 3]
    xj = all_status[, 4]
    lincomb_status = cbind(Ai*xi+Aj*xj, Ai*(1-xi)+Aj*(1-xj), (1-Ai)+(1-Aj))
    
    # network types: 
    # 1: both healthy regardless of selfishness
    # 2: one sick altruistic and one healthy selfish
    # 3: one sick selfish and one healthy altruistic
    # 4: both sick and altruistic
    # 5: both sick and one of them selfish
    # 6: both sick and both selfish
    unique_status = unique(lincomb_status)
    status_cases = list()
    # Indices of each possible combination of health status and selfishness
    rate_type = sapply(1: nrow(all_status), function(i){
      which(colSums(abs(t(unique_status) - lincomb_status[i, ])) == 0)
    })
    
    
    # summarization stats to be recorded
    time = NULL
    event.type = NULL
    preval = NULL
    dens = NULL
    # add involved individual labels 
    per1 = NULL
    per2 = NULL
    
    # more stats to record for inference
    tot.event = 0
    
    
    
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=6))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=6))
    b.pr = a.pr/avg.pr
    
    Bayes.tab <- MLE.tab <- NULL
    
    vars = c("beta","gamma",
             paste0("alpha", 1: nrow(unique_status)),
             paste0("omega", 1: nrow(unique_status)))
    params = c(bet, gam, alpha.r, alpha.d)
    
    
    # initialize ER(p) graph if the initial netwok is not provided
    if(is.null(init.net)){
      init.net = sample_gnp(n = N, p = init.p, directed = F)
    }
    # convert to adjacency matrix
    adjmat = get.adjacency(init.net, type=c("both"), names=F)
    # adjmat = as(adjmat,"dgTMatrix")
    adjmat = as.matrix(adjmat)
    # save the initial network structure (adj.mat)
    init.adj = adjmat
    
    ## check whether adjmat is symmetric and nnzero() returns an even number
    # if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
    #   stop(paste("The initial adjmat isn't symmetric, where nonzero =", 
    #              nnzero(adjmat),"\n"))
    # }
    
    # initialize epidemic
    infected = sample(N, init.infec)
    epid = rep(0,N); epid[infected] = 1
    # save initial cases
    init.infected = infected
    
    # initialize data objects to keep track of
    susceptible = which(epid==0); infected = which(epid==1); recovered = NULL;
    altruistic = c(1: N)[!c(1: N) %in% selfish]
    
    healthy = union(susceptible,recovered)
    # number of infected, susceptible, and recovered
    N.I = init.infec; N.S = N - N.I; N.R = sum(epid==-1)
    N.H = N.S + N.R # number of healthy
    N.IA = length(intersect(infected, altruistic))
    N.HS = length(intersect(healthy, selfish))
    N.IS = length(intersect(infected, selfish))
    N.HA = length(intersect(healthy, altruistic))
    
    N.c = numeric(nrow(unique_status)); N.d = numeric(nrow(unique_status))
    N.c[1] = nnzero(adjmat[healthy, healthy])/2
    N.c[2] = nnzero(adjmat[intersect(infected, altruistic), intersect(healthy, selfish)])
    N.c[3] = nnzero(adjmat[intersect(infected, selfish), intersect(healthy, altruistic)])
    N.c[4] = nnzero(adjmat[intersect(infected, altruistic), intersect(infected, altruistic)])/2
    N.c[5] = nnzero(adjmat[intersect(infected, altruistic), intersect(infected, selfish)])
    N.c[6] = nnzero(adjmat[intersect(infected, selfish), intersect(infected, selfish)])/2
    
    NMax = c(N.H*(N.H-1)/2, N.IA*N.HS, N.IS*N.HA, N.IA*(N.IA-1)/2, N.IA*N.IS, N.IS*(N.IS-1)/2)
    N.d = NMax - N.c
    N.SI = nnzero(adjmat[susceptible,infected])
    
    health_sick = list(healthy, infected)
    altru_self = list(altruistic, selfish)
    # iterations
    big.sums <- rep(0, length(params))
    SIR_counts = c()
    t.cur = 0
    while (t.cur < tmax) {
      # NOTES #
      # a. check to make sure N.S also > 0
      # b. handle edge cases (N.S = 1 and/or N.I = 1)
      
      # # check whether all N's are non-negative
      # Ns = c(N.SI, N.I, N.d, N.c)
      # if(any(Ns < 0) | all(Ns == 0)){
      #   msg = paste("Weird stuff happened! All the Ns =", paste(Ns,collapse = ", "), "\n")
      #   stop(msg)
      # }
      
      # 1. calculate sum of risks and sample the next event time
      Lambda = bet * N.SI + N.I * gam + sum(alpha.r * N.d) + sum(alpha.d * N.c)
      t.next = t.cur + rexp(1, rate = Lambda)
      
      
      # 2. determine the event type
      ratios = c(bet * N.SI, N.I * gam, alpha.r * N.d, alpha.d * N.c) / Lambda
      (z = sample(2 + 2*6, 1, prob = ratios))
      
      # 2.1 deal with network change
      
      if (z %in% 3: (3+length(N.d)-1)){
        connect_type = z - 2
        involved_status = matrix(all_status[rate_type == connect_type, ], ncol = 4, byrow = T)
        all_pairs = c()
        for (ii in 1: nrow(involved_status)) {
          # ii = 1
          temp_status = involved_status[ii, ] + 1
          all_pairs = rbind(all_pairs, 
                            expand.grid(intersect(health_sick[[temp_status[1]]], altru_self[[temp_status[2]]]),
                                        intersect(health_sick[[temp_status[3]]], altru_self[[temp_status[4]]])) %>% as.matrix() %>% rowSort()
          )
        }
        all_pairs = all_pairs[all_pairs[, 1] != all_pairs[, 2], ]
        all_pairs_tmp = matrix(all_pairs[which(adjmat[unique(all_pairs)] == 0), ], ncol = 2)
        selected_pair = as.integer(all_pairs[sample(1: nrow(all_pairs_tmp), 1), ])
        adjmat[selected_pair[1], selected_pair[2]] = 1
        adjmat[selected_pair[2], selected_pair[1]] = 1
      }
      
      if (z %in% 9: (9+length(N.c)-1)){
        connect_type = z - 2 - length(N.d)
        involved_status = matrix(all_status[rate_type == connect_type, ], ncol = 4, byrow = T)
        all_pairs = c()
        for (ii in 1: nrow(involved_status)) {
          # ii = 1
          temp_status = involved_status[ii, ] + 1
          all_pairs = rbind(all_pairs, 
                            expand.grid(intersect(health_sick[[temp_status[1]]], altru_self[[temp_status[2]]]),
                                        intersect(health_sick[[temp_status[3]]], altru_self[[temp_status[4]]])) %>% as.matrix() %>% rowSort()
          )
        }
        all_pairs = all_pairs[all_pairs[, 1] != all_pairs[, 2], ]
        all_pairs_tmp = matrix(all_pairs[which(adjmat[unique(all_pairs)] == 1), ], ncol = 2)
        selected_pair = as.integer(all_pairs[sample(1: nrow(all_pairs_tmp), 1), ])
        adjmat[selected_pair[1], selected_pair[2]] = 0
        adjmat[selected_pair[2], selected_pair[1]] = 0
      }
      
      
      #2.2 deal with epidemic progress
      if(z==2){
        # recover an infected person
        rec.ind = sample(N.I,1)
        recover = infected[rec.ind]
        infected = infected[-rec.ind]
        N.I = N.I - 1 
        
        epid[recover] = -1
        
        recovered = c(recovered, recover)
        N.R = N.R + 1
        health_sick[[1]] = c(health_sick[[1]], recover)
        health_sick[[2]] = health_sick[[2]][health_sick[[2]] != recover]
      }
      
      
      if(z==1){
        # select a susceptible person to infect
        if(N.S == 1){
          # only one person to infect
          new.infec = susceptible
          susceptible = NULL
        }else{
          # multiple susceptible inviduals to choose from
          # if(N.I == 1){
          #   # only one infected person
          #   num.nei.i = sum(adjmat[infected, ])
          #   expand.grid(susceptible, infected)
          # }else{
          #   num.nei.i = rowSums(adjmat[infected, ])
          # }
          SI_pair = expand.grid(susceptible, infected)
          new.infec = SI_pair[sample(1: nrow(SI_pair), 1), 1]
          # infec.ind = sample(N.S,1,prob = num.nei.i)
          # new.infec = susceptible[infec.ind]
          susceptible = susceptible[susceptible != new.infec]
          
        }
        epid[new.infec] = 1
        N.S = N.S - 1
        infected = c(infected, new.infec); 
        # add the newly infected person into the list of sick persons
        health_sick[[2]] = c(health_sick[[2]], new.infec)
        # remove the newly infected person into the list of healthy persons
        health_sick[[1]] = health_sick[[1]][health_sick[[1]] != new.infec]
        N.I = N.I + 1
        # N.SI = N.SI + sum(which(as.logical(adjmat[new.infec, ])) %in% susceptible)
      }
      if (any(diag(adjmat) != 0)){
        print(paste0("adjmat diagonal not all zero"))
        break
      }
      healthy = health_sick[[1]]
      
      N.H = N.S + N.R # number of healthy
      N.IA = length(intersect(infected, altruistic))
      N.HS = length(intersect(healthy, selfish))
      N.IS = length(intersect(infected, selfish))
      N.HA = length(intersect(healthy, altruistic))
      
      N.c[1] = nnzero(adjmat[healthy, healthy])/2
      N.c[2] = nnzero(adjmat[intersect(infected, altruistic), intersect(healthy, selfish)])
      N.c[3] = nnzero(adjmat[intersect(infected, selfish), intersect(healthy, altruistic)])
      N.c[4] = nnzero(adjmat[intersect(infected, altruistic), intersect(infected, altruistic)])/2
      N.c[5] = nnzero(adjmat[intersect(infected, altruistic), intersect(infected, selfish)])
      N.c[6] = nnzero(adjmat[intersect(infected, selfish), intersect(infected, selfish)])/2
      
      NMax = c(N.H*(N.H-1)/2, N.IA*N.HS, N.IS*N.HA, N.IA*(N.IA-1)/2, N.IA*N.IS, N.IS*(N.IS-1)/2)
      N.d = NMax - N.c
      N.SI = nnzero(adjmat[susceptible, infected])
      SIR_counts = rbind(SIR_counts, 
                         c(N.SI, N.S, N.I, N.R))
      # }
      
      # A sanity check: if N.I==0, stop the simulation
      if(N.I==0){
        # if(verbose){
        #   cat(paste0("At time ", t.next, ", "),msg,"\n")
        # }
        # still update inference statistics
        del.t = t.next - t.cur
        z = as.integer(z)
        tot.event = tot.event + 1
        
        big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
        
        # still have to record the event
        time = c(time, t.next)
        event.type = c(event.type, z)
        # preval = c(preval, N.I)
        # dens = c(dens, nnzero(adjmat))
        # here: a recovery must be the latest event
        per1 = c(per1, recover)
        per2 = c(per2, NA)
        
        # cat("Disease extinct at time", t.next, ", simulation stops.\n")
        break
      }
      
      # 3. update inference statistics and do inference occasionally
      # 3.1 update stats
      del.t = t.next - t.cur
      z = as.integer(z)
      tot.event = tot.event + 1
      
      big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
      
      
      
      # 4. report the event and summarize and set t.cur
      if (!is.null(lag_period)){
        if (as.integer(t.next/lag_period) == as.integer(t.cur/lag_period) + 1){
          alpha.d = alpha.d*c(1+N.I/nrow(adjmat), 1+0.5*N.I/nrow(adjmat), 1+0.5*N.I/nrow(adjmat), 1, 1, 1)
          alpha.r = alpha.r*c(1-N.I/nrow(adjmat), 1-0.5*N.I/nrow(adjmat), 1-0.5*N.I/nrow(adjmat), 1, 1, 1)
        }
      }
      t.cur = t.next
      
      # record the event
      time = c(time, t.cur)
      event.type = c(event.type, z)
      # preval = c(preval, N.I)
      # dens = c(dens, nnzero(adjmat))
      # record labels of individual(s) involved
      if (z >= 3){
        per1 = c(per1, selected_pair[1])
        per2 = c(per2, selected_pair[2])
      }else if (z==2){
        per1 = c(per1, recover)
        per2 = c(per2, NA)
      }else{
        per1 = c(per1, new.infec)
        per2 = c(per2, NA)
      }
      
      
    }
    event.log = data.frame(time = time, event = event.type,
                           per1 = per1, per2 = per2)
    if (!is.null(lag_period)){
      # return results/summary
      event.counts_tmp = table(event.type)
      event.counts = vector("integer", length = length(big.sums))
      event.counts[names(event.counts_tmp) %>% as.integer()] = event.counts_tmp
      
      # inference at the end
      
      # Bayesian inference: Gamma sampling
      a.post = a.pr + event.counts
      b.post = b.pr + big.sums
      m.post = a.post/b.post
      lb.post = qgamma(.025, a.post, b.post)
      ub.post = qgamma(.975, a.post, b.post)
      this.chunk = cbind(m.post, lb.post, ub.post, event.counts, tot.event)
      Bayes.tab = rbind(Bayes.tab, this.chunk)
      
      # MLE
      MLEs = event.counts/big.sums
      MLE.tab = rbind.data.frame(event.counts, MLEs)
      names(MLE.tab) = vars
      
      return(list(G0 = init.adj, I0 = init.infected,
                  events = event.log, infer=Bayes.tab,
                  infer2 = MLE.tab, SIR_record = SIR_counts,
                  truth = params))
    }else{
      return(list(G0 = init.adj, I0 = init.infected,
                  events = event.log,SIR_record = SIR_counts,
                  truth = params))
    }
    
    # }else{
    #   return(list(G0 = init.adj, I0 = init.infected,
    #               events = event.log,truth = params))
    # }
    
  }

stochastic_coevolve_infer6 <- 
  function(N=75, 
           tmax=100, event.max = 2e5,
           priors,
           bet=0.1, gam=0.15, 
           init.infec = 1,
           alpha.r=c(0.015, 0.0095, 0.0125, 0.004, 0.007, 0.0100), 
           alpha.d=c(0.10, 0.15, 0.07, 0.20, 0.12, 0.04),
           init.net = NULL, selfish = c(1: N)[as.logical(rbinom(N, 1, 0.5))], lag_period = 7,
           init.p = 0.1){
    # set.seed(10001)
    # N=500; tmax=100; bet=0.1; gam=0.15; event.max = 10; init.infec = 1; init.net = NULL; init.p = 0.1
    # selfish = c(1: N)[as.logical(rbinom(N, 1, 0.5))]; init.net = NULL; lag_period = NULL
    # 
    # priors = data.frame(count = rep(1,4), avg = c(0.075, 0.1, 0.005, 0.05))
    # alpha.r = c(0.015, 0.0095, 0.0125, 0.004, 0.007, 0.0100) # corresponds to c(0.005, 0.002, 0.0075)
    # alpha.d = c(0.10, 0.25, 0.07, 0.20, 0.12, 0.04) # corresponds to c(0.02, 0.1, 0.05)
    
    all_status = expand.grid(lapply(1: 4, function(i){c(0, 1)})) %>% as.matrix()
    colnames(all_status) = c("Ai", "xi", "Aj", "xj")
    Ai = all_status[, 1]
    xi = all_status[, 2]
    Aj = all_status[, 3]
    xj = all_status[, 4]
    lincomb_status = cbind(Ai*xi+Aj*xj, Ai*(1-xi)+Aj*(1-xj), (1-Ai)+(1-Aj))
    
    # network types: 
    # 1: both healthy regardless of selfishness
    # 2: one sick altruistic and one healthy selfish
    # 3: one sick selfish and one healthy altruistic
    # 4: both sick and altruistic
    # 5: both sick and one of them selfish
    # 6: both sick and both selfish
    unique_status = unique(lincomb_status)
    status_cases = list()
    # Indices of each possible combination of health status and selfishness
    rate_type = sapply(1: nrow(all_status), function(i){
      which(colSums(abs(t(unique_status) - lincomb_status[i, ])) == 0)
    })
    
    
    # summarization stats to be recorded
    time = NULL
    event.type = NULL
    preval = NULL
    dens = NULL
    # add involved individual labels 
    per1 = NULL
    per2 = NULL
    
    # more stats to record for inference
    tot.event = 0
    
    
    
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=6))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=6))
    b.pr = a.pr/avg.pr
    
    Bayes.tab <- MLE.tab <- NULL
    
    vars = c("beta","gamma",
             paste0("alpha", 1: nrow(unique_status)),
             paste0("omega", 1: nrow(unique_status)))
    params = c(bet, gam, alpha.r, alpha.d)
    
    
    # initialize ER(p) graph if the initial netwok is not provided
    if(is.null(init.net)){
      init.net = sample_gnp(n = N, p = init.p, directed = F)
    }
    # convert to adjacency matrix
    adjmat = get.adjacency(init.net, type=c("both"), names=F)
    # adjmat = as(adjmat,"dgTMatrix")
    adjmat = as.matrix(adjmat)
    # save the initial network structure (adj.mat)
    init.adj = adjmat
    
    ## check whether adjmat is symmetric and nnzero() returns an even number
    # if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
    #   stop(paste("The initial adjmat isn't symmetric, where nonzero =", 
    #              nnzero(adjmat),"\n"))
    # }
    
    # initialize epidemic
    infected = sample(N, init.infec)
    epid = rep(0,N); epid[infected] = 1
    # save initial cases
    init.infected = infected
    
    # initialize data objects to keep track of
    susceptible = which(epid==0); infected = which(epid==1); recovered = NULL;
    altruistic = c(1: N)[!c(1: N) %in% selfish]
    
    healthy = union(susceptible,recovered)
    # number of infected, susceptible, and recovered
    N.I = init.infec; N.S = N - N.I; N.R = sum(epid==-1)
    N.H = N.S + N.R # number of healthy
    N.IA = length(intersect(infected, altruistic))
    N.HS = length(intersect(healthy, selfish))
    N.IS = length(intersect(infected, selfish))
    N.HA = length(intersect(healthy, altruistic))
    
    N.c = numeric(nrow(unique_status)); N.d = numeric(nrow(unique_status))
    N.c[1] = nnzero(adjmat[healthy, healthy])/2
    N.c[2] = nnzero(adjmat[intersect(infected, altruistic), intersect(healthy, selfish)])
    N.c[3] = nnzero(adjmat[intersect(infected, selfish), intersect(healthy, altruistic)])
    N.c[4] = nnzero(adjmat[intersect(infected, altruistic), intersect(infected, altruistic)])/2
    N.c[5] = nnzero(adjmat[intersect(infected, altruistic), intersect(infected, selfish)])
    N.c[6] = nnzero(adjmat[intersect(infected, selfish), intersect(infected, selfish)])/2
    
    NMax = c(N.H*(N.H-1)/2, N.IA*N.HS, N.IS*N.HA, N.IA*(N.IA-1)/2, N.IA*N.IS, N.IS*(N.IS-1)/2)
    N.d = NMax - N.c
    N.SI = nnzero(adjmat[susceptible,infected])
    
    health_sick = list(healthy, infected)
    altru_self = list(altruistic, selfish)
    # iterations
    big.sums <- rep(0, length(params))
    SIR_counts = c()
    t.cur = 0
    
    disconn_factor <- conn_factor <- rep(1, length(N.d))
    while (t.cur < tmax) {
      # profvis({
      if (tot.event %% 1000 == 0) print(paste0("Now running the ", tot.event, "-th event"))
      if (tot.event >= event.max){break}
      # NOTES #
      # a. check to make sure N.S also > 0
      # b. handle edge cases (N.S = 1 and/or N.I = 1)
      
      # # check whether all N's are non-negative
      # Ns = c(N.SI, N.I, N.d, N.c)
      # if(any(Ns < 0) | all(Ns == 0)){
      #   msg = paste("Weird stuff happened! All the Ns =", paste(Ns,collapse = ", "), "\n")
      #   stop(msg)
      # }
      
      # 1. calculate sum of risks and sample the next event time
      Lambda = bet * N.SI + N.I * gam + sum(alpha.r * N.d) + sum(alpha.d * N.c)
      t.next = t.cur + rexp(1, rate = Lambda)
      
      
      # 2. determine the event type
      ratios = c(bet * N.SI, N.I * gam, alpha.r * N.d, alpha.d * N.c) / Lambda
      (z = sample(2 + 2*6, 1, prob = ratios))
      
      # 2.1 deal with network change
      
      if (z %in% 3: (3+length(N.d)-1)){
        connect_type = z - 2
        
        if (connect_type == 1){# both healthy regardless of selfishness
          all_pairs = comboGrid(rep(list(health_sick[[1]]), 2), repetition = F)
        }
        if (connect_type == 4){# both sick and altruistic
          all_pairs = comboGrid(rep(list(intersect(health_sick[[2]], altru_self[[1]])), 2), repetition = F)
        }
        if (connect_type == 6){# both sick and both selfish
          all_pairs = comboGrid(rep(list(intersect(health_sick[[2]], altru_self[[2]])), 2), repetition = F)
        }
        if (connect_type == 2){# one sick altruistic and one healthy selfish
          all_pairs = comboGrid(list(intersect(health_sick[[2]], altru_self[[1]]), 
                                     intersect(health_sick[[1]], altru_self[[2]])), repetition = F)
        }
        if (connect_type == 3){# one sick selfish and one healthy altruistic
          all_pairs = comboGrid(list(intersect(health_sick[[2]], altru_self[[2]]), 
                                     intersect(health_sick[[1]], altru_self[[1]])), repetition = F)
        }
        if (connect_type == 5){# both sick and one of them selfish
          all_pairs = comboGrid(list(intersect(health_sick[[2]], altru_self[[2]]), 
                                     intersect(health_sick[[2]], altru_self[[1]])), repetition = F)
        }
        
        
        # all_pairs = c()
        # for (ii in 1: nrow(involved_status)) {
        #   # ii = 1
        #   temp_status = involved_status[ii, ] + 1
        #   all_pairs = rbind(all_pairs, 
        #                     expand.grid(intersect(health_sick[[temp_status[1]]], altru_self[[temp_status[2]]]),
        #                                 intersect(health_sick[[temp_status[3]]], altru_self[[temp_status[4]]])) %>% as.matrix() %>% rowSort()
        #   )
        # }
        # all_pairs = all_pairs[all_pairs[, 1] != all_pairs[, 2], ]
        selected_pair = all_pairs[sample(which(adjmat[all_pairs] == 0), 1), ]
        # all_pairs_tmp = matrix(all_pairs[which(adjmat[unique(all_pairs)] == 0), ], ncol = 2)
        # selected_pair = as.integer(all_pairs[sample(1: nrow(all_pairs_tmp), 1), ])
        adjmat[selected_pair[1], selected_pair[2]] = 1
        adjmat[selected_pair[2], selected_pair[1]] = 1
      }
      
      if (z %in% 9: (9+length(N.c)-1)){
        connect_type = z - 2 - length(N.d)
        if (connect_type == 1){# both healthy regardless of selfishness
          all_pairs = comboGrid(rep(list(health_sick[[1]]), 2), repetition = F)
        }
        if (connect_type == 4){# both sick and altruistic
          all_pairs = comboGrid(rep(list(intersect(health_sick[[2]], altru_self[[1]])), 2), repetition = F)
        }
        if (connect_type == 6){# both sick and both selfish
          all_pairs = comboGrid(rep(list(intersect(health_sick[[2]], altru_self[[2]])), 2), repetition = F)
        }
        if (connect_type == 2){# one sick altruistic and one healthy selfish
          all_pairs = comboGrid(list(intersect(health_sick[[2]], altru_self[[1]]), 
                                     intersect(health_sick[[1]], altru_self[[2]])), repetition = F)
        }
        if (connect_type == 3){# one sick selfish and one healthy altruistic
          all_pairs = comboGrid(list(intersect(health_sick[[2]], altru_self[[2]]), 
                                     intersect(health_sick[[1]], altru_self[[1]])), repetition = F)
        }
        if (connect_type == 5){# both sick and one of them selfish
          all_pairs = comboGrid(list(intersect(health_sick[[2]], altru_self[[2]]), 
                                     intersect(health_sick[[2]], altru_self[[1]])), repetition = F)
        }
        
        # all_pairs = c()
        # for (ii in 1: nrow(involved_status)) {
        #   # ii = 1
        #   temp_status = involved_status[ii, ] + 1
        #   all_pairs = rbind(all_pairs, 
        #                     expand.grid(intersect(health_sick[[temp_status[1]]], altru_self[[temp_status[2]]]),
        #                                 intersect(health_sick[[temp_status[3]]], altru_self[[temp_status[4]]])) %>% as.matrix() %>% rowSort()
        #   )
        # }
        # all_pairs = all_pairs[all_pairs[, 1] != all_pairs[, 2], ]
        
        # all_pairs = all_pairs[all_pairs[, 1] != all_pairs[, 2], ]
        selected_pair = all_pairs[sample(which(adjmat[all_pairs] == 1), 1), ]
        # all_pairs_tmp = matrix(all_pairs[which(adjmat[unique(all_pairs)] == 1), ], ncol = 2)
        # selected_pair = as.integer(all_pairs[sample(1: nrow(all_pairs_tmp), 1), ])
        adjmat[selected_pair[1], selected_pair[2]] = 0
        adjmat[selected_pair[2], selected_pair[1]] = 0
      }
      
      
      #2.2 deal with epidemic progress
      if(z==2){
        # recover an infected person
        rec.ind = sample(N.I,1)
        recover = infected[rec.ind]
        infected = infected[-rec.ind]
        N.I = N.I - 1 
        
        epid[recover] = -1
        
        recovered = c(recovered, recover)
        N.R = N.R + 1
        health_sick[[1]] = c(health_sick[[1]], recover)
        health_sick[[2]] = health_sick[[2]][health_sick[[2]] != recover]
      }
      
      
      if(z==1){
        # select a susceptible person to infect
        if(N.S == 1){
          # only one person to infect
          new.infec = susceptible
          susceptible = NULL
        }else{
          # multiple susceptible inviduals to choose from
          # if(N.I == 1){
          #   # only one infected person
          #   num.nei.i = sum(adjmat[infected, ])
          #   expand.grid(susceptible, infected)
          # }else{
          #   num.nei.i = rowSums(adjmat[infected, ])
          # }
          SI_pair = expand.grid(susceptible, infected)
          new.infec = SI_pair[sample(1: nrow(SI_pair), 1), 1]
          # infec.ind = sample(N.S,1,prob = num.nei.i)
          # new.infec = susceptible[infec.ind]
          susceptible = susceptible[susceptible != new.infec]
          
        }
        epid[new.infec] = 1
        N.S = N.S - 1
        infected = c(infected, new.infec); 
        # add the newly infected person into the list of sick persons
        health_sick[[2]] = c(health_sick[[2]], new.infec)
        # remove the newly infected person into the list of healthy persons
        health_sick[[1]] = health_sick[[1]][health_sick[[1]] != new.infec]
        N.I = N.I + 1
        # N.SI = N.SI + sum(which(as.logical(adjmat[new.infec, ])) %in% susceptible)
      }
      if (any(diag(adjmat) != 0)){
        print(paste0("adjmat diagonal not all zero"))
        break
      }
      healthy = health_sick[[1]]
      
      N.H = N.S + N.R # number of healthy
      N.IA = length(intersect(infected, altruistic))
      N.HS = length(intersect(healthy, selfish))
      N.IS = length(intersect(infected, selfish))
      N.HA = length(intersect(healthy, altruistic))
      
      N.c[1] = nnzero(adjmat[healthy, healthy])/2
      N.c[2] = nnzero(adjmat[intersect(infected, altruistic), intersect(healthy, selfish)])
      N.c[3] = nnzero(adjmat[intersect(infected, selfish), intersect(healthy, altruistic)])
      N.c[4] = nnzero(adjmat[intersect(infected, altruistic), intersect(infected, altruistic)])/2
      N.c[5] = nnzero(adjmat[intersect(infected, altruistic), intersect(infected, selfish)])
      N.c[6] = nnzero(adjmat[intersect(infected, selfish), intersect(infected, selfish)])/2
      
      NMax = c(N.H*(N.H-1)/2, N.IA*N.HS, N.IS*N.HA, N.IA*(N.IA-1)/2, N.IA*N.IS, N.IS*(N.IS-1)/2)
      N.d = NMax - N.c
      N.SI = nnzero(adjmat[susceptible, infected])
      SIR_counts = rbind(SIR_counts, 
                         c(N.SI, N.S, N.I, N.R))
      # }
      
      # A sanity check: if N.I==0, stop the simulation
      if(N.I==0){
        # if(verbose){
        #   cat(paste0("At time ", t.next, ", "),msg,"\n")
        # }
        # still update inference statistics
        del.t = t.next - t.cur
        z = as.integer(z)
        tot.event = tot.event + 1
        
        big.sums = big.sums + c(N.SI, N.I, N.d*conn_factor, N.c*disconn_factor) * del.t
        
        # still have to record the event
        time = c(time, t.next)
        event.type = c(event.type, z)
        # preval = c(preval, N.I)
        # dens = c(dens, nnzero(adjmat))
        # here: a recovery must be the latest event
        per1 = c(per1, recover)
        per2 = c(per2, NA)
        
        # cat("Disease extinct at time", t.next, ", simulation stops.\n")
        break
      }
      
      # 3. update inference statistics and do inference occasionally
      # 3.1 update stats
      del.t = t.next - t.cur
      z = as.integer(z)
      tot.event = tot.event + 1
      
      big.sums = big.sums + c(N.SI, N.I, N.d*conn_factor, N.c*disconn_factor) * del.t
      
      
      
      # 4. report the event and summarize and set t.cur
      if (!is.null(lag_period)){
        if (as.integer(t.next/lag_period) == as.integer(t.cur/lag_period) + 1){
          disconn_factor = c(1+0.5*N.I/nrow(adjmat), 1+N.I/nrow(adjmat), 1+N.I/nrow(adjmat), 1, 1, 1)
          conn_factor = c(1-0.5*N.I/nrow(adjmat), 1-N.I/nrow(adjmat), 1-N.I/nrow(adjmat), 1, 1, 1)
          alpha.d = alpha.d*disconn_factor
          alpha.r = alpha.r*conn_factor
        }
      }
      t.cur = t.next
      
      # record the event
      time = c(time, t.cur)
      event.type = c(event.type, z)
      # preval = c(preval, N.I)
      # dens = c(dens, nnzero(adjmat))
      # record labels of individual(s) involved
      if (z >= 3){
        per1 = c(per1, selected_pair[1])
        per2 = c(per2, selected_pair[2])
      }else if (z==2){
        per1 = c(per1, recover)
        per2 = c(per2, NA)
      }else{
        per1 = c(per1, new.infec)
        per2 = c(per2, NA)
      }
      
      # })
    }
    
    event.log = data.frame(time = time, event = event.type,
                           per1 = per1, per2 = per2)
    
    # return results/summary
    event.counts_tmp = table(event.type)
    event.counts = vector("integer", length = length(big.sums))
    event.counts[names(event.counts_tmp) %>% as.integer()] = event.counts_tmp
    
    # inference at the end
    
    # Bayesian inference: Gamma sampling
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    m.post = a.post/b.post
    lb.post = qgamma(.025, a.post, b.post)
    ub.post = qgamma(.975, a.post, b.post)
    this.chunk = cbind(m.post, lb.post, ub.post, event.counts, tot.event)
    Bayes.tab = rbind(Bayes.tab, this.chunk)
    
    # MLE
    MLEs = event.counts/big.sums
    MLE.tab = rbind.data.frame(event.counts, MLEs)
    names(MLE.tab) = vars
    
    return(list(G0 = init.adj, I0 = init.infected,
                events = event.log, infer=Bayes.tab,
                infer2 = MLE.tab, SIR_record = SIR_counts,
                truth = params))
    
    # }else{
    #   return(list(G0 = init.adj, I0 = init.infected,
    #               events = event.log,truth = params))
    # }
    
  }

# 10/13/2024 Arguement Nevent_max added, which will also stop the program if Nevent_max or tmax whichever comes earlier 
stochastic_coevolve_infer7 = function(N=50, tmax=100, Nevent_max = 50000,
                                      bet=0.05, gam=0.2, init.infec = 1, model = "SIS",
                                      alpha.r=0.005, alpha.d=0.05, quarantine=F,
                                      init.net = NULL, init.p = 0.1,
                                      verbose = F, plot = F,
                                      infer.interval = 100, 
                                      Bayes = T, MLE = T, priors = NULL, samples = 500,
                                      infer.verbose = F, infer.plot = F, return.infer = F,
                                      net.deg.plot = F, net.interval = 100,
                                      demo.sleep = F){
  ## alpha.* = c(SSrate, SIrate, IIrate) ##
  #set.seed(seed)
  
  # check edge rates validity
  if(quarantine & (length(alpha.r) < 3 | length(alpha.d) < 3)){
    stop("Under quarantine, need to specify 3 values for edge connection/disconnection rates respectively!")
  }
  # summarization stats to be recorded
  time = NULL
  event.type = NULL
  preval = NULL
  dens = NULL
  # add involved individual labels 
  per1 = NULL
  per2 = NULL
  
  # more stats to record for inference
  if(quarantine){
    event.counts = numeric(8)
    big.sums = numeric(8)
  }else{
    event.counts = numeric(4)
    big.sums = numeric(4)
  }
  tot.event = 0
  
  if(Bayes){
    if(quarantine & nrow(priors)==4){
      a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
      avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
      b.pr = a.pr/avg.pr
    }else{
      a.pr = priors$count
      b.pr = a.pr/priors$avg
    }
    if(quarantine){
      vars = c("beta","gamma",
               "alpha.SS","alpha.SI","alpha.II",
               "omega.SS","omega.SI","omega.II")
    }else{
      vars = c("beta","gamma","alpha","omega")
    }
    Bayes.tab = NULL
  }
  if (MLE){
    MLE.tab = NULL
  }
  params = c(bet, gam, alpha.r, alpha.d)
  
  
  # initialize ER(p) graph if the initial netwok is not provided
  if(is.null(init.net)){
    init.net = sample_gnp(n = N, p = init.p, directed = F)
  }
  # convert to adjacency matrix
  adjmat = get.adjacency(init.net, type=c("both"), names=F)
  adjmat = as(adjmat,"dgTMatrix")
  # save the initial network structure (adj.mat)
  init.adj = adjmat
  
  
  ## check whether adjmat is symmetric and nnzero() returns an even number
  if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
    stop(paste("The initial adjmat isn't symmetric, where nonzero =", 
               nnzero(adjmat),"\n"))
  }
  
  # initialize epidemic
  infected = sample(N, init.infec)
  epid = rep(0,N); epid[infected] = 1
  cat("Initially infected: ",infected,"\n")
  # save initial cases
  init.infected = infected
  
  # initialize data objects to keep track of
  susceptible = which(epid==0)
  N.I = init.infec; N.S = N - N.I
  t.cur = 0
  if(quarantine){
    N.c = numeric(3); N.d = numeric(3)
    adjmat.ii = adjmat[infected,infected]
    N.c[3] = nnzero(adjmat.ii)/2
    N.d[3] = N.I * (N.I - 1)/2 - N.c[3]
    # both need this:
    adjmat.si = adjmat[susceptible,infected]
    if(model=="SIR"){
      recovered = NULL; N.R = 0
      adjmat.rs = adjmat[union(susceptible,recovered), union(susceptible,recovered)]
      N.c[1] = nnzero(adjmat.rs)/2
      N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
      adjmat.ri = adjmat.si
      N.c[2] = nnzero(adjmat.ri)
      N.d[2] = N.S * N.I - N.c[2]
    }else{
      adjmat.ss = adjmat[susceptible,susceptible]
      N.c[1] = nnzero(adjmat.ss)/2
      N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
      
      N.c[2] = nnzero(adjmat.si)
      N.d[2] = N.S * N.I - N.c[2]
    }
    # at the start: N.R = 0, so this works
    N.SI = N.c[2]
  }else{
    N.c = nnzero(adjmat)/2
    N.d = N * (N-1)/2 - N.c
    
    adjmat.si = adjmat[susceptible,infected]
    N.SI = nnzero(adjmat.si)
  }
  
  
  # iterations
  counter = 10
  while ((t.cur < tmax) & (tot.event <= Nevent_max)) {
    
    # NOTES #
    # a. check to make sure N.S also > 0
    # b. handle edge cases (N.S = 1 and/or N.I = 1)
    
    # # check whether all N's are non-negative
    # Ns = c(N.SI, N.I, N.d, N.c)
    # if(any(Ns < 0) | all(Ns == 0)){
    #   msg = paste("Weird stuff happened! All the Ns =", paste(Ns,collapse = ", "), "\n")
    #   stop(msg)
    # }
    
    # 1. calculate sum of risks and sample the next event time
    Lambda = bet * N.SI + N.I * gam + sum(alpha.r * N.d) + sum(alpha.d * N.c)
    t.next = t.cur + rexp(1, rate = Lambda)
    
    if ((t.cur < counter) & (t.next > counter)) {
      print(paste0(counter, " Days have passed"))
      counter = counter + 10
    }
    
    # 2. determine the event type
    ratios = c(bet * N.SI, N.I * gam, alpha.r * N.d, alpha.d * N.c)
    # 2.1 deal with network change
    if(quarantine){
      
      z = sample(8, 1, prob = ratios)
      
      if(z==3){
        # reconnect a random S-S/S-R/R-R pair
        if(model == "SIS"){
          num.dis = (N.S - 1) - rowSums(adjmat.ss)
          from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
          to.cands = which(adjmat.ss[from.ind,]==0)
        }else{
          num.dis = (N.S + N.R - 1) - rowSums(adjmat.rs)
          from.ind = sample(N.S + N.R, 1, prob = num.dis); from = union(susceptible, recovered)[from.ind]
          to.cands = which(adjmat.rs[from.ind,]==0)
        }
        to.cands = to.cands[to.cands!=from.ind]
        if(length(to.cands) == 1){
          to.ind = to.cands
        }else{
          to.ind = sample(to.cands,1)
        }
        if(model == "SIS"){
          to = susceptible[to.ind]
          adjmat.ss[as.integer(from.ind),as.integer(to.ind)] = 1
          adjmat.ss[as.integer(to.ind),as.integer(from.ind)] = 1
          msg = paste0("reconnect an S-S pair, (",from,",",to,")")
        }else{
          to = union(susceptible, recovered)[to.ind]
          adjmat.rs[as.integer(from.ind),as.integer(to.ind)] = 1
          adjmat.rs[as.integer(to.ind),as.integer(from.ind)] = 1
          msg = paste0("reconnect an S-S/S-R/R-R pair, (",from,",",to,")")
        }
        adjmat[as.integer(from),as.integer(to)] = 1
        adjmat[as.integer(to),as.integer(from)] = 1
        N.c[1] = N.c[1] + 1; N.d[1] = N.d[1] - 1
      }
      if(z==4){
        # reconnect a random S-I/R-I pair
        if(model == "SIS"){
          if (N.S == 1){
            from = susceptible
            to.ind = sample(which(adjmat.si==0),1); to = infected[to.ind]
          }else if (N.I == 1){
            num.dis = N.I - adjmat.si
            from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
            to = infected
          }else{
            num.dis = N.I - rowSums(adjmat.si)
            from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
            to.ind = sample(which(adjmat.si[from.ind,]==0),1); to = infected[to.ind]
          }
          msg = paste0("reconnect an S-I pair, (",from,",",to,").")
        }else{
          #N.SR = N.S + N.R; 
          # 03/27/2020: debug
          sus.recov = union(susceptible, recovered)
          N.SR = length(sus.recov)
          if (N.SR == 1){
            from = sus.recov
            # 03/27/2020: debug
            num.dis = N.SR - as.vector(adjmat.ri)
            to.ind = sample(N.I, 1, prob = num.dis); to = infected[to.ind]
            #to.ind = sample(which(adjmat.ri==0),1); to = infected[to.ind]
          }else if (N.I == 1){
            num.dis = N.I - as.vector(adjmat.ri)
            from.ind = sample(N.SR,1,prob = num.dis); from = sus.recov[from.ind]
            to = infected
          }else{
            num.dis = N.I - rowSums(adjmat.ri)
            # debugging...
            #cat(N.S, N.R, N.I, '\n')
            #cat(length(susceptible), length(recovered), length(infected),'\n')
            #print(dim(adjmat.ri))
            from.ind = sample(N.SR,1,prob = num.dis); from = sus.recov[from.ind]
            
            # 03/28/2020: debug
            vec.dis = 1 - as.vector(adjmat.ri[from.ind,])
            to.ind = sample(N.I,1,prob = vec.dis); to = infected[to.ind]
            #to.ind = sample(which(adjmat.ri[from.ind,]==0),1); to = infected[to.ind]
          }
          msg = paste0("reconnect an S-I/R-I pair, (",from,",",to,").")
        }
        adjmat[as.integer(from),as.integer(to)] = 1
        adjmat[as.integer(to),as.integer(from)] = 1
        # 03/27/2020
        # try to fix a bug
        if(length(susceptible)>0){
          adjmat.si = adjmat[susceptible,infected]
          N.SI = nnzero(adjmat.si)
        }else{
          N.SI = 0
        }
        
        if(model == "SIR"){ adjmat.ri = adjmat[sus.recov,infected] }
        
        N.c[2] = N.c[2] + 1; N.d[2] = N.d[2] - 1
        #N.SI = N.c[2]
      }
      if(z==5){
        # reconnect a random I-I pair
        num.dis = (N.I - 1) - rowSums(adjmat.ii)
        from.ind = sample(N.I,1,prob = num.dis); from = infected[from.ind]
        to.cands = which(adjmat.ii[from.ind,]==0)
        to.cands = to.cands[to.cands!=from.ind]
        if(length(to.cands) == 1){
          to.ind = to.cands
        }else{
          to.ind = sample(to.cands,1)
        }
        to = infected[to.ind]
        adjmat.ii[as.integer(from.ind),as.integer(to.ind)] = 1
        adjmat.ii[as.integer(to.ind),as.integer(from.ind)] = 1
        adjmat[as.integer(from),as.integer(to)] = 1
        adjmat[as.integer(to),as.integer(from)] = 1
        N.c[3] = N.c[3] + 1; N.d[3] = N.d[3] - 1
        
        msg = paste0("reconnect an I-I pair, (",from,",",to,").")
      }
      if(z==6){
        # disconnect a random S-S/R-S/R-R pair
        dis.ind = sample(N.c[1], 1) 
        if(model == "SIS"){
          from.ind = adjmat.ss@i[dis.ind] + 1; to.ind = adjmat.ss@j[dis.ind] + 1
          adjmat.ss[as.integer(from.ind),as.integer(to.ind)] = 0
          adjmat.ss[as.integer(to.ind),as.integer(from.ind)] = 0
          
          from = susceptible[from.ind]; to = susceptible[to.ind]
          
          msg = paste0("disconnect an S-S pair, (",from,",",to,").")
        }else{
          sus.recov = union(susceptible, recovered)
          from.ind = adjmat.rs@i[dis.ind] + 1; to.ind = adjmat.rs@j[dis.ind] + 1
          adjmat.rs[as.integer(from.ind),as.integer(to.ind)] = 0
          adjmat.rs[as.integer(to.ind),as.integer(from.ind)] = 0
          
          from = sus.recov[from.ind]; to = sus.recov[to.ind]
          
          msg = paste0("disconnect an S-S/S-R/R-R pair, (",from,",",to,").")
        }
        adjmat[as.integer(from),as.integer(to)] = 0
        adjmat[as.integer(to),as.integer(from)] = 0
        
        N.c[1] = N.c[1] - 1; N.d[1] = N.d[1] + 1
        
        # ## sometimes from/to is just empty... check it out
        # if(length(from)!=1){
        #   stop("from.ind = ", from.ind, 
        #        " but susceptible length is ", length(susceptible),
        #        " N.c[1] =", N.c[1], " dis.ind =", dis.ind, 
        #        " while adjmat.ss@i has length ", length(adjmat.ss@i))
        # }
        # ##
      }
      if(z==7){
        # disconnect a random S-I/R-I pair
        # dis.ind = sample(N.c[2], 1)
        if(model == "SIS"){
          if (N.S == 1){
            from = susceptible
            to.ind = sample(which(adjmat.si==1),1); to = infected[to.ind]
            adjmat[as.integer(from),as.integer(to)] = 0
            adjmat[as.integer(to),as.integer(from)] = 0
            adjmat.si = adjmat[susceptible,infected]
          }else if (N.I == 1){
            num.con = as.vector(adjmat.si)
            from.ind = sample(N.S,1,prob = num.con); from = susceptible[from.ind]
            to = infected
            adjmat[as.integer(from),as.integer(to)] = 0
            adjmat[as.integer(to),as.integer(from)] = 0
            adjmat.si = adjmat[susceptible,infected]
          }else{
            dis.ind = sample(N.c[2], 1)
            from = adjmat.si@i[dis.ind] + 1; to = adjmat.si@j[dis.ind] + 1
            adjmat.si[as.integer(from),as.integer(to)] = 0
            
            from = susceptible[from]; to = infected[to]
            adjmat[as.integer(from),as.integer(to)] = 0
            adjmat[as.integer(to),as.integer(from)] = 0
          }
          msg = paste0("disconnect an S-I pair, (",from,",",to,").")
        }else{
          N.SR = N.S + N.R; sus.recov = union(susceptible, recovered)
          if (N.SR == 1){
            from = sus.recov
            to.ind = sample(which(adjmat.ri==1),1); to = infected[to.ind]
            adjmat[as.integer(from),as.integer(to)] = 0
            adjmat[as.integer(to),as.integer(from)] = 0
            adjmat.ri = adjmat[sus.recov,infected]
          }else if (N.I == 1){
            num.con = adjmat.ri
            from.ind = sample(N.SR,1,prob = num.con); from = sus.recov[from.ind]
            to = infected
            adjmat[as.integer(from),as.integer(to)] = 0
            adjmat[as.integer(to),as.integer(from)] = 0
            adjmat.ri = adjmat[sus.recov,infected]
          }else{
            dis.ind = sample(N.c[2], 1)
            from = adjmat.ri@i[dis.ind] + 1; to = adjmat.ri@j[dis.ind] + 1
            # debugging...
            #cat('dis.ind=',dis.ind,'\n')
            #cat('from=',from,', to=',to, '\n')
            #cat(dim(adjmat.ri), '\n')
            #cat(length(adjmat.ri@i),  length(adjmat.ri@j), N.c[2], '\n')
            
            if(length(adjmat.ri@i) != N.c[2]){
              cat('Discrepancy when disconnecting H-I...\n')
              stop(paste0("At time ", t.next, ", ",msg,"\n"))
            }
            
            #adjmat.ri[as.integer(from),as.integer(to)] = 0
            
            from = sus.recov[from]; to = infected[to]
            adjmat[as.integer(from),as.integer(to)] = 0
            adjmat[as.integer(to),as.integer(from)] = 0
            adjmat.ri = adjmat[sus.recov,infected]
          }
          if(length(susceptible)>0){
            adjmat.si = adjmat[susceptible,infected]
          }else{
            adjmat.si = NULL
          }
          msg = paste0("disconnect an S-I/R-I pair, (",from,",",to,").")
        }
        
        N.c[2] = N.c[2] - 1; N.d[2] = N.d[2] + 1
        N.SI = nnzero(adjmat.si)
        
      }
      if(z==8){
        # disconnect a random I-I pair
        dis.ind = sample(N.c[3], 1) 
        from = adjmat.ii@i[dis.ind] + 1; to = adjmat.ii@j[dis.ind] + 1
        adjmat.ii[as.integer(from),as.integer(to)] = 0
        adjmat.ii[as.integer(to),as.integer(from)] = 0
        
        from = infected[from]; to = infected[to]
        adjmat[as.integer(from),as.integer(to)] = 0
        adjmat[as.integer(to),as.integer(from)] = 0
        
        N.c[3] = N.c[3] - 1; N.d[3] = N.d[3] + 1
        
        msg = paste0("disconnect an I-I pair, (",from,",",to,").")
      }
    }else{
      z = sample(4, 1, prob = ratios)
      if(z==3){
        # reconnect a random pair
        num.dis = (N - 1) - rowSums(adjmat)
        from = sample(N,1,prob = num.dis)
        to.cands = which(adjmat[from,]==0)
        to.cands = to.cands[to.cands!=from]
        if(length(to.cands) == 1){
          to = to.cands
        }else{
          to= sample(to.cands,1)
        }
        adjmat[as.integer(from), as.integer(to)] = 1
        adjmat[as.integer(to),as.integer(from)] = 1
        
        N.c = N.c + 1; N.d = N.d - 1
        
        # UPDATE: if any of S or I is empty now
        if(length(susceptible) == 0 | length(infected) == 0){
          adjmat.si = NULL; N.SI = 0
        }else{
          adjmat.si = adjmat[susceptible,infected]
          N.SI = nnzero(adjmat.si)
        }
        
        msg = paste0("reconnect a pair, (",from,",",to,").")
      }
      if(z==4){
        # disconnect a random pair
        dis.ind = sample(N.c, 1)
        from = adjmat@i[dis.ind] + 1; to = adjmat@j[dis.ind] + 1
        adjmat[from,to] = 0; adjmat[to,from] = 0
        N.c = N.c - 1; N.d = N.d + 1
        
        # UPDATE: if any of S or I is empty now
        if(length(susceptible)==0 | length(infected)==0){
          adjmat.si = NULL; N.SI = 0
        }else{
          adjmat.si = adjmat[susceptible,infected]
          N.SI = nnzero(adjmat.si)
        }
        
        msg = paste0("disconnect a pair, (",from,",",to,").")
      }
    }
    
    #2.2 deal with epidemic progress
    if(z==2){
      # recover an infected person
      rec.ind = sample(N.I,1)
      recover = infected[rec.ind]
      infected = infected[-rec.ind]
      N.I = N.I - 1 
      if(model == "SIS"){
        susceptible = c(susceptible,recover)
        epid[recover] = 0
        N.S = N.S + 1
      }else{
        epid[recover] = -1
        if(quarantine){
          recovered = c(recovered, recover)
          N.R = N.R + 1
        }
      }
      
      msg = paste0("individual ", recover,
                   " is recovered. Disease prevalence = ",
                   N.I/N, ".")
    }
    if(z==1){
      # select a susceptible person to infect
      if(N.S == 1){
        # only one person to infect
        new.infec = susceptible
        susceptible = NULL
      }else{
        # multiple susceptible inviduals to choose from
        if(N.I == 1){
          # only one infected person
          num.nei.i = adjmat.si
        }else{
          num.nei.i = rowSums(adjmat.si)
        }
        infec.ind = sample(N.S,1,prob = num.nei.i)
        new.infec = susceptible[infec.ind]
        susceptible = susceptible[-infec.ind]
      }
      epid[new.infec] = 1
      N.S = N.S - 1
      infected = c(infected, new.infec); N.I = N.I + 1
      
      msg = paste0("individual ", new.infec,
                   " is infected. Disease prevalence = ",
                   N.I/N, ".")
    }
    
    # A sanity check: if N.I==0, stop the simulation
    if(N.I==0){
      if(verbose){
        cat(paste0("At time ", t.next, ", "),msg,"\n")
      }
      # still update inference statistics
      del.t = t.next - t.cur
      z = as.integer(z)
      event.counts[z] = event.counts[z] + 1
      tot.event = tot.event + 1
      
      big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
      
      # still have to record the event
      time = c(time, t.next)
      event.type = c(event.type, z)
      preval = c(preval, N.I)
      dens = c(dens, nnzero(adjmat))
      # here: a recovery must be the latest event
      per1 = c(per1, recover)
      per2 = c(per2, NA)
      
      cat("Disease extinct at time", t.next, ", simulation stops.\n")
      break
    }
    
    
    # # check whether adjmat is still symmetric--things have gone wrong somewhere
    # if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
    #   stop(paste("The adjmat isn't symmetric any more, where nonzero =", 
    #              nnzero(adjmat), "and that's after z =", z,
    #              "from.ind =", from.ind, "to.ind = ", to.ind, 
    #              "to.cands =", paste(to.cands,collapse = ","),
    #              "from =", from, "to =", to,"\n"))
    # }
    
    
    # 2.2.* make updates about S-I, S-S, I-I connections if z==1 or 2
    # while handling edge cases of N.S <= 1 and/or N.I == 1
    if(z %in% c(1,2)){
      if(N.S == 0){
        adjmat.si = NULL
        N.SI = 0
        if(quarantine && model == "SIR" && N.R == 0){
          adjmat.ri = NULL
          N.c[2] = 0
        }
        # 03/27/2020: debug
        if(quarantine && model == "SIR" && N.R > 0){
          adjmat.ri = adjmat[recovered,infected]
          N.c[2] = nnzero(adjmat.ri)
        }
      }else{
        adjmat.si = adjmat[susceptible,infected]
        N.SI = nnzero(adjmat.si)
        if(quarantine && model == "SIR"){
          adjmat.ri = adjmat[union(susceptible,recovered),infected]
          N.c[2] = nnzero(adjmat.ri)
        }
      }
      if(quarantine){
        if(model == "SIS"){
          if(N.S <= 1){
            adjmat.ss = 0
            N.c[1] = 0
            N.d[1] = 0
          }else{
            adjmat.ss = adjmat[susceptible,susceptible]
            N.c[1] = nnzero(adjmat.ss)/2
            N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
          }
          N.c[2] = N.SI
          N.d[2] = N.S * N.I - N.c[2]
        }else{
          #N.SR = N.S + N.R; 
          sus.recov = union(susceptible,recovered)
          # 03/27/2020: debug
          N.SR = length(sus.recov)
          
          if(N.SR <= 1){
            adjmat.rs = 0
            N.c[1] = 0
            N.d[1] = 0
          }else{
            adjmat.rs = adjmat[sus.recov,sus.recov]
            N.c[1] = nnzero(adjmat.rs)/2
            N.d[1] = N.SR * (N.SR - 1)/2 - N.c[1]
          }
          N.d[2] = N.SR * N.I - N.c[2]
        }
        if(N.I == 1){
          adjmat.ii = 0
          N.c[3] = 0
          N.d[3] = 0
        }else{
          adjmat.ii = adjmat[infected,infected]
          N.c[3] = nnzero(adjmat.ii)/2
          N.d[3] = N.I * (N.I - 1)/2 - N.c[3]
        }
      }
    }
    
    # 3. update inference statistics and do inference occasionally
    # 3.1 update stats
    del.t = t.next - t.cur
    z = as.integer(z)
    event.counts[z] = event.counts[z] + 1
    tot.event = tot.event + 1
    
    big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
    
    # 3.2 inference if ...
    if(tot.event %% infer.interval == 0){
      if(Bayes){
        # Bayesian inference: Gamma sampling
        a.post = a.pr + event.counts
        b.post = b.pr + big.sums
        if(infer.verbose){
          cat("With", tot.event,"events, posterior means and variances:\n",
              a.post/b.post, "\n",
              a.post/(b.post^2), "\n")
        }
        # record posterior stats
        if(return.infer){
          m.post = a.post/b.post
          lb.post = qgamma(.025, a.post, b.post)
          ub.post = qgamma(.975, a.post, b.post)
          this.chunk = cbind(m.post, lb.post, ub.post, event.counts, tot.event)
          Bayes.tab = rbind(Bayes.tab, this.chunk)
        }
        # plot
        if(infer.plot){
          par(mfrow=c(2,2))
          for(ix in 1:length(params)){
            v = vars[ix]; p = params[ix]
            sams = rgamma(samples, shape=a.post[ix], rate = b.post[ix])
            xl = range(sams, p)
            hist(sams, main=paste("Posterior samples for",vars[ix]),
                 xlab = vars[ix], xlim=xl, col="lightblue")
            abline(v=p,col="red")
            if(demo.sleep & ix%%4 == 0){
              Sys.sleep(0.2)
            }
          }
        }
      }
      if (MLE){
        # MLE
        MLEs = event.counts/big.sums
        if(infer.verbose){
          cat("With",tot.event, "events, MLEs:","\n", MLEs, "\n")
        }
        MLE.tab = rbind(MLE.tab, c(tot.event, MLEs, event.counts))
      }
    }
    
    # updated 02/04/2020
    # 3.3 examine network degree distribution if
    if(net.deg.plot & (tot.event %% net.interval == 0)){
      degs = rowSums(adjmat)
      avg.deg = mean(degs)
      poi.degs = rpois(N, avg.deg)
      
      deg.dat = data.frame(Degrees = degs, DegPoi = poi.degs)
      print(
        ggplot(deg.dat) + geom_density(aes(x = Degrees, color="Actual degree")) +
          geom_density(aes(x = DegPoi, color="Poisson")) +
          labs(color="Distribution", 
               caption = paste("At t =", round(t.cur, digits = 2))) + 
          theme_bw(base_size = 14) %+replace%
          theme(legend.position = c(0.8,0.8))
      )
    }
    
    
    
    # 4. report the event and summarize and set t.cur
    
    t.cur = t.next
    
    # record the event
    time = c(time, t.cur)
    event.type = c(event.type, z)
    preval = c(preval, N.I)
    dens = c(dens, nnzero(adjmat))
    # record labels of individual(s) involved
    if (z >= 3){
      per1 = c(per1, from)
      per2 = c(per2, to)
    }else if (z==2){
      per1 = c(per1, recover)
      per2 = c(per2, NA)
    }else{
      per1 = c(per1, new.infec)
      per2 = c(per2, NA)
    }
    
    
    if(verbose){
      cat(paste0("At time ", t.cur, ", "),msg,"\n")
    }
    
  }
  
  # return results/summary
  cat("Simulation finished. At the end: \n", 
      "Disease prevalence =", N.I/N, "\n",
      "Network density = ", nnzero(adjmat)/(N*(N-1)),"\n")
  event.log = data.frame(time = time, event = event.type,
                         per1 = per1, per2 = per2,
                         preval = preval/N, dens = dens/(N*(N-1)))
  
  # inference at the end
  if(Bayes){
    # Bayesian inference: Gamma sampling
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    cat("Posterior means and variances with", tot.event, "events in total:\n",
        a.post/b.post, "\n",
        a.post/(b.post^2), "\n")
    cat("Event type counts:\n", event.counts,"\n")
    # record posterior stats
    if(return.infer){
      m.post = a.post/b.post
      lb.post = qgamma(.025, a.post, b.post)
      ub.post = qgamma(.975, a.post, b.post)
      this.chunk = cbind(m.post, lb.post, ub.post, event.counts, tot.event)
      Bayes.tab = rbind(Bayes.tab, this.chunk)
    }
    # plot
    if(infer.plot){
      par(mfrow=c(2,2))
      for(ix in 1:length(params)){
        v = vars[ix]; p = params[ix]
        sams = rgamma(samples, shape=a.post[ix], rate = b.post[ix])
        xl = range(sams, p)
        hist(sams, main=paste("Posterior samples for",vars[ix]),
             xlab = vars[ix], xlim=xl, col="lightblue")
        abline(v=p,col="red")
        if(demo.sleep & ix%%4 == 0){
          Sys.sleep(0.2)
        }
      }
    }
  }
  if(MLE){
    # MLE
    MLEs = event.counts/big.sums
    cat("MLEs with", tot.event, "events in total:","\n", MLEs, "\n")
    cat("Event type counts:\n", event.counts,"\n")
    
    MLE.tab = rbind(MLE.tab, c(tot.event, MLEs, event.counts))
    MLE.tab = as.data.frame(MLE.tab)
    if(quarantine){
      names(MLE.tab) = c("event.num","beta","gamma",
                         "alpha.SS","alpha.SI","alpha.II",
                         "omega.SS","omega.SI","omega.II",
                         "infection","recovery",
                         "SS.connect","SI.connect","II.connect",
                         "SS.discon","SI.discon","II.discon")
    }else{
      names(MLE.tab) = c("event.num","beta","gamma","alpha","omega",
                         "infection","recovery","connect","discon")
    }
  }
  
  if(plot){
    subt = paste(model," beta =",bet," gamma =", gam, 
                 " alpha.r =", alpha.r, " alpha.d =", alpha.d)
    plot(preval ~ time, data = event.log, 
         xlab="Time", ylab="Disease prevalence", 
         main = "Disease Prevalence vs. Time", sub = subt,
         ylim = c(0,1), type="l",lwd=2)
    plot(dens ~ time, data = event.log,
         xlab="Time", ylab="Edge density", 
         main = "Network Edge Density vs. Time", sub = subt,
         type="l",lwd=2)
  }
  
  if(infer.plot){
    #params = c(bet, gam, alpha.r, alpha.d)
    if(Bayes){
      par(mfrow=c(1,1))
    }else{
      # a function to calculate the frequentist 95% CI for MLEs
      MLE.CI <- function(x){
        l.hat = x[1]; n.obs = x[2]
        if (is.na(l.hat) || l.hat == 0){
          lb = 0; ub = 0
        }else if (n.obs >= 20){
          lb = l.hat * (1 + qnorm(.025)/sqrt(n.obs))
          ub = l.hat * (1 + qnorm(.975)/sqrt(n.obs))
        }else{
          lb = l.hat * qchisq(.025, 2*n.obs) / (2*n.obs)
          ub = l.hat * qchisq(.975, 2*n.obs) / (2*n.obs)
        }
        return(c(lb,ub))
      }
      
      lcol = length(names(MLE.tab)[-1])
      vars = names(MLE.tab)[-1][1:(lcol/2)]
      counts = names(MLE.tab)[-1][(lcol/2 + 1):lcol]
      e.num = MLE.tab$event.num
      if(quarantine){
        par(mfrow=c(2,2))
      }else{
        par(mfrow=c(2,2))
      }
      for(ix in 1:length(vars)){
        v = vars[ix]
        if(all(is.na(MLE.tab[[v]]))){
          plot.new()
        }else{
          tab.v = MLE.tab[,c(v,counts[ix])]
          bounds.v = t(apply(tab.v, 1, MLE.CI))
          bounds.v = as.data.frame(bounds.v)
          names(bounds.v) = c("lb","ub")
          #s.m = mean(MLE.tab[[v]], na.rm=T)
          #s.sd = sd(MLE.tab[[v]], na.rm=T)
          yl = range(min(bounds.v$lb), max(bounds.v$ub), params[ix])
          plot(MLE.tab[[v]]~e.num, xlab="# events", ylab=v,
               main = paste("MLE for", v, "v.s. # events"), 
               type="l", ylim=yl)
          lines(bounds.v$lb~e.num, lty=2, col="gray")
          lines(bounds.v$ub~e.num, lty=2, col="gray")
        }
        abline(h=params[ix], col="red")
      }
      par(mfrow=c(1,1))
    }
  }
  
  # add: print out big.sums if Bayesian
  if(Bayes){
    cat("big.sums: \n", big.sums,"\n")
  }
  
  # 02/04/2020: degree distribution at the end
  if(net.deg.plot){
    degs = rowSums(adjmat)
    avg.deg = mean(degs)
    poi.degs = rpois(N, avg.deg)
    
    deg.dat = data.frame(Degrees = degs, DegPoi = poi.degs)
    print(
      ggplot(deg.dat) + geom_density(aes(x = Degrees, color="Actual degree")) +
        geom_density(aes(x = DegPoi, color="Poisson")) +
        labs(color="Distribution", 
             caption = paste("At the end, t =", round(t.cur, digits = 2))) + 
        theme_bw(base_size = 14) %+replace%
        theme(legend.position = c(0.8,0.8))
    )
  }
  
  
  # return event data
  # MODIFICATION: return true parametes too
  if(return.infer){
    if(Bayes){
      Bayes.tab = as.data.frame(Bayes.tab)
      Bayes.tab$var = vars
      infer.tab = Bayes.tab
    }else{
      infer.tab = MLE.tab
    }
    return(list(G0 = init.adj, I0 = init.infected,
                events = event.log, infer=infer.tab,
                infer2 = MLE.tab,
                truth = params))
  }else{
    return(list(G0 = init.adj, I0 = init.infected,
                events = event.log,truth = params))
  }
  
}
# 05/17/2019
# repeated simulations
# 05/21/2019
# Modified: use the infer2 function instead
# and filter out results with too few total events

# RIGHT NOW: it only supports Bayesian inference plotting
# 06/03/2019: 
# shrink the Bayesian plot sizes
# and only retain the legend for the 1 st plot
rep_stochastic_coevolve_infer <- function(..., Bayesian = F, seed = 42, rep=10,
                                          savedat = F, fname="exdat",
                                          event.thres = 800, pdfpath = NULL){
  set.seed(seed)
  
  # go through repetitions
  bigdat = NULL
  i = 0
  while(i < rep){
    res = stochastic_coevolve_infer2(..., Bayes = Bayesian)
    # only keep the data when tot.event >= event.thres
    if(nrow(res$events) >= event.thres){
      i = i+1
      if(savedat){
        fpath = paste0(fname,"_",i,".rds")
        saveRDS(res, fpath)
      }
      infer.dat = res$infer
      infer.dat$rep = i
      bigdat = rbind(bigdat,infer.dat)
      cat("\n\n Repetion No.",i,"done.\n\n")
    }
  }
  
  # make plots
  params = res$truth
  if(Bayesian){
    vars = unique(bigdat$var)
    names(params) = vars
    bigdat$rep = as.factor(bigdat$rep)
    if(!is.null(pdfpath)){
      pdf(paste0(pdfpath,".pdf"), width = 5, height = 2.8)
    }
    for(v in vars){
      vardat = bigdat %>% filter(var==v)
      pp = ggplot(vardat, aes(x=tot.event)) + 
        geom_hline(yintercept = params[v], size=2, color="black") +
        geom_ribbon(aes(ymin = lb.post, ymax = ub.post, fill = rep), 
                    alpha=.2) +
        geom_line(aes(y=m.post, color=rep), size = 1) +
        labs(x="# total events", y="value", title=v, 
             color="repetition") +
        theme_classic() 
      if(which(vars==v) > 1){
        print(
          pp + guides(fill=FALSE, color=FALSE)
        )
      }else{
        print(
          pp + guides(fill=FALSE)
        )
      }
    }
    if(!is.null(pdfpath)){
      dev.off()
    }
  }else{
    lcol = length(names(bigdat)[-1])
    vars = names(bigdat)[-1][1:(lcol/2)]
    counts = names(bigdat)[-1][(lcol/2 + 1):lcol]
    par(mfrow=c(2,2))
    if(!is.null(pdfpath)){
      pdf(paste0(pdfpath,".pdf"), width = 8, height = 5)
    }
    for(ix in 1:length(vars)){
      v = vars[ix]
      yl = range(bigdat[[v]], params[ix], na.rm = T)
      xl = range(bigdat$event.num)
      for(i in 1:rep){
        tab.v = bigdat[bigdat$rep == i, c("event.num",v,counts[ix])]
        if(i == 1){
          plot(tab.v[[v]] ~ tab.v$event.num, xlab="# events", ylab=v,
               main = paste("MLEs for", v, "v.s. # events"), 
               type="l", ylim=yl, xlim=xl)
        }else{
          lines(tab.v[[v]] ~ tab.v$event.num)
        }
      }
      abline(h=params[ix], col="red", lwd=2)
    }
    if(!is.null(pdfpath)){
      dev.off()
    }
  }
  
  return(bigdat)
}
