parTab <- read.csv("~/Documents/GitHub/ct_inference/pars/parTab_test_seeirr.csv")
pars1 <- parTab$values
names(pars1) <- parTab$names

pars2 <- c(beta = pars1["R0"]/pars1["infectious"], sigma = 1/pars1["latent"],
           alpha = 1/pars1["incubation"], gamma = 1/pars1["infectious"],
           omega = 1/pars1["recovery"])
names(pars2) <- c("beta","sigma","alpha","gamma","omega")

init <- c(1-pars1["I0"],0,0,pars1["I0"],0,0,0)

ts <- seq(0,365,by=1)

y <- solveSEEIRRModel_rlsoda(ts, init, pars2, TRUE)
colnames(y) <- c("t","S","E1","E2","I","R1","R2","inc")
y <- as_tibble(y)
y$Detectable <- y$E2 + y$I + y$R1
y1 <- y %>% as_tibble %>% pivot_longer(-t)

ggplot(y1) + geom_line(aes(x=t,y=value,col=name))

create_post_func_seeirr <- function(parTab, data, ts, PRIOR_FUNC=NULL,ver="likelihood"){
  par_names <- parTab$names
  test_times <- unique(data$date)
  observed_prev <- data$POS
  N_obs <- data$POS + data$NEG

  f <- function(pars){
    names(pars) <- par_names
    prev <- detectable_SEEIRRModel(pars, ts)
    if(ver == "likelihood"){
      prev <- prev[which(ts %in% test_times)]
      lik <- sum(dbinom(observed_prev, N_obs, prev,log=TRUE))
      if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
      return(lik)
    } else {
      return(prev)
    }
  }
  f
}

prior_func_me <- function(pars){
  names(pars) <- parTab$names
  p1 <- dnorm(pars["R0"], 0, 100,log=TRUE)
  p2 <- dnorm(pars["latent"],2,0.5,log=TRUE)
  p3 <- dnorm(pars["incubation"],2,0.5,log=TRUE)
  p4 <- dnorm(pars["infectious"],10,2,log=TRUE)
  p5 <- dnorm(pars["recovery"],5, 1, log=TRUE)
  return(sum(p1, p2,p3,p4,p5))
}


model_func <- create_post_func_seeirr(parTab, ts=ts,NULL, NULL,"model")
dat <- f(parTab$values)

N <- 1000
obs_times <- c(50,150)
dat1 <- rbinom(length(dat), N,dat)
dat1 <- tibble(date=ts[1:length(dat)],POS=dat1,NEG=N -dat1)
dat1 <- dat1 %>% filter(date %in% obs_times)


f1 <- create_post_func_seeirr(parTab, dat1,ts=ts, PRIOR_FUNC=prior_func_me, "likelihood")
f1(parTab$values)

mcmcPars <- c("iterations"=50000,"popt"=0.234,"opt_freq"=1000,
               "thin"=1,"adaptive_period"=20000,"save_block"=1000)
## Get random starting values
parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_times)

startTab <- generate_start_tab(parTab)
covMat <- diag(nrow(startTab))
mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)

output <- run_MCMC(parTab=parTab, data=dat1, mcmcPars=mcmcPars,
                   filename="~/Documents/viral_load_model_test/chains/seeirr",
                   CREATE_POSTERIOR_FUNC = create_post_func_seeirr, mvrPars = mvrPars,
                   PRIOR_FUNC=prior_func_me, ts=ts,ver="likelihood")
chain <- read.csv(output$file)
chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]
plot(coda::as.mcmc(chain))
par(ask=FALSE,mfrow=c(1,1))

samps <- sample(unique(chain$sampno), 1000)
trajs <- matrix(0, nrow=1000, ncol=365)

for(i in 1:length(samps)){
  samp <- samps[i]
  tmp_pars <- get_index_pars(chain, samp)
  trajs[i,] <- model_func(tmp_pars)[1:365]
}
preds <- t(apply(trajs,2, function(x) quantile(x, c(0.025,0.5,0.975))))
plot(preds[,2],type='l',ylim=c(0,1))
lines(preds[,1])
lines(preds[,3])
lines(dat,col="blue")
