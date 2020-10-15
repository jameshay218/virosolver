chain <- read.csv("test1_univariate_chain.csv")
chain <- chain[chain$sampno > 10000,]
samps <- sample(unique(chain$sampno),1000)
i <- 10
convert_names <- which(colnames(chain) %in% c("prob",paste0("prob.",1:(length(times)))))-1
tmp_pars <- get_index_pars(chain, samps[i])
names(tmp_pars)[convert_names] <- "prob"
omg <- gaussian_process_model(tmp_pars, times)
plot(omg,type='l')
lines(prob_infection,col="red")

preds <- matrix(0, nrow=length(samps),ncol=length(omg))

for(i in seq_along(samps)){

  tmp_pars <- get_index_pars(chain, samps[i])
  names(tmp_pars)[convert_names] <- "prob"
  omg <- gaussian_process_model(tmp_pars, times)
  preds[i,] <- omg
  #lines(omg,type='l',col=i)
}

lines(prob_infection,col="red",lwd=2)

quants <- apply(preds, 2, function(x) quantile(x, c(0.025,0.5,0.975)))
plot(quants[3,],type='l',ylim=c(0,0.01))
lines(quants[2,])
lines(quants[1,])
lines(prob_infection,col="red",lwd=2)
