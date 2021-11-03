read_mcmc_chains_tab <- function(chain_file,mcmc_1xs,input){
    ## If no chain file available, then return NULL
    if(is.null(chain_file))
        return(NULL)
    
    mcmc_thin <- input$mcmc_thin
    mcmc_burnin <- input$mcmc_adaptive
    mcmc_multi <- input$mcmc_multi
    
    ## Load in MCMC chains but only estimated parameters
    mcmc_chains_unfixed <- lazymcmc::load_mcmc_chains(chain_file$datapath,parTab=mcmc_1xs$seir_pars,
                                                      unfixed=TRUE,thin=mcmc_thin,burnin=mcmc_burnin,multi=mcmc_multi,chainNo=TRUE)
    mcmc_chains_unfixed <- as.data.frame(mcmc_chains_unfixed[[2]])
    
    ## Reshape this data frame into long format
    mcmc_chains_unfixed$chain <- as.factor(mcmc_chains_unfixed$chain)
    mcmc_chains_unfixed <- mcmc_chains_unfixed %>% tidyr::pivot_longer(-c("sampno","chain"), names_to="variable",values_to='value') 
    
    ## Load in MCMC chains exactly as saved each iteration
    mcmc_chains_fixed <- lazymcmc::load_mcmc_chains(chain_file$datapath,parTab=mcmc_1xs$seir_pars,
                                                    unfixed=FALSE,thin=mcmc_thin,burnin=mcmc_burnin,multi=mcmc_multi,chainNo=TRUE)
    mcmc_chains_fixed <- as.data.frame(mcmc_chains_fixed[[2]])
    
    return(list(mcmc_chains_unfixed, mcmc_chains_fixed))
}