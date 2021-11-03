## Assumes user has simulated data 
plot_vk <- function(param_table, ct_data=sim_cts) {
  pars <- param_table$values
  names(pars) <- param_table$names
  
  ## Solve the Ct model over a range of times since infection (referred to as "ages")
  test_ages <- seq(1,50,by=1) #FIXME: user sets "ages" (OR read from input)!
  
  ## This gives the modal Ct value
  cts <- viral_load_func(pars, test_ages)
  
  if (is.null(ct_data)) {
    p_ct_model <- ggplot(data.frame(ct=c(40,cts),t=c(0,test_ages))) + 
      geom_line(aes(x=t,y=ct)) + 
      scale_y_continuous(trans="reverse",
                         limits=c(40,10)) +
      theme_bw() +
      ylab("Modal Ct value") +
      xlab("Days since infection")
  }
  else {
    p_ct_model <- ggplot(ct_data %>% filter(ct < 40)) +
      geom_point(aes(x=age,y=ct),alpha=0.10) +
      geom_violin(data=(ct_data %>% filter(ct <40) %>% filter(age %% 5 == 0 | age == 1)),
                  aes(x=age,group=age,y=ct),scale="width",fill="blue",
                  draw_quantiles=c(0.025,0.5,0.975),
                  alpha=0.5) + 
      geom_line(data=data.frame(ct=c(40,cts),t=c(0,test_ages)), 
                aes(x=t,y=ct)) + 
      scale_y_continuous(trans="reverse",
                         limits=c(40,10)) +
      theme_bw() +
      ylab("Modal Ct value") +
      xlab("Days since infection")
  }
  ggplotly(
    p <- p_ct_model,
    tooltip="all",
    dynamicTicks=FALSE
  )
}