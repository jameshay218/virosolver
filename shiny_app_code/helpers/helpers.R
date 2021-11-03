library(tidyverse)
library(lubridate)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(devtools)


test_func <- function(){
  x <- c(1,2,3,4,5)
  y <- c(2,3,4,5,6)
  plot(x,y)
}

long_dat <- function (data, district="Delhi") {
  data <- data %>% filter(District == "Delhi")
  delhi_data_long <- data %>% pivot_longer(-c(Date,State,District))
  delhi_data_long <- delhi_data_long %>% rename(Type=name,cases=value) %>%
    group_by(Type, State, District) %>%
    mutate(new_cases=cases-lag(cases,1)) %>%
    ## Smooth incidence into 7 day rolling average
    mutate(cases_7day=zoo::rollmean(new_cases, k=7,fill=0))
  
  return(delhi_data_long)
}

plot_cases <- function(data, region_name) {
  data_long <- long_dat(data)
  #browser()
  ggplot(data_long %>% filter(Type %in% c("Confirmed","Deceased"))) +
    geom_line(aes(x=Date,y=cases_7day,col=Type)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_color_manual(values=c("Confirmed"="blue","Deceased"="orange")) +
    ylab("New cases") +
    xlab("Reporting date") +
    ggtitle("New cases over time in district of Delhi, India") +
    facet_wrap(~Type,scales="free")
}

long_ct_dat <- function(ct_dat) {
  colnames(ct_dat) <- c("id","date","age_gender","remarks","E","N")
  ct_dat$date <- mdy(ct_dat$date)
  head(ct_dat)
  
  ## Only look at Cts <= 35
  ct_dat_long <- ct_dat %>% pivot_longer(-c(id,date,age_gender,remarks)) %>% 
    rename(gene=name,Ct=value) %>% 
    filter(Ct <= 35) %>%
    mutate(Ct_round=round(Ct, 0))
  
  return(ct_dat_long)
}

plot_cts <- function(ct_dat) {
  ct_dat_long <- long_ct_dat(ct_dat)
  ## Plot N gene all Cts over time
  p_raw <- ct_dat_long %>% ggplot() + 
    geom_point(aes(x=date,y=Ct_round),size=0.25,alpha=0.25) + 
    geom_smooth(aes(x=date,y=Ct_round)) + 
    scale_y_continuous(trans="reverse") +
    scale_x_date(breaks="1 month")+
    ggtitle("Ct values over time") +
    xlab("Date") +
    ylab("Ct") +
    theme_classic() +
    facet_wrap(~gene, nrow=2)

  return(p_raw)
}

plot_daily_mean <- function(ct_dat) { 
  ct_dat_long <- long_ct_dat(ct_dat)
  ct_dat_summary <- ct_dat_long %>% 
    group_by(date, gene) %>%
    summarize(skew_ct=moments::skewness(Ct),
              mean_ct=mean(Ct),
              N=n())
  
  p_daily_mean <- ct_dat_summary %>% 
    ggplot() + 
    geom_point(aes(x=date,y=mean_ct),size=1,alpha=1) + 
    geom_smooth(aes(x=date,y=mean_ct),span=0.2) + 
    scale_y_continuous(trans="reverse")+
    scale_x_date(breaks="1 month")+
    ggtitle("Daily mean and smoothed mean Cts") +
    xlab("Date") +
    ylab("Mean Ct") +
    theme_classic()  +
    facet_wrap(~gene, nrow=2)
  
  return(p_daily_mean)
}

plot_daily_skew <- function(ct_dat) {
  ct_dat_long <- long_ct_dat(ct_dat)
  
  ct_dat_summary <- ct_dat_long %>% 
    group_by(date, gene) %>%
    summarize(skew_ct=moments::skewness(Ct),
              mean_ct=mean(Ct),
              N=n())
  
  p_daily_skew <- ct_dat_summary %>% 
    ggplot() + 
    geom_point(aes(x=date,y=skew_ct),size=1,alpha=1) + 
    geom_smooth(aes(x=date,y=skew_ct),span=0.2) + 
    scale_x_date(breaks="1 month")+
    scale_y_continuous(trans="reverse")+
    ggtitle("Daily skew and smoothed skew of Cts") +
    xlab("Date") +
    ylab("Skewness of Ct values") +
    theme_classic()  +
    facet_wrap(~gene, nrow=2)
  
  return(p_daily_skew)
}