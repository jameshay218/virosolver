###### Data Visualization Function Definitions ##########
## emake_long()
## Convert epidemic trajectory data into long format
## FIXME: this function assumes all numeric columns pertain to case counts. 
## FIXME: certain plots will only be interpretable if filtering is implemented
## for the epi tab. 
emake_long <- function (data=sample_epiDat, filters=hash()) 
  {
  filter_choices <- c()
  if (!is.empty(filters)) {
    for (filter in ls(filters)) {
      f_col <- str_split(filter,'_')[[1]][1] # FIXME: add safeguards here! 
      f_val <- filters[[filter]]
      filter_choices <- c(filter_choices, f_val)
      if (f_val != "All (No Filtering)") {
        data <- data %>% filter(data[[paste(f_col)]] == paste(f_val))
      }
    }
  }
  date_col <- colnames(data[, sapply(data, class) %in% c('Date')]) #requires exactly 1 date column
  ## if not formatted as date, requires date column to be called "date"
  if(is.null(date_col)) { date_col <- colnames(data)[tolower(colnames(data)) == "date"]} 
  names(data)[tolower(names(data)) == 'date'] <- 'date'
  date_col <- 'date'
  data[[date_col]] <- as.Date(data[[date_col]])
  
  browser()
  pivot_cols <- colnames(data[, sapply(data, class) %in% c('character', 'factor','Date')])
  gb_cols <- colnames(data[, sapply(data, class) %in% c('character', 'factor')])
  data_long <- data %>% pivot_longer(-pivot_cols) #pivot_longer(-c(Date,State,District))
  data_long <- data_long %>% rename(Type=name,cases=value) %>% 
    group_by_at(gb_cols) %>%
    mutate(new_cases=cases-lag(cases,1)) %>%
    ## Smooth incidence into 7 day rolling average
    mutate(cases_7day=zoo::rollmean(new_cases, k=7,fill=0))
  
  return(list(data_long,filter_choices))
}


## get_unique_subs()
## This function takes in a dataframe and column name
## and returns a UI selectInput object and the object's 
## name for referencing.
##
## Options for filtering are selected from unique data values
## from the given column. 
## This function is called from within the get_options function. 
get_unique_subs <- function(data, col,ns, numeric=FALSE) {
  unique_items <- unique(data[[col]])
  unique_items <- sort(unique_items[!is.na(unique_items)])
  obj_name <- paste0(col, "_options") ##FIXME: not checking for duplicate column names 
  if (numeric) {
    if (is.null(unique_items)) {
      return()
    }
    start = unique_items[1]
    stop = unique_items[length(unique_items)]
    choice_obj <- tags$div(
      'id'=ns(paste(obj_name)),
      sliderInput(
        inputId = ns(paste(obj_name)), 
        label = paste("Select", col, " range: "), 
        value= c(start,stop),
        min = start,
        max = stop,
        step = (stop - start) / 10
      )
    )
  }
  else {
    choice_obj <- tags$div(
      'id'=ns(paste(obj_name)),
      pickerInput(
        inputId = ns(paste(obj_name)), 
        label = paste("Select ", col, " value for filtering: "), 
        choices =  c("",unique_items), options = list(`actions-box` = TRUE), 
        selected="",
        multiple = TRUE
      )
    )
  }
  
  choice_cmpd<-list(obj_name, choice_obj)
  return(choice_cmpd)
}

## get_options()
## This function takes in a dataframe and a list
## of columns to filter on and
## returns a list of UI selectInput objects 
## corresponding to each character column
## present in the dataframe. 
get_options <- function(data,cols,ns) {
  col_names <- cols
  choices <- list()
  obj_names <- c()
  for (col in col_names) {
    if (is.character(data[[col]])) {
      choices <- list(choices, get_unique_subs(data, col,ns=ns))
      obj_names <- c(obj_names, paste0(col,"_options"))
    }
    if (is.numeric(data[[col]])) {
      choices <- list(choices, get_unique_subs(data, col,ns=ns, numeric=TRUE))
      obj_names <- c(obj_names, paste0(col,"_options"))
    }
  }
  return(list(choices,obj_names))
} 


## emake_grs()
## Calculates growth rate as a log ratio of 
## cases tomorrow relative to today
emake_grs <- function(data_long) { 
  grs_data <- data_long %>% 
    ## Get growth rate as log ratio of cases tomorrow relative to today
    mutate(gr=log(lead(cases_7day,1)/cases_7day)) %>% # FIXME: NaNs produced. 
    ## If NA or NaN growth rate, then set to 0
    mutate(gr = ifelse(is.na(gr),0, gr))
  return(grs_data)
}

## plot_cases()
## Plot incident cases per reporting date
plot_cases <- function(epi_dat_long,filters=NULL, 
                       type_filters=c("Confirmed","Deceased")) {
  data_long <- epi_dat_long
  filter_choices <- filters
  ggplot(data_long %>% filter(Type %in% type_filters)) +
    geom_line(aes(x=date,y=cases_7day,col=Type)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=type_filters) +  #FIXME: hardcoding -- India-data specific 
    ylab("New cases") +
    xlab("Reporting date") +
    ggtitle(paste0("New cases over time filtered by ", paste(filter_choices,collapse=" & "))) +
    facet_wrap(~Type,scales="free")
}

### Ct Data

## ctmake_long()
## Requires EXACTLY ONE date column
ctmake_long <- function(ct_dat, filters=hash(),ct_cols=c('ct','ct_round'), ct_thresh=35) { # Remarks
  colnames(ct_dat) <-  tolower(colnames(ct_dat))
  date_col <- colnames(ct_dat[, sapply(ct_dat, class) %in% c('Date')]) #requires exactly 1 date column
  colnames(ct_dat)[colnames(ct_dat) == date_col] <- 'date'
  
  ##NOTE: if date column is not formatted as date, requires date column to be called "date"
  if(is.null(date_col)) { date_col <- colnames(ct_dat)[tolower(colnames(ct_dat)) == "date"]} 
  colnames(ct_dat)[colnames(ct_dat) == date_col] <- 'date'
  date_col <- 'date'
  ct_dat[[date_col]] <- as.Date(ct_dat[[date_col]])
  stry_cols <- colnames(ct_dat)[!colnames(ct_dat) %in% tolower(ct_cols)]
  ct_dat_long <- ct_dat %>% pivot_longer(-tolower(stry_cols))

  for (filter in ls(filters)) {
    if (!is.null(filter)) {
      val <- filters[[filter]]
      filter <-  tolower(str_split(filter,'_options')[[1]][1])
      if (is.character(val)) { ct_dat_long <- ct_dat_long %>% filter(!!rlang::sym(filter) == rlang::sym(val)) }
      else if (is.vector(val) & (filter %in% tolower(ct_cols))) { ct_dat_long <- ct_dat_long %>% filter( name == rlang::sym(filter) &
                                                         value >= val[1] & 
                                                         value <= val[2]) } ## FIXME: currently only filters one at a time! 
      else if (is.vector(val)) { ct_dat_long <- ct_dat_long %>% filter(!!rlang::sym(filter) >= val[1] & 
                                                                    !!rlang::sym(filter) >= val[2]) }
    }
  }
  ## Only look at Cts <= thresh
  ct_dat_long <- ct_dat_long %>%
    rename(Ct=value) %>% 
    filter(Ct <= ct_thresh) %>%
    mutate(Ct_round=round(Ct, 0))
  
  return(ct_dat_long)
}

## Plots Cts for each gene over time
plot_ct_raw <- function(ct_dat_long) {
  p_raw <- ct_dat_long %>% ggplot() + 
    geom_point(aes(x=date,y=Ct_round),size=0.25,alpha=0.25) + 
    geom_smooth(aes(x=date,y=Ct_round)) + 
    scale_y_continuous(trans="reverse") +
    scale_x_date(breaks="1 month")+
    ggtitle("Ct values over time") +
    xlab("Date") +
    ylab("Ct") +
    theme_classic()
}

summarize_ct <- function(ct_dat_long) {
  ct_dat_summary <- ct_dat_long %>% 
    group_by(date) %>%
    summarize(skew_ct=moments::skewness(Ct),
              mean_ct=mean(Ct),
              N=n())
  return(ct_dat_summary)
}

plot_ct_mean <- function(ct_dat_long, ct_dat_summary) {
  p_daily_mean <- ct_dat_summary %>% 
    ggplot() + 
    geom_point(aes(x=date,y=mean_ct),size=1,alpha=1) + 
    geom_smooth(aes(x=date,y=mean_ct),span=0.2) + 
    scale_y_continuous(trans="reverse")+
    scale_x_date(breaks="1 month")+
    ggtitle("Daily mean and smoothed mean Cts") +
    xlab("Date") +
    ylab("Mean Ct") +
    theme_classic() 
}

plot_ct_skew <- function(ct_dat_long, ct_dat_summary) {
  p_daily_skew <- ct_dat_summary %>% 
    ggplot() + 
    geom_point(aes(x=date,y=skew_ct),size=1,alpha=1) + 
    geom_smooth(aes(x=date,y=skew_ct),span=0.2) + 
    scale_x_date(breaks="1 month")+
    scale_y_continuous(trans="reverse")+
    ggtitle("Daily skew and smoothed skew of Cts") +
    xlab("Date") +
    ylab("Skewness of Ct values") +
    theme_classic()
}

combine_vis_dat <- function(ct_dat_long, ct_dat_summary, grs_dat, e_dat_long, type_filters=c("Confirmed","Deceased")) {
  grs_dat <- grs_dat 
  epi_data_long <- e_dat_long 
  comb_dat <- ct_dat_summary %>% left_join(grs_dat)
  comb_dat <- comb_dat %>% filter(Type %in% type_filters) %>% filter(N >= 10)
  return(comb_dat)
}

mean_gr_scatter <- function(comb_dat) {
  p_mean_gr_scatter <- comb_dat %>% ggplot(aes(x=mean_ct,y=gr)) +
    geom_hline(yintercept=0,col="black",linetype="dashed",size=0.75) +
    geom_point(size=0.75) +
    stat_cor(method="pearson") +
    scale_x_continuous(trans="reverse") +
    geom_smooth(method="lm",col="black") +
    xlab("Mean Ct") +
    ylab("Growth rate of\n 7-day average cases") +
    theme_classic()
  p_skew_gr_scatter <- comb_dat %>% ggplot(aes(x=skew_ct,y=gr)) +
    geom_hline(yintercept=0,col="black",linetype="dashed",size=0.75) +
    geom_point(size=0.75) +
    stat_cor(method="pearson") +
    geom_smooth(method="lm",col="black") +
    xlab("Skewness of Ct values") +
    ylab("Growth rate of\n 7-day average cases") +
    theme_classic()
}


p_both_scatter <- function(comb_dat) {
  comb_dat %>% ggplot(aes(x=skew_ct,y=mean_ct,col=gr)) +
    geom_point(size=1) +
    stat_cor(method="pearson",label.y=31) +
    scale_color_viridis_c(name="Growth rate") +
    ylab("Mean Ct") +
    xlab("Skewness of Ct values") +
    theme_classic()
}

p_mean_time <- function(comb_dat) {
  dates <- sort(comb_dat$date)
  if (!is.null(dates)) { 
    start_date <- dates[1]
    end_date <- tail(dates,1)
  }
  my_x_axis <- scale_x_date(limits=as.Date(c(start_date-10,end_date+10)),breaks = "14 days")
  p_mean_time <- comb_dat %>%  #%>% filter(gene=="E")
    ggplot(aes(x=date,y=mean_ct)) + 
    geom_point() +
    geom_smooth(span=0.2) +
    scale_y_continuous(trans="reverse",limits=c(31,20),breaks=seq(20,31,by=2)) +
    ylab("Mean N Ct") +
    xlab("Date") +
    theme_classic() +
    my_x_axis
}

p_skew_time <- function(comb_dat) {
  dates <- sort(comb_dat$date)
  if (!is.null(dates)) { 
    start_date <- dates[1]
    end_date <- tail(dates,1)
  }
  my_x_axis <- scale_x_date(limits=as.Date(c(start_date-10,end_date+10)),breaks = "14 days")
  plots <- list()
  iter <- 1
  skew_plot <- ggplot(comb_dat,aes(x=date,y=skew_ct)) + 
  geom_point() +
  geom_smooth(span=0.2) +
  scale_y_continuous(limits=c(-1.5,1.5),breaks=seq(-1.5,1.5,by=0.5)) +
  ylab(paste0("Skewness of Ct"))+
  xlab("Date") +
  theme_classic() +
  my_x_axis +
  theme(axis.text.x=element_text(angle=45,hjust=1))
  
  plots[[iter]] <- skew_plot
  iter <- iter + 1

  plots
}

p_cases_confirmed <- function(data_long, comb_dat) {
  data_long <- data_long
  dates <- sort(comb_dat$date)
  if (!is.null(dates)) { 
    start_date <- dates[1]
    end_date <- tail(dates,1)
  }
  my_x_axis <- scale_x_date(limits=as.Date(c(start_date-10,end_date+10)),breaks = "14 days")
  ggplot(data_long %>% 
           filter(Type %in% c("Confirmed"))%>% ## FIXME: hardcoding "Confirmed". 
           filter(date >= min(comb_dat$date))) +
    geom_line(aes(x=date,y=cases_7day)) +
    theme_classic() +
    ylab("New cases") +
    xlab("Reporting date") +
    my_x_axis
}


p_gr_confirmed <- function(data_grs, comb_dat) {
  dates <- sort(comb_dat$date)
  # Creates x-axis dates using max and min from combined dataset. 
  if (!is.null(dates)) { 
    start_date <- dates[1]
    end_date <- tail(dates,1)
  }
  my_x_axis <- scale_x_date(limits=as.Date(c(start_date-10,end_date+10)),breaks = "14 days")
  ggplot(data_grs %>% 
           filter(date >= min(comb_dat$date)),
         aes(x=date,y=gr)) +
    geom_hline(yintercept=0,linetype="dashed",col="black",size=0.75) +
    geom_line() +
    geom_smooth(span=0.2) +
    theme_classic() +
    ylab("Growth rate") +
    xlab("Reporting date") +
    my_x_axis
}

combine_plots <- function(p_cases_confirmed,p_gr_confirmed,p_mean_time,p_skew_time) {
  p_main1 <- p_cases_confirmed/p_gr_confirmed/p_mean_time/p_skew_time
}

violin_plots <- function(comb_dat, ct_dat_long, data_grs) {
  dates <- sort(ct_dat_long$date)
  if (!is.null(dates)) { 
    start_date <- dates[1]
    end_date <- tail(dates,1)
  }
  make_calendar <- seq(as.Date(start_date-10),as.Date(end_date+10),by="1 day")
  made_calendar <- as_tibble(cbind(make_calendar,MMWRweek(make_calendar))) %>% rename(date=make_calendar)
  epi_cal <- made_calendar %>% group_by(MMWRyear,MMWRweek) %>% mutate(min_date=min(date))
  
  ## Get Epi calendar dates
  epi_cal$EpiWeek <- paste0(epi_cal$MMWRweek,"/",epi_cal$MMWRyear)
  epi_cal$EpiWeek <- factor(epi_cal$EpiWeek,levels=unique(epi_cal$EpiWeek))
  epi_cal <- epi_cal %>% group_by(EpiWeek) %>% mutate(min_date=min(date))
  epi_week_dat <- left_join(ct_dat_long,epi_cal)
  
  plots <- list()

  iter <- 1
  p_violins <- ggplot(epi_week_dat) +
    geom_violin(aes(x=min_date,y=Ct_round,group=min_date),scale="width",width=5,fill="grey70",draw_quantiles=c(0.025,0.5,0.975),alpha=0.5) +
    geom_smooth(aes(x=min_date,y=Ct_round),span=0.2) +
    scale_y_continuous(trans="reverse") +
    scale_x_date(breaks="1 month",limits=as.Date(c(start_date,end_date))) +
    ylab("Ct distribution") +
    xlab("Date grouped by epi week") +
    theme_classic()
  plots[[iter]] <- p_violins #FIXME: list structure is a relic from plotting by gene 
  iter <- iter + 1

  plots
}