#source("/Users/kotah/Desktop/virosolver/ShinySolver/ShinySolver/R/helpers/vk_helpers.R")

#set shorter names
pt <- example_gp_partab 
row.names(pt) <- c(example_gp_partab$names[1:19],c(1:(384-19)))

vk_content <- function() {
  content <- fluidRow(
       sidebarPanel(
         fileInput("vk_data","Upload Ct Value Data (CSV)", accept=".csv"),
         sliderInput(
           inputId="overall_prob",
           label="Overall Probability",
           min=pt["overall_prob",]$lower_bound,
           max=pt["overall_prob",]$upper_bound,
           value=pt["overall_prob",]$values
         ),
         bsPopover("overall_prob", "The overall probability of contracting the disease (FIXME: double check this w/ vignettes)", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="tshift",
           label="Tshift",
           min=pt["tshift",]$lower_bound,
           max=pt["tshift",]$upper_bound,
           value=pt["tshift",]$values
         ),
         sliderInput(
           inputId="desired_mode",
           label="desired mode",
           min=pt["desired_mode",]$lower_bound,
           max=pt["desired_mode",]$upper_bound,
           value=pt["desired_mode",]$values
         ),
         sliderInput(
           inputId="viral_peak",
           label="viral peak",
           min=pt["viral_peak",]$lower_bound,
           max=pt["viral_peak",]$upper_bound,
           value=pt["viral_peak",]$values
         ),
         sliderInput(
           inputId="obs_sd",
           label="observed standard deviation",
           min=pt["obs_sd",]$lower_bound,
           max=pt["obs_sd",]$upper_bound,
           value=pt["obs_sd",]$values
         ),
         sliderInput(
           inputId="sd_mod",
           label="model standard deviation",
           min=pt["sd_mod",]$lower_bound,
           max=pt["sd_mod",]$upper_bound,
           value=pt["sd_mod",]$values
         ),
         sliderInput(
           inputId="sd_mod_wane",
           label="waning model standard deviation",
           min=pt["sd_mod_wane",]$lower_bound,
           max=pt["sd_mod_wane",]$upper_bound,
           value=pt["sd_mod_wane",]$values
         ),
         sliderInput(
           inputId="true_0",
           label="true Ct at start of infection",
           min=pt["true_0",]$lower_bound,
           max=pt["true_0",]$upper_bound,
           value=pt["true_0",]$values
         ),
         sliderInput(
           inputId="intercept",
           label="Ct captured by assay at start of infection",
           min=pt["intercept",]$lower_bound,
           max=pt["intercept",]$upper_bound,
           value=pt["intercept",]$values
         ),
         sliderInput(
           inputId="LOD",
           label="limit of detection",
           min=pt["LOD",]$lower_bound,
           max=pt["LOD",]$upper_bound,
           value=pt["LOD",]$values
         ),
         sliderInput(
           inputId="t_switch",
           label="time between peak viral load and beginning of waning (plateau) period",
           min=pt["t_switch",]$lower_bound,
           max=pt["t_switch",]$upper_bound,
           value=pt["t_switch",]$values
         ),
         sliderInput(
           inputId="level_switch",
           label="level switch --- ",
           min=pt["level_switch",]$lower_bound,
           max=pt["level_switch",]$upper_bound,
           value=pt["level_switch",]$values
         ),
         sliderInput(
           inputId="wane_rate2",
           label="wane rate",
           min=pt["wane_rate2",]$lower_bound,
           max=pt["wane_rate2",]$upper_bound,
           value=pt["wane_rate2",]$values
         ),
         sliderInput(
           inputId="prob_detect",
           label="probability detection waning period",
           min=pt["prob_detect",]$lower_bound,
           max=pt["prob_detect",]$upper_bound,
           value=pt["prob_detect",]$values
         )#,
         #width=2
       ),
       mainPanel(
         plotlyOutput("vk_iplot")
         #width=10
       )
  )
}

vk_tab <- tabPanel("Viral Kinetics", vk_content())
