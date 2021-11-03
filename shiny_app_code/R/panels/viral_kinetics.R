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
         bsPopover("overall_prob", "The overall probability of contracting the virus", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="tshift",
           label="Tshift",
           min=pt["tshift",]$lower_bound,
           max=pt["tshift",]$upper_bound,
           value=pt["tshift",]$values
         ),
         bsPopover("tshift", "Time Shift", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="desired_mode",
           label="Mode",
           min=pt["desired_mode",]$lower_bound,
           max=pt["desired_mode",]$upper_bound,
           value=pt["desired_mode",]$values
         ),
         bsPopover("desired_mode", "Mode", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="viral_peak",
           label="Viral Peak",
           min=pt["viral_peak",]$lower_bound,
           max=pt["viral_peak",]$upper_bound,
           value=pt["viral_peak",]$values
         ),
         bsPopover("viral_peak", "Viral Peak", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="obs_sd",
           label="SD (Observed)",
           min=pt["obs_sd",]$lower_bound,
           max=pt["obs_sd",]$upper_bound,
           value=pt["obs_sd",]$values
         ),
         bsPopover("obs_sd", "Observed Standard Deviation", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="sd_mod",
           label="SD (Model)",
           min=pt["sd_mod",]$lower_bound,
           max=pt["sd_mod",]$upper_bound,
           value=pt["sd_mod",]$values
         ),
         bsPopover("sd_mod", "Model Standard Deviation", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="sd_mod_wane",
           label="SD (Waning)",
           min=pt["sd_mod_wane",]$lower_bound,
           max=pt["sd_mod_wane",]$upper_bound,
           value=pt["sd_mod_wane",]$values
         ),
         bsPopover("sd_mod_wane", "Waning model standard deviation", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="true_0",
           label="Ct (True)",
           min=pt["true_0",]$lower_bound,
           max=pt["true_0",]$upper_bound,
           value=pt["true_0",]$values
         ),
         bsPopover("true_0", "True Ct at start of infection.", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="intercept",
           label="Ct (Observed)",
           min=pt["intercept",]$lower_bound,
           max=pt["intercept",]$upper_bound,
           value=pt["intercept",]$values
         ),
         bsPopover("intercept", "Ct measured at start of infection", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="LOD",
           label="LOD",
           min=pt["LOD",]$lower_bound,
           max=pt["LOD",]$upper_bound,
           value=pt["LOD",]$values
         ),
         bsPopover("LOD", "Limit of detection", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="t_switch",
           label="t (Switch)",
           min=pt["t_switch",]$lower_bound,
           max=pt["t_switch",]$upper_bound,
           value=pt["t_switch",]$values
         ),
         bsPopover("t_switch", "Time between peak viral load and beginning of waning (plateau) period", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="level_switch",
           label="Level Switch",
           min=pt["level_switch",]$lower_bound,
           max=pt["level_switch",]$upper_bound,
           value=pt["level_switch",]$values
         ),
         bsPopover("level_switch", "Level switch", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="wane_rate2",
           label="Wane Rate",
           min=pt["wane_rate2",]$lower_bound,
           max=pt["wane_rate2",]$upper_bound,
           value=pt["wane_rate2",]$values
         ),
         bsPopover("wane_rate2", "Wane rate", placement = "right", trigger = "hover",
                   options = list(container="body")),
         sliderInput(
           inputId="prob_detect",
           label="Probability of Detection",
           min=pt["prob_detect",]$lower_bound,
           max=pt["prob_detect",]$upper_bound,
           value=pt["prob_detect",]$values
         ),
         bsPopover("prob_detect", "Probability of detection during the waning period", placement = "right", trigger = "hover",
                   options = list(container="body"))
       ),
       mainPanel(
         plotlyOutput("vk_iplot")
       )
  )
}

vk_tab <- tabPanel("Viral Kinetics", vk_content())
