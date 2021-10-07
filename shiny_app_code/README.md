This README.md currently functions as a notebook recording the quirks, areas for improvement, and unimplemented features of the virosolver shiny web app. 

server.R 
---------------------------

plots.R
---------------------------

global.R
---------------------------

UI.R
---------------------------

helpers/plots.R 
---------------------------

helpers/helpers.R
--------------------------

helpers/mcmc_helpers.R
--------------------------

helpers/vk_helpers.R
--------------------------

panels/data_visualization.R
--------------------------
1. Epi tab filtering is not implemented. It should be simple enough to follow the same process for filtering as used for  the Ct values tab. 
2. Plot generation is exceedingly slow; it takes at least 1-2 minutes for plots to load. 
3. Filtering in Ct tab does not have any safegaurds. It is very possible to fillter too much and break the plots due to lack of data. 
4. Vi_plot and ps_time objects in dv_plots() function are of type list, this is why they are parsed in for loops prior to inclusion in the overall plots object.5. In vis_content, the empty div titled placeholder is necessary to allow for filtering input options to be added to / removed from the DOM later on.
6. Vis_content is currently the only modularized code for a panel. Others can follow this format in order to de-clutter the server.R file.  

panels/run_model.R
---------------------------

panels/viral_kinetics.R
---------------------------

panels/primary_output.R
---------------------------


