rootdir <- dirname(rstudioapi::getSourceEditorContext()$path)

source(paste0(rootdir,"/UI.R"))
source(paste0(rootdir,"/server.R"))

shinyApp(
  ui=viro_UI,
  server=viro_server
)