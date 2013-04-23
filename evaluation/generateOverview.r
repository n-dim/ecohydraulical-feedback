generateOverview <- function(path){
  grDevices::pdf.options(useDingbats = FALSE); 
  require(knitr)
  oldwd <- getwd()
  setwd('./knitr/overview/')
  file <- knit2pdf('./overview.Rnw', encoding='UTF-8')
  setwd(oldwd)
  file.copy(file.path('./knitr/overview/', file), path)
  file.copy(file.path('./knitr/overview/figure'), recursive=T, path)
}


filename <- 'overview.tex'
generateOverview(paste0('../../TestOutput/'))