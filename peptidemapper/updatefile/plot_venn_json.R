library(VennDiagram)
library(rjson)
library(gplots)


#inputfilename="Katherine_basic_search_venn.json"

plotVennDiagram_json <- function(inputfilename,output_prot,output_pep){
  m1<-fromJSON(file=inputfilename)
  m=m1$prodataseries
  venn_plot_prot=venn.diagram(  m, alpha = 0.3, filename = NULL,
                         fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(m)],
                         cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(m)],
                         cat.cex = 1, cat.dist =c(0.2,0.2,0.2,0.2,0.23)[1:length(m)],
                         margin = 0.1)
  m=m1$pepseqdataseries
  venn_plot_pep=venn.diagram(  m, alpha = 0.3, filename = NULL,
                           fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(m)],
                           cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(m)],
                           cat.cex = 1, cat.dist =c(0.2,0.2,0.2,0.2,0.23)[1:length(m)],
                           margin = 0.1)
  png(output_prot,res = 175, bg = "#EEEEEE",units = "px",width = 1000,height = 1000)
  grid.draw(venn_plot_prot)
  dev.off()
  png(output_pep,res = 175, bg = "#EEEEEE",units = "px",width = 1000,height = 1000)
  grid.draw(venn_plot_pep)
  dev.off()
  
}