args <- commandArgs(TRUE)
inputfilename=args[1]
output_prot=args[2]
output_pep=args[3]
source("json_plot_venn_to_pallab.R")
plotVennDiagram_json(inputfilename,output_prot,output_pep)
