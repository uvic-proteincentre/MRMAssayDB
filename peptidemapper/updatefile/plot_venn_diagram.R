
library(VennDiagram)
# library(gplots)
plotVennDiagram <- function(inputfilename,venn_pep_fname="venn_pep.png",venn_prot_fname="venn_prot.png"){
  m<-read.csv(inputfilename,check.names=F,header = T,sep="\t")
  # dim(m)
  # colnames(m)
  #m1=m[,c("PeptideTracker ID","Passel ID","SRMAtlas ID","Cptac ID","Panoramaweb ID")]
  m1=m[,c("PeptideTracker ID","Passel ID","SRMAtlas ID","Cptac ID","Panoramaweb ID","UniProtKB Accession")]
  #write.table(m1,"ReportBook_mother_file_v2.csv",sep = "\t",row.names = F)
  colnames(m1)=c("PeptideTracker","PASSEL","SRMAtlas","CPTAC","PanoramaWeb","uniprot")
  m1=as.matrix(m1)
  m1[,"uniprot"]=sapply(m1[,"uniprot"],function(x) {strsplit(x = x,split = "-",fixed = T)[[1]][1]})
  
  venn_pep <- draw.quintuple.venn(
    area1 = length(which(!is.na(m1[,1]))),
    area2 = length(which(!is.na(m1[,2]))),
    area3 = length(which(!is.na(m1[,3]))),
    area4 = length(which(!is.na(m1[,4]))),
    area5 = length(which(!is.na(m1[,5]))),
    n12 = length(which(!is.na(m1[,1])&!is.na(m1[,2]))),
    n13 = length(which(!is.na(m1[,1])&!is.na(m1[,3]))),
    n14 = length(which(!is.na(m1[,1])&!is.na(m1[,4]))),
    n15 = length(which(!is.na(m1[,1])&!is.na(m1[,5]))),
    n23 = length(which(!is.na(m1[,2])&!is.na(m1[,3]))),
    n24 = length(which(!is.na(m1[,2])&!is.na(m1[,4]))),
    n25 = length(which(!is.na(m1[,2])&!is.na(m1[,5]))),
    n34 = length(which(!is.na(m1[,3])&!is.na(m1[,4]))),
    n35 = length(which(!is.na(m1[,3])&!is.na(m1[,5]))),
    n45 = length(which(!is.na(m1[,4])&!is.na(m1[,5]))),
    n123 = length(which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,3]))),
    n124 = length(which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,4]))),
    n125 = length(which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,5]))),
    n134 = length(which(!is.na(m1[,1])&!is.na(m1[,3])&!is.na(m1[,4]))),
    n135 = length(which(!is.na(m1[,1])&!is.na(m1[,3])&!is.na(m1[,5]))),
    n145 = length(which(!is.na(m1[,1])&!is.na(m1[,4])&!is.na(m1[,5]))),
    n234 = length(which(!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,4]))),
    n235 = length(which(!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,5]))),
    n245 = length(which(!is.na(m1[,2])&!is.na(m1[,4])&!is.na(m1[,5]))),
    n345 = length(which(!is.na(m1[,3])&!is.na(m1[,4])&!is.na(m1[,5]))),
    
    n1234 = length(which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,4]))),
    n1235 = length(which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,5]))),
    n1245 = length(which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,4])&!is.na(m1[,5]))),
    n1345 = length(which(!is.na(m1[,1])&!is.na(m1[,3])&!is.na(m1[,4])&!is.na(m1[,5]))),
    n2345 = length(which(!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,4])&!is.na(m1[,5]))),
    
    n12345 = length(which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,4])&!is.na(m1[,5]))),
    category = colnames(m1)[1:5],
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.cex = 1, cat.dist =c(0.2,0.2,0.2,0.2,0.22),
    cat.pos = c(0, 320, 215, 145, 40),
    margin = 0.1,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
            1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    ind = TRUE
  );
  
  
  venn_prot <- draw.quintuple.venn(
    area1 = length(unique(m1[which(!is.na(m1[,1])),"uniprot"])),
    area2 = length(unique(m1[which(!is.na(m1[,2])),"uniprot"])),
    area3 = length(unique(m1[which(!is.na(m1[,3])),"uniprot"])),
    area4 = length(unique(m1[which(!is.na(m1[,4])),"uniprot"])),
    area5 = length(unique(m1[which(!is.na(m1[,5])),"uniprot"])),
    n12 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,2])),"uniprot"])),
    n13 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,3])),"uniprot"])),
    n14 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,4])),"uniprot"])),
    n15 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,5])),"uniprot"])),
    n23 = length(unique(m1[which(!is.na(m1[,2])&!is.na(m1[,3])),"uniprot"])),
    n24 = length(unique(m1[which(!is.na(m1[,2])&!is.na(m1[,4])),"uniprot"])),
    n25 = length(unique(m1[which(!is.na(m1[,2])&!is.na(m1[,5])),"uniprot"])),
    n34 = length(unique(m1[which(!is.na(m1[,3])&!is.na(m1[,4])),"uniprot"])),
    n35 = length(unique(m1[which(!is.na(m1[,3])&!is.na(m1[,5])),"uniprot"])),
    n45 = length(unique(m1[which(!is.na(m1[,4])&!is.na(m1[,5])),"uniprot"])),
    n123 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,3])),"uniprot"])),
    n124 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,4])),"uniprot"])),
    n125 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,5])),"uniprot"])),
    n134 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,3])&!is.na(m1[,4])),"uniprot"])),
    n135 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,3])&!is.na(m1[,5])),"uniprot"])),
    n145 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,4])&!is.na(m1[,5])),"uniprot"])),
    n234 = length(unique(m1[which(!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,4])),"uniprot"])),
    n235 = length(unique(m1[which(!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,5])),"uniprot"])),
    n245 = length(unique(m1[which(!is.na(m1[,2])&!is.na(m1[,4])&!is.na(m1[,5])),"uniprot"])),
    n345 = length(unique(m1[which(!is.na(m1[,3])&!is.na(m1[,4])&!is.na(m1[,5])),"uniprot"])),
    
    n1234 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,4])),"uniprot"])),
    n1235 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,5])),"uniprot"])),
    n1245 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,4])&!is.na(m1[,5])),"uniprot"])),
    n1345 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,3])&!is.na(m1[,4])&!is.na(m1[,5])),"uniprot"])),
    n2345 = length(unique(m1[which(!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,4])&!is.na(m1[,5])),"uniprot"])),
    
    n12345 = length(unique(m1[which(!is.na(m1[,1])&!is.na(m1[,2])&!is.na(m1[,3])&!is.na(m1[,4])&!is.na(m1[,5])),"uniprot"])),
    category = colnames(m1)[1:5],
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.cex = 1, cat.dist =c(0.2,0.2,0.2,0.2,0.22),
    cat.pos = c(0, 320, 215, 145, 40),
    margin = 0.1,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
            1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    ind = TRUE
  );
  
  png(venn_pep_fname,res = 175, bg = "#EEEEEE",units = "px",width = 1000,height = 1000)
  grid.draw(venn_pep)
  dev.off()
  
  png(venn_prot_fname,res = 175, bg = "#EEEEEE",units = "px",width = 1000,height = 1000)
  grid.draw(venn_prot)
  dev.off()
  
}
