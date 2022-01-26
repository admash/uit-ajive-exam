# Build CpG <-> Gene dataset in tidy format


## Gene names is listed in the UCSC_RefGene_Name column of the file
## CpGs are listed in Name
#BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19" , force=TRUE)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#data("IlluminaHumanMethylation450kanno.ilmn12.hg19")

annot450k = as.data.frame(minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)) %>%
  dplyr::select(cgName=Name, geneName=UCSC_RefGene_Name, regionType=UCSC_RefGene_Group) 

afile = file("annot.csv", "w")
    cat('"cgName","geneName","regionType"', file=afile)
for(i in 1:nrow(annot450k)){
  if(i %% 1000 == 0) message(i)
  geneNames = strsplit(annot450k[i,"geneName"],";")[[1]]
  regionTypes = strsplit(annot450k[i,"regionType"],";")[[1]]
  if(length(geneNames) == 0) geneNames = ""
  if(length(regionTypes) == 0) regionTypes = ""
  for(j in 1:length(geneNames)){
    cat('"',paste0(c(annot450k[i,"cgName"],geneNames[j], regionTypes[j]), collapse='","'), '"\n', file=afile, sep="")
  }
}
close(afile)

system("sort annot.csv | uniq > cpg-450k-tidy-annot.csv")