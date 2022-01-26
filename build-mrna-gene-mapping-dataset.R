
# information of mRNA genes
# gene names in 3rd column ("Symbol")
#library(lumi)
#OLD#mRNAmapping <- lumi::nuID2RefSeqID(mRNA.highvar.names, lib.mapping='lumiHumanIDMapping', returnAllInfo=TRUE)
 mRNAmapping <- lumi::nuID2RefSeqID(NULL, lib.mapping='lumiHumanIDMapping', returnAllInfo=TRUE) %>% 
   mutate(mRNAname = rownames(.), geneName=Symbol) %>% 
   `rownames<-`(NULL) %>%
   filter(geneName != "")
saveRDS(mRNAmapping[c("mRNAname", "geneName")],"mRNA-gene-mapping.RDS")
