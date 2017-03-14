BH_correction_Pvalue<-function(P_values,FDR=0.10)
{
  #Benjamini Hochberg
  # One good technique for controlling the false discovery rate was briefly 
  # mentioned by Simes (1986) and developed in detail by Benjamini and Hochberg (1995). 
  # Put the individual P values in order, from smallest to largest. The smallest P value 
  # has a rank of i=1, then next smallest has i=2, etc. Compare each individual P value 
  # to its Benjamini-Hochberg critical value, (i/m)Q, where i is the rank, m is the total 
  # number of tests, and Q is the false discovery rate you choose. The largest P value 
  # that has P<(i/m)Q is significant, and all of the P values smaller than it are also 
  # significant, even the ones that aren't less than their Benjamini-Hochberg critical value.
  
  Pvals_sorted<-sort(P_values)
  ranks<-1:length(Pvals_sorted)
  total_test<-length(Pvals_sorted)
  Q<-FDR
  BH_critical_val<-(ranks/total_test)*Q
  BH_crit_signif<-which(Pvals_sorted<BH_critical_val)
  
  if(length(BH_crit_signif)==0){warning("No P-values significant under Benjamini-Hochberg FDR correction")}
  
  if(length(BH_crit_signif)!=0)
  {
    largest_BH_crit_value<-BH_crit_signif[length(BH_crit_signif)]
    BH_crit_Pval<-Pvals_sorted[which(BH_crit_signif==largest_BH_crit_value)]
    names(BH_crit_Pval)<-NULL
    message(paste("The P-value under which you can assume significance at FDR=",Q,"is",BH_crit_Pval))
    return(BH_crit_Pval)
  }
  
}

Enrichment_Test_GO<-function(Gene_list,Gene_Background,GO_Annot_table,geneid_col=1,annot_col=1,FDR=0.05,statistic="HyperGeometric",output_file="GO_Analysis")
{
  #Test list setup
  GO_slim_table<-GO_Annot_table #variable renaming convenience
  #------ This code section unfortunately only picks single matches, 
  #------ a gene may map to more than one ontology term, REMEMBER!  
  
  #list_ixmatch<-match(Gene_list,GO_slim_table[,geneid_col]) #find all indices in gene list which match the GO slim table, you need it to print summaries
  #no_hits<-which(is.na(list_ixmatch)) # indices in the list where hits did not occur
  #hits<-list_ixmatch[!is.na(list_ixmatch)] # indices of the rows where it matched the GO table
  #message(paste(paste(Gene_list[no_hits],sep = ",",collapse = ","),"from the test list had no Gene Ontology mappings"))
  hits<-unlist(sapply(Gene_list,function(x){grep(pattern = x,x = GO_slim_table[,geneid_col],fixed = TRUE)}))
  hit_terms<-GO_slim_table[hits,annot_col] #retrieving the annotations for the successful GO hits 
  hit_term_counts<-table(hit_terms) # creating a frequency summary of annotations
  
  #------ This code section ALSO unfortunately only picks single matches, 
  #------ a gene may map to more than one ontology term, REMEMBER!  
  #Background setup
  #bg_ixmatch<-match(Gene_Background,GO_slim_table[,geneid_col]) # indices of GOslimtable where Gene_Background mapped successfully
  #bg_no_hits<-which(is.na(bg_ixmatch)) # retrieve those indices for Gene Background where no matches were found in the GO table list
  #bg_hits<-which(!is.na(bg_ixmatch)) #retrieve those indices for the Gene Background which match the GOtable
  #message(paste(paste(Gene_Background[bg_no_hits],sep = ",",collapse = ","),"from the background had no Gene Ontology mappings"))
  #bg_ixmatch<-bg_ixmatch[!is.na(bg_ixmatch)] #remove the NA's from the GOTable match indices
  bg_hits_ix<-unlist(sapply(Gene_Background,function(x){grep(pattern=x,x=GO_slim_table[,geneid_col],fixed=TRUE)}))
  bg_hit_terms<-GO_slim_table[bg_hits_ix,annot_col] # subset the GO annotations for the hit from the background list
  bg_term_counts<-table(bg_hit_terms) # count the number of times the annotation comes up  
  
  #Determine for the test list, how many are in the category and how many aren't
  #Do the same for the background list, how many in that category and how many aren't
  #construct a contingency matrix
  #do the test sequentially for each category
  #return a list for each category tested
  relevant_bg_term_counts<-bg_term_counts[match(names(hit_term_counts),names(bg_term_counts))]
  
  
  if(statistic=="Fisher")
  {
    #Fisher Exact Test
    message("Applying Fisher Exact test for Enrichment Analysis")
    in_gp_gl<-hit_term_counts
    nin_gp_gl<-sapply(hit_term_counts,function(x){sum(hit_term_counts)-x})
    in_gp_bg<-relevant_bg_term_counts
    nin_gp_bg<-sapply(relevant_bg_term_counts,function(x){sum(relevant_bg_term_counts)-x})
    contingency_vector_Fisher<-rbind(in_gp_gl,nin_gp_gl,in_gp_bg,nin_gp_bg)
    Fisher_test_list<-apply(contingency_vector_Fisher,2,function(x){fisher.test(matrix(x,nrow=2),alternative = "greater")})
    Fisher_Pvals<-unlist(lapply(Fisher_test_list,function(x){x$p.value}))
    BH_cut_Pval<-BH_correction_Pvalue(Fisher_Pvals,FDR=FDR)
    
    if(length(BH_cut_Pval)==0)
    {
      # If none of the P-values are significant after correction, there is no point
      # in making enrichment plots. So just return the list of the result of statistical
      # tests.
      invisible(Fisher_test_list)
      warning("No P-values found significant after correction")
    }
    
    Cleared_Pval_Categories_index<-which(Fisher_Pvals<=BH_cut_Pval)
    
    if(length(Cleared_Pval_Categories_index)!=0)
    {
      Stat_signif_GO_cat_contingency<-contingency_vector_Fisher[,Cleared_Pval_Categories_index]
      pdf(paste("GO_enrichment",output_file,".pdf",sep="_",collapse=""),width = 5)
      library(RColorBrewer)
      cols<-brewer.pal(10,name = "Set3")
      if(is.vector(Stat_signif_GO_cat_contingency))
      {
        Stat_signif_GO_cat_contingency<-matrix(Stat_signif_GO_cat_contingency,4,1,byrow=TRUE)
        colnames(Stat_signif_GO_cat_contingency)<-names(Cleared_Pval_Categories_index)
      }
      bp<-barplot(Stat_signif_GO_cat_contingency[c(1,3),],beside = FALSE,col=c(cols[4],cols[7]),log="x",horiz = TRUE,main="GO Enrichment",las=1,cex.names = 0.01)
      text(x=Stat_signif_GO_cat_contingency[c(1),], y=bp, labels=colnames(Stat_signif_GO_cat_contingency), pos=4, xpd=1,cex=0.4,offset=0.1)
      dev.off()
    }else(warning("No GO category cleared the corrected P-value cut-off"))
  }
  
  if(statistic=="Binomial")
  {
    #Binomial Test
    message("Applying Binomial test for Enrichment Analysis")
    in_gp_gl<-hit_term_counts
    total_in_gp_gl<-rep(sum(hit_term_counts),length(in_gp_gl))
    in_gp_bg<-relevant_bg_term_counts
    total_in_gp_bg<-rep(sum(relevant_bg_term_counts),length(in_gp_bg))
    contingency_vector_Binomial<-rbind(in_gp_gl,total_in_gp_gl,in_gp_bg,total_in_gp_bg)
    Binomial_test_list<-apply(contingency_vector_Binomial,2,function(x){binom.test(x = x[1],n = x[2],p = x[3]/x[4],alternative = "greater")})  
    Binomial_Pvals<-unlist(lapply(Binomial_test_list,function(x){x$p.value}))
    BH_cut_Pval<-BH_correction_Pvalue(Binomial_Pvals,FDR=FDR)
    
    if(length(BH_cut_Pval)==0)
    {
      # If none of the P-values are significant after correction, there is no point
      # in making enrichment plots. So just return the list of the result of statistical
      # tests.
      invisible(Binomial_test_list)
      warning("No P-values found significant after correction")
    }
    
    Cleared_Pval_Categories_index<-which(Binomial_Pvals<=BH_cut_Pval)
    
    if(length(Cleared_Pval_Categories_index)!=0)
    {
      Stat_signif_GO_cat_contingency<-contingency_vector_Binomial[,Cleared_Pval_Categories_index]
      pdf(paste("GO_enrichment",output_file,".pdf",sep="_",collapse=""),width = 5)
      library(RColorBrewer)
      cols<-brewer.pal(10,name = "Set3")
      if(is.vector(Stat_signif_GO_cat_contingency))
      {
        Stat_signif_GO_cat_contingency<-matrix(Stat_signif_GO_cat_contingency,4,1,byrow=TRUE)
        colnames(Stat_signif_GO_cat_contingency)<-names(Cleared_Pval_Categories_index)
      }
      bp<-barplot(Stat_signif_GO_cat_contingency[c(1,3),],beside = FALSE,col=c(cols[4],cols[7]),log="x",horiz = TRUE,main="GO Enrichment",las=1,cex.names = 0.01)
      text(x=Stat_signif_GO_cat_contingency[c(1),], y=bp, labels=colnames(Stat_signif_GO_cat_contingency), pos=4, xpd=1,cex=0.4,offset=0.1)
      dev.off()
    }else(warning("No GO category cleared the corrected P-value cut-off"))
  }
  
  if(statistic=="HyperGeometric")
  {
    #Hypergeometric test
    message("Applying Hypergeometric test for Enrichment Analysis")
    in_gp_gl<-hit_term_counts
    in_gp_bg<-relevant_bg_term_counts
    total_in_gp_bg<-rep(sum(relevant_bg_term_counts),length(in_gp_bg))
    total_in_gp_gl<-rep(sum(hit_term_counts),length(in_gp_gl))
    contingency_vector_HyperGM<-rbind(in_gp_gl,in_gp_bg,total_in_gp_bg,total_in_gp_gl)
    HyperGeoTest_list<-apply(contingency_vector_HyperGM,2,function(x){min(1-cumsum(dhyper(0:(x[1]-1),x[2],x[3],x[4])))})  
    HyperGeo_Pvals<-HyperGeoTest_list
    BH_cut_Pval<-BH_correction_Pvalue(HyperGeo_Pvals,FDR=FDR)
    
    if(length(BH_cut_Pval)==0)
    {
      # If none of the P-values are significant after correction, there is no point
      # in making enrichment plots. So just return the list of the result of statistical
      # tests.
      invisible(HyperGeoTest_list)
      warning("No P-values found significant after correction")
    }
    
    Cleared_Pval_Categories_index<-which(HyperGeo_Pvals<=BH_cut_Pval)
    
    if(length(Cleared_Pval_Categories_index)!=0)
    {
      Stat_signif_GO_cat_contingency<-contingency_vector_HyperGM[,Cleared_Pval_Categories_index]
      pdf(paste("GO_enrichment",output_file,".pdf",sep="_",collapse=""),width = 5)
      library(RColorBrewer)
      cols<-brewer.pal(10,name = "Set3")
      if(is.vector(Stat_signif_GO_cat_contingency))
      {
        Stat_signif_GO_cat_contingency<-matrix(Stat_signif_GO_cat_contingency,4,1,byrow=TRUE)
        colnames(Stat_signif_GO_cat_contingency)<-names(Cleared_Pval_Categories_index)
      }
      bp<-barplot(Stat_signif_GO_cat_contingency[c(1,2),],beside = FALSE,col=c(cols[4],cols[7]),log="x",horiz = TRUE,main="GO Enrichment",las=1,cex.names = 0.01)
      text(x=Stat_signif_GO_cat_contingency[c(1),], y=bp, labels=colnames(Stat_signif_GO_cat_contingency), pos=4, xpd=1,cex=0.4,offset=0.1)
      dev.off()
    }else(warning("No GO category cleared the corrected P-value cut-off"))
  }
  
}