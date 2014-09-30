#Dr. Susanne Weller
#22/09/2014
#Association calculations


association.pvalue <- function(a,b) {
  
if(a>=b){
  return(NA)
}  
else{
  if (class(genclinv7.ass[[a]]) == "factor" & class(genclinv7.ass[[b]])=="factor"){
    
    sum_a <- summary(genclinv7.ass[[a]])
    observedfreq_a <- sum_a[2]/length(which(!is.na(genclinv7.ass[[a]])))
    sum_b <-  summary(genclinv7.ass[[b]])
    observedfreq_b<- sum_b[2]/length(which(!is.na(genclinv7.ass[[b]])))
    expectedfreq_1 <- observedfreq_a*observedfreq_b
    expectedfreq_0 <- 1-expectedfreq_1
    table_a_b <- table(genclinv7.ass[[a]],genclinv7.ass[[b]])
    observedcount_1 <- table_a_b[4]
    observedcount_0 <- sum(table_a_b) - table_a_b[4]
    binom <- binom.test(observedcount_1, length(which(!is.na(genclinv7.ass[[a]]))),expectedfreq_1)    
    pvalue <- binom$p.value
    return(pvalue)
}
  if (class(genclinv7.ass[[a]]) == "factor" & class(genclinv7.ass[[b]])=="numeric"){
    
    #Wilcoxon test for two classes
    if (length(levels(genclinv7.ass[[a]]))==2 ){
    wil <-wilcox.test(genclinv7.ass[[b]] ~ genclinv7.ass[[a]], data=genclinv7.ass)
    pvalue <- wil$p.value
    return(pvalue)
    }
    #Kruskal Wallis rank sum test for more than one class
    else{
    krus <- kruskal.test(genclinv7.ass[[b]] ~ genclinv7.ass[[a]], data=genclinv7.ass)
    pvalue <- krus$p.value
    return(pvalue)
    }
  }
  if (class(genclinv7.ass[[a]]) == "numeric" & class(genclinv7.ass[[b]])=="factor"){
    #Wilcoxon test for two classes
    if (length(levels(genclinv7.ass[[b]]))==2 ){
      wil <-wilcox.test(genclinv7.ass[[a]] ~ genclinv7.ass[[b]], data=genclinv7.ass)
      pvalue <- wil$p.value
      return(pvalue)
    }
    #Kruskal Wallis rank sum test for more than one class
    else{
      krus <- kruskal.test(genclinv7.ass[[a]] ~ genclinv7.ass[[b]], data=genclinv7.ass)
      pvalue <- krus$p.value
      return(pvalue)
    }
    
  }
  if (class(genclinv7.ass[[a]]) == "numeric" & class(genclinv7.ass[[b]])=="numeric"){
    #Spearman correlation for continuous variables
    cor <-cor.test(genclinv7.ass[[a]],genclinv7.ass[[b]], data=genclinv7.ass, method="spearman")
    pvalue <- cor$p.value
    return(pvalue)
  }
}
}