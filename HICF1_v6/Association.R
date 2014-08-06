#Dr. Susanne Weller
#23/07/2014
#Association calculations


association.pvalue <- function(a,b) {
  
if(a>=b){
  return(NA)
}  
else{
  if (class(genclinv6[[a]]) == "factor" & class(genclinv6[[b]])=="factor"){
    
    sum_a <- summary(genclinv6[[a]])
    observedfreq_a <- sum_a[2]/length(which(!is.na(genclinv6[[a]])))
    sum_b <-  summary(genclinv6[[b]])
    observedfreq_b<- sum_b[2]/length(which(!is.na(genclinv6[[b]])))
    expectedfreq_1 <- observedfreq_a*observedfreq_b
    expectedfreq_0 <- 1-expectedfreq_1
    table_a_b <- table(genclinv6[[a]],genclinv6[[b]])
    observedcount_1 <- table_a_b[4]
    observedcount_0 <- sum(table_a_b) - table_a_b[4]
    binom <- binom.test(observedcount_1, length(which(!is.na(genclinv6[[a]]))),expectedfreq_1)    
    pvalue <- binom$p.value
    return(pvalue)
}
  if (class(genclinv6[[a]]) == "factor" & class(genclinv6[[b]])=="numeric"){
    
    #Wilcoxon test for two classes
    if (length(levels(genclinv6[[a]]))==2 ){
    wil <-wilcox.test(genclinv6[[b]] ~ genclinv6[[a]], data=genclinv6)
    pvalue <- wil$p.value
    return(pvalue)
    }
    #Kruskal Wallis rank sum test for more than one class
    else{
    krus <- kruskal.test(genclinv6[[b]] ~ genclinv6[[a]], data=genclinv6)
    pvalue <- krus$p.value
    return(pvalue)
    }
  }
  if (class(genclinv6[[a]]) == "numeric" & class(genclinv6[[b]])=="factor"){
    #Wilcoxon test for two classes
    if (length(levels(genclinv6[[b]]))==2 ){
      wil <-wilcox.test(genclinv6[[a]] ~ genclinv6[[b]], data=genclinv6)
      pvalue <- wil$p.value
      return(pvalue)
    }
    #Kruskal Wallis rank sum test for more than one class
    else{
      krus <- kruskal.test(genclinv6[[a]] ~ genclinv6[[b]], data=genclinv6)
      pvalue <- krus$p.value
      return(pvalue)
    }
    
  }
  if (class(genclinv6[[a]]) == "numeric" & class(genclinv6[[b]])=="numeric"){
    #Spearman correlation for continuous variables
    cor <-cor.test(genclinv6[[a]],genclinv6[[b]], data=genclinv6, method="spearman")
    pvalue <- cor$p.value
    return(pvalue)
  }
}
}