# https://github.com/pmartinezarbizu/pairwiseAdonis
#
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=999)
{
    
    co <- combn(unique(as.character(factors)),2)
    pairs <- c()
    Df <- c()
    SumsOfSqs <- c()
    F.Model <- c()
    R2 <- c()
    p.value <- c()
    
    
    for(elem in 1:ncol(co)){
        if(inherits(x, 'dist')){
            x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                            factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
        }
        
        else  (
            if (sim.function == 'daisy'){
                x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
            } 
            else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
        )
        
        x2 = data.frame(Fac = factors[factors %in% c(co[1,elem],co[2,elem])])
        
        ad <- adonis2(x1 ~ Fac, data = x2,
                      permutations = perm);
        pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
        Df <- c(Df,ad$Df[1])
        SumsOfSqs <- c(SumsOfSqs,ad$SumOfSqs[1])
        F.Model <- c(F.Model,ad$F[1]);
        R2 <- c(R2,ad$R2[1]);
        p.value <- c(p.value,ad$`Pr(>F)`[1])
    }
    p.adjusted <- p.adjust(p.value,method=p.adjust.m)
    
    sig = c(rep('',length(p.adjusted)))
    sig[p.adjusted <= 0.05] <-'.'
    sig[p.adjusted <= 0.01] <-'*'
    sig[p.adjusted <= 0.001] <-'**'
    sig[p.adjusted <= 0.0001] <-'***'
    pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
    
    if(!is.null(reduce)){
        pairw.res <- subset (pairw.res, grepl(reduce,pairs))
        pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
        
        sig = c(rep('',length(pairw.res$p.adjusted)))
        sig[pairw.res$p.adjusted <= 0.1] <-'.'
        sig[pairw.res$p.adjusted <= 0.05] <-'*'
        sig[pairw.res$p.adjusted <= 0.01] <-'**'
        sig[pairw.res$p.adjusted <= 0.001] <-'***'
        pairw.res <- data.frame(pairw.res[,1:7],sig)
    }
    class(pairw.res) <- c("pwadonis", "data.frame")
    return(pairw.res)
} 


### Method summary
summary.pwadonis = function(object, ...) {
    cat("Result of pairwise.adonis:\n")
    cat("\n")
    print(object, ...)
    cat("\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}