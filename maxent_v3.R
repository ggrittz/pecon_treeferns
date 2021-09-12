library(ecospat)
library(dismo)
library(blockCV)
library(dplyr)
library(ENMeval)
library(ROCR)

####Climate database to be modeled
cl='chelsa' #chelsa or 'wc' for worldclim

####Dataframe with coordinates
data = read.csv('banco de dados/data_ref.csv',sep=';')

sp1=unique(data$sp) #species names to loop

{
  res.eval=list()
  for (j in sp1) {
    
    print(paste('Modelando',j))
    
    ####Loading the registers of the species to be modeled
    occ=data.frame(filter(data,sp==j)[,c('x','y')],sp=1)
    occ=na.omit(unique(data.frame(occ,extract(env$current[[1]],occ[,c('x','y')])))) [,c('x','y','sp')]
    
    ####Spatial thinning
    if(nrow(occ)>15){
      occ=ecospat.occ.desaggregation(occ,0.2)
    }
    ####Pseudoabsences
    bg=data.frame(randomPoints(env$current[[1]],10000),sp=0)
    
    
    ####Spatially partitioning the registres in blocks
    fold=ENMeval::get.block(occ,bg)
    occ=data.frame(occ,fold=fold$occ.grp)
    bg=data.frame(bg,fold=fold$bg.grp)
    
    ####Final dataframe
    md=rbind(occ,bg)
    md=data.frame(md,extract(env$current,md[,c('x','y')]))
    
    ####Hypertunning parameters grid
    grid=expand.grid(beta=seq(0.5,5,by=.5),
                     fc=c('L','H', 'LQ','LQP','LQH','LQHP','LQHPT'),
                     AIC=NA,auc=NA,boyce=NA)
    
   ####Running MaxEnt 
    res=list()
    for (i in 1:nrow(grid)) {
      for (k in 1:5) {
        
        train=md[which(md$fold!=k),]
        test=md[which(md$fold==k),]
        
        md1=maxent(train[,5:ncol(train)],
                   train$sp,
                   args=prepPara(betamultiplier = grid[i,'beta'],
                                 userfeatures =grid[i,'fc']))
        
        if(k!=5){
          
        pred=data.frame('pred'=predict(md1,test,args=c("outputformat=cloglog")),'sp'=na.omit(test)$sp)
        res[k]=list(pred)
        
        } else {
          
          r=predict(md1,env$current,args=c("outputformat=raw"))
        }
      }
      
      ####Evaluation
      pred=do.call('rbind',res)
      AIC=calc.aicc(get.params(md1),occ[,1:2],r)[1,1]
      
      pred2=prediction(pred$pred,pred$sp)
      auc=performance(pred2,'auc')@y.values[[1]]
      boyce=ecospat.boyce(pred[,'pred'],
                          pred[which(pred$sp==1),'pred'],PEplot = FALSE)$Spearman.cor
      
      grid[i,'AIC']=AIC
      grid[i,'auc']=auc
      grid[i,'boyce']=boyce
      print(i)
    }
    
    #Selecting best fit model
    best=grid %>%
      arrange((AIC)) %>%
      head(1)
      
    res.eval[j]=list(best)
    
    ####Running again, but only the best fit model new

      md.final=maxent(md[,5:ncol(md)],
                      md$sp,
                      args=prepPara(betamultiplier = best[1,'beta'],
                                    userfeatures =best[1,'fc']))
    
    
    
    ####Projecting the model in each future environmental layers (RCP26, RCP85)
      for (c in names(env)) {
        proj=env[[c]]
        names(proj)=paste('bio',bios,sep = '')
        proj=predict(md.final,proj,args=c("outputformat=cloglog"))
        writeRaster(proj, paste('mapas/v3/',cl,'/',j,'_',c,'.tif',sep=''),overwrite=T)
      }
  }
  ####Saving results
   result= do.call(rbind,res.eval)
   write.csv(result,paste('mapas/v3/',cl,'/results.csv',sep=''))
}
