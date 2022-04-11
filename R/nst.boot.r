nst.boot<-function(nst.result,group=NULL,rand=999,trace=TRUE,
                   two.tail=FALSE,out.detail=FALSE,between.group=FALSE,
                   nworker=1, SES=FALSE, RC=FALSE)
{
  if(is.null(nst.result$details)){stop("Bootstrapping need detailed output of NST.")}
  rand.mean=nst.result$details$rand.mean
  Dmax=nst.result$details$Dmax
  obs3=nst.result$details$obs3
  if(is.null(group)){group=nst.result$details$group}else{
    group.old=nst.result$details$group
    group=group[match(rownames(group.old),rownames(group)),,drop=FALSE]
    if(sum(is.na(group[,1]))>0)
    {
      warning("some old samples are not in the group file for NST.PANOVA.")
      group=group[!is.na(group[,1]),,drop=FALSE]
    }
  }
  #meta.group=nst.result$details$meta.group
  dist.ran=nst.result$details$dist.ran
  rm(nst.result)
  
  grp.lev=unique(group[,1])
  ik1=lapply(1:length(grp.lev),function(i){c(i,i)})
  if(between.group & (length(grp.lev)>1))
  {
    cbn=utils::combn(length(grp.lev),2)
    ik2=lapply(1:ncol(cbn),function(i){c(cbn[1,i],cbn[2,i])})
  }
  
  
  nbt<-function(ik, group, grp.lev, obs3, dist.ran, Dmax, rand, trace, between.group, SES, RC)
  {
    # functions
    id.2col<-function(samps,names2)
    {
      combs=utils::combn(length(samps),2)
      name1=paste0(names2[,1],"___",names2[,2])
      name2=paste0(names2[,2],"___",names2[,1])
      namei=paste0(samps[combs[1,]],"___",samps[combs[2,]])
      id1=match(namei,name1)
      id2=match(namei,name2)
      id1[is.na(id1)]=id2[is.na(id1)]
      id1[!is.na(id1)]
    }
    
    if(between.group)
    {
      id.2colb<-function(samps1,samps2,names2)
      {
        combs=expand.grid(1:length(samps1),1:length(samps2))
        name1=paste0(names2[,1],"___",names2[,2])
        name2=paste0(names2[,2],"___",names2[,1])
        namei=paste0(samps1[combs[,1]],"___",samps2[combs[,2]])
        id1=match(namei,name1)
        id2=match(namei,name2)
        id1[is.na(id1)]=id2[is.na(id1)]
        id1[!is.na(id1)]
      }
    }
    
    
    stnst<-function(obs3b,dist.ranb,Dmax,SES,RC)
    {
      obs = as.numeric(obs3b[,3])
      rand.mean = apply(dist.ranb, 1, mean)
      
      if(SES)
      {
        rand.sd=apply(dist.ranb,1,stats::sd,na.rm=TRUE)
        sesb=(obs-rand.mean)/rand.sd
        sesb[which(obs==rand.mean)]=0
      }
      if(RC)
      {
        obsmm=matrix(obs,nrow=nrow(dist.ranb),ncol=ncol(dist.ranb))
        a1=rowSums(dist.ranb>obsmm,na.rm = TRUE)/ncol(dist.ranb)
        a2=rowSums(dist.ranb==obsmm,na.rm = TRUE)/ncol(dist.ranb)
        rcb=(0.5-(a1+(a2/2)))*2
      }
      
      GD=rand.mean/cbind(dist.ranb,obs) # G / D and G / Gk
      EC=(Dmax-rand.mean)/(Dmax-cbind(dist.ranb,obs)) # E / C and E / Ek
      EC[is.nan(EC)]=1
      ECGD=EC
      ECGD[EC>1]=GD[EC>1]
      
      Eij=(Dmax-rand.mean)/Dmax
      Cij = (Dmax - obs)/Dmax
      Dsij = obs/Dmax
      Gsij = rand.mean/Dmax
      
      ####
      ECijx=Eij/Cij; CEijx=Cij/Eij
      ECijx[is.nan(ECijx)]=1
      CEijx[is.nan(CEijx)]=1
      DGijx=Dsij/Gsij;GDijx=Gsij/Dsij
      DGijx[is.nan(DGijx)]=1
      GDijx[is.nan(GDijx)]=1
      
      MSTij = (ECijx) * (DGijx)
      MSTij[which(MSTij > 1)] = ((CEijx) * (GDijx))[which(MSTij > 1)]
      #####
      #MSTij = (Eij/Cij) * (Dsij/Gsij)
      #MSTij[which(MSTij > 1)] = ((Cij/Eij) * (Gsij/Dsij))[which(MSTij > 1)]
      MSTi=mean(MSTij,na.rm = TRUE)
      
      STmin=Eij
      STmin[which(obs>rand.mean)]=(1-Eij)[which(obs>rand.mean)]
      
      ECGD.mi=colMeans(ECGD)
      STi.max=max(ECGD.mi)
      STi.min=mean(STmin)
      STi=ECGD.mi[[length(ECGD.mi)]]
      if(STi==STi.min){NSTi=0}else{NSTi=(STi-STi.min)/(STi.max-STi.min)}
      
      outxx=c(ST=STi,NST=NSTi,MST=MSTi)
      if(SES){outxx=c(outxx,SES=mean(sesb,na.rm = TRUE))}
      if(RC){outxx=c(outxx,RC=mean(rcb,na.rm = TRUE))}
      outxx
    }
    
    # calculation
    i=ik[1]
    k=ik[2]
    out=list()
    if(trace){message("Now bootstrapping group i=",i," k=",k," in ",length(grp.lev),". ",date())}
    
    sampi=rownames(group)[which(group[,1]==grp.lev[i])]
    if(i==k)
    {
      idi=which((obs3[,1] %in% sampi)&(obs3[,2] %in% sampi))
    }else{
      sampk=rownames(group)[which(group[,1]==grp.lev[k])]
      idi=which(((obs3[,1] %in% sampi)&(obs3[,2] %in% sampk))|((obs3[,1] %in% sampk)&(obs3[,2] %in% sampi)))
    }
    
    stnsti=stnst(obs3b = obs3[idi,,drop=FALSE],dist.ranb = dist.ran[idi,,drop=FALSE],Dmax=Dmax,SES=SES,RC=RC)
    if(i==k)
    {
      stnstij=t(sapply(1:rand,
                       function(j)
                       {
                         sampij=sample(sampi,size = length(sampi),replace = TRUE)
                         idij=id.2col(sampij,obs3[,1:2,drop=FALSE])
                         t=1
                         while((length(idij)<1)&(t<100))
                         {
                           sampij=sample(sampi,size = length(sampi),replace = TRUE)
                           idij=id.2col(sampij,obs3[,1:2,drop=FALSE])
                           t=t+1
                         }
                         stnst(obs3b = obs3[idij,,drop=FALSE],dist.ranb = dist.ran[idij,,drop=FALSE],
                               Dmax=Dmax,SES=SES,RC=RC)
                       }))
    }else{
      stnstij=t(sapply(1:rand,
                       function(j)
                       {
                         sampij=sample(sampi,size = length(sampi),replace = TRUE)
                         sampkj=sample(sampk,size = length(sampk),replace = TRUE)
                         idikj=id.2colb(sampij,sampkj,obs3[,1:2,drop=FALSE])
                         stnst(obs3b = obs3[idikj,,drop=FALSE],dist.ranb = dist.ran[idikj,,drop=FALSE],
                               Dmax=Dmax,SES = SES,RC=RC)
                       }))
    }
    stnstij=rbind(obs=stnsti,stnstij)
    STqtij=stats::quantile(stnstij[,1])
    NSTqtij=stats::quantile(stnstij[,2])
    MSTqtij=stats::quantile(stnstij[,3])
    if(SES){SESqtij=stats::quantile(stnstij[,"SES"]);names(SESqtij)<-c("Min","Quantile25","Median","Quantile75","Max")}
    if(RC){RCqtij=stats::quantile(stnstij[,"RC"]);names(RCqtij)<-c("Min","Quantile25","Median","Quantile75","Max")}
    names(STqtij)<-names(NSTqtij)<-names(MSTqtij)<-c("Min","Quantile25","Median","Quantile75","Max")
    
    
    STbp=grDevices::boxplot.stats(stnstij[,1],do.conf = FALSE,do.out = TRUE)
    NSTbp=grDevices::boxplot.stats(stnstij[,2],do.conf = FALSE,do.out = TRUE)
    MSTbp=grDevices::boxplot.stats(stnstij[,3],do.conf = FALSE,do.out = TRUE)
    if(SES){SESbp=grDevices::boxplot.stats(stnstij[,"SES"],do.conf = FALSE,do.out = TRUE)}
    if(RC){RCbp=grDevices::boxplot.stats(stnstij[,"RC"],do.conf = FALSE,do.out = TRUE)}
    
    STboxp=STbp$stats
    NSTboxp=NSTbp$stats
    MSTboxp=MSTbp$stats
    if(SES){SESboxp=SESbp$stats;names(SESboxp)<-c("LowerWhisker","LowerHinge","Median","HigherHinge","HigherWhisker")}
    if(RC){RCboxp=RCbp$stats;names(RCboxp)<-c("LowerWhisker","LowerHinge","Median","HigherHinge","HigherWhisker")}
    names(STboxp)<-names(NSTboxp)<-names(MSTboxp)<-c("LowerWhisker","LowerHinge","Median","HigherHinge","HigherWhisker")
    
    out$ST.stats=c(obs=stnsti[[1]],mean=mean(stnstij[,1],na.rm = TRUE),stdev=stats::sd(stnstij[,1],na.rm = TRUE),STqtij,STboxp)
    out$NST.stats=c(obs=stnsti[[2]],mean=mean(stnstij[,2],na.rm = TRUE),stdev=stats::sd(stnstij[,2],na.rm = TRUE),NSTqtij,NSTboxp)
    out$MST.stats=c(obs=stnsti[[3]],mean=mean(stnstij[,3],na.rm = TRUE),stdev=stats::sd(stnstij[,3],na.rm = TRUE),MSTqtij,MSTboxp)
    if(SES){out$SES.stats=c(obs=stnsti[["SES"]],mean=mean(stnstij[,"SES"],na.rm = TRUE),stdev=stats::sd(stnstij[,"SES"],na.rm = TRUE),SESqtij,SESboxp)}
    if(RC){out$RC.stats=c(obs=stnsti[["RC"]],mean=mean(stnstij[,"RC"],na.rm = TRUE),stdev=stats::sd(stnstij[,"RC"],na.rm = TRUE),RCqtij,RCboxp)}
    out$ST.out=STbp$out
    out$NST.out=NSTbp$out
    out$MST.out=MSTbp$out
    if(SES){out$SES.out=SESbp$out}
    if(RC){out$SES.out=SESbp$out}
    out$ST.boot=stnstij[,1]
    out$NST.boot=stnstij[,2]
    out$MST.boot=stnstij[,3]
    if(SES){out$SES.boot=stnstij[,"SES"]}
    if(RC){out$RC.boot=stnstij[,"RC"]}
    out$ST.noout=out$ST.boot[which(!(out$ST.boot %in% out$ST.out))]
    out$NST.noout=out$NST.boot[which(!(out$NST.boot %in% out$NST.out))]
    out$MST.noout=out$MST.boot[which(!(out$MST.boot %in% out$MST.out))]
    if(SES){out$SES.noout=out$SES.boot[which(!(out$SES.boot %in% out$SES.out))]}
    if(RC){out$RC.noout=out$RC.boot[which(!(out$RC.boot %in% out$RC.out))]}
    
    if(i==k)
    {
      out$id=grp.lev[i]
    }else{
      out$id=paste0(grp.lev[i],".vs.",grp.lev[k])
    }
    out
  }
  
  if(nworker==1)
  {
    bt1=lapply(1:length(ik1),
               function(i)
               {
                 nbt(ik=ik1[[i]], group=group, grp.lev=grp.lev, obs3=obs3,
                     dist.ran=dist.ran, Dmax=Dmax, rand=rand, trace=trace,
                     between.group=between.group, SES = SES, RC = RC)
               })
    if(between.group & (length(grp.lev)>1))
    {
      bt2=lapply(1:length(ik2),
                 function(i)
                 {
                   nbt(ik=ik2[[i]], group=group, grp.lev=grp.lev, obs3=obs3,
                       dist.ran=dist.ran, Dmax=Dmax, rand=rand, trace=trace,
                       between.group=between.group, SES = SES, RC = RC)
                 })
    }
  }else{
    c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
    if(class(c1)[1]=='try-error'){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
    if(class(c1)[1]=='try-error'){c1 <- try(parallel::makeCluster(nworker, setup_strategy = "sequential"))}
    message("Now parallel computing within groups. begin at ", date(), ". Please wait...")
    bt1 <- parallel::parLapply(c1, ik1, nbt, group, grp.lev, obs3, dist.ran, Dmax, rand, trace, between.group, SES, RC)
    if(between.group & (length(grp.lev)>1))
    {
      message("Now parallel computing between groups. begin at ", date(), ". Please wait...")
      bt2 <- parallel::parLapply(c1, ik2, nbt, group, grp.lev, obs3, dist.ran, Dmax, rand, trace, between.group, SES, RC)
    }
    parallel::stopCluster(c1)
    gc()
  }
  
  ST.stats=lapply(1:length(bt1),function(i){bt1[[i]]$ST.stats})
  NST.stats=lapply(1:length(bt1),function(i){bt1[[i]]$NST.stats})
  MST.stats=lapply(1:length(bt1),function(i){bt1[[i]]$MST.stats})
  if(SES){SES.stats=lapply(1:length(bt1),function(i){bt1[[i]]$SES.stats})}
  if(RC){RC.stats=lapply(1:length(bt1),function(i){bt1[[i]]$RC.stats})}
  
  ST.out=lapply(1:length(bt1),function(i){bt1[[i]]$ST.out})
  NST.out=lapply(1:length(bt1),function(i){bt1[[i]]$NST.out})
  MST.out=lapply(1:length(bt1),function(i){bt1[[i]]$MST.out})
  if(SES){SES.out=lapply(1:length(bt1),function(i){bt1[[i]]$SES.out})}
  if(RC){RC.out=lapply(1:length(bt1),function(i){bt1[[i]]$RC.out})}
  
  ST.boot=lapply(1:length(bt1),function(i){bt1[[i]]$ST.boot})
  NST.boot=lapply(1:length(bt1),function(i){bt1[[i]]$NST.boot})
  MST.boot=lapply(1:length(bt1),function(i){bt1[[i]]$MST.boot})
  if(SES){SES.boot=lapply(1:length(bt1),function(i){bt1[[i]]$SES.boot})}
  if(RC){RC.boot=lapply(1:length(bt1),function(i){bt1[[i]]$RC.boot})}
  
  ST.noout=lapply(1:length(bt1),function(i){bt1[[i]]$ST.noout})
  NST.noout=lapply(1:length(bt1),function(i){bt1[[i]]$NST.noout})
  MST.noout=lapply(1:length(bt1),function(i){bt1[[i]]$MST.noout})
  if(SES){SES.noout=lapply(1:length(bt1),function(i){bt1[[i]]$SES.noout})}
  if(RC){RC.noout=lapply(1:length(bt1),function(i){bt1[[i]]$RC.noout})}
  
  if(between.group & (length(grp.lev)>1))
  {
    STb.stats=lapply(1:length(bt2),function(i){bt2[[i]]$ST.stats})
    NSTb.stats=lapply(1:length(bt2),function(i){bt2[[i]]$NST.stats})
    MSTb.stats=lapply(1:length(bt2),function(i){bt2[[i]]$MST.stats})
    if(SES){SESb.stats=lapply(1:length(bt2),function(i){bt2[[i]]$SES.stats})}
    if(RC){RCb.stats=lapply(1:length(bt2),function(i){bt2[[i]]$RC.stats})}
    
    STb.out=lapply(1:length(bt2),function(i){bt2[[i]]$ST.out})
    NSTb.out=lapply(1:length(bt2),function(i){bt2[[i]]$NST.out})
    MSTb.out=lapply(1:length(bt2),function(i){bt2[[i]]$MST.out})
    if(SES){SESb.out=lapply(1:length(bt2),function(i){bt2[[i]]$SES.out})}
    if(RC){RCb.out=lapply(1:length(bt2),function(i){bt2[[i]]$RC.out})}
    
    STb.boot=lapply(1:length(bt2),function(i){bt2[[i]]$ST.boot})
    NSTb.boot=lapply(1:length(bt2),function(i){bt2[[i]]$NST.boot})
    MSTb.boot=lapply(1:length(bt2),function(i){bt2[[i]]$MST.boot})
    if(SES){SESb.boot=lapply(1:length(bt2),function(i){bt2[[i]]$SES.boot})}
    if(RC){RCb.boot=lapply(1:length(bt2),function(i){bt2[[i]]$RC.boot})}
    
    STb.noout=lapply(1:length(bt2),function(i){bt2[[i]]$ST.noout})
    NSTb.noout=lapply(1:length(bt2),function(i){bt2[[i]]$NST.noout})
    MSTb.noout=lapply(1:length(bt2),function(i){bt2[[i]]$MST.noout})
    if(SES){SESb.noout=lapply(1:length(bt2),function(i){bt2[[i]]$SES.noout})}
    if(RC){RCb.noout=lapply(1:length(bt2),function(i){bt2[[i]]$RC.noout})}
    
    bgn=lapply(1:length(bt2),function(i){bt2[[i]]$id})
  }
  
  ltom<-function(l.stats,l.out)
  {
    lenm=max(sapply(l.out,length))
    if(lenm>0)
    {
      outm=t(sapply(1:length(l.out),
               function(i)
               {
                 c(l.stats[[i]],l.out[[i]],rep(NA,lenm-length(l.out[[i]])))
               }))
      colnames(outm)=c(names(l.stats[[1]]),paste0("Outlier",1:lenm))
    }else{
      if(length(l.stats)==1){outm=t(l.stats[[1]])}else{outm=Reduce(rbind,l.stats)}
    }
    outm
  }
  
  ST.sum=ltom(ST.stats, ST.out)
  NST.sum=ltom(NST.stats, NST.out)
  MST.sum=ltom(MST.stats, MST.out)
  if(SES){SES.sum=ltom(SES.stats, SES.out);rownames(SES.sum)<-grp.lev}
  if(RC){RC.sum=ltom(RC.stats, RC.out);rownames(RC.sum)<-grp.lev}
  
  rownames(ST.sum)<-rownames(NST.sum)<-rownames(MST.sum)<-grp.lev
  
  if(between.group & (length(grp.lev)>1))
  {
    STb.sum=ltom(STb.stats, STb.out)
    NSTb.sum=ltom(NSTb.stats, NSTb.out)
    MSTb.sum=ltom(MSTb.stats, MSTb.out)
    if(SES){SESb.sum=ltom(SESb.stats, SESb.out);rownames(SESb.sum)<-unlist(bgn)}
    if(RC){RCb.sum=ltom(RCb.stats, RCb.out);rownames(RCb.sum)<-unlist(bgn)}
    rownames(STb.sum)<-rownames(NSTb.sum)<-rownames(MSTb.sum)<-unlist(bgn)
  }
  
  p.count<-function(x,y,obs.dxy,two.tail=FALSE)
  {
    xm=matrix(x,nrow=length(x),ncol=length(y))
    ym=matrix(y,nrow=length(x),ncol=length(y),byrow = TRUE)
    ns=sum(!is.na(x))*sum(!is.na(y))
    EPS <- sqrt(.Machine$double.eps)
    if(obs.dxy>=0)
    {
      p=sum(xm<=(ym+EPS),na.rm = TRUE)/ns
    }else{
      p=sum(xm>=(ym-EPS),na.rm = TRUE)/ns
    }
    if(two.tail){p=2*p}
    p
  }
  
  st.comp<-nst.comp<-mst.comp<-list();k=1
  if(SES){ses.comp<-list()}
  if(RC){rc.comp<-list()}
  if(length(grp.lev)>1)
  {
    for(i in 1:(length(grp.lev)-1))
    {
      for(j in (i+1):(length(grp.lev)))
      {
        if(trace){message("Now comparing between group i=",i," and group j=",j,". ",date())}
        ST.wtij=stats::wilcox.test(ST.boot[[i]],ST.boot[[j]])
        STno.wtij=stats::wilcox.test(ST.noout[[i]],ST.noout[[j]])
        ST.obsdij=ST.stats[[i]][[1]]-ST.stats[[j]][[1]]
        ST.cpij=p.count(ST.boot[[i]],ST.boot[[j]],ST.obsdij,two.tail = two.tail)
        STno.cpij=p.count(ST.noout[[i]],ST.noout[[j]],ST.obsdij,two.tail = two.tail)
        
        NST.wtij=stats::wilcox.test(NST.boot[[i]],NST.boot[[j]])
        NSTno.wtij=stats::wilcox.test(NST.noout[[i]],NST.noout[[j]])
        NST.obsdij=NST.stats[[i]][[1]]-NST.stats[[j]][[1]]
        NST.cpij=p.count(NST.boot[[i]],NST.boot[[j]],NST.obsdij,two.tail = two.tail)
        NSTno.cpij=p.count(NST.noout[[i]],NST.noout[[j]],NST.obsdij,two.tail = two.tail)
        
        MST.wtij=stats::wilcox.test(MST.boot[[i]],MST.boot[[j]])
        MSTno.wtij=stats::wilcox.test(MST.noout[[i]],MST.noout[[j]])
        MST.obsdij=MST.stats[[i]][[1]]-MST.stats[[j]][[1]]
        MST.cpij=p.count(MST.boot[[i]],MST.boot[[j]],MST.obsdij,two.tail = two.tail)
        MSTno.cpij=p.count(MST.noout[[i]],MST.noout[[j]],MST.obsdij,two.tail = two.tail)
        
        st.comp[[k]]=c(group1=grp.lev[i],group2=grp.lev[j],
                       group1.obs=ST.stats[[i]][[1]],group2.obs=ST.stats[[j]][[1]],
                       w.value=ST.wtij$statistic,p.wtest=ifelse(two.tail,ST.wtij$p.value,(ST.wtij$p.value)/2),
                       p.count=ST.cpij,w.value.noOut=STno.wtij$statistic,
                       p.wtest.noOut=ifelse(two.tail,STno.wtij$p.value,(STno.wtij$p.value)/2),
                       p.count.noOut=STno.cpij)
        nst.comp[[k]]=c(group1=grp.lev[i],group2=grp.lev[j],
                        group1.obs=NST.stats[[i]][[1]],group2.obs=NST.stats[[j]][[1]],
                        w.value=NST.wtij$statistic,p.wtest=ifelse(two.tail,NST.wtij$p.value,(NST.wtij$p.value)/2),
                        p.count=NST.cpij,w.value.noOut=NSTno.wtij$statistic,
                        p.wtest.noOut=ifelse(two.tail,NSTno.wtij$p.value,(NSTno.wtij$p.value)/2),
                        p.count.noOut=NSTno.cpij)
        mst.comp[[k]]=c(group1=grp.lev[i],group2=grp.lev[j],
                        group1.obs=MST.stats[[i]][[1]],group2.obs=MST.stats[[j]][[1]],
                        w.value=MST.wtij$statistic,p.wtest=ifelse(two.tail,MST.wtij$p.value,(MST.wtij$p.value)/2),
                        p.count=MST.cpij,w.value.noOut=MSTno.wtij$statistic,
                        p.wtest.noOut=ifelse(two.tail,MSTno.wtij$p.value,(MSTno.wtij$p.value)/2),
                        p.count.noOut=MSTno.cpij)
        
        if(SES)
        {
          SES.wtij=stats::wilcox.test(SES.boot[[i]],SES.boot[[j]])
          SESno.wtij=stats::wilcox.test(SES.noout[[i]],SES.noout[[j]])
          SES.obsdij=SES.stats[[i]][[1]]-SES.stats[[j]][[1]]
          SES.cpij=p.count(SES.boot[[i]],SES.boot[[j]],SES.obsdij,two.tail = two.tail)
          SESno.cpij=p.count(SES.noout[[i]],SES.noout[[j]],SES.obsdij,two.tail = two.tail)
          ses.comp[[k]]=c(group1=grp.lev[i],group2=grp.lev[j],
                          group1.obs=SES.stats[[i]][[1]],group2.obs=SES.stats[[j]][[1]],
                          w.value=SES.wtij$statistic,p.wtest=ifelse(two.tail,SES.wtij$p.value,(SES.wtij$p.value)/2),
                          p.count=SES.cpij,w.value.noOut=SESno.wtij$statistic,
                          p.wtest.noOut=ifelse(two.tail,SESno.wtij$p.value,(SESno.wtij$p.value)/2),
                          p.count.noOut=SESno.cpij)
        }
        if(RC)
        {
          RC.wtij=stats::wilcox.test(RC.boot[[i]],RC.boot[[j]])
          RCno.wtij=stats::wilcox.test(RC.noout[[i]],RC.noout[[j]])
          RC.obsdij=RC.stats[[i]][[1]]-RC.stats[[j]][[1]]
          RC.cpij=p.count(RC.boot[[i]],RC.boot[[j]],RC.obsdij,two.tail = two.tail)
          RCno.cpij=p.count(RC.noout[[i]],RC.noout[[j]],RC.obsdij,two.tail = two.tail)
          rc.comp[[k]]=c(group1=grp.lev[i],group2=grp.lev[j],
                          group1.obs=RC.stats[[i]][[1]],group2.obs=RC.stats[[j]][[1]],
                          w.value=RC.wtij$statistic,p.wtest=ifelse(two.tail,RC.wtij$p.value,(RC.wtij$p.value)/2),
                          p.count=RC.cpij,w.value.noOut=RCno.wtij$statistic,
                          p.wtest.noOut=ifelse(two.tail,RCno.wtij$p.value,(RCno.wtij$p.value)/2),
                          p.count.noOut=RCno.cpij)
        }
        k=k+1
      }
    }
    if(length(st.comp)==1)
    {
      st.comps=t(st.comp[[1]])
      nst.comps=t(nst.comp[[1]])
      mst.comps=t(mst.comp[[1]])
      if(SES){ses.comps=t(ses.comp[[1]])}
      if(RC){rc.comps=t(rc.comp[[1]])}
    }else{
      st.comps=Reduce(rbind,st.comp)
      nst.comps=Reduce(rbind,nst.comp)
      mst.comps=Reduce(rbind,mst.comp)
      rownames(st.comps)<-rownames(nst.comps)<-rownames(mst.comps)<-c()
      if(SES){ses.comps=Reduce(rbind,ses.comp);rownames(ses.comps)<-c()}
      if(RC){rc.comps=Reduce(rbind,rc.comp);rownames(rc.comps)<-c()}
    }
    outcp=data.frame(Index=c(rep("ST",nrow(st.comps)),
                             rep("NST",nrow(nst.comps)),
                             rep("MST",nrow(mst.comps))),
                     rbind(st.comps,nst.comps,mst.comps),
                     stringsAsFactors = FALSE)
    if(SES){outcp=rbind(outcp,data.frame(Index=rep("SES",nrow(ses.comps)),
                                         ses.comps,stringsAsFactors = FALSE))}
    if(RC){outcp=rbind(outcp,data.frame(Index=rep("RC",nrow(rc.comps)),
                                        rc.comps,stringsAsFactors = FALSE))}
  }else{
    outcp=NULL
  }
  
  if(between.group)
  {
    ml=list(ST.sum,STb.sum,NST.sum,NSTb.sum,MST.sum,MSTb.sum)
    ml.indexs=c(rep("ST",nrow(ml[[1]])+nrow(ml[[2]])),
                rep("NST",nrow(ml[[3]])+nrow(ml[[4]])),
                rep("MST",nrow(ml[[5]])+nrow(ml[[6]])))
    if(SES){ml=c(ml,list(SES.sum,SESb.sum));ml.indexs=c(ml.indexs,rep("SES",nrow(ml[[length(ml)-1]])+nrow(ml[[length(ml)]])))}
    if(RC){ml=c(ml,list(RC.sum,RCb.sum));ml.indexs=c(ml.indexs,rep("RC",nrow(ml[[length(ml)-1]])+nrow(ml[[length(ml)]])))}
    
  }else{
    ml=list(ST.sum,NST.sum,MST.sum)
    ml.indexs=c(rep("ST",nrow(ml[[1]])),rep("NST",nrow(ml[[2]])),rep("MST",nrow(ml[[3]])))
    if(SES){ml=c(ml,list(SES.sum));ml.indexs=c(ml.indexs,rep("SES",nrow(ml[[length(ml)]])))}
    if(RC){ml=c(ml,list(RC.sum));ml.indexs=c(ml.indexs,rep("RC",nrow(ml[[length(ml)]])))}
  }
  ml.grps=unlist(lapply(ml,function(m){rownames(m)}))
  ml.len=sapply(ml,ncol)
  outml=Reduce(rbind,lapply(ml,function(m){outm=cbind(m,matrix(NA,nrow=nrow(m),ncol=max(ml.len)-ncol(m)));colnames(outm)=c();outm}))
  colnames(outml)=colnames(ml[[which.max(ml.len)]])
  outlr=data.frame(Index=ml.indexs,Group=ml.grps, outml,stringsAsFactors = FALSE)
  rownames(outlr)=c()
  
  output=list(summary=outlr,compare=outcp)
  
  if(out.detail)
  {
    names(ST.boot)<-names(NST.boot)<-names(MST.boot)<-names(ST.noout)<-names(NST.noout)<-names(MST.noout)<-grp.lev
    if(SES){names(SES.boot)<-names(SES.noout)<-grp.lev}
    if(RC){names(RC.boot)<-names(RC.noout)<-grp.lev}
    
    detail=list(ST.boot=ST.boot,NST.boot=NST.boot,MST.boot=MST.boot,
                ST.boot.rmout=ST.noout,NST.boot.rmout=NST.noout,MST.boot.rmout=MST.noout)
    
    if(SES){detail=c(detail,list(SES.boot=SES.boot,SES.boot.rmout=SES.noout))}
    if(RC){detail=c(detail,list(RC.boot=RC.boot,RC.boot.rmout=RC.noout))}
    
    if(between.group)
    {
      names(STb.boot)<-names(NSTb.boot)<-names(MSTb.boot)<-names(STb.noout)<-names(NSTb.noout)<-names(MSTb.noout)<-unlist(bgn)
      if(SES){names(SESb.boot)<-names(SESb.noout)<-unlist(bgn)}
      if(RC){names(RCb.boot)<-names(RCb.noout)<-unlist(bgn)}
      detail=c(detail,list(STb.boot=STb.boot,NSTb.boot=NSTb.boot,MSTb.boot=MSTb.boot,
                           STb.boot.rmout=STb.noout,NSTb.boot.rmout=NSTb.noout,MSTb.boot.rmout=MSTb.noout))
      if(SES){detail=c(detail,list(SESb.boot=SESb.boot,SESb.boot.rmout=SESb.noout))}
      if(RC){detail=c(detail,list(RCb.boot=RCb.boot,RCb.boot.rmout=RCb.noout))}
    }
    output=c(output,list(detail=detail))
  }
  output
}

