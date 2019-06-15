nst.boot<-function(nst.result,group=NULL,rand=999,trace=TRUE,
                   two.tail=FALSE,out.detail=FALSE,between.group=FALSE,
                   nworker=1)
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
  
  
  nbt<-function(ik, group, grp.lev, obs3, dist.ran, Dmax, rand, trace, between.group)
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
    
    
    stnst<-function(obs3b,dist.ranb,Dmax)
    {
      obs = as.numeric(obs3b[,3])
      rand.mean = apply(dist.ranb, 1, mean)
      GD=rand.mean/cbind(dist.ranb,obs) # G / D and G / Gk
      EC=(Dmax-rand.mean)/(Dmax-cbind(dist.ranb,obs)) # E / C and E / Ek
      EC[is.nan(EC)]=1
      ECGD=EC
      ECGD[EC>1]=GD[EC>1]
      
      Eij=(Dmax-rand.mean)/Dmax
      STmin=Eij
      STmin[which(obs>rand.mean)]=(1-Eij)[which(obs>rand.mean)]
      
      ECGD.mi=colMeans(ECGD)
      STi.max=max(ECGD.mi)
      STi.min=mean(STmin)
      STi=ECGD.mi[[length(ECGD.mi)]]
      if(STi==STi.min){NSTi=0}else{NSTi=(STi-STi.min)/(STi.max-STi.min)}
      
      c(ST=STi,NST=NSTi)
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
    
    stnsti=stnst(obs3b = obs3[idi,,drop=FALSE],dist.ranb = dist.ran[idi,,drop=FALSE],Dmax=Dmax)
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
                         stnst(obs3b = obs3[idij,,drop=FALSE],dist.ranb = dist.ran[idij,,drop=FALSE],Dmax=Dmax)
                       }))
    }else{
      stnstij=t(sapply(1:rand,
                       function(j)
                       {
                         sampij=sample(sampi,size = length(sampi),replace = TRUE)
                         sampkj=sample(sampk,size = length(sampk),replace = TRUE)
                         idikj=id.2colb(sampij,sampkj,obs3[,1:2,drop=FALSE])
                         stnst(obs3b = obs3[idikj,,drop=FALSE],dist.ranb = dist.ran[idikj,,drop=FALSE],Dmax=Dmax)
                       }))
    }
    stnstij=rbind(stnsti,stnstij)
    STqtij=stats::quantile(stnstij[,1])
    NSTqtij=stats::quantile(stnstij[,2])
    names(STqtij)<-names(NSTqtij)<-c("Min","Quantile25","Median","Quantile75","Max")
    STbp=grDevices::boxplot.stats(stnstij[,1],do.conf = FALSE,do.out = TRUE)
    NSTbp=grDevices::boxplot.stats(stnstij[,2],do.conf = FALSE,do.out = TRUE)
    STboxp=STbp$stats
    NSTboxp=NSTbp$stats
    names(STboxp)<-names(NSTboxp)<-c("LowerWhisker","LowerHinge","Median","HigherHinge","HigherWhisker")
    
    out$ST.stats=c(obs=stnsti[1],mean=mean(stnstij[,1]),stdev=stats::sd(stnstij[,1]),STqtij,STboxp)
    out$NST.stats=c(obs=stnsti[2],mean=mean(stnstij[,2]),stdev=stats::sd(stnstij[,2]),NSTqtij,NSTboxp)
    out$ST.out=STbp$out
    out$NST.out=NSTbp$out
    out$ST.boot=stnstij[,1]
    out$NST.boot=stnstij[,2]
    out$ST.noout=out$ST.boot[which(!(out$ST.boot %in% out$ST.out))]
    out$NST.noout=out$NST.boot[which(!(out$NST.boot %in% out$NST.out))]
    
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
                     dist.ran=dist.ran, Dmax=Dmax, rand=rand, trace=trace, between.group=between.group)
               })
    if(between.group & (length(grp.lev)>1))
    {
      bt2=lapply(1:length(ik2),
                 function(i)
                 {
                   nbt(ik=ik2[[i]], group=group, grp.lev=grp.lev, obs3=obs3,
                       dist.ran=dist.ran, Dmax=Dmax, rand=rand, trace=trace, between.group=between.group)
                 })
    }
  }else{
    c1 <- parallel::makeCluster(nworker, type = "PSOCK")
    message("Now parallel computing within groups. begin at ", date(), ". Please wait...")
    bt1 <- parallel::parLapply(c1, ik1, nbt, group, grp.lev, obs3, dist.ran, Dmax, rand, trace, between.group)
    if(between.group & (length(grp.lev)>1))
    {
      message("Now parallel computing between groups. begin at ", date(), ". Please wait...")
      bt2 <- parallel::parLapply(c1, ik2, nbt, group, grp.lev, obs3, dist.ran, Dmax, rand, trace, between.group)
    }
    parallel::stopCluster(c1)
    gc()
  }
  
  ST.stats=lapply(1:length(bt1),function(i){bt1[[i]]$ST.stats})
  NST.stats=lapply(1:length(bt1),function(i){bt1[[i]]$NST.stats})
  ST.out=lapply(1:length(bt1),function(i){bt1[[i]]$ST.out})
  NST.out=lapply(1:length(bt1),function(i){bt1[[i]]$NST.out})
  ST.boot=lapply(1:length(bt1),function(i){bt1[[i]]$ST.boot})
  NST.boot=lapply(1:length(bt1),function(i){bt1[[i]]$NST.boot})
  ST.noout=lapply(1:length(bt1),function(i){bt1[[i]]$ST.noout})
  NST.noout=lapply(1:length(bt1),function(i){bt1[[i]]$NST.noout})
  if(between.group & (length(grp.lev)>1))
  {
    STb.stats=lapply(1:length(bt2),function(i){bt2[[i]]$ST.stats})
    NSTb.stats=lapply(1:length(bt2),function(i){bt2[[i]]$NST.stats})
    STb.out=lapply(1:length(bt2),function(i){bt2[[i]]$ST.out})
    NSTb.out=lapply(1:length(bt2),function(i){bt2[[i]]$NST.out})
    STb.boot=lapply(1:length(bt2),function(i){bt2[[i]]$ST.boot})
    NSTb.boot=lapply(1:length(bt2),function(i){bt2[[i]]$NST.boot})
    STb.noout=lapply(1:length(bt2),function(i){bt2[[i]]$ST.noout})
    NSTb.noout=lapply(1:length(bt2),function(i){bt2[[i]]$NST.noout})
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
  rownames(ST.sum)<-rownames(NST.sum)<-grp.lev
  
  if(between.group & (length(grp.lev)>1))
  {
    STb.sum=ltom(STb.stats, STb.out)
    NSTb.sum=ltom(NSTb.stats, NSTb.out)
    rownames(STb.sum)<-rownames(NSTb.sum)<-unlist(bgn)
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
  
  st.comp<-nst.comp<-list();k=1
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
        st.comp[[k]]=c(group1=grp.lev[i],group2=grp.lev[j],
                       ST.group1=ST.stats[[i]][1],ST.group2=ST.stats[[j]][1],
                       w.value=ST.wtij$statistic,p.wtest=ifelse(two.tail,ST.wtij$p.value,(ST.wtij$p.value)/2),
                       p.count=ST.cpij,w.value.noOut=STno.wtij$statistic,
                       p.wtest.noOut=ifelse(two.tail,STno.wtij$p.value,(STno.wtij$p.value)/2),
                       p.count.noOut=STno.cpij)
        nst.comp[[k]]=c(group1=grp.lev[i],group2=grp.lev[j],
                       NST.group1=NST.stats[[i]][1],NST.group2=NST.stats[[j]][1],
                       w.value=NST.wtij$statistic,p.wtest=ifelse(two.tail,NST.wtij$p.value,(NST.wtij$p.value)/2),
                       p.count=NST.cpij,w.value.noOut=NSTno.wtij$statistic,
                       p.wtest.noOut=ifelse(two.tail,NSTno.wtij$p.value,(NSTno.wtij$p.value)/2),
                       p.count.noOut=NSTno.cpij)
        k=k+1
      }
    }
    st.comps=Reduce(rbind,st.comp)
    nst.comps=Reduce(rbind,nst.comp)
    rownames(st.comps)<-rownames(nst.comps)<-c()
  }else{
    st.comps=NULL
    nst.comps=NULL
  }
  
  output=list(ST.summary=ST.sum,NST.summary=NST.sum,ST.compare=st.comps,NST.compare=nst.comps)
  if(between.group)
  {
    output=c(output,list(ST.between.summary=STb.sum,NST.between.summary=NSTb.sum))
  }
  if(out.detail)
  {
    names(ST.boot)<-names(NST.boot)<-names(ST.noout)<-names(NST.noout)<-grp.lev
    detail=list(ST.boot=ST.boot,NST.boot=NST.boot,ST.boot.rmout=ST.noout,NST.boot.rmout=NST.noout)
    if(between.group)
    {
      names(STb.boot)<-names(NSTb.boot)<-names(STb.noout)<-names(NSTb.noout)<-unlist(bgn)
      detail=c(detail,list(STb.boot=STb.boot,NSTb.boot=NSTb.boot,STb.boot.rmout=STb.noout,NSTb.boot.rmout=NSTb.noout))
    }
    output=c(output,list(detail=detail))
  }
  output
}

