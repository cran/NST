nst.panova<-function(nst.result,group=NULL,rand=999,trace=TRUE,SES=FALSE,RC=FALSE)
{
  if(is.null(nst.result$details)){stop("PANOVA need detailed output of NST.")}
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
  obs=as.numeric(obs3[,3])
  
  if(SES)
  {
    rand.sd=apply(dist.ran,1,stats::sd,na.rm=TRUE)
    sesb=(obs-rand.mean)/rand.sd
    sesb[which(obs==rand.mean)]=0
  }else{sesb=NULL}
  if(RC)
  {
    obsmm=matrix(obs,nrow=nrow(dist.ran),ncol=ncol(dist.ran))
    a1=rowSums(dist.ran>obsmm,na.rm = TRUE)/ncol(dist.ran)
    a2=rowSums(dist.ran==obsmm,na.rm = TRUE)/ncol(dist.ran)
    rcb=(0.5-(a1+(a2/2)))*2
  }else{rcb=NULL}
  
  GD=rand.mean/cbind(dist.ran,obs) # G / D and G / Gk
  EC=(Dmax-rand.mean)/(Dmax-cbind(dist.ran,obs)) # E / C and E / Ek
  EC[is.nan(EC)]=1
  ECGD=EC
  ECGD[EC>1]=GD[EC>1]
  
  Cij=(Dmax-obs)/Dmax
  Eij=(Dmax-rand.mean)/Dmax
  Dsij = obs/Dmax
  Gsij = rand.mean/Dmax
  
  ####
  ECijx=Eij/Cij; CEijx=Cij/Eij
  ECijx[is.nan(ECijx)]=1
  CEijx[is.nan(CEijx)]=1
  DGijx=Dsij/Gsij;GDijx=Gsij/Dsij
  DGijx[is.nan(DGijx)]=1
  GDijx[is.nan(GDijx)]=1
  
  MSTsij = (ECijx) * (DGijx)
  MSTsij[which(MSTsij > 1)] = ((CEijx) * (GDijx))[which(MSTsij > 1)]
  #####
  #MSTsij = (Eij/Cij) * (Dsij/Gsij)
  #MSTsij[which(MSTsij > 1)] = ((Cij/Eij) * (Gsij/Dsij))[which(MSTsij > 1)]
  
  STmin=Eij
  STmin[which(obs>rand.mean)]=(1-Eij)[which(obs>rand.mean)]
  
  grp.lev=unique(group[,1])
  perms<-rands<-Rname<-groups<-list()
  k=1
  sampn=nrow(group)
  permk=rbind(1:sampn,as.matrix(permute::shuffleSet(sampn,rand)))
  perms[[k]]=t(sapply(1:nrow(permk),function(i){rownames(group)[permk[i,]]}))
  groups[[k]]=group
  rands[[k]] = nrow(perms[[k]])
  Rname[[k]]=c("all.grp","all.grp")
  for(i in 1:(length(grp.lev)-1))
  {
    sampi=rownames(group)[group[,1]==grp.lev[i]]
    for(j in (i+1):length(grp.lev))
    {
      sampj=rownames(group)[group[,1]==grp.lev[j]]
      k=k+1
      sampij=c(sampi,sampj)
      permk=rbind(1:length(sampij),as.matrix(permute::shuffleSet(length(sampij),rand)))
      perms[[k]]=t(sapply(1:nrow(permk),function(i){sampij[permk[i,]]}))
      groups[[k]]=group[match(sampij,rownames(group)),,drop=FALSE]
      rands[[k]]=nrow(perms[[k]])
      Rname[[k]]=c(grp.lev[i],grp.lev[j])
    }
  }
  
  eps=.Machine$double.eps
  if(trace){trace.seq=seq(from=1,to=rands[[1]],by=200)}else{trace.seq=NULL}
  
  stnst<-function(u,perms,rands,groups,obs3,ECGD,STmin,trace.seq,MSTsij,SES,RC,sesb,rcb)
  {
    if(!is.null(trace.seq)){if(u %in% trace.seq){message("Now randomization u=",u," in ",rands[[1]],". ",date())}}
    STs<-NSTs<-MSTs<-grps<-STij<-NSTij<-MSTij<-grpij<-list()
    if(SES){SESs<-SESij<-list()}
    if(RC){RCs<-RCij<-list()}
    for(k in 1:length(perms))
    {
      STs[[k]]<-NSTs[[k]]<-MSTs[[k]]<-grps[[k]]<-STij[[k]]<-NSTij[[k]]<-MSTij[[k]]<-grpij[[k]]<-list()
      if(SES){SESs[[k]]<-SESij[[k]]<-list()}
      if(RC){RCs[[k]]<-RCij[[k]]<-list()}
      if(u>nrow(perms[[k]]))
      {
        STs[[k]]<-NSTs[[k]]<-MSTs[[k]]<-grps[[k]]<-STij[[k]]<-NSTij[[k]]<-MSTij[[k]]<-grpij[[k]]<-NA
        if(SES){SESs[[k]]<-SESij[[k]]<-NA}
        if(RC){RCs[[k]]<-RCij[[k]]<-NA}
      }else{
        groupk=groups[[k]]
        rownames(groupk)=perms[[k]][u,]
        grpk.lev=unique(groupk[,1])
        for(i in 1:length(grpk.lev))
        {
          sampi=rownames(groupk)[which(groupk[,1]==grpk.lev[i])]
          ni=length(sampi)
          idi=which((obs3[,1] %in% sampi)&(obs3[,2] %in% sampi))
          ECGDi=ECGD[idi,,drop=FALSE]
          ECGD.mi=colMeans(ECGDi)
          STi.max=max(ECGD.mi)
          STi.min=mean(STmin[idi])
          STi=ECGD.mi[[length(ECGD.mi)]]
          if(STi==STi.min){NSTi=0}else{NSTi=(STi-STi.min)/(STi.max-STi.min)}
          
          STij[[k]][[i]]=ECGD[idi,ncol(ECGD)]
          NSTij[[k]][[i]]=(STij[[k]][[i]]-STi.min)/(STi.max-STi.min)
          MSTij[[k]][[i]]=MSTsij[idi]
          #NSTij[[k]][[i]]=(STij[[k]][[i]]-STmin[idi])/(STi.max-STmin[idi]) # less reasonable
          grpij[[k]][[i]]=rep(grpk.lev[i],length(idi))
          STs[[k]][[i]]=STi
          NSTs[[k]][[i]]=NSTi
          MSTs[[k]][[i]]=mean(MSTsij[idi],na.rm = TRUE)
          grps[[k]][[i]]=grpk.lev[i]
          
          if(SES){SESij[[k]][[i]]=sesb[idi];SESs[[k]][[i]]=mean(sesb[idi],na.rm = TRUE)}
          if(RC){RCij[[k]][[i]]=rcb[idi];RCs[[k]][[i]]=mean(rcb[idi],na.rm = TRUE)}
        }
      }
    }
    dFP<-function(Vij,grpij,Vi,grpi)
    {
      Fobs<-Pobs<-Rname<-Vobs<-dobs<-list()
      for(k in 1:length(Vij))
      {
        if(sum(is.na(Vij[[k]]))>0)
        {
          Fobs[[k]]<-Pobs[[k]]<-dobs[[k]]<-NA
          Rname[[k]]<-Vobs[[k]]<-c(NA,NA)
        }else{
          aovs=try(stats::aov(unlist(Vij[[k]])~unlist(grpij[[k]])),silent = TRUE)
          if(class(aovs)[1]=="try-error")
          {
            Fobs[[k]]=NA
            Pobs[[k]]=NA
          }else{
            aovss=summary(aovs)
            Fobs[[k]]=aovss[[1]]$`F value`[1]
            Pobs[[k]]=aovss[[1]]$`Pr(>F)`[1]
          }
          if(k==1){dobs[[k]]=NA}else{dobs[[k]]=Vi[[k]][[1]]-Vi[[k]][[2]]}
          if(k==1){Rname[[k]]=c("all.grp","all.grp")}else{Rname[[k]]=c(grpi[[k]][[1]],grpi[[k]][[2]])}
          if(k==1){Vobs[[k]]=c(NA,NA)}else{Vobs[[k]]=c(Vi[[k]][[1]],Vi[[k]][[2]])}
        }
      }
      Fobs=unlist(Fobs)
      Pobs=unlist(Pobs)
      dobs=unlist(dobs)
      Vout=data.frame(Reduce(rbind,Rname),Reduce(rbind,Vobs),stringsAsFactors = FALSE)
      rownames(Vout)=c();colnames(Vout)=c("group1","group2",paste0("Index.",c("group1","group2")))
      list(Vobs=Vout,Fvalue=Fobs,Pvalue=Pobs,dif=dobs)
    }
    outx=list(ST.dfp=dFP(Vij=STij,grpij=grpij,Vi=STs,grpi=grps),
         NST.dfp=dFP(Vij=NSTij,grpij=grpij,Vi=NSTs,grpi=grps),
         MST.dfp=dFP(Vij=MSTij,grpij=grpij,Vi=MSTs,grpi=grps))
    if(SES){outx=c(outx,list(SES.dfp=dFP(Vij=SESij,grpij=grpij,Vi=SESs,grpi=grps)))}
    if(RC){outx=c(outx,list(RC.dfp=dFP(Vij=RCij,grpij=grpij,Vi=RCs,grpi=grps)))}
    outx
  }
  
  rres=lapply(1:rands[[1]],stnst,perms,rands,groups,obs3,ECGD,STmin,trace.seq,MSTsij,SES,RC,sesb,rcb)
  ind.names=c("ST","NST","MST")
  if(SES){ind.names=c(ind.names,"SES")}
  if(RC){ind.names=c(ind.names,"RC")}
  
  output=list()
  for(j in 1:length(rres[[1]]))
  {
    fsj=sapply(1:length(rres),function(i){rres[[i]][[j]]$Fvalue})
    fobsj=fsj[,1]
    pobsj=rres[[1]][[j]]$Pvalue
    dsj=matrix(sapply(1:length(rres),function(i){rres[[i]][[j]]$dif}),nrow=length(rres[[1]][[j]]$dif),ncol=length(rres))
    dobsj=dsj[,1]
    paovpj=rowSums(fsj>=matrix((fobsj-eps),nrow=length(fobsj),ncol=ncol(fsj)),na.rm = TRUE)/rowSums(!is.na(fsj))
    pdj.m=dsj-matrix(dobsj,nrow=nrow(dsj),ncol=ncol(dsj))
    pdj=(rowSums(pdj.m>eps,na.rm = TRUE)+(rowSums(pdj.m<=eps,na.rm = TRUE)*0.5))/rowSums(!is.na(pdj.m))
    pdj[which(pdj>0.5)]=1-pdj[which(pdj>0.5)]
    
    output[[j]]=data.frame(index=ind.names[j],rres[[1]][[j]]$Vobs,Difference=dobsj,
                           F.obs=fobsj,P.anova=pobsj,P.panova=paovpj,P.perm=pdj,
                           stringsAsFactors = FALSE)
  }
  out=Reduce(rbind,output)
  rownames(out)=c()
  out
}