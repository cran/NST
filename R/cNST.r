cNST<-function(beta.obs, beta.rand.list, group, Dmax=1,
               between.group=FALSE,SES=FALSE,RC=FALSE,
               output.detail=FALSE)
{
  beta.obs=as.matrix(beta.obs)
  beta.r1=as.matrix(beta.rand.list[[1]])
  if(sum(rownames(beta.obs)!=rownames(beta.r1))>0){stop("sample IDs in beta.obs and beta.rand.list do not match.")}
  obs3=dist.3col(beta.obs)
  obs=as.numeric(obs3[,3])
  dist.ran=sapply(1:length(beta.rand.list),function(i){as.numeric(dist.3col(beta.rand.list[[i]])[,3])})
  
  if(min(min(obs,na.rm = TRUE),min(dist.ran,na.rm = TRUE))<0){stop("beta diversity value should not be negative.")}
  if(!is.na(Dmax)){if(Dmax<(max(max(obs,na.rm = TRUE),max(dist.ran,na.rm = TRUE)))){stop("Some beta diversity values are larger than Dmax. Please check.")}}
  
  rand.mean=apply(dist.ran, 1, mean,na.rm=TRUE)
  if(SES)
  {
    rand.sd=apply(dist.ran,1,stats::sd,na.rm=TRUE)
    ses=obs3
    ses[,3]=(obs-rand.mean)/rand.sd
    ses[which(obs==rand.mean),3]=0
    colnames(ses)[3]=paste0("SES.",colnames(obs3)[3])
  }
  if(RC)
  {
    obsmm=matrix(obs,nrow=nrow(dist.ran),ncol=ncol(dist.ran))
    a1=rowSums(dist.ran>obsmm,na.rm = TRUE)/ncol(dist.ran)
    a2=rowSums(dist.ran==obsmm,na.rm = TRUE)/ncol(dist.ran)
    rc=obs3
    rc[,3]=(0.5-(a1+(a2/2)))*2
    colnames(rc)[3]=paste0("RC.",colnames(obs3)[3])
  }
  
  if(is.na(Dmax))
  {
    rand.amax=apply(dist.ran, 1, max,na.rm=TRUE)
    rand.bmax=apply(dist.ran, 1, function(v){max(stats::density(v,na.rm = TRUE)$x)})
    d.up1=(sum((obs<=1)*(rand.amax<=1))==length(obs))*(rand.bmax>1)
    rand.max=d.up1+((1-d.up1)*rand.bmax)
    Dmax=max(obs,rand.max)
  }
  
  GD=rand.mean/cbind(dist.ran,obs) # G / D and G / Gk
  EC=(Dmax-rand.mean)/(Dmax-cbind(dist.ran,obs)) # E / C and E / Ek
  EC[is.nan(EC)]=1
  ECGD=EC
  ECGD[which(EC>1)]=GD[which(EC>1)]
  
  Cij=(Dmax-obs)/Dmax
  Eij=(Dmax-rand.mean)/Dmax
  Dsij=obs/Dmax
  Gsij=rand.mean/Dmax
  
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
  #MSTij=(Eij/Cij)*(Dsij/Gsij)
  #MSTij[which(MSTij>1)]=((Cij/Eij)*(Gsij/Dsij))[which(MSTij>1)]
  indexs=data.frame(obs3[,1:2],D.ij=obs,G.ij=rand.mean,Ds.ij=Dsij,
                    Gs.ij=Gsij,C.ij=Cij, E.ij=Eij, ST.ij=ECGD[,ncol(ECGD)],MST.ij=MSTij)
  colnames(indexs)[3:ncol(indexs)]=paste0(colnames(indexs)[3:ncol(indexs)],".",colnames(obs3)[3])
  STmin=Eij
  STmin[which(obs>rand.mean)]=(1-Eij)[which(obs>rand.mean)]
  
  grp.lev=unique(group[,1])
  outi<-outij<-list()
  if(between.group){outik<-outikj<-list();bgi=1}
  
  for(i in 1:length(grp.lev))
  {
    sampi=rownames(group)[which(group[,1]==grp.lev[i])]
    ni=length(sampi)
    idi=which((obs3[,1] %in% sampi)&(obs3[,2] %in% sampi))
    ECGDi=ECGD[idi,,drop=FALSE]
    ECGD.mi=colMeans(ECGDi,na.rm = TRUE)
    STi.max=max(ECGD.mi,na.rm = TRUE)
    STi.min=mean(STmin[idi])
    STi=ECGD.mi[[length(ECGD.mi)]]
    if(STi==STi.min){NSTi=0}else{NSTi=(STi-STi.min)/(STi.max-STi.min)}
    MSTi=mean(MSTij[idi],na.rm=TRUE)
    outi[[i]]=data.frame(group=grp.lev[i],size=length(idi),ST.i=STi,NST.i=NSTi,MST.i=MSTi,stringsAsFactors = FALSE)
    
    STij=ECGD[idi,ncol(ECGD)]
    NSTij=(STij-STi.min)/(STi.max-STi.min)
    MSTiji=MSTij[idi]
    outij[[i]]=data.frame(group=rep(grp.lev[i],length(idi)),
                          obs3[idi,1:2,drop=FALSE],C.ij=Cij[idi],E.ij=Eij[idi],
                          ST.ij=STij,NST.ij=NSTij,MST.ij=MSTiji,stringsAsFactors = FALSE)
    
    if(SES)
    {
      ses.i=mean(as.numeric(ses[idi,3]),na.rm = TRUE)
      outi[[i]]=data.frame(outi[[i]],SES.i=ses.i,stringsAsFactors = FALSE)
      outij[[i]]=data.frame(outij[[i]],SES.ij=as.numeric(ses[idi,3]),stringsAsFactors = FALSE)
    }
    if(RC)
    {
      rc.i=mean(as.numeric(rc[idi,3]),na.rm = TRUE)
      outi[[i]]=data.frame(outi[[i]],RC.i=rc.i,stringsAsFactors = FALSE)
      outij[[i]]=data.frame(outij[[i]],RC.ij=as.numeric(rc[idi,3]),stringsAsFactors = FALSE)
    }
    
    if(between.group)
    {
      if(i<length(grp.lev))
      {
        for(k in (i+1):length(grp.lev))
        {
          sampk=rownames(group)[which(group[,1]==grp.lev[k])]
          nk=length(sampk)
          idik=which(((obs3[,1] %in% sampi)&(obs3[,2] %in% sampk))|((obs3[,1] %in% sampk)&(obs3[,2] %in% sampi)))
          
          ECGDik=ECGD[idik,,drop=FALSE]
          ECGD.mik=colMeans(ECGDik)
          STik.max=max(ECGD.mik)
          STik.min=mean(STmin[idik])
          STik=ECGD.mik[[length(ECGD.mik)]]
          if(STik==STik.min){NSTik=0}else{NSTik=(STik-STik.min)/(STik.max-STik.min)}
          MSTik=mean(MSTij[idik])
          outik[[bgi]]=data.frame(group=paste0(grp.lev[i],".vs.",grp.lev[k]),size=length(idik),
                                  ST.i=STik,NST.i=NSTik,MST.i=MSTik,stringsAsFactors = FALSE)
          
          STikj=ECGD[idik,ncol(ECGD)]
          NSTikj=(STikj-STik.min)/(STik.max-STik.min)
          MSTikj=MSTij[idik]
          outikj[[bgi]]=data.frame(group=rep(paste0(grp.lev[i],".vs.",grp.lev[k]),length(idik)),
                                   obs3[idik,1:2,drop=FALSE],C.ij=Cij[idik],E.ij=Eij[idik],
                                   ST.ij=STikj,NST.ij=NSTikj,MST.ij=MSTikj,stringsAsFactors = FALSE)
          
          if(SES)
          {
            sesik=mean(as.numeric(ses[idik,3]),na.rm = TRUE)
            outik[[bgi]]=data.frame(outik[[bgi]],SES.i=sesik,stringsAsFactors = FALSE)
            outikj[[bgi]]=data.frame(outikj[[bgi]],SES.ij=as.numeric(ses[idik,3]),stringsAsFactors = FALSE)
          }
          if(RC)
          {
            rcik=mean(as.numeric(rc[idik,3]),na.rm = TRUE)
            outik[[bgi]]=data.frame(outik[[bgi]],RC.i=rcik,stringsAsFactors = FALSE)
            outikj[[bgi]]=data.frame(outikj[[bgi]],RC.ij=as.numeric(rc[idik,3]),stringsAsFactors = FALSE)
          }
          
          bgi=bgi+1
        }
      }
    }
  }
  outi=Reduce(rbind,outi)
  outij=Reduce(rbind,outij)
  colnames(outi)[3:ncol(outi)]=paste0(colnames(outi)[3:ncol(outi)],".",colnames(obs3)[3])
  colnames(outij)[4:ncol(outij)]=paste0(colnames(outij)[4:ncol(outij)],".",colnames(obs3)[3])
  output=list(index.pair=indexs,index.grp=outi,index.pair.grp=outij,Dmax=Dmax)
  if(between.group)
  {
    outik=Reduce(rbind,outik)
    outikj=Reduce(rbind,outikj)
    colnames(outik)[3:ncol(outik)]=paste0(colnames(outik)[3:ncol(outik)],".",colnames(obs3)[3])
    colnames(outikj)[4:ncol(outikj)]=paste0(colnames(outikj)[4:ncol(outikj)],".",colnames(obs3)[3])
    output=c(output,list(index.between=outik,index.pair.between=outikj))
  }
  if(SES){output$index.pair=data.frame(output$index.pair,ses[,3,drop=FALSE],stringsAsFactors = FALSE)}
  if(RC){output$index.pair=data.frame(output$index.pair,rc[,3,drop=FALSE],stringsAsFactors = FALSE)}
  if(output.detail)
  {
    details=list(rand.mean=rand.mean,Dmax=Dmax,obs3=obs3,dist.ran=dist.ran,group=group)
    output=c(output,list(details=details))
  }
  output
}