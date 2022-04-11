pNST<-function(comm, tree=NULL, pd=NULL,pd.desc=NULL,pd.wd=NULL,pd.spname=NULL,
               group, meta.group=NULL, abundance.weighted=TRUE, rand=1000,
               output.rand=FALSE, taxo.null.model=NULL, phylo.shuffle=TRUE,
               exclude.conspecifics=FALSE, nworker=4, LB=FALSE,
               between.group=FALSE, SES=FALSE, RC=FALSE, dirichlet=FALSE)
{
  requireNamespace("iCAMP")
  if(sum(is.na(comm))>0)
  {
    comm[is.na(comm)]=0
    warning("NA is not allowed in comm. automatically filled zero.")
  }
  if(sum(comm<0,na.rm = TRUE)>0)
  {
    stop("Negative value is not allowed in comm. data transformation is not applicable to betaMNTD.")
  }
  if(max(rowSums(comm,na.rm = TRUE))<=1 & (!dirichlet))
  {
    warning("The values in comm are less than 1, thus considered as proportional data, Dirichlet distribution is used to assign abundance in null model.")
    dirichlet=TRUE
  }
  
  groupck<-function(group)
  {
    grp.tab=table(as.vector(group[,1]))
    if(sum(grp.tab<2)>0){stop("some group(s) has only one sample. impossible to perform beta diversity analysis.")}
    if(sum(grp.tab<6)>0){warning("some groups have less than 6 samples, for which NST can be calculated but not recommened.")}
    invisible()
  }
  
  aslist<-function(a){if(is.null(a)){NULL}else{out=list(a);names(out)=deparse(substitute(a));out}}
  sampc=NST::match.name(rn.list = c(aslist(comm),aslist(group),aslist(meta.group)))
  comm=sampc$comm
  group=sampc$group
  groupck(group)
  if(!is.null(meta.group)){meta.group=sampc$meta.group}
  
  if(is.null(pd.desc))
  {
    if(is.null(pd)){
      if(is.null(pd.wd))
      {
        if(is.null(tree)){stop("Tree or pd, at least one should be provided.")}
        warning("Since pd.wd is not specified, a new folder is created in current working directory.")
        time.code=format(Sys.time(),"%y%m%d%H%M")
        pd.wd=paste0(getwd(),"/pdbig.",time.code)
        if(dir.exists(pd.wd)){stop("Newly named pd.wd happens to exist. better specify a pd.wd.")}else{dir.create(pd.wd)}
      }
      
      if(file.exists(paste0(pd.wd,"/pd.desc")))
      {
        warning("Attention: the pd.wd already has a file named pd.desc, which is directly used. Please double check whether this is the phylogenetic distance file you need!")
        pd.big=list()
        pd.big$tip.label=utils::read.csv(paste0(pd.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
        pd.big$pd.wd=pd.wd
        pd.big$pd.file="pd.desc"
        pd.big$pd.name.file="pd.taxon.name.csv"
      }else{
        pd.big=iCAMP::pdist.big(tree = tree, wd=pd.wd, nworker = nworker)
      }
      pd.desc=pd.big$pd.file
      pd.spname=pd.big$tip.label
    }else{
      spc=NST::match.name(cn.list=list(comm=comm),both.list = list(pd=pd))
      comm=spc$comm
      pd=spc$pd
      big=FALSE
    }
  }
  
  if(!is.null(pd.desc))
  {
    spc=NST::match.name(name.check = pd.spname,cn.list=list(comm=comm))
    comm=spc$comm
    big=TRUE
    if(length(pd.spname) != ncol(comm)){stop("pd.spname has some OTUs not in community matrix.")}
  }
  
  if(abundance.weighted){wtn="wt"}else{wtn="uw"}
  if(big)
  {
    obsm=as.matrix(NST::bmntd.big(comm = comm,pd.desc = pd.desc,pd.spname = pd.spname,pd.wd = pd.wd,
                                   abundance.weighted = abundance.weighted,exclude.conspecifics = exclude.conspecifics))
  }else{
    obsm=as.matrix(iCAMP::bmntd(comm = comm, pd = pd, abundance.weighted = abundance.weighted,
                                exclude.conspecifics = exclude.conspecifics))
  }
  obs3=NST::dist.3col(obsm)
  obs=as.numeric(obs3[,3])
  
  permx<-function(comx,rand)
  {
    nnz=which(colSums(comx)>0)
    permat=lapply(1:rand,
                    function(i)
                    {
                      out=1:ncol(comx)
                      out[nnz]=nnz[sample(1:length(nnz),size=length(nnz))]
                      out
                    })
    permat
  }
  
  if(is.null(meta.group))
  {
    permat=permx(comx=comm,rand=rand)
  }else{
    meta.lev=unique(meta.group[,1])
    permat=list()
    for(j in 1:length(meta.lev))
    {
      sampj=rownames(meta.group)[which(meta.group[,1]==meta.lev[j])]
      permat[[j]]=permx(comx=comm[which(rownames(comm) %in% sampj),,drop=FALSE],rand=rand)
    }
  }
  
  if(!is.null(taxo.null.model))
  {
    # Null models definition
    null.models=data.frame(sp.freq=c("equip", "equip", "equip",
                                     "prop", "prop", "prop",
                                     "fix", "fix", "fix",
                                     "fix","fix", "fix", "fix"),
                           samp.rich=c("equip", "prop", "fix",
                                       "equip", "prop", "fix",
                                       "equip", "prop", "fix",
                                       "fix", "fix", "fix", "fix"),
                           swap.method=c("not", "not", "not",
                                         "not", "not", "not",
                                         "not", "not", "swap",
                                         "swap", "tswap", "quasiswap", "backtrack"),
                           stringsAsFactors = FALSE)
    rownames(null.models)=c("EE", "EP", "EF",
                            "PE", "PP", "PF",
                            "FE", "FP", "FF",
                            "FF.swap", "FF.tswap",
                            "FF.quasiswap", "FF.backtrack")
    
    sp.freq=null.models[taxo.null.model,"sp.freq"]
    samp.rich=null.models[taxo.null.model,"samp.rich"]
    swap.method=null.models[taxo.null.model,"swap.method"]
  }else{
    sp.freq<-samp.rich<-swap.method<-NULL
    phylo.shuffle=TRUE
  }
  
  
  dist.rand<-function(i,permat,comm,meta.group=NULL,pd=NULL,
                      pd.desc=NULL,pd.wd=NULL,pd.spname=NULL,
                      abundance.weighted=TRUE,big=TRUE,
                      sp.freq=NULL,samp.rich=NULL,swap.method=NULL,
                      phylo.shuffle=TRUE,exclude.conspecifics=FALSE,dirichlet=FALSE)
  {
    message("Now randomizing i=",i,". ",date())
    requireNamespace("NST")
    if(abundance.weighted){abundance.null="region"}else{abundance.null="not"}
    if(is.null(meta.group))
    {
      if(is.null(sp.freq))
      {
        comr=comm[,permat[[i]],drop=FALSE]
      }else{
        comr=taxo.null(comm = comm,sp.freq = sp.freq,samp.rich = samp.rich,
                       swap.method = swap.method,burnin = i,
                       abundance = abundance.null,region.meta = NULL,dirichlet=dirichlet)
        
        if(phylo.shuffle)
        {
          comr=comr[,permat[[i]],drop=FALSE]
        }
      }
    }else{
      meta.lev=unique(meta.group[,1])
      comr=comm
      for(j in 1:length(meta.lev))
      {
        idj=which(rownames(comm) %in% (rownames(meta.group)[which(meta.group[,1]==meta.lev[j])]))
        if(is.null(sp.freq))
        {
          comr[idj,]=comm[idj,permat[[j]][[i]]]
        }else{
          comrj=taxo.null(comm = comm[idj,,drop=FALSE],sp.freq = sp.freq,samp.rich = samp.rich,
                          swap.method = swap.method,burnin = i,
                          abundance = abundance.null,region.meta = NULL,dirichlet=dirichlet)
          if(phylo.shuffle)
          {
            comr[idj,]=comrj[,permat[[j]][[i]],drop=FALSE]
          }else{
            comr[idj,]=comrj
          }
        }
      }
    }
    rownames(comr)=rownames(comm)
    colnames(comr)=colnames(comm)
    
    if(big)
    {
      disrandm=as.matrix(NST::bmntd.big(comm = comr,pd.desc = pd.desc,pd.spname = pd.spname,pd.wd = pd.wd,spname.check = FALSE,
                                         abundance.weighted = abundance.weighted,exclude.conspecifics = exclude.conspecifics))
    }else{
      disrandm=as.matrix(iCAMP::bmntd(comm=comr,pd=pd,abundance.weighted = abundance.weighted,
                                      exclude.conspecifics = exclude.conspecifics))
    }
    disrand=NST::dist.3col(disrandm)
    as.numeric(disrand[,3])
  }
  
  if(nworker==1)
  {
    dist.ran=sapply(1:rand,dist.rand,permat=permat,comm=comm,meta.group=meta.group,
                    pd=pd,pd.desc=pd.desc,pd.wd=pd.wd,pd.spname=pd.spname,
                    abundance.weighted=abundance.weighted,big=big,
                    sp.freq=sp.freq,samp.rich=samp.rich,swap.method=swap.method,
                    phylo.shuffle=phylo.shuffle,exclude.conspecifics=exclude.conspecifics,
                    dirichlet=dirichlet)
  }else{
    c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
    if(class(c1)[1]=='try-error'){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
    if(class(c1)[1]=='try-error'){c1 <- try(parallel::makeCluster(nworker, setup_strategy = "sequential"))}
    message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
    if(LB)
    {
      dist.ran=parallel::parSapplyLB(c1,1:rand,dist.rand,
                                     permat=permat,comm=comm,meta.group=meta.group,
                                     pd=pd,pd.desc=pd.desc,pd.wd=pd.wd,pd.spname=pd.spname,
                                     abundance.weighted=abundance.weighted,big=big,
                                     sp.freq=sp.freq,samp.rich=samp.rich,swap.method=swap.method,
                                     phylo.shuffle=phylo.shuffle,exclude.conspecifics=exclude.conspecifics,
                                     dirichlet=dirichlet)
    }else{
      dist.ran=parallel::parSapply(c1,1:rand,dist.rand,
                                   permat=permat,comm=comm,meta.group=meta.group,
                                   pd=pd,pd.desc=pd.desc,pd.wd=pd.wd,pd.spname=pd.spname,
                                   abundance.weighted=abundance.weighted,big=big,
                                   sp.freq=sp.freq,samp.rich=samp.rich,swap.method=swap.method,
                                   phylo.shuffle=phylo.shuffle,exclude.conspecifics=exclude.conspecifics,
                                   dirichlet=dirichlet)
    }
    parallel::stopCluster(c1)
  }
  rand.mean=apply(dist.ran, 1, mean)
  if(SES)
  {
    rand.sd=apply(dist.ran,1,stats::sd)
    ses=obs3
    ses[,3]=(obs-rand.mean)/rand.sd
    ses[which(obs==rand.mean),3]=0
    colnames(ses)[3]=paste0("bNTI.",wtn)
  }
  if(RC)
  {
    obsmm=matrix(obs,nrow=nrow(dist.ran),ncol=ncol(dist.ran))
    a1=rowSums(dist.ran>obsmm)/ncol(dist.ran)
    a2=rowSums(dist.ran==obsmm)/ncol(dist.ran)
    rc=obs3
    rc[,3]=(0.5-(a1+(a2/2)))*2
    colnames(rc)[3]=paste0("RC.bMNTD.",wtn)
  }
  rand.amax=apply(dist.ran, 1, max)
  rand.bmax=apply(dist.ran, 1, function(v){max(stats::density(v)$x)})
  d.up1=(sum((obs<=1)*(rand.amax<=1))==length(obs))*(rand.bmax>1)
  rand.max=d.up1+((1-d.up1)*rand.bmax)
  Dmax=max(obs,rand.max)
  
  GD=rand.mean/cbind(dist.ran,obs) # G / D and G / Gk
  EC=(Dmax-rand.mean)/(Dmax-cbind(dist.ran,obs)) # E / C and E / Ek
  EC[is.nan(EC)]=1
  ECGD=EC
  ECGD[EC>1]=GD[EC>1]
  
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
                    Gs.ij=Gsij,C.ij=Cij, E.ij=Eij, ST.ij=ECGD[,ncol(ECGD)], MST.ij=MSTij)
  colnames(indexs)[3:ncol(indexs)]=paste0(colnames(indexs)[3:ncol(indexs)],".bMNTD")
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
    ECGD.mi=colMeans(ECGDi)
    STi.max=max(ECGD.mi)
    STi.min=mean(STmin[idi])
    STi=ECGD.mi[[length(ECGD.mi)]]
    if(STi==STi.min){NSTi=0}else{NSTi=(STi-STi.min)/(STi.max-STi.min)}
    MSTi=mean(MSTij[idi],na.rm=TRUE)
    outi[[i]]=data.frame(group=grp.lev[i],size=ni,ST.i=STi,NST.i=NSTi,MST.i=MSTi,stringsAsFactors = FALSE)
    
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
  colnames(outi)[3:ncol(outi)]=paste0(colnames(outi)[3:ncol(outi)],".bMNTD")
  colnames(outij)[4:ncol(outij)]=paste0(colnames(outij)[4:ncol(outij)],".bMNTD")
  
  output=list(index.pair=indexs,index.grp=outi,index.pair.grp=outij,Dmax=Dmax,dist.method="bMNTD")
  if(between.group)
  {
    outik=Reduce(rbind,outik)
    outikj=Reduce(rbind,outikj)
    colnames(outik)[3:ncol(outik)]=paste0(colnames(outik)[3:ncol(outik)],".bMNTD")
    colnames(outikj)[4:ncol(outikj)]=paste0(colnames(outikj)[4:ncol(outikj)],".bMNTD")
    output=c(output,list(index.between=outik,index.pair.between=outikj))
  }
  if(SES){output$index.pair=data.frame(output$index.pair,ses[,3,drop=FALSE],stringsAsFactors = FALSE)}
  if(RC){output$index.pair=data.frame(output$index.pair,rc[,3,drop=FALSE],stringsAsFactors = FALSE)}
  if(output.rand)
  {
    details=list(rand.mean=rand.mean,Dmax=Dmax,obs3=obs3,dist.ran=dist.ran,group=group,meta.group=meta.group)
    output=c(output,list(details=details))
  }
  output
}