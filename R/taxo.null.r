taxo.null<-function(comm,sp.freq=c("not","equip","prop","prop.ab","fix"),
                    samp.rich=c("not","equip","prop","fix"),
                    swap.method=c("not","swap","tswap","quasiswap","backtrack"),burnin=0,
                    abundance=c("not","shuffle","local","region"),region.meta=NULL)
{
  # universe of null model
  # 2016.3.30 revised by Daliang. Include an option to count frequency as regional abundance sum when sp.freq="prop.ab"
  # originated by Daliang Ning on 2015.9.27
  comm.b=comm
  comm.b[comm.b>0]=1
  #if(abundance!="not"){source(paste(code.wd,"/ab.assign.r",sep = ""))}
  sp.freq=sp.freq[1];samp.rich=samp.rich[1];abundance=abundance[1];swap.method=swap.method[1]
  
  if(sp.freq=="not"|samp.rich=="not")
  {
    if(abundance=="not")
    {
      out=comm.b
    }else if(abundance=="shuffle"){
      sp.num=ncol(comm)
      out=t(sapply(1:nrow(comm),
             function(i)
             {
               res=rep(0,sp.num)
               xx=comm[i,comm[i,]>0]
               res[comm[i,]>0]=xx[sample.int(length(xx))]
               res
             }))
      rownames(out)<-rownames(comm.b)
      colnames(out)<-colnames(comm.b)
    }else if(abundance=="local"){
      prob.ab=comm
      samp.ab=rowSums(comm)
      out=ab.assign(comm.b = comm.b,samp.ab = samp.ab,prob.ab = prob.ab)
    }else if(abundance=="region"){
      prob.ab=matrix(colSums(comm),nrow = nrow(comm),ncol = ncol(comm),byrow = TRUE)
      samp.ab=rowSums(comm)
      out=ab.assign(comm.b = comm.b,samp.ab = samp.ab,prob.ab = prob.ab)
    }else{stop("Abundance assignment method input error.")}
  }else{
    # all the 9 combinations
    if(sp.freq=="equip")
    {
      if(samp.rich=="equip")
      {
        comr.b=matrix(sample(comm.b),nrow=nrow(comm.b),ncol(comm.b))
      }else if(samp.rich=="prop"){
        samp.richness=rowSums(comm.b)
        size.rand=sum(comm.b)
        samp.num=nrow(comm)
        sp.num=ncol(comm)
        prob.sp=as.vector(matrix(samp.richness,nrow = samp.num,ncol=sp.num))
        comr.b=matrix(0,nrow = samp.num,ncol=sp.num)
        id.out=sample(1:(samp.num*sp.num),size =size.rand,replace = FALSE,prob = prob.sp)
        comr.b[id.out]=1
      }else if(samp.rich=="fix"){
        samp.richness=rowSums(comm.b)
        sp.num=ncol(comm)
        comr.b=t(sapply(1:nrow(comm), function(i){res=rep(0,sp.num);id=sample(1:sp.num,samp.richness[i]);res[id]=1;res}))
      }else{stop("Sample richness randomization method input error.")}
    }else if(sp.freq=="prop"|sp.freq=="prop.ab"){
      if(sp.freq=="prop")
      {
        sp.freqence=colSums(comm.b)
      }else{
        sp.freqence=colSums(comm) # add this option 2016.3.30
      }
      
      if(samp.rich=="equip")
      {
        size.rand=sum(comm.b)
        samp.num=nrow(comm)
        sp.num=ncol(comm)
        prob.sp=as.vector(matrix(sp.freqence,nrow = samp.num,ncol=sp.num,byrow = TRUE))
        comr.b=matrix(0,nrow = samp.num,ncol=sp.num)
        id.out=sample(1:(samp.num*sp.num),size =size.rand,replace = FALSE,prob = prob.sp)
        comr.b[id.out]=1
      }else if(samp.rich=="prop"){
        samp.richness=rowSums(comm.b)
        size.rand=sum(comm.b)
        samp.num=nrow(comm)
        sp.num=ncol(comm)
        prob.sp1=matrix(sp.freqence,nrow = samp.num,ncol=sp.num,byrow = TRUE)
        prob.sp2=matrix(samp.richness,nrow = samp.num,ncol=sp.num,byrow = FALSE)
        prob.sp=as.vector(prob.sp1*prob.sp2)
        comr.b=matrix(0,nrow = samp.num,ncol=sp.num)
        id.out=sample(1:(samp.num*sp.num),size =size.rand,replace = FALSE,prob = prob.sp)
        comr.b[id.out]=1
      }else if(samp.rich=="fix"){
        samp.richness=rowSums(comm.b)
        samp.num=nrow(comm)
        sp.num=ncol(comm)
        comr.b=t(matrix(sapply(1:samp.num,
                        function(i)
                        {
                          res=rep(0,sp.num)
                          id=sample(1:sp.num,samp.richness[i],replace = FALSE,prob = sp.freqence)
                          res[id]=1
                          res
                        }),ncol=samp.num))
      }else{stop("Sample richness randomization method input error.")}
    }else if(sp.freq=="fix"){
      if(samp.rich=="equip")
      {
        sp.freqence=colSums(comm.b)
        sp.num=ncol(comm)
        samp.num=nrow(comm)
        comr.b=sapply(1:sp.num, function(i){res=rep(0,samp.num);id=sample(1:samp.num,sp.freqence[i]);res[id]=1;res})
      }else if(samp.rich=="prop"){
        sp.freqence=colSums(comm.b)
        samp.richness=rowSums(comm.b)
        sp.num=ncol(comm)
        samp.num=nrow(comm)
        comr.b=sapply(1:sp.num, function(i){res=rep(0,samp.num);id=sample(1:samp.num,sp.freqence[i],replace = FALSE,prob = samp.richness);res[id]=1;res})
      }else if(samp.rich=="fix"){
        if(!(swap.method %in% c("swap","tswap","quasiswap","backtrack")))
        {stop("Swap method input error. Swap method is necessary for fixed sp.freq and samp.rich.")}
        if(burnin==0){stop("burnin must be input for Swap method, indicating which round of randomizaiton is running.")}
        #library(vegan)
        nm<-vegan::nullmodel(comm.b,swap.method)
        sm<-stats::simulate(nm,burnin=burnin,nsim=1,thin=5)
        comr.b<-matrix(sm,nrow=nrow(comm.b),ncol=ncol(comm.b))
      }else{stop("Sample richness randomization method input error.")}
    }else{stop("Species frequency randomization method input error.")}
    
    # abundance assignment
    if(!(abundance %in% c("not","shuffle","local","region"))){stop("Abundance assignment method input error.")}
    if(abundance=="shuffle")
    {
      if(samp.rich=="fix")
      {
        out<-t(sapply(1:nrow(comr.b), function(i){res=rep(0,ncol(comr.b));xx=comm[i,comm[i,]>0];res[comr.b[i,]>0]=xx[sample.int(length(xx))];res}))
        rownames(out)<-rownames(comm.b)
        colnames(out)<-colnames(comm.b)
      }else{
        warning("The abundance assignment cannot be 'shuffle', thus change to per local.")
        abundance="local"
      }
    }
    
    if(abundance=="not")
    {
      out=comr.b
    }else if(abundance=="local"){
      samp.ab=rowSums(comm)
      samp.richness=rowSums(comr.b)
      samp.richness.old=rowSums(comm.b)
      prob.ab=t(sapply(1:nrow(comm),
                     function(i)
                     {
                       res=rep(0,ncol(comm))
                       ab.old=comm[i,which(comm[i,]>0)]
                       ab.rank=ab.old[order(ab.old,decreasing = TRUE)]
                       id.new=which(comr.b[i,]>0)
                       if(samp.richness[i]<=samp.richness.old[i])
                       {
                         res[id.new[sample.int(length(id.new))]]=ab.rank[1:length(id.new)]
                       }else{
                         res[id.new[sample.int(length(id.new))]]=c(ab.rank,rep(min(ab.rank),(length(id.new)-length(ab.rank))))
                       }
                       res
                     }))
      out=ab.assign(comm.b = comr.b,samp.ab = samp.ab,prob.ab = prob.ab)
    }else if(abundance=="region"){
      if(is.null(region.meta))
      {
        prob.ab=matrix(colSums(comm),nrow = nrow(comm),ncol = ncol(comm),byrow = TRUE)
      }else{
        if(is.null(names(region.meta)))
        {
          if(length(region.meta)!=ncol(comm)){stop("region.meta setting is wrong.")}else{meta.abs=region.meta}
        }else{
          meta.abs=rep(0,ncol(comm))
          region.meta=region.meta[which(names(region.meta) %in% colnames(comm))]
          meta.abs[match(names(region.meta),colnames(comm))]=region.meta
        }
        prob.ab=matrix(meta.abs,nrow = nrow(comm),ncol = ncol(comm),byrow = TRUE)
      }
      samp.ab=rowSums(comm)
      out=ab.assign(comm.b = comr.b,samp.ab = samp.ab,prob.ab = prob.ab)
    }
  }
  rownames(out)<-rownames(comm)
  colnames(out)<-colnames(comm)
  out
}