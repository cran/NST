beta.g<-function(comm,dist.method="bray",abundance.weighted=TRUE,
                 as.3col=FALSE,out.list=TRUE, transform.method=NULL, logbase=2)
{
  #library(vegan)
  if(length(abundance.weighted)<length(dist.method))
  {
    abundance.weighted=c(abundance.weighted,rep(abundance.weighted[length(abundance.weighted)],(length(dist.method)-length(abundance.weighted))))
  }else{
    abundance.weighted=abundance.weighted[1:length(dist.method)]
  }
  
  id.ab=which(abundance.weighted)
  dist.methodcr=dist.method
  dist.method[which(dist.method=="sorensen")]="bray"
  dist.method[which(dist.method=="ruzicka")]="jaccard"
  dist.methodcr[id.ab[which(dist.methodcr[id.ab]=="sorensen")]]="bray"
  dist.methodcr[id.ab[which(dist.methodcr[id.ab]=="jaccard")]]="ruzicka"
  id.uw=which(!abundance.weighted)
  dist.methodcr[id.uw[which(dist.methodcr[id.uw]=="bray")]]="sorensen"
  dist.methodcr[id.uw[which(dist.methodcr[id.uw]=="ruzicka")]]="jaccard"
  
  # comm data transform
  if(!is.null(transform.method))
  {
    if(inherits(transform.method,"function"))
    {
      comm=transform.method(comm)
    }else{
      requireNamespace("vegan")
      comm=vegan::decostand(comm,method = transform.method, logbase=logbase, na.rm = TRUE)
    }
  }
  
  # distance
  if(length(dist.method)==1)
  {
    if(dist.method=="mGower"){
      coms=comm
      coms[coms==0]=0.1
      coms=log10(coms)+1
      res=vegan::vegdist(coms,method = "altGower",binary = (!abundance.weighted),diag=TRUE,upper = TRUE,na.rm=TRUE)  
    }else if(dist.method=="mEuclidean"){
      ABJ=vegan::designdist(comm,method = "A+B-J",terms = "binary")
      res=vegan::vegdist(comm,method = "euclidean",binary = (!abundance.weighted),diag=TRUE,upper = TRUE,na.rm=TRUE)
      res=res/ABJ
    }else if(dist.method=="mManhattan"){
      ABJ=vegan::designdist(comm,method = "A+B-J",terms = "binary")
      res=vegan::vegdist(comm,method = "manhattan",binary = (!abundance.weighted),diag=TRUE,upper = TRUE,na.rm=TRUE)
      res=res/ABJ
    }else if(dist.method=="chao.jaccard"){
      res=chaojaccard(comm)
      abundance.weighted=TRUE
    }else if(dist.method=="chao.sorensen"){
      res=chaosorensen(comm)
      abundance.weighted=TRUE
    }else{
      res=vegan::vegdist(comm,method=dist.method,binary = (!abundance.weighted),diag = TRUE,upper = TRUE,na.rm = TRUE)
    }
    if(as.3col){res=dist.3col(res);colnames(res)[3]=dist.methodcr}
  }else{
    res=list()
    for(i in 1:length(dist.method))
    {
      if(dist.method[i]=="mGower"){
        coms=comm
        coms[coms==0]=0.1
        coms=log10(coms)+1
        res[[i]]=vegan::vegdist(coms,method = "altGower",binary = (!abundance.weighted[i]),diag=TRUE,upper = TRUE,na.rm=TRUE)  
      }else if(dist.method[i]=="mEuclidean"){
        ABJ=vegan::designdist(comm,method = "A+B-J",terms = "binary")
        res[[i]]=vegan::vegdist(comm,method = "euclidean",binary = (!abundance.weighted[i]),diag=TRUE,upper = TRUE,na.rm=TRUE)
        res[[i]]=res[[i]]/ABJ
      }else if(dist.method[i]=="mManhattan"){
        ABJ=vegan::designdist(comm,method = "A+B-J",terms = "binary")
        res[[i]]=vegan::vegdist(comm,method = "manhattan",binary = (!abundance.weighted[i]),diag=TRUE,upper = TRUE,na.rm=TRUE)
        res[[i]]=res[[i]]/ABJ
      }else if(dist.method[i]=="chao.jaccard"){
        res[[i]]=chaojaccard(comm,to.dist=FALSE)
        abundance.weighted[i]=TRUE
      }else if(dist.method[i]=="chao.sorensen"){
        res[[i]]=chaosorensen(comm,to.dist=FALSE)
        abundance.weighted[i]=TRUE
      }else{
        res[[i]]=vegan::vegdist(comm,method=dist.method[i],binary = (!abundance.weighted[i]),diag = TRUE,upper = TRUE,na.rm = TRUE)
      }
      if(as.3col|(!out.list)){res[[i]]=dist.3col(res[[i]])}
    }
    wtn=c("uw","wt")[abundance.weighted+1]
    names(res)=paste0(dist.methodcr,".",wtn)
    if(!out.list)
    {
      res=cbind(res[[1]][,1:2,drop=FALSE],sapply(1:length(res),function(i){res[[i]][,3,drop=FALSE]}))
      colnames(res)[3:ncol(res)]=paste0(dist.methodcr,".",wtn)
    }
  }
  res
}

# chaosorensen
chaosorensen<-function(comm,dissimilarity=TRUE,to.dist=TRUE)
{
  # function cj was copied from chao.sorensen in R package "fassil"
  # Vavrek, Matthew J. 2011. fossil: palaeoecological and palaeogeographical analysis tools. Palaeontologia Electronica, 14:1T. http://palaeo-electronica.org/2011_1/238/index.html
  cs<-function (x, y) 
  {
    n <- sum(x)
    m <- sum(y)
    if (length(x[y > 0 & x == 2]) == 0) 
      f2plus <- 1
    else f2plus <- length(x[y > 0 & x == 2])
    p1 <- sum(x[y > 0]/n)
    p2 <- ((m - 1)/m) * (length(x[y > 0 & x == 1])/(2 * f2plus))
    p3 <- sum(x[y == 0]/n)
    u <- p1 + p2 * p3
    if (u > 1) 
      u <- 1
    if (length(y[x > 0 & y == 2]) == 0) 
      fplus2 <- 1
    else fplus2 <- length(y[x > 0 & y == 2])
    q1 <- sum(y[x > 0]/m)
    q2 <- ((n - 1)/n) * (length(y[x > 0 & y == 1])/(2 * fplus2))
    q3 <- sum(y[x == 0]/m)
    v <- q1 + q2 * q3
    if (v > 1) 
      v <- 1
    if (u == 0 & v == 0) 
      c.s <- 0
    else c.s <- (2 * u * v)/(u + v)
    return(c.s)
  }
  
  res=matrix(NA,nrow=nrow(comm),ncol=nrow(comm))
  for(i in 1:(nrow(comm)-1))
  {
    for(j in (i+1):nrow(comm))
    {
      res[i,j]<-res[j,i]<-cs(comm[i,],comm[j,])
    }
  }
  diag(res)=1
  rownames(res)<-colnames(res)<-rownames(comm)
  if(dissimilarity) res=1-res
  if(to.dist) res=stats::as.dist(res)
  res
}

# chaojaccard
chaojaccard<-function(comm,dissimilarity=TRUE,to.dist=TRUE)
{
  # function cj was copied from chao.jaccard in R package "fassil"
  # Vavrek, Matthew J. 2011. fossil: palaeoecological and palaeogeographical analysis tools. Palaeontologia Electronica, 14:1T. http://palaeo-electronica.org/2011_1/238/index.html
  cj<-function (x, y) 
  {
    n <- sum(x)
    m <- sum(y)
    if (length(x[y > 0 & x == 2]) == 0) 
      f2plus <- 1
    else f2plus <- length(x[y > 0 & x == 2])
    p1 <- sum(x[y > 0]/n)
    p2 <- ((m - 1)/m) * (length(x[y > 0 & x == 1])/(2 * f2plus))
    p3 <- sum(x[y == 0]/n)
    u <- p1 + p2 * p3
    if (u > 1) 
      u <- 1
    if (length(y[x > 0 & y == 2]) == 0) 
      fplus2 <- 1
    else fplus2 <- length(y[x > 0 & y == 2])
    q1 <- sum(y[x > 0]/m)
    q2 <- ((n - 1)/n) * (length(y[x > 0 & y == 1])/(2 * fplus2))
    q3 <- sum(y[x == 0]/m)
    v <- q1 + q2 * q3
    if (v > 1) 
      v <- 1
    if (u == 0 & v == 0) 
      c.j <- 0
    else c.j <- (u * v)/(u + v - (u * v))
    return(c.j)
  }
  
  res=matrix(NA,nrow=nrow(comm),ncol=nrow(comm))
  for(i in 1:(nrow(comm)-1))
  {
    for(j in (i+1):nrow(comm))
    {
      res[i,j]<-res[j,i]<-cj(comm[i,],comm[j,])
    }
  }
  diag(res)=1
  rownames(res)<-colnames(res)<-rownames(comm)
  if(dissimilarity) res=1-res
  if(to.dist) res=stats::as.dist(res)
  res
}