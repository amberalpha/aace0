#' top-level derive method for covariance table
#'
#' loads data to global if not exists; loops through date estimating co, convert to data.table, save to rd
#'
#' Details (body)
#'
#' @param su securityuniverse object
#' @param da date of class Date
#' @param verbose flag for printing date
#' @param nfac number of factors
#' @param ... passed to cewrap
#' @family top level
#' @examples
#' \dontrun{
#' deraaco()
#' }
#' @export
#' @import data.table
#' @import zoo
deraaco <- function(su=getrdatv(type="su"),si=getrdatv(type="si"),verbose=TRUE,...) {
  if(! exists("prem.g",envir=globalenv())) getbdhgl()
  nfac <- getsi("nfac")
  da <- su[,sort(unique(date))]
  x <- vector("list",length(da))
  for(i in seq_along(da)) {
    if(verbose) print(da[i])
    x[[i]] <- data.table(cewrap(pa=getpa(su=su,da=da[i]),nfac=nfac,da=da[i],...))
  }
  putrdatv(rbindlist(x),type="co")
}

#' get global zoo panels
#'
#' loads timeseries required for co estimation
#'
#' @param field named list of derived panels to load; they are assigned in the global environment, with the name in the list
#' @family internal
#' @examples
#' \dontrun{
#' getbdhgl()
#' exists('prem.g') #true
#' }
#' @export
getbdhgl <- function(field=list(
                              list(prem.g="x0700redoto"),
                              list(mcap.g="x0702mcap"),
                              list(dipr.g="x0702dipr"),
                              list(erpr.g="x0702erpr"),
                              list(bopr.g="x0702bopr"),
                              list(best.g="x0702best"),
                              list(vix.g="x0502vix")),
                      myclass='zoo') {
  x <- lapply(field,function(field) {assign(x=names(field),value=getstep(unlist(field),n='001',myclass=myclass),envir=globalenv())})
}


#' return windowed zoo panel from global panel
#'
#' Description (body)
#'
#' Details (body)
#'
#' @param su securityuniverse object
#' @param da date of class Date, forms the datum for the lags in win
#' @param bui identifiers for inclusion
#' @param field name of global zoo panel from which to extract
#' @param win window normally negative (prior) lags referred to da
#' @family utility
#' @examples
#' \dontrun{
#' getpa()
#' }
#' @export
getpa <- function(su=su.g,da=su[,max(date)],bui=su[date==su[date<=da,max(date)],bui],field="prem.g",win=-(1:230)) {
  if(0==length(bui)*length(da)) {
    return(NULL)
  } else {
    get(field,envir=globalenv())[as.Date(offda(da,win)),bui]
  }
}

#' derive ce
#'
#' high-level covariance estimation
#'
#' takes panel of returns, applies normalisations, estimates ce, adds benchmark data, returns as dataframe
#'
#' @param pa zoo panel (normally returns)
#' @param normalise flag to redistribute - see aautil::zoonorm()
#' @param nfac # factors in PCA
#' @param applyvix flag to scale historical returns to 'today's vol'
#' @param mixvix bayes adjustment to vix rescaling, strength of shrinkage to mean
#' @param bench type of benchmark weights
#' @param ... passed to fms2
#' @examples
#' \dontrun{
#' cewrap()
#' }
#' @export
cewrap <- function (pa=getpa(), normalise="NONE", nfac=20, applyvix=TRUE, mixvix=.5, bench="equal",da=max(index(pa)),...) {
  if (applyvix) {
     vix <- na.locf(vix.g )
     vixinverse <- 1/(mixvix * vix + (1 - mixvix) *
                        rollapply(vix,6 * 52, mean, align = "right", na.rm = T, fill = 'extend'))
     pa <- sweep(pa, MARGIN = 1, FUN = "*",
                 STAT = vixinverse[index(pa),]/as.numeric(vixinverse[as.Date(da), ]))
  }
  if(normalise!="NONE") { pa <- zoonorm(pa,dimension=normalise) }
  fm <- fms2(
          pa,
          range.factors = c(nfac,nfac),
          ...)
  fm1 <- addbench(fm,mcap = mcap.g[match(as.Date(da), index(mcap.g)),,drop = FALSE], bench = bench)
  dfrce(fm1,da=da)
}


#' extract object of class 'ce' from a co table
#'
#'
#'
#' largely for historical reasons ce objects are converted to a dataframe - this function converts the dataframe back to a 'ce'
#'
#' @param cetab a data.frame or data.table
#' @param dat the date to extract
#' @family conversion
#' @examples
#' \dontrun{
#' require(aace)
#' vcvce(dtce(getrdatv("jo","co"))) #covariance matrices
#' }
#' @export
dtce <- function (cetab, dat = cetab[,max(date)])
{
  stopifnot(length(dat) == 1)
  #stopifnot(valda(dat) && length(dat) == 1)
  cetab <- data.frame(cetab[date==dat])
  jbui <- grep("bui",colnames(cetab))
  jloadings <- grep("loadings",colnames(cetab))
  jfmp <- grep("fmp",colnames(cetab))
  jhpl <- grep("hpl",colnames(cetab))
  jmethod <- grep("method",colnames(cetab))
  jfull <- grep("full",colnames(cetab))
  juniqueness <- grep("uniqueness",colnames(cetab))
  jsdev <- grep("sdev",colnames(cetab))
  jqua <- grep("qua",colnames(cetab))
  jmvp <- match("mvp",colnames(cetab))
  jmcp <- match("mcp",colnames(cetab))
  jevp <- match("evp",colnames(cetab))
  jbeta <- grep("beta",colnames(cetab))
  jpoco <- grep("poco",colnames(cetab))
  bui <- cetab[,jbui,drop=TRUE]
  attributes(bui) <- NULL
  ifull <- which(cetab[,jfull]=="1")
  factors <- 1:ncol(cetab[,jloadings,drop=FALSE])
  ans <- list(
    loadings=as.matrix(cetab[,jloadings,drop=FALSE]),
    fmp=as.matrix(cetab[,jfmp,drop=FALSE]),
    hpl=as.matrix(cetab[,jhpl,drop=FALSE]),
    method=as.matrix(cetab[,jmethod,drop=FALSE]),
    full=as.matrix(cetab[,jfull,drop=FALSE]),
    uniqueness=as.matrix(cetab[,juniqueness,drop=FALSE]),
    sdev=as.matrix(cetab[,jsdev,drop=FALSE]),
    qua=as.matrix(cetab[,jqua,drop=FALSE]),
    evp=as.matrix(cetab[,jevp,drop=FALSE]),
    mcp=as.matrix(cetab[,jmcp,drop=FALSE]),
    mvp=as.matrix(cetab[,jmvp,drop=FALSE]),
    beta=as.matrix(cetab[,jbeta,drop=FALSE]),
    poco=as.matrix(cetab[,jpoco,drop=FALSE]),
    weight=NULL,
    call="restored"
  )
  #change mode for logical
  mode(ans$full) <- "numeric"
  mode(ans$full) <- "logical"
  #label
  dimnames(ans$loadings) <- list(bui,psz("loadings",factors))
  dimnames(ans$fmp) <- list(bui,psz("fmp",factors))
  dimnames(ans$hpl) <- list(bui,psz("hpl",factors))
  dimnames(ans$method) <- list(bui,psz("method",factors))
  dimnames(ans$full) <- list(bui,"full")
  dimnames(ans$uniqueness) <- list(bui,"uniqueness")
  dimnames(ans$sdev) <- list(bui,"sdev")
  dimnames(ans$evp) <- list(bui,"evp")
  dimnames(ans$mcp) <- list(bui,"mcp")
  dimnames(ans$mvp) <- list(bui,"mvp")
  dimnames(ans$beta) <- list(bui,c("evp","mcp","mvp"))
  dimnames(ans$poco) <- list(bui,"poco")
  class(ans) <- "ce"
  #stopifnot(valce(ans))
  ans
}

#dfrce - convert ce to a dataframe for storage as co
#'
#'
#'
#' largely for historical reasons ce objects are converted to a dataframe co
#'
#' @param x an object of class ce; see package aace
#' @param dat the date to assign it
#' @family conversion
#' @examples
#' \dontrun{
#' require(aace)
#' dfrce(dtce(getrdatv("jo","co"))) #a round trip: get a stored co object, convert to ce for one date, convert back to co
#' }
#' @export
dfrce <- function(x,dat=max(index(x)))
{
  colnames(x$loadings) <- 1:ncol(x$loadings)  #fix colnames because data.frame prepends object name
  colnames(x$fmp) <- 1:ncol(x$fmp)
  colnames(x$hpl) <- 1:ncol(x$fmp)
  colnames(x$method) <- 1:ncol(x$method)
  res <- data.frame(
    date=dat,
    bui=as.matrix(rownames(x$loadings)),
    sdev=x$sdev,
    uniqueness=x$uniqueness,
    full=as.character(as.numeric(x$full)),
    loadings=x$loadings,
    fmp=x$fmp,
    hpl=x$hpl,
    method=metce(x),
    qua=x$qua,
    evp=x$evp,
    mcp=x$mcp,
    mvp=x$mvp,
    poco=x$poco,
    betaevp=x$beta[,"evp"],
    betamcp=x$beta[,"mcp"],
    betamvp=x$beta[,"mvp"]
  )
  names(res) <- gsub("\\.","",names(res))
  res
}



#internals----------------------------------------
#' @export
addbench <- function(
  ce,
  mcap=t(ce$loadings[,1,drop=FALSE]*NA),
  bench=c("equal","minvar","mcap")
)
{
  bench <- match.arg(bench)
  stopifnot(all(buice(ce)%in%colnames(mcap)))
  bui <- buice(ce)
  n <- length(bui)
  mcap <- coredata(mcap)[,bui,drop=TRUE] #nb no lagging is done to mcap, should already be lagged 1 period (ie opening mcap)
  if(mean(!is.na(mcap))<0.1) mcap <- mcap*NA      #insufficient data for cap-weighted
  if(mean(!is.na(mcap))>0) mcap[is.na(mcap)] <- min(mcap,na.rm=TRUE) #na assigned minimum weight
  mcap[mcap==0] <- min(mcap[mcap!=0],na.rm=TRUE) #0 assigned minimum weight
  vv <- vcvce(ce)$T
  stopifnot(all(colnames(vv)==colnames(mcap)))
  mvp <- solve(vv,rep(1,ncol(vv)))
  evp <- rep(1/ncol(vv),ncol(vv))
  mcp <- as.numeric(mcap/sum(mcap,na.rm=TRUE))
  mvp <- mvp/sum(mvp)
  evp <- evp/sum(abs(evp))
  ce$mvp <- matrix(mvp,n,1,dimn=list(bui,NULL))
  ce$evp <- matrix(evp,n,1,dimn=list(bui,NULL))
  ce$mcp <- matrix(mcp,n,1,dimn=list(bui,NULL))
  bpf <- switch(bench,
                equal=ce$evp,
                minvar=ce$mvp,
                mcap=ce$mcp
  )
  ce$beta <- cbind(
    vv%*%evp*1./as.numeric(((t(evp)%*%vv)%*%evp)),
    vv%*%mcp*1./as.numeric(((t(mcp)%*%vv)%*%mcp)),
    vv%*%mvp*1./as.numeric(((t(mvp)%*%vv)%*%mvp))
  )
  dimnames(ce$beta) <- list(bui,c("evp","mcp","mvp"))
  pol <- as.numeric(sign(t(ldgce(ce))%*%bpf))
  pol[which(pol==0)] <- 1
  ce$bench <- bench
  ce$fmp <- sweep(ce$fmp,STAT=pol,MAR=2,FUN="*")
  ce$hpl <- sweep(ce$hpl,STAT=pol,MAR=2,FUN="*")
  ce$loadings <- sweep(ce$loadings,STAT=pol,MAR=2,FUN="*")
  ce$poco <- matrix(0,n,1,dimn=list(bui,NULL))
  ce
}


#--------na interpolation section

# getzco - covariance-based (PCA) interpolator
#
# @export
# getzco <- function(co=getrd(100),part=c('z','psi')) {
#  part <- match.arg(part)
#  ce <- dtce(co)
#  if(part=='z') {
#  list(
#    M=fmpce(ce)[,1,drop=FALSE]%*%t(ldgce(ce)[,1,drop=FALSE]),
#    S=fmpce(ce)[,-1,drop=FALSE]%*%t(ldgce(ce)[,-1,drop=FALSE]),
#    T=fmpce(ce)%*%t(ldgce(ce))
#  )
#  } else {
#    list(
#      M=fmpce(ce)[,1,drop=FALSE],
#      S=fmpce(ce)[,-1,drop=FALSE],
#      T=fmpce(ce)
#    )
#  }
# }


#both the interp functions could/should iterate which would be easy with a loop on penultimate line and no paz

#' interpolates using pca
#'
#' single step interpolation of M,S components using co, starting from 0
#' @param co covariance object
#' @param pa panel with NA
#' @param comp text flag M/S/T for component to use
#' @export
interpce <- function(co=getrdatv("jo","co"),pa=getpa(su=getrdatv("jo","su")),comp=c('T','M','S')) {
  comp <- match.arg(comp)
  ce <- dtce(co)
  stopifnot(all(sort(colnames(pa))==sort(buice(ce))))
  i <- is.na(coredata(pa))
  paz <- pa
  coredata(paz)[i] <- 0 #this does not affect the result if pa has same na as estimation window
  coredata(pa)[i] <- coredata(mscecomp(x=ce,ret=paz)[[comp]])[i]
  pa
}


#' interpolates using regression on te
#'
#' single step interpolation of M,S components using co, starting from 0
#' @param pa panel
#' @param lo a 'loadings' object, named list of z, p, l (components postmult, phi, lambda)
#' @export
interpte <- function(pa=getpa(su=getrdatv("jo","su")),lo=getlote()$lo) {
  i <- is.na(coredata(pa))
  coredata(pa)[i] <- 0
  coredata(pa)[i] <- coredata(pa%*%lo$z)[i]
  pa
}

#' generate covariance-based (PCA) interpolator
#'
#' returns three 'loadings' objects named lo, l, p
#' @param co covariance object
#' @export
getloco <- function(co=getrd(100)) {
  ce <- dtce(co)
  fmp <- fmpce(ce)
  ldg <- ldgce(ce)
  list(
    lo=list(
        z=fmp%*%t(ldg),
        l=ldg,
        p=fmp
      ),
    M=list(
      z=fmp[,1,drop=FALSE]%*%t(ldg[,1,drop=FALSE]),
      l=ldg[,1,drop=FALSE],
      p=fmp[,1,drop=FALSE]
    ),
    S=list(
      z=fmp[,-1,drop=FALSE]%*%t(ldg[,-1,drop=FALSE]),
      l=ldg[,-1,drop=FALSE],
      p=fmp[,-1,drop=FALSE]
    )
  )
}

# getzte - industry-based interpolator
#
# the te object is read from rd, and has normally been read back from directories
# this works but is cumbersome/not needed so has been replaced by getzte that uses autoprune, see below
# getzte <- function(te=getrd(103),su=getrdatv("jo","su",2),da=su[,max(date)],loocv=FALSE) {
#   buix <- su[date==da,unique(bui)]
#   te <- setkey(te[buix,list(bui,bcode,BICS_REVENUE_PERC_LEVEL_ASSIGNED)],bui)
#   setkey(te[,BRPLA:=sum(BICS_REVENUE_PERC_LEVEL_ASSIGNED),list(bui,bcode)],bui,bcode)[,BICS_REVENUE_PERC_LEVEL_ASSIGNED:=NULL]
#   te <- unique(te)
#   x <- tabtomat(data.frame(te))
#   x[is.na(x)] <- 0
#   stopifnot(all(99<apply(x,1,sum)) & all(apply(x,1,sum)<101))
#   sol <- list(T=x%*%solve(crossprod(x))%*%t(x))
#   if(loocv) {
#     for( i in 1:nrow(x) ) {
#       sol[[i+1]] <- x[-i,]%*%solve(crossprod(x[-i,]))%*%t(x)
#     }
#     names(sol)[2:(1+nrow(x))] <- paste0('x',1:nrow(x))
#   }
#   sol
# }


#' integer industries prune
#'
#' returns data.table with cols: bui bcode BRPLA
#' @param su security universe
#' @param da date
#' @param nmin min # per node
#' @export
pruneztei <- function(su=getrdatv("jo","su",2),da=su[,max(date)],nmin=3) {
  buix <- su[date==da,unique(bui)]
  te <- getbdp()[buix][,list(bui,bcode=BICS_LEVEL_CODE_ASSIGNED,BRPLA=1)]
  tesums <- setcolorder(setkey(te[,list(agg=sum(BRPLA)),bcode],bcode)[,startcode:=bcode][,bagg:=agg],c('startcode','agg','bcode','bagg'))
  clen <- tesums[,max(nchar(bcode))]
  while(0<clen) {
    tesums[,acode:=bcode][((bagg<nmin)&(nchar(bcode)==clen)),acode:=substr(bcode,1,nchar(bcode)-2)]
    tesums[,bcode:=acode][,acode:=NULL][,bagg:=NULL][,bagg:=sum(agg),bcode]
    clen <- clen-2
  }
  te <- setkey(te,bcode)[setkey(tesums,startcode)][,bcode:=i.bcode][,i.bcode:=NULL][,agg:=NULL][,bagg:=NULL][bcode=="",bcode:="00"]
  te[is.na(bcode),bcode:='1']
  te
}

#' fractional industries prune
#'
#' returns data.table with cols: bui bcode BRPLA
#' currently silently drops bui that have no fractional info - this needs a rethink
#' @param su security universe
#' @param da date
#' @param nmin min # per node
#' @param nmin min aggregate weight per node
#' @export
pruneztef <- function(su=getrdatv("jo","su",2),da=su[,max(date)],wmin=2,nmin=4) {
  buix <- su[date==da,unique(bui)]
  #next line has allow=TRUE which copes with empty files
  te <-  getbdp(mnem=data.table(data.frame(field = c("BICS_REVENUE_PERC_LEVEL_ASSIGNED"))))[buix,allow=TRUE][,list(bui,bcode=BICS_LEVEL_CODE_ASSIGNED,BRPLA=BICS_REVENUE_PERC_LEVEL_ASSIGNED/100)][!is.na(BRPLA)]
  tesums <- setcolorder(setkey(te[,list(agg=sum(BRPLA),count=.N),bcode],bcode)[,startcode:=bcode][,bagg:=agg][,bcount:=count],c('startcode','agg','count','bcode','bagg','bcount'))
  clen <- tesums[,max(nchar(bcode))]
  while(0<clen) {
    tesums[,acode:=bcode][(  ((bagg<wmin)|(bcount<nmin))  &(nchar(bcode)==clen)),acode:=substr(bcode,1,nchar(bcode)-2)]
    tesums[,bcode:=acode][,acode:=NULL][,bagg:=NULL][,bagg:=sum(agg),bcode][,bcount:=sum(count),bcode]
    clen <- clen-2
  }
  te <- setkey(te,bcode)[setkey(tesums,startcode)][,bcode:=i.bcode][,i.bcode:=NULL][,agg:=NULL][,bagg:=NULL][,count:=NULL][,bcount:=NULL][bcode=="",bcode:="00"]
  setkey(te[,list(BRPLA=sum(BRPLA)),list(bui,bcode)],bui,bcode)[]
}

###REPLACED WITH getlote but different o/p structure
# z for integer industries pruned to nmin/node
#
# for part=psi, returns the postmultiplication matrix for factor scores; for part=z it is the matrix for 'fit'
# getztei <- function(su=getrdatv("jo","su",2),da=su[,max(date)],loocv=FALSE,nmin=3,wmin=3,type=c('i','f'),part=c('z','psi')) {
#   type <- match.arg(type)
#   part <- match.arg(part)
#   if(type=='i') {
#     te <- pruneztei(su=su,da=da,nmin=nmin)
#   } else {
#     te <- pruneztef(su=su,da=da,nmin=nmin,wmin=wmin)
#   }
#   x <- tabtomat(data.frame(te))
#   x[is.na(x)] <- 0
#   if(part=='z') {
#     postx <- t(x)
#   } else {
#     postx <- diag(ncol(x))
#   }
#   sol <- list(T=x%*%solve(crossprod(x))%*%postx)
#   if(loocv) {
#     for( i in 1:nrow(x) ) {
#       sol[[i+1]] <- x[-i,]%*%solve(crossprod(x[-i,]))%*%postx
#     }
#     names(sol)[2:(1+nrow(x))] <- paste0('x',1:nrow(x))
#   }
#   sol
# }

#' loadings object from industry object
#'
#' @param te industries
#' @param loocv flac to drop
#' @export
getlote <- function(te=pruneztei(),loocv=FALSE) {
  x <- tabtomat(data.frame(te))
  x[is.na(x)] <- 0
  mylo <- function(x,i=NULL) { #local function to drop identifiers from x
    if(!is.null(i)) { xl <- x[-i,,drop=FALSE] } else { xl <- x }
    list(
      z = xl%*%solve(crossprod(xl))%*%t(x),
      p = xl%*%solve(crossprod(xl)),
      l = x
    )
  }
  sol <- list(lo=mylo(x))
  if(loocv) {
    for( i in 1:nrow(x) ) {
      sol[[i+1]] <- mylo(x,i)
    }
    names(sol)[2:(1+nrow(x))] <- paste0('lo',1:nrow(x))
  }
  sol
}

# #' returns fitted and loocv fit on pa, using pca
# loocvi <- function(pa=getpa(su),...) {
# #  mhatijT <- mhatijS <- mhatijM <- mhatiT <- mhatiS <- mhatiM <- mhatT <- mhatS <- mhatM <- m <- pa
#   mhatijT <- mhatijS <- mhatijM <- mhatiT <- mhatiS <- mhatiM <- m <- pa
#   for(i in 1:nrow(m)) {
#     mi <- m[-i,]
#     ce <- dtce(data.table(cewrap(pa=mi,...)))
#     #this section identical to mhati, but harder to adapt for ij
# #     mscec <- mscecomp(ce,m[i,,drop=FALSE])
# #     mhatM[i,] <- mscec$M
# #     mhatS[i,] <- mscec$S
# #     mhatT[i,] <- mscec$T
#     #
#     fmp <- ce$fmp/as.numeric(ce$sdev)
#     ldg <- ce$loadings*as.numeric(ce$sdev)
#     sco <- mz(coredata(pa)[i, fulce(ce), drop = FALSE] %*% fmp[fulce(ce),, drop = FALSE])
#     mhatiM[i,] <- mz(sco[, 1, drop = FALSE] %*% t(ldg[, 1, drop = FALSE]))
#     mhatiS[i,] <- mz(sco[,-1, drop = FALSE] %*% t(ldg[,-1, drop = FALSE]))
#     mhatiT[i,] <- mz(sco[,  , drop = FALSE] %*% t(ldg[,  , drop = FALSE]))
#     zdrop <- function(fmp,ldg,j=c("M","S","T")) {
#       jj <- switch(match.arg(j),M=1,S=2:ncol(fmp),T=1:ncol(fmp))
#       z <- fmp[,jj,drop=F]%*%t(ldg[,jj,drop=F])
#       dd <- diag(z)
#       iinv <- which((1e-10<dd) & (dd<1))
#       scal <- rep(1,length(dd))
#       scal[iinv] <- 1/(1-dd[iinv])
#       diag(z) <- 0
#       sweep(z,MAR=2,FUN="*",STAT=scal)
#     }
#     #zd <- zdrop(fmp,ldg,"S")
#     mhatijM[i,] <- m[i, fulce(ce), drop = FALSE] %*% (zdrop(fmp,ldg,"M")[fulce(ce),])
#     mhatijS[i,] <- m[i, fulce(ce), drop = FALSE] %*% (zdrop(fmp,ldg,"S")[fulce(ce),])
#     mhatijT[i,] <- m[i, fulce(ce), drop = FALSE] %*% (zdrop(fmp,ldg,"T")[fulce(ce),])
#   }
# #   mhatlist <- lapply(list(M=mhatM,S=mhatS,T=mhatT),as.numeric)
#   mhatilist <- lapply(list(M=mhatiM,S=mhatiS,T=mhatiT),as.numeric)
#   mhatijlist <- lapply(list(M=mhatijM,S=mhatijS,T=mhatijT),as.numeric)
#   mfitlist <- lapply(mscecomp(dtce(data.table(cewrap(pa=m,...))),m),as.numeric)
# #   names(mhatlist) <- paste0('mhat',names(mscec))
#   names(mhatilist) <- paste0('mtwiddlei',c('M','S','T'))
#   names(mhatijlist) <- paste0('mtwiddleij',c('M','S','T'))
#   names(mfitlist) <- paste0('mhat',c('M','S','T'))
# # c(list(act=as.numeric(m)),mhatlist,mhatilist,mhatijlist,mfitlist)
#   c(list(act=as.numeric(m)),mhatilist,mhatijlist,mfitlist)
# }



#' @export
fmswrap <- function(x,
                    lambda=getp(sn='co',pn='lambda'),
                    k=getp(sn='co',pn='k'),
                    shrinkb=getp(sn='co',pn='shrinkb')) {
  fms2(x=x,lambda=lambda,range.factors=c(k,k),shrinkb=shrinkb,weight=seq(1,1,length=nrow(x)))
}





#' returns fitted and loocv fit on pa, using pca
#'
#' drops each row in turn and estimates new co; fits the row - this only makes sense for ce and is the only version that makes sense for ce
#' @param pa panel
#' @param ... passed to cewrap and thence to fms2
#' @export
loocviFun <- function(zd=gett('zd'),...) {
  m <- zd
  mhatijMS <- mhatijS <- mhatijM <- mhatiMS <- mhatiS <- mhatiM <- m*NA
  for(i in 1:nrow(m)) {
    mi <- m[-i,]
    ce <- fmswrap(x=mi,...)
    fmp <- ce$fmp/as.numeric(ce$sdev)
    ldg <- ce$loadings*as.numeric(ce$sdev)
    sco <- mz(coredata(zd)[i, fulce(ce), drop = FALSE] %*% fmp[fulce(ce),, drop = FALSE])
    mhatiM[i,] <- mz(sco[, 1, drop = FALSE] %*% t(ldg[, 1, drop = FALSE]))
    mhatiS[i,] <- mz(sco[,-1, drop = FALSE] %*% t(ldg[,-1, drop = FALSE]))
    mhatiMS[i,] <- mz(sco[,  , drop = FALSE] %*% t(ldg[,  , drop = FALSE]))
    mhatijM[i,] <- m[i, fulce(ce), drop = FALSE] %*% (zdrop(fmp,ldg,"M")[fulce(ce),])
    mhatijS[i,] <- m[i, fulce(ce), drop = FALSE] %*% (zdrop(fmp,ldg,"S")[fulce(ce),])
    mhatijMS[i,] <- m[i, fulce(ce), drop = FALSE] %*% (zdrop(fmp,ldg,"MS")[fulce(ce),])
  }
  ce <- fmswrap(x=mi,...)
  #cen <- fmswrap(x=rollmeanr(mi,k=5),...)
  mfit  <- lapply(lapply(lapply(mscecomp(ce ,m),coredata),data.table,keep.rownames=T),melt,id='rn')
  #mfitn <- lapply(lapply(lapply(mscecomp(cen,m),coredata),data.table,keep.rownames=T),melt,id='rn')
  ms <- setkey(rbind(
    melt(data=data.table(coredata(m),keep.rownames=T),id='rn')[,comp:='T'][,xv:='no'][,nper:=1],
    melt(data=data.table(coredata(mhatiM),keep.rownames=T),id='rn')[,comp:='M'][,xv:='i'][,nper:=1],
    melt(data=data.table(coredata(mhatiS),keep.rownames=T),id='rn')[,comp:='S'][,xv:='i'][,nper:=1],
    melt(data=data.table(coredata(mhatiMS),keep.rownames=T),id='rn')[,comp:='MS'][,xv:='i'][,nper:=1],
    melt(data=data.table(coredata(mhatijM),keep.rownames=T),id='rn')[,comp:='M'][,xv:='ij'][,nper:=1],
    melt(data=data.table(coredata(mhatijS),keep.rownames=T),id='rn')[,comp:='S'][,xv:='ij'][,nper:=1],
    melt(data=data.table(coredata(mhatijMS),keep.rownames=T),id='rn')[,comp:='MS'][,xv:='ij'][,nper:=1],
    mfit$S[,comp:='S'][,xv:='no'][,nper:=1],
    mfit$M[,comp:='M'][,xv:='no'][,nper:=1],
    mfit$T[,comp:='MS'][,xv:='no'][,nper:=1]#,
    #mfitn$S[,comp:='S'][,xv:='no'][,nper:=5],
    #mfitn$M[,comp:='M'][,xv:='no'][,nper:=5],
    #mfitn$T[,comp:='MS'][,xv:='no'][,nper:=5]
  ),rn,variable)
  totl <- melt(data=data.table(coredata(m),keep.rownames=T),id='rn')[,list(rn,variable,Total=value)]
  setkey(ms,rn,variable)
  setkey(totl,rn,variable)
  loocvid <- setnames(totl[ms],old=c('rn','variable'),new=c('date','ticker'))
  putt(loocvid)
}

#' @export
zdrop <- function(fmp,ldg,j=c("M","S","MS"),dodrop=T) {
  jj <- switch(match.arg(j),M=1,S=2:ncol(fmp),MS=1:ncol(fmp))
  z <- fmp[,jj,drop=F]%*%t(ldg[,jj,drop=F])
  if(dodrop) {
    dd <- diag(z)
    iinv <- which((1e-10<dd) & (dd<1))
    scal <- rep(1,length(dd))
    scal[iinv] <- 1/(1-dd[iinv])
    diag(z) <- 0
    z <- sweep(z,MAR=2,FUN="*",STAT=scal)
  }
  z
}

#calc the z matrix (fit operator) and optionally substitute out the LHS from the RHS
#' @export
zdropFun <- function(zd=gett('zd'),comp=c('M','S','MS'),dodrop=T,nper=1,...) {
  nper <- as.numeric(nper)
  comp <- match.arg(comp)
  if(nper>1) { zd <- rollmeanr(zd,k=nper) }
  ce <- fmswrap(x=zd,...)
  #ce <- fms2(x=zd,...)
  fmp <- ce$fmp/as.numeric(ce$sdev)
  ldg <- ce$loadings*as.numeric(ce$sdev)
  zdropd <- zdrop(fmp=fmp,ldg=ldg,j=comp,dodrop=dodrop)
  putt(zdropd)
}

#' returns fitted and loocv fit on pa, using te
#'
#' drops each column in turn of pa, ie each identifier bui - note there is NO iteration, no substitutio of the 'left out', it is just fitted which is more correct
#' @param pa panel
#' @param lo loadings object
#' @export
loocvj <- function(pa=getpa(su),lo=getlote(loocv=TRUE)) {
  #stopifnot(all(unlist(lapply(lapply(lapply(data.frame(t(data.frame(lapply(z,colnames)))),duplicated),"[",-1),all))))
  stopifnot(all.equal(colnames(pa),colnames(lo$lo$z)))
  mhat <- m <- coredata(pa)
  for(j in 1:ncol(m)) { #drop each column (identifier bui) in turn and postmultiply by the (x'x)-1.X
    mj <- m[,-j]
    xxj <- lo[[paste0('lo',j)]][['z']]
    mhat[,j] <- mj%*%(xxj[,j,drop=FALSE])
  }
  mfit <- m%*%lo[['lo']][['z']]
  list(act=as.numeric(m),mhat=as.numeric(mhat),mfit=as.numeric(mfit))
}

#' returns fitted and loocv fit on pa, using te
#'
#' iterloocv - drops each observation in turn and replaces it with an iterated
#' @param pa panel
#' @param lo loadings object
#' @param niter # iterations
#' @param myfun text name of function to apply to actual, fit, hat
#' @export
iterloocv <- function(pa=getpa(su),lo=getloco(),niter=2,myfun=c("as.numeric","mattotab","as.matrix")) {
  myfun <- get(match.arg(myfun))
  #stopifnot(all(unlist(lapply(lapply(lapply(data.frame(t(data.frame(lapply(z,colnames)))),duplicated),"[",-1),all))))
  stopifnot(all.equal(colnames(pa),colnames(lo$lo$z)))
  mhat <- m <- coredata(pa)
  comp <- names(lo)
  mfitlist <- mhatlist <- NULL
  for(i in 1:nrow(m)) {
    mi <- mean(m[i,])
    for(j in 1:ncol(m)) {
      mrow <- m[i,,drop=FALSE]
      mrow[,j] <- mi
      for(l in 1:niter) {mrow[1,j] <- (mrow%*%lo$lo$z)[1,j]}
      mhat[i,j] <- mrow[,j]
    }
  }
  for(k in seq_along(comp)) {
    mhatlist[[comp[k]]] <- mhat%*%lo[[comp[k]]]$z
    mfitlist[[comp[k]]] <- m%*%lo[[comp[k]]]$z
  }
  names(mhatlist) <- paste0('mhat',comp)
  names(mfitlist) <- paste0('mfit',comp)
  c(list(act=lapply(list(m),myfun)[[1]]),lapply(mhatlist,myfun),lapply(mfitlist,myfun))
}

# ilcvsumm <- function(x=iterloocv(...),...) {
#   s1 <- summary(lm(act ~ fit,data=data.frame(x)))
#   s2 <- summary(lm(act ~ loocv,data=data.frame(x)))
#   tab <- matrix(NA,5,2,dimnames=list(c('fit','loocv','diff','loocvM','loocvS'),c('Rsq','coef')))
#   tab[1,1] <- s1$r.squared
#   tab[1,2] <- s1$coefficients[2,1]
#   tab[2,1] <- s2$r.squared
#   tab[2,2] <- s2$coefficients[2,1]
#   tab[3,] <- tab[1,]-tab[2,]
#   if('loocvM'%in%names(x)) {
#     s3 <- summary(lm(act ~ loocvM,data=data.frame(x)))
#     s4 <- summary(lm(act ~ loocvS,data=data.frame(x)))
#     tab[4,1] <- s3$r.squared
#     tab[4,2] <- s3$coefficients[2,1]
#     tab[5,1] <- s4$r.squared
#     tab[5,2] <- s4$coefficients[2,1]
#   }
#   tab
# }


#THE FOLLOWING REPLACED WITH aasdlpkg::lcv4
#' summary of loadings object
#'
# ilcvsFun <- function(loocvid=gett('loocvid'),...) {
#   co <- setkey(loocvid[comp!='T',list(coef=summary(lm(Total~value))$coefficients[2,1]),'comp,xv,nper'],comp,xv,nper)
#   r2 <- setkey(loocvid[comp!='T',list(r2=summary(lm(Total~value))$r.squared),'comp,xv,nper'],comp,xv,nper)
#   ilcvsd <- co[r2]
#   putt(ilcvsd)
# }

# #summarises correlation of 'start' with final iteration
# iterate0 <- function(n=10,niter=5,FUN=mean,...) {
#   suppressWarnings(do.call(FUN,list(unlist(lapply(lapply(lapply(1:n,FUN=iterate1,niter=niter,...),FUN=cor),FUN=`[`,1,niter)))))
# }
#
# #wrapper to iterate2, applies 'drop' ie sets NA a specified fraction of data
# iterate1 <- function(seed=1,pa=getbdh(su),z=getzco()$T,idropfraction=0,initial=mean(pa,na.rm=TRUE),niter=5) {
#   m <- coredata(pa)
#   if(0<idropfraction) {
#     mask <- m-m+1
#     idrop <- round(length(mask)*idropfraction)
#     mask[1:idrop] <- NA
#     for(i in 1:nrow(mask)) mask[i,] <- mask[i,sample(1:ncol(m),size=ncol(m))]
#     m <- m*mask
#   }
#   ina <- which(is.na(m))
#   if(length(ina)==0) {
#     return()
#   } else {
#     x <- iterate2(m,z,initial,niter)
#   }
#   if(0<idropfraction) x[1:idrop,'start'] <- coredata(pa)[ina]
#   x
# }
#
# #iterate2 - inner function applying z
# iterate2 <- function(m,z,initial,niter=5) {
#   ina <- which(is.na(m))
#   res <- matrix(0*NA,length(ina),niter+2,dimnames=list(rownames(m)[ina],c('start',c(as.character(0:niter)))))
#   res[,"0"] <- initial
#   m[ina] <- initial
#   for(i in 1:niter) {
#     m[ina] <- (m%*%z)[ina]
#     res[,as.character(i)] <- m[ina]
#   }
#   res
# }

#' script in a function
#'
#' @export
runco <- function(nmin=2:14,wmin=2:10) {
  require(aapa)
  require(aate)
  aatopselect('test')
  getbdhgl()
  su <- getrd(99)
  da <- su[,max(date)]
  co <- getrd(150)
  pa <- getpa(su)
  x <- vector("list")

  x[[length(x)+1]] <- ilcvsumm(iterloocv(pa,lo=getloco()))[,nmin:=NA][,wmin:=NA][,iterated:="iterated"][,model:="PCA"]
  x[[length(x)+1]] <- ilcvsumm(loocvi(pa))[,nmin:=NA][,wmin:=NA][,iterated:="non-iterated"][,model:="PCA"]

  #fractional - nmin
  for(i in seq_along(nmin)) {
    print(i)
    loit <- getlote(te=pruneztef(su=su,da=da,nmin=nmin[i]))
    lonit <- getlote(te=pruneztef(nmin=nmin[i],wmin=min(wmin)),loocv=T)
    x[[length(x)+1]] <- ilcvsumm(iterloocv(pa,lo=loit))[,nmin:=nmin[i]][,wmin:=min(wmin)][,iterated:="iterated"][,model:="fractional"]
    x[[length(x)+1]] <- ilcvsumm(loocvj(pa,lo=lonit))[,nmin:=nmin[i]][,wmin:=min(wmin)][,iterated:="non-iterated"][,model:="fractional"]
  }
  #fractional - wmin
  for(i in seq_along(wmin)) {
    print(i)
    loit <- getlote(te=pruneztef(su=su,da=da,nmin=min(nmin),wmin=wmin[i]))
    lonit <- getlote(te=pruneztef(nmin=min(nmin),wmin=wmin[i]),loocv=T)
    x[[length(x)+1]] <- ilcvsumm(iterloocv(pa,lo=loit))[,nmin:=min(nmin)][,wmin:=wmin[i]][,iterated:="iterated"][,model:="fractional"]
    x[[length(x)+1]] <- ilcvsumm(loocvj(pa,lo=lonit))[,nmin:=min(nmin)][,wmin:=wmin[i]][,iterated:="non-iterated"][,model:="fractional"]
  }
  #integer
  for(i in seq_along(nmin)) {
    print(i)
    loit <- getlote(te=pruneztei(su=su,da=da,nmin=nmin[i]))
    lonit <- getlote(te=pruneztef(nmin=nmin[i]),loocv=T)
    x[[length(x)+1]] <- ilcvsumm(iterloocv(pa,lo=loit))[,nmin:=nmin[i]][,wmin:=NA][,iterated:="iterated"][,model:="integer"]
    x[[length(x)+1]] <- ilcvsumm(loocvj(pa,lo=lonit))[,nmin:=nmin[i]][,wmin:=NA][,iterated:="non-iterated"][,model:="integer"]
  }

  putrdatv(body(runco),app='ilcv',type=length(x))
  putrdatv(x,app='ilcv',type='nmin')
}

#' accessor/filter for 'co' results
#'
#' @param x named list defining output filter
#' @export
getco <- function(x=list(model='fractional',iterated='non-iterated',nmin=8,component='mhat')) {
  co <- setkeyv(rbindlist(getrdatv(app='ilcv',type='nmin')),names(x))
  co[x]
}

#fitplot - in a full panel that contains NA, drop size points from the non-NA, and crossplot fitted and actual
#' @export
fitplot <- function(size=min(1000,0.1*length(iok)),vbl=list(
  book.price=zoonorm(log(1+getpa(su=su,fi='bopr.g')^-1)),
  upside=zoonorm(log(getpa(su=su,fi='best.g')))),
  comp=c("M","S","T"),
  su=getrdatv("jo","su",v=2),
  ...) {
  #par(mfrow=c(1,1))
  for(i in seq_along(vbl)) {
    for(j in seq_along(comp)) {
      set.seed(123)
      print(i)
      pa <- vbl[[i]]
      iok <- which(!is.na(coredata(pa)))
      izap <- sample(iok,siz=size,rep=F)
      act <- coredata(pa)[izap]
      coredata(pa)[izap]<-NA
      xint <- interpce(pa=pa,comp=comp[j]) #needs co... not fixed
      int <- coredata(xint)[izap]
      stopifnot(all(abs(act - int)>1e-10))
      lmx <- summary(lm(act~int))
      sub <- paste0('R-squared ',round(lmx$r.squared,2),' coeff: ',round(lmx$coefficients[2,1],2),' tstat: ',round(lmx$coefficients[2,3],1))
      plot(int,act,pch=20,col='grey',ylab='actual',xlab='fitted',main=paste0(names(vbl)[i]," (",comp[j],")"),sub=sub,...)
      grid()
      if(comp[j]=="T") {abline(0,1)}
    }
  }
}

#'plot stocks to illustrate yield fit via PCA
#'
#' @export
examplebui <- function(co=getrd(100),iselect=4) {
  x0 <- getrd(100)[order(uniqueness),list(bui,uniqueness,rank(uniqueness))]
  y1 <- getrd(100)[order(uniqueness),list(bui,uniqueness,rank(uniqueness))]
  x1 <- y1[,bui]
  y2 <- getbdp(mnem=data.table(data.frame(field = c("BICS_REVENUE_PERC_LEVEL_ASSIGNED"))))[x0[,bui]][,.N,bui][order(-N)]
  x2 <- y2[,bui]
  #up : large magnitude, small residual
  pa1 <- pa0 <- zoonorm(log(getpa(su=su,fi='best.g')))
  for(j in 1:ncol(pa)) {
    pax <- pa0
    pax[,j] <- NA
    pa1[,j] <- interpce(pa=pax)[,j] #needs co... not fixed
  }
  y3 <- sort(apply(abs(pa0-pa1),2,mean,na.rm=T))
  x3 <- names(y3) #this excludes many
  y4 <- sort(apply(-pa0,2,mean,na.rm=T))
  x4 <- names(y4)
  #f1 loading large
  y5 <- getrd(100)[order(-loadings1),list(bui,loadings1,rank(loadings1))]
  x5 <- y5[,bui]
  dfx <- data.frame(x2,x3,x4,x5)
  #bui <- Reduce(intersect,dfx[1:40,1:3])
  dfy1 <- setkey(y1[,list(bui,uniqueness)],bui)
  dfy2 <- setkey(y2,bui)
  dfy3 <- setkey(data.table(bui=names(y3),upresid=y3),bui)
  dfy4 <- setkey(data.table(bui=names(y4),up=-y4),bui)
  dfy5 <- setkey(y5[,list(bui,loadings1)],bui)
  bui <- Reduce(intersect,dfx[1:62,1:3])
  bdp <- getbdp()
  rep <- setkey(dfy1[dfy2[dfy3[dfy4[dfy5]]]][bui],bui)
  toptab <- bdp[,list(bui,NAME)][rep]
  par(mfrow=c(2,2))
  for(i in 1:4) {
    bui0 <- toptab[i,bui]
    plot(cbind(pa0[,bui0],pa1[,bui0]),scr=1,col=1:2,main=bdp[bui0,NAME],ylim=c(-.1,.3),ylab='upside',xlab="")
    grid()
  }
  par(mfrow=c(1,1))
  bui0 <- toptab[iselect,bui]
  setnames(setkeyv(getbdp()[,list(bui,NAME)][setkey(data.table(fmpce(dtce(co))%*%t(ldgce(dtce(co))[bui0,,drop=F]),bui=buice(dtce(co))),bui)],bui0)[],old=bui0,new=getbdp()[bui0,NAME])[]
}


colm <- function(pa,z,tau=-2:2) {
  nona <- !is.na(coredata(pa))
  x <- iterloocv(pa=pa,z=z,myfun='as.m')
  x1 <- x2 <- pa
  for(i in seq_along(tau)) {
    x1 <- x1+lag(x2,k=tau[i],na.pad=TRUE)
  }
  x3 <- data.table(mattotab(x1))[!is.na(field)]

#get valid M,S,R,bui,t and same for t+tau
#regress
#T on lagged MS, delta MS, lagged R
#MS on lagged MS and delta MS
#R on lagged R
}

#' interpolation of missing
#'
#' @param pa panel
#' @param co covariance
#' @param te industry
#' @param phis ar1 on systematic components
#' @param phir ar1 on residual
#' @param type i/f/p string for integer, fractional, pca
#' @export
arco <- function(pa=getpa(),co=getrdatv("jo","co"),te=pruneztei(),phis=0,phir=0,type=c('i','f','p')) {
  type <- match.arg(type)
  if(type=='i') {
    lo <- getlote(te=te)$lo
    pai <- interpte(pa=pa,lo=lo)
  } else if(type=='f') {
    lo <- getlote(te=te)$lo
    pai <- interpte(pa=pa,lo=lo)
  } else if(type=='p') {
    lo <- getloco(co=co)$lo
    pai <- interpce(co=co,pa=pa)
  }
  fit00 <- blend <- fit <- pa*NA
  lo$l[is.na(lo$l)] <- 0
  sbar <- apply(pai%*%lo$p,2,mean) #score means
  rbar <- apply(pai-pai%*%lo$z,2,mean) #residual means

  y00 <- coredata(pai[1,,drop=FALSE]) #was *0
  s00 <- y00%*%lo$p
  r00 <- y00*0

  for(i in 1:nrow(pa)) {
    s01 <- sbar + (s00-sbar)*phis #last period (0) updated score using no this-period data
    r01 <- rbar + (r00-rbar)*phir #last period (0) updated residual using no this-period data
    y01 <- s01%*%t(lo$l)+r01 #updated one step ahead (1) total forecast

    y11 <- coredata(pa[i,])
    j11 <- is.na(y11)
    y11[j11] <- y01[j11] #substitute NA with updated total forecast

    s00 <- y11%*%lo$p #this period score
    r00 <- y11-s00%*%t(lo$l) #this period residual
    fit[i,] <- as.numeric(y01) #forecast
    fit00[i,] <- blend[i,] <- as.numeric(y11) #observed where observed; otherwise forecast
    fit00[i,j11] <- (s00%*%t(lo$l) + r00)[j11] #observed scores where not observed data
  }
  list(
    act=pa,
    fit=zoo(fit,index(pa)),
    blend=zoo(blend,index(pa)),
    res=pa-fit,
    fit00=zoo(fit00,index(pa)),
    MSE=sum((pa-fit)^2,na.rm=TRUE)
  )
}



myscreen <- function(da='2014-12-31',j=1:100,thresh=c(valid=.9,size=2000),prev="EQ0000000000047098",nreq=20) {
  stopifnot((thresh[1]<=1))
  da1 <- max(ca[as.character(date)<=da,date])
  x1 <- zoo(rollapplyr(data=!is.na(prem.g[,j]),width=52,mean,fill=NA),index(prem.g))[as.Date(da1)]
  j1 <- coredata(x1)>=thresh[1]
  stopifnot(nreq<=sum(j1,na.rm=TRUE))
  b1 <- colnames(j1)[j1]
  x2 <- mcap.g[as.Date(da1),b1]
  x2[,is.na(x2)] <- 0
  j2 <- thresh[2]<x2
  if(nreq>=sum(j2,na.rm=TRUE))   browser()
  b2 <- colnames(x2[,j2])
  setkey(data.table(bui=b2,inprev=-as.numeric(b2%in%prev),mcap=-as.numeric(coredata(x2)[,b2])),inprev,mcap)[1:nreq,]
}

myrepscr <- function(da=seq.Date(from=as.Date('2006-12-31'),by=365.25,len=8),...) {
  prev <- NULL
  ll <- NULL
  for(i in seq_along(da)) {
    ll[[i]] <- myscreen(prev=prev,da=da[i],...)
    prev <- ll[[i]][,bui]
  }
  ll
}


#' security universe scoring high on multiple criteria (in last period, no na, high mcap)
#'
#' @param years universe is reviewed each yearend
#' @param nsel number to select
#' @param x panel for testing
#' @export
maxdensu <- function(years=2009:2014,nsel=1000,x=best.g) {
  patt <- paste0("^",years,".?")
  ll <- NULL
  for( i in seq_along(years)) {
    test <- x[grep(patt=patt[i],rownames(x)),]
    x1 <- sort(apply(is.na(test),2,sum))[1:nsel]
    ll[[i]] <- data.table(date=max(index(test)),bui=names(x1),nas=as.numeric(x1))
  }
  rbindlist(ll)
}

#' access cod as a named list of ce
#'
#' @param cod of class co
#' @export
getcod <- function(cod=gett('cod')) {
  da <- cod[,unique(date)]
  ll <- structure(vector('list',length(da)),names=da)
  for(i in seq_along(ll)) {ll[[i]]<-aaco::dtce(cod[date==da[i]])}
  ll
}


#returns a list co, tei, tef
# getz <- function(ico=100,isu=99,da=su[,max(date)],nmin=3,wmin=3) {
#   co <- getrd(ico)
#   list(co=list(
#           z=getzco(co=co,part='z'),
#           p=getzco(co=co,part='psi'),
#           l=getzco(co=co,part='lambda') #does not exist
#           ),
#        #tei=getztei(su,da,loocv,type,part)
#        tei=list(
#           z=getztei(isu=isu,da=da,loocv,type,part='z'),
#           p=getztei(isu=isu,da=da,loocv,type,part='psi'),
#           l=getztei(isu=isu,da=da,loocv,type,part='lambda') #not exists
#           ),
#        tef=
#      )
# }


#' @export
mscecomp <- function(x,ret) {
  sco <- scoce(x, ret)
  ldg <- ldgce(x)
  stopifnot(all.equal(colnames(ret),rownames(ldg)))
  list(M=mz(sco[,1,drop=FALSE] %*% t(ldg[,1,drop=FALSE])),S=mz(sco[,-1,drop=FALSE] %*% t(ldg[,-1,drop=FALSE])),T=mz(sco %*% t(ldg)))
}


