readcube <- function(fname, drpcat=NULL) {
    require(FITSio)
    ff <- file(fname, "rb")
    hd0 <- parseHdr(readFITSheader(ff))
    mangaid <- hd0[grep("MANGAID", hd0)+1]
    plateifu <- hd0[grep("PLATEIFU", hd0)+1]
    ra <- as.numeric(hd0[grep("IFURA", hd0)+1])
    dec <- as.numeric(hd0[grep("IFUDEC", hd0)+1])
    l <- as.numeric(hd0[grep("IFUGLON", hd0)+1])
    b <- as.numeric(hd0[grep("IFUGLAT", hd0)+1])    
    ebv <- as.numeric(hd0[grep("EBVGAL", hd0)+1])
    hd0 <- parseHdr(readFITSheader(ff))
    cra <- as.numeric(hd0[grep("CRVAL1", hd0)+1])    
    cdec <- as.numeric(hd0[grep("CRVAL2", hd0)+1])    
    cpix1 <- as.numeric(hd0[grep("CRPIX1", hd0)+1])    
    cpix2 <- as.numeric(hd0[grep("CRPIX2", hd0)+1])    
    cd1 <- as.numeric(hd0[grep("CD1_1", hd0)+1])    
    cd2 <- as.numeric(hd0[grep("CD2_2", hd0)+1])    
    close(ff)
    
    flux <- readFITS(fname, hdu=1)$imDat
    ivar <- readFITS(fname, hdu=2)$imDat
    mask <- readFITS(fname, hdu=3)$imDat
    lambda <- readFITS(fname, hdu=4)$imDat
    ivar[ivar <= 0] <- NA
    ivar[mask >= 1024] <- NA
    flux[is.na(ivar)] <- NA
    extc <- elaw(lambda, ebv)
    nx <- dim(flux)[1]
    ny <- dim(flux)[2]
    for (ix in 1:nx) {
        for (iy in 1:ny) {
            flux[ix,iy,] <- flux[ix,iy,]*extc
            ivar[ix,iy,] <- ivar[ix,iy,]/extc^2
        }
    }
    snr <- apply(sqrt(pmax(flux,0)^2*ivar), c(1,2), median, na.rm=TRUE)
    xpos <- (1:nx-cpix1)*cd1
    ypos <- (1:ny-cpix2)*cd2
    phi <- function(x,y) atan2(-y, x)
    rho <- function(x,y) sqrt(x^2 + y^2)
    long <- outer(xpos, ypos, phi)
    lat <- atan(180/pi/outer(xpos, ypos, rho))
    dec.f <- 180/pi*asin(sin(lat)*sin(cdec*pi/180) - cos(lat)*cos(long)*cos(cdec*pi/180))
    ra.f <- cra + 180/pi*atan2(cos(lat)*sin(long),
                               sin(lat)*cos(cdec*pi/180) + cos(lat)*cos(long)*sin(cdec*pi/180))
    gimg <- readFITS(fname, hdu=8)$imDat*elaw(4640.4,ebv)
    rimg <- readFITS(fname, hdu=9)$imDat*elaw(6122.3,ebv)
    iimg <- readFITS(fname, hdu=10)$imDat*elaw(7439.5,ebv)
    zimg <- readFITS(fname, hdu=11)$imDat*elaw(8897.1,ebv)
    if (!is.null(drpcat)) {
        ind <- which(drpcat$plateifu == plateifu)
        z <- drpcat$nsa_z[ind]
        zdist <- drpcat$nsa_zdist[ind]
    } else {
        z <- zdist <- NA
    }
    list(meta=list(mangaid=mangaid, plateifu=plateifu, ra=ra, dec=dec, l=l, b=b, 
                   ebv=ebv, z=z, zdist=zdist, cpix=c(cpix1, cpix2)),
        xpos=xpos, ypos=ypos,
        lambda=lambda, ra.f=ra.f, dec.f=dec.f, flux=flux, ivar=ivar, snr=snr,
         gimg=gimg, rimg=rimg, iimg=iimg, zimg=zimg)
}

readrss <- function(fname, drpcat=NULL, ndither=3) {
    require(FITSio)
    ff <- file(fname, "rb")
    hd0 <- parseHdr(readFITSheader(ff))
    mangaid <- hd0[grep("MANGAID", hd0)+1]
    plateifu <- hd0[grep("PLATEIFU", hd0)+1]
    nexp <- as.numeric(hd0[grep("NEXP", hd0)+1])
    ra <- as.numeric(hd0[grep("IFURA", hd0)+1])
    dec <- as.numeric(hd0[grep("IFUDEC", hd0)+1])
    l <- as.numeric(hd0[grep("IFUGLON", hd0)+1])
    b <- as.numeric(hd0[grep("IFUGLAT", hd0)+1])    
    ebv <- as.numeric(hd0[grep("EBVGAL", hd0)+1])
    close(ff)
    
    flux <- readFITS(fname, hdu=1)$imDat
    ivar <- readFITS(fname, hdu=2)$imDat
    mask <- readFITS(fname, hdu=3)$imDat
    lambda <- as.vector(readFITS(fname, hdu=5)$imDat)
    xpos <- readFITS(fname, hdu=9)$imDat
    ypos <- readFITS(fname, hdu=10)$imDat
    xpos <- apply(xpos, 2, mean)
    ypos <- apply(ypos, 2, mean)
    ivar[ivar <= 0] <- NA
    ivar[mask >= 1024] <- NA
    flux[is.na(ivar)] <- NA
    extc <- elaw(lambda, ebv)
    flux <- flux*extc
    ivar <- ivar/extc^2
    nl <- nrow(flux)
    ns <- nexp/ndither
    npos <- ncol(flux)/ns
    xy <- data.frame(x=xpos, y=ypos)
    txy <- SearchTrees::createTree(xy)
    neighbors <- as.vector(SearchTrees::knnLookup(txy, newdat=xy[1:npos,], k=ns))
    flux <- flux[,neighbors]
    ivar <- ivar[,neighbors]
    snr <- apply(sqrt(pmax(flux,0)^2*ivar), 2, median, na.rm=TRUE)
    flux <- array(t(flux), dim=c(npos, ns, nl))
    ivar <- array(t(ivar), dim=c(npos, ns, nl))
    snr <- matrix(snr, npos, ns)
    xpos <- matrix(xpos[neighbors], npos, ns)
    ypos <- matrix(ypos[neighbors], npos, ns)
    if (!is.null(drpcat)) {
        ind <- which(drpcat$plateifu == plateifu)
        z <- drpcat$nsa_z[ind]
        zdist <- drpcat$nsa_zdist[ind]
    } else {
        z <- zdist <- NA
    }
    list(meta=list(mangaid=mangaid, plateifu=plateifu, nexp=nexp, ra=ra, dec=dec, l=l, b=b, 
                   ebv=ebv, z=z, zdist=zdist),
         lambda=lambda, flux=flux, ivar=ivar, snr=snr,
         xpos=xpos, ypos=ypos
        )
}

stackrss <- function(gdat, dz=NULL) {
    meta <- gdat$meta
    ivar <- apply(gdat$ivar, c(1,3), sum, na.rm=TRUE)
    ivar[ivar <= 0] <- NA
    flux <- apply(gdat$flux*gdat$ivar, c(1,3), sum, na.rm=TRUE)/ivar
    snr <- apply(sqrt(pmax(flux,0)^2*ivar), 1, median, na.rm=TRUE)
    nr <- nrow(flux)
    nc <- ncol(flux)
    dim(flux) <- c(nr, 1, nc)
    dim(ivar) <- c(nr, 1, nc)
    snr <- matrix(snr, nr, 1)
    if (!is.null(dz)) {
        dz <- rowMeans(dz, na.rm=TRUE)
    }
    xpos <- rowMeans(gdat$xpos)
    ypos <- rowMeans(gdat$ypos)
    dec.f <- meta$dec + ypos/3600
    ra.f <- meta$ra - xpos/3600/cos(dec.f*pi/180)
    list(meta=meta, lambda=gdat$lambda, flux=flux, ivar=ivar, snr=snr, 
         xpos=xpos, ypos=ypos, ra.f=ra.f, dec.f=dec.f, dz=dz)
}

## Galactic extinction correction from Fitzpatrick (1998): http://arxiv.org/abs/astro-ph/9809387v1
## This is spline fit portion valid from near-UV to near-IR and R=3.1

elaw <- function(lambda, ebv) {
  il <- c(0,0.377,0.820,1.667,1.828,2.141,2.433,3.704,3.846)
  al <- c(0,0.265,0.829,2.688,3.055,3.806,4.315,6.265,6.591)
  fai <- splinefun(il,al)
  10^(0.4*fai(10000/lambda)*ebv)
}

