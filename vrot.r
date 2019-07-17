## get redshift offsets from a data cube or stacked rss file

regrid <- function(lambda.out, lib) {
    lambda.in <- lib$lambda
    nc <- ncol(lib)
    lib.out <- matrix(NA, length(lambda.out), nc)
    lib.out[,1] <- lambda.out
    for (i in 2:nc) {
        lib.out[,i] <- approx(lambda.in, lib[,i], xout=lambda.out)$y
    }
    colnames(lib.out) <- colnames(lib)
    data.frame(lib.out)
}


getdz <- function(gdat, lib, snrthresh=5, nlthresh=2000, dzlim=0.003, searchinterval=1e-4) {
  dims <- dim(gdat$flux)
  nr <- dims[1]
  nc <- dims[2]
  dz <- matrix(NA, nr, nc)
  dz.err <- matrix(NA, nr, nc)
  z <- gdat$meta$z
  lambda <- gdat$lambda
  fmla <- as.formula(paste("f ~", paste(names(lib)[-1], collapse="+")))
  fn <- function(x) {
    lambda.out <- lambda/(1+z+x)
    lib.out <- regrid(lambda.out, lib)
    lfit <- lm(fmla, data=lib.out, weights=iv)
    deviance(lfit)
  }
  fn_v <- Vectorize(fn)
  z_search <- seq(-dzlim, dzlim, by=searchinterval)
  
  for (i in 1:nr) {
    for (j in 1:nc) {
      if (is.na(gdat$snr[i,j]) || gdat$snr[i,j]<=snrthresh) next
        f <- gdat$flux[i, j, ]
        if (length(which(!is.na(f))) < nlthresh) next
          iv <- gdat$ivar[i, j, ]
          dev_grid <- fn_v(z_search)
          z0 <- z_search[which.min(dev_grid)]
          bestz <- Rsolnp::solnp(pars=z0, fun=fn, 
                                 LB=z0-searchinterval, UB=z0+searchinterval, 
                                 control=list(trace=0))
          dz[i,j] <- bestz$pars
          dz.err[i,j] <- sqrt(2./bestz$hessian)
    }
  }
  list(dz=dz, dz.err=dz.err)
}

## peculiar velocity

vrel <- function(z, zbar) {
  299792.458*(z-zbar)/(1+zbar)
}

## simple 2D spatial interpolation

fillpoly <- function(ra, dec, zvals, dxy=0.5) {
    ny <- max(round((max(dec)-min(dec))*3600/dxy), 100)
    nx <- ny
    allok <- complete.cases(ra, dec, zvals)
    zsurf <- akima::interp(x=ra[allok], y=dec[allok], z=zvals[allok], 
                           linear=FALSE, nx=nx, ny=ny)
    list(ra=zsurf$x, dec=zsurf$y, z=zsurf$z)
}

## an image function using ggplot2
    
ggimage <- function(zmat, x=NULL, y=NULL, col=viridis::viridis(256),
                    xlab=NULL, ylab=NULL, legend=NULL, title=NULL,
                    xrev=FALSE, yrev=FALSE, asp=1,
                    addContour=FALSE, binwidth=NULL
                   ) {
    require(ggplot2)
    if (is.null(x)) x <- 1:nrow(zmat)
    if (is.null(y)) y <- 1:ncol(zmat)
    if (!is.null(dim(x))) {
      x <- colMeans(x)
      y <- rowMeans(y)
    }
    xy <- expand.grid(x,y)
    df <- data.frame(cbind(xy, as.vector(zmat)))
    names(df) <- c("x","y","z")
    g1 <- ggplot(df, aes(x=x, y=y, z=z)) + 
        geom_raster(aes(fill=z)) + 
        scale_fill_gradientn(colors=col, na.value="#FFFFFF00") +
        coord_fixed(ratio = asp)
    if (addContour) {
        if (!is.null(binwidth)) {
            g1 <- g1 + geom_contour(binwidth=binwidth, na.rm=TRUE)
        } else {
            g1 <- g1 + geom_contour(na.rm=TRUE)
        }
    }
    if (xrev) g1 <- g1 + scale_x_reverse()
    if (yrev) g1 <- g1 + scale_y_reverse()
    if (!is.null(xlab)) {
        g1 <- g1 + xlab(xlab)
    }
    if (!is.null(ylab)) {
        g1 <- g1 + ylab(ylab)
    }
    if (!is.null(legend)) {
        g1 <- g1 + labs(fill=legend)
    }
    if (!is.null(title)) {
        g1 <- g1 + ggtitle(title)
    }
    g1
}

rainbow <- c("blue", "cyan", "green", "yellow", "red")

## disk galaxy rotation curve model

vrot <- function(gdat.stack, dz.stack, phi.guess, sd.phi.guess=10,
                      ci.guess=0.7, sd.ci.guess=0.1,
                      sd.kc0=0.5,
                      order=3, N_r=100,
                      r_eff=1, remax=1.5, v_norm=100, PLOTS=TRUE,
                      smodel=NULL,
                      stanfile="vrot.stan", stanfiledir=".", 
                      seed=220755, iter=2000, warmup=1000,
                      chains=4, cores=4,
                      control=list(adapt_delta=0.975, max_treedepth=15)) {
  require(rstan)
  require(ggplot2)
  v <- vrel(dz.stack$dz+gdat.stack$meta$z, gdat.stack$meta$z)
  v_sample <- v[!is.na(v)]
  dv <- (dz.stack$dz.err/dz.stack$dz)*v
  df <- data.frame(x=gdat.stack$xpos, y=gdat.stack$ypos, v=v, dv=dv)
  df <- df[complete.cases(df),]
  rmax <- max(sqrt(gdat.stack$xpos^2+gdat.stack$ypos^2)[!is.na(dz.stack$dz)])
  if (rmax %% 2 != 0) rmax = rmax+1
  r_norm <- r_eff*remax
  r_post <- seq(0, 0.99, length=N_r)*rmax/r_norm
  
  vlos_data <- list(order=order, N=nrow(df), x=df$x/r_norm, y=df$y/r_norm, 
                    v=df$v/v_norm, dv=df$dv/v_norm, 
                    phi0=cos(phi.guess*pi/180), sd_phi0=sin(phi.guess*pi/180)*sd.phi.guess*pi/180, 
                    si0 = sqrt(1-ci.guess^2), sd_si0=sd.ci.guess,
                    sd_kc0 = sd.kc0,
                    r_norm=r_norm, v_norm=v_norm,
                    N_r=N_r, r_post=r_post)
  
  
  inits <- function(chain_id) {
    p <- phi.guess*pi/180
    si <- sqrt(1-ci.guess^2)
    list(phi=p, si=si, x_c=0., y_c=0., v_sys=0.)
  }
  if (!is.null(smodel)) {
    stanfit <- sampling(smodel, data=vlos_data, init=inits, seed=seed,
                        iter=iter, warmup=warmup, control=control, chains=chains, cores=cores)
  } else {
    stanfit <- stan(file=file.path(stanfiledir, stanfile), data=vlos_data, init=inits, 
                  seed=seed, iter=iter, warmup=warmup, control=control, chains=chains, cores=cores)
  }
  post <- extract(stanfit)
  vf_obs <- fillpoly(gdat.stack$ra.f, gdat.stack$dec.f, v)
  v_model <- rep(NA, length(v))
  v_model[!is.na(v)] <- colMeans(post$v_model)*v_norm
  vf_model <- fillpoly(gdat.stack$ra.f, gdat.stack$dec.f, v_model)
  r_q <- apply(post$r*remax, 2, quantile, probs=c(0.025, 0.5, 0.975))
  rownames(r_q) <- paste("r", rownames(r_q), sep="_")
  v_r_q <- apply(post$v_rot*v_norm, 2, quantile, probs=c(0.025, 0.5, 0.975))
  rownames(v_r_q) <- paste("v_rot", rownames(v_r_q), sep="_")
  v_e_q <- apply(post$v_exp*v_norm, 2, quantile, probs=c(0.025, 0.5, 0.975))
  rownames(v_e_q) <- paste("v_exp", rownames(v_e_q), sep="_")
  df_v <- data.frame(cbind(t(r_q), t(v_r_q), t(v_e_q)))
  graphs <- list()
  i <- 1
  df <- data.frame(x=gdat.stack$xpos, y=gdat.stack$ypos, v=v)
  graphs[[i]] <- ggplot(df) + geom_point(aes(x=x, y=y, color=v), size=7) +
                              scale_colour_gradientn(colors=rainbow, na.value="gray95") +
                              coord_fixed(ratio=1) +
                              ggtitle(paste("mangaid", gdat.stack$meta$mangaid, ": observed fiber velocities"))
  i <- i+1                              
  graphs[[i]] <- ggimage(vf_obs$z, vf_obs$ra, vf_obs$dec, col=rainbow, asp=1/cos(pi*gdat.stack$meta$dec/180), addContour=T) + 
                          xlim(rev(range(gdat.stack$ra.f))) + ylim(range(gdat.stack$dec.f)) +
                          xlab(expression(alpha)) + ylab(expression(delta)) +
                          ggtitle(paste("mangaid", gdat.stack$meta$mangaid, ": observed velocity field"))
  i <- i+1
  graphs[[i]] <- ggimage(vf_model$z, vf_model$ra, vf_model$dec, col=rainbow, asp=1/cos(pi*gdat.stack$meta$dec/180), addContour=T) + 
                          xlim(rev(range(gdat.stack$ra.f))) + ylim(range(gdat.stack$dec.f)) +
                          xlab(expression(alpha)) + ylab(expression(delta)) +
                          ggtitle(paste("mangaid", gdat.stack$meta$mangaid, ": model velocity field"))
  if (exists("v_res", post)) {
    v_res <- rep(NA, length(v))
    v_res[!is.na(v)] <- colMeans(post$v_res)*v_norm
    vf_res <- fillpoly(gdat.stack$ra.f, gdat.stack$dec.f, v_res)
  } else {
    vf_res <- list(z=vf_obs$z-vf_model$z, ra=vf_obs$ra, dec=vf_obs$dec)
  }
  i <- i+1
  graphs[[i]] <- ggimage(vf_res$z, vf_res$ra, vf_res$dec, asp=1/cos(pi*gdat.stack$meta$dec/180), addContour=T) + 
                        xlim(rev(range(gdat.stack$ra.f))) + ylim(range(gdat.stack$dec.f)) +
                        xlab(expression(alpha)) + ylab(expression(delta)) +
                        ggtitle(paste("mangaid", gdat.stack$meta$mangaid, ": residual velocity field"))
  df1 <- data.frame(r=as.vector(abs(post$r))*remax, v_rot=as.vector(post$v_rot)*v_norm, v_exp=as.vector(post$v_exp)*v_norm)
  i <- i+1
  graphs[[i]] <- ggplot(df_v) + geom_point(aes(x=r_50., y=v_rot_50.)) +
          geom_errorbar(aes(x=r_50., ymin=v_rot_2.5., ymax=v_rot_97.5.)) +
          geom_errorbarh(aes(x=r_50., y=v_rot_50., xmin=r_2.5., xmax=r_97.5.)) +
          stat_density_2d(aes(x=r, y=v_rot, fill=..level..), data=df1, geom="polygon", alpha=0.5, show.legend=FALSE) + 
          ggtitle(paste("mangaid", gdat.stack$meta$mangaid, ": joint distribution of r, v_rot")) +
          xlab(expression(r/r[eff]))
  i <- i+1
  graphs[[i]] <- ggplot(df_v) + geom_point(aes(x=r_50., y=v_exp_50.)) +
          geom_errorbar(aes(x=r_50., ymin=v_exp_2.5., ymax=v_exp_97.5.)) +
          geom_errorbarh(aes(x=r_50., y=v_exp_50., xmin=r_2.5., xmax=r_97.5.)) +
          stat_density_2d(aes(x=r, y=v_exp, fill=..level..), data=df1, geom="polygon", alpha=0.5, show.legend=FALSE) + 
          ggtitle(paste("mangaid", gdat.stack$meta$mangaid, ": joint distribution of r, v_exp")) +
          xlab(expression(r/r[eff]))
  if(exists("vrot_post", post)) {
  	  vq <- t(apply(v_norm * post$vrot_post, 2, quantile, probs=c(0.025, 0.5, 0.975)))
  	  df2 <- data.frame(r=r_post*remax, ymin=vq[,"2.5%"], ymed=vq[,"50%"], ymax=vq[,"97.5%"])
          i <- i+1
  	  graphs[[i]] <- ggplot(df2) + geom_line(aes(x=r, y=ymed)) + geom_ribbon(aes(x=r, ymin=ymin, ymax=ymax), alpha=0.5) +
          xlab(expression(r/r[eff]))
  } else {
  	  df2 <- NULL
  }
  if (PLOTS) {
  	  for (i in seq_along(graphs)) {
  	  	  x11()
  	  	  plot(graphs[[i]])
  	  }
  }
  list(stanfit=stanfit, r_eff=r_eff, remax=remax, v_norm=v_norm, vlos_data=vlos_data, 
       v_sample=v_sample, vf_obs=vf_obs, vf_model=vf_model, vf_res=vf_res, df_v=df_v, df1=df1, df2=df2,
       graphs=graphs
  )
}

## disk galaxy rotation curve model -- gaussian process version
  
vrot_gp <- function(gdat.stack, dz.stack, phi.guess, sd.phi.guess=10,
                      ci.guess=0.7, sd.ci.guess=0.1,
                      sd.kc0=0.5,
                      order=3, N_xy=25^2,
                      r_eff=1, remax=1.5, v_norm=100, PLOTS=TRUE,
                      smodel=NULL,
                      stanfile="vrot_gp.stan", stanfiledir="~/vrot_stanmodels", 
                      seed=220755, iter=1000, warmup=500,
                      chains=4, cores=4,
                      control=list(adapt_delta=0.9, max_treedepth=10)) {
  require(rstan)
  require(spmutils)
  require(ggplot2)
  require(akima)
  v <- vrel(dz.stack$dz+gdat.stack$meta$z, gdat.stack$meta$z)
  v_sample <- v[!is.na(v)]
  dv <- (dz.stack$dz.err/dz.stack$dz)*v
  df <- data.frame(x=gdat.stack$xpos, y=gdat.stack$ypos, v=v, dv=dv)
  allok <- complete.cases(df)
  df <- df[allok,]
  
  n_r <- round(sqrt(N_xy))
  N_xy <- n_r^2
  rho <- seq(1/n_r, 1, length=n_r)
  theta <- seq(0, (n_r-1)*2*pi/n_r, length=n_r)
  rt <- expand.grid(theta, rho)
  x_pred <- rt[,2]*cos(rt[,1])
  y_pred <- rt[,2]*sin(rt[,1])
  r_norm <- r_eff*remax
  
  vlos_data <- list(order=order, N=nrow(df), x=df$x/r_norm, y=df$y/r_norm, 
                    v=df$v/v_norm, dv=df$dv/v_norm, 
                    phi0=phi.guess*pi/180, sd_phi0=sd.phi.guess*pi/180, 
                    ci0 = ci.guess, sd_ci0=sd.ci.guess,
                    sd_kc0 = sd.kc0,
                    r_norm=r_norm, v_norm=v_norm,
                    N_xy=N_xy, N_r=n_r, x_pred=x_pred, y_pred=y_pred
                    )
  
  
  inits <- function(chain_id) {
    p <- phi.guess*pi/180
    ci <- ci.guess
    list(phi=p, cos_i=ci, x_c=0., y_c=0., v_sys=0.)
  }
  if (!is.null(smodel)) {
    stanfit <- sampling(smodel, data=vlos_data, init=inits, seed=seed,
                        iter=iter, warmup=warmup, control=control, chains=chains, cores=cores)
  } else {
    stanfit <- stan(file=file.path(stanfiledir, stanfile), data=vlos_data, init=inits, 
                  seed=seed, iter=iter, warmup=warmup, control=control, chains=chains, cores=cores)
  }
  post <- extract(stanfit)
  graphs <- list()
  vf_model <- rep(NA, length(gdat.stack$xpos))
  vf_res <- rep(NA, length(gdat.stack$ypos))
  vf_model[allok] <- colMeans(post$vf_model)*v_norm
  vf_res[allok] <- colMeans(post$vf_res)*v_norm
  df <- data.frame(x=gdat.stack$xpos, y=gdat.stack$ypos, v=v, vf_model=vf_model, vf_res=vf_res)
  g <- 1
  graphs[[g]] <- ggplot(df) + geom_point(aes(x=x, y=y, color=v), size=7) +
                              scale_colour_gradientn(colors=rainbow, na.value="gray95") +
                              coord_fixed(ratio=1) +
                              ggtitle(paste("mangaid", gdat.stack$meta$mangaid, ": observed fiber velocities"))
  g <- g+1
  df <- df[allok,]
  df1 <- akima::interp(x=df$x, y=df$y, z=df$vf_model, nx=200, ny=200)
  graphs[[g]] <- ggimage(df1$z, df1$x, df1$y, col=rainbow, xlab="x", ylab="y", legend="vf_model", addContour=TRUE)
  g <- g+1
  df1 <- akima::interp(x=df$x, y=df$y, z=df$vf_res, nx=200, ny=200)
  graphs[[g]] <- ggimage(df1$z, df1$x, df1$y, xlab="x", ylab="y", legend="vf_res", addContour=TRUE)
  
  v_pred <- colMeans(post$v_pred)
  sd_v_pred <- apply(post$v_pred, 2, sd)
  v_model <- colMeans(post$v_model)
  v_res <- colMeans(post$v_res)
  x <- x_pred*r_norm
  y <- y_pred*r_norm
  rho <- c(0, rho)*remax
  df <- akima::interp(x=x, y=y, z=v_pred, nx=200, ny=200)
  g <- g+1
  graphs[[g]] <- ggimage(df$z, df$x, df$y, col=rainbow, xlab="x", ylab="y", legend="v_pred", addContour=TRUE)
  df <- akima::interp(x=x, y=y, z=v_model, nx=200, ny=200)
  g <- g+1
  graphs[[g]] <- ggimage(df$z, df$x, df$y, col=rainbow, xlab="x", ylab="y", legend="v_model", addContour=TRUE)
  df <- akima::interp(x=x, y=y, z=v_res, nx=200, ny=200)
  g <- g+1
  graphs[[g]] <- ggimage(df$z, df$x, df$y, xlab="x", ylab="y", legend="v_res", addContour=TRUE)
  df <- akima::interp(x=x, y=y, z=sd_v_pred, nx=200, ny=200)
  g <- g+1
  graphs[[g]] <- ggimage(df$z, df$x, df$y, xlab="x", ylab="y", legend="sd_v_pred", addContour=TRUE)
  vrot_pred <- post$vrot_pred
  vexp_pred <- post$vexp_pred
  v_r_q <- data.frame(cbind(rho, t(apply(vrot_pred, 2, quantile, probs=c(.025, 0.5, 0.975)))))
  names(v_r_q) <- c("rho", "vr_025", "vr_500", "vr_975")
  g <- g+1
  graphs[[g]] <- ggplot(v_r_q) + geom_line(aes(x=rho, y= vr_500)) + geom_ribbon(aes(x=rho, ymin=vr_025, ymax=vr_975), alpha=0.5) +
                  xlab(expression(r/r[eff]))+ylab(expression(V[rot]))
  v_e_q <- data.frame(cbind(rho, t(apply(vexp_pred, 2, quantile, probs=c(.025, 0.5, 0.975)))))
  names(v_e_q) <- c("rho", "ve_025", "ve_500", "ve_975")
  g <- g+1
  graphs[[g]] <- ggplot(v_e_q) + geom_line(aes(x=rho, y= ve_500)) + geom_ribbon(aes(x=rho, ymin=ve_025, ymax=ve_975), alpha=0.5) +
                  xlab(expression(r/r[eff]))+ylab(expression(V[exp]))

  r_q <- apply(post$r*remax, 2, quantile, probs=c(0.025, 0.5, 0.975))
  rownames(r_q) <- paste("r", rownames(r_q), sep="_")
  v_r_q <- apply(post$v_rot*v_norm, 2, quantile, probs=c(0.025, 0.5, 0.975))
  rownames(v_r_q) <- paste("v_rot", rownames(v_r_q), sep="_")
  df_v <- data.frame(cbind(t(r_q), t(v_r_q)))
  df <- data.frame(r=as.vector(abs(post$r))*remax, v_rot=as.vector(post$v_rot)*v_norm)

  g <- g+1
  graphs[[g]] <- ggplot(df_v) + geom_point(aes(x=r_50., y=v_rot_50.)) +
          geom_errorbar(aes(x=r_50., ymin=v_rot_2.5., ymax=v_rot_97.5.)) +
          geom_errorbarh(aes(x=r_50., y=v_rot_50., xmin=r_2.5., xmax=r_97.5.)) +
          ggtitle(paste("mangaid", gdat.stack$meta$mangaid, ": joint distribution of r, v_rot")) +
          xlab(expression(r/r[eff])) + ylab("v (km/sec)")
  if (PLOTS) {
    for (i in seq_along(graphs)) {
      x11()
      plot(graphs[[i]])
    }
  }

  list(stanfit=stanfit, r_eff=r_eff, remax=remax, v_norm=v_norm, vlos_data=vlos_data, 
       v_sample=v_sample, x_pred=x_pred, y_pred=y_pred,
       graphs=graphs
  )
}
  
