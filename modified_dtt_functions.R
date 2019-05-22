
modified_disparity <- 
  function (phy = NULL, data, index = c("avg.sq", "avg.manhattan", "discrete", "binary")) 
  {
    if (!is.null(phy)) {
      if (!"phylo" %in% class(phy)) 
        stop("Supply 'phy' as a 'phylo' object")
      td <- treedata(phy, data)
      phy = td$phy
      data = td$data
      desc = geiger:::.cache.descendants(phy)$tips
      nb.tip <- length(td$phy$tip.label)
      nb.node <- td$phy$Nnode
      result <- numeric()
      for (i in 1:nb.node) {
        l <- desc[[nb.tip + i]]
        d <- td$data[phy$tip.label[l], ]
        result[i] <- modified.disparity(d, index)
      }
      names(result) = nb.tip + 1:nb.node
      return(result)
    }
    else {
      return(modified.disparity(data, index))
    }
  }

modified.disparity <- function (data, index = c("avg.sq", "avg.manhattan", "discrete", "binary")) 
{
  disp = match.arg(index, c("avg.sq", "avg.manhattan", "discrete", "binary"))
  if (disp == "avg.sq") {
    d <- dist(data, method = "euclidean")^2
    r <- mean(d)
  }
  if (disp == "avg.manhattan") {
    d <- dist(data, method = "manhattan")
    r <- mean(d)
  }
  if (disp == "binary") {
    # variance of a binomial variable
     data <- ifelse(grepl("^N", data), yes = 0, no = 1)
     bin_var <- function(n, p){n * p * (1 - p)}
     r <- bin_var(n = 1, 
                  p = sum(data)/length(data))
    
    #d <- daisy(data_df, metric = "gower")
    #r <- mean(d)
    # number of states
    #r <- length(unique(data[!is.na(data)]))
  }
  if (disp == "discrete"){
    # number of states
    r <- length(unique(data[!is.na(data)]))
  }
  return(r)
}




modified.dtt <- 
  function (phy, data, disp = c("avg.sq", "avg.manhattan", "discrete", "binary")) 
  {
    disp = match.arg(disp, c("avg.sq", "avg.manhattan", "discrete", "binary"))
    phy$node.label <- NULL
    td <- treedata(phy, data)
    phy2 <- td$phy
    phy <- new2old.phylo(td$phy)
    result <- numeric()
    node.depth <- branching.times(phy2)
    stem.depth <- numeric()
    stem.depth[1] <- node.depth[1]
    for (i in 2:phy2$Nnode) {
      anc <- which(as.numeric(phy$edge[, 2]) == -i)
      stem.depth[i] <- node.depth[names(node.depth) == phy2$edge[anc, 
                                                                 1]]
    }
    ltt <- sort(node.depth, decreasing = TRUE)
    node.depth <- node.depth/max(ltt)
    stem.depth <- stem.depth/max(ltt)
    ltt <- ltt/max(ltt)
    if (length(dim(td$data)) == 2) {
      d <- modified_disparity(phy2, td$data, index = disp)
      result[1] <- d[1]
      for (i in 2:length(ltt)) {
        x <- d[stem.depth >= ltt[i - 1] & node.depth < ltt[i - 
                                                             1]]
        if (length(x) == 0) 
          result[i] = 0
        else result[i] <- mean(x)
      }
      result[length(ltt) + 1] <- 0
      if (result[1] > 0) 
        result <- result/result[1]
    }
    if (length(dim(td$data)) > 2) {
      for (i in 1:dim(td$data)[3]) {
        pp <- as.matrix(td$data[, , i])
        d <- modified_disparity(phy2, pp, index = disp)
        y <- numeric()
        y[1] <- d[1]
        for (j in 2:length(ltt)) {
          x <- d[stem.depth >= ltt[j - 1] & node.depth < 
                   ltt[j - 1]]
          if (length(x) == 0) 
            y[j] = 0
          else y[j] <- mean(x)
        }
        y[length(ltt) + 1] <- 0
        if (y[1] > 0) 
          y <- y/y[1]
        result <- cbind(result, y)
      }
    }
    
    return(result)
  }

modified_dtt <- function (phy, data, index = c("avg.sq", "avg.manhattan", "discrete", "binary"), 
                          mdi.range = c(0, 1), nsim = 999, CI = 0.95, plot = FALSE, calculateMDIp = T) 
{
  disp = match.arg(index, c("avg.sq", "avg.manhattan", "discrete", "binary"))
  td <- treedata(phy, data)
  dtt.data <- modified.dtt(td$phy, td$data, disp = disp)
  ltt <- sort(branching.times(td$phy), decreasing = TRUE)
  ltt <- c(0, (max(ltt) - ltt)/max(ltt))
  dtt.sims = NULL
  MDI = NULL
  ylim = c(range(pretty(dtt.data)))
  if (is.numeric(nsim) & is.numeric(td$data)) {
    if (nsim > 0) {
      s <- ratematrix(td$phy, td$data)
      sims <- sim.char(td$phy, s, nsim)
      dtt.sims <- modified.dtt(td$phy, sims)
      mean.sims <- apply(dtt.sims, 1, mean)
      median.sims <- apply(dtt.sims, 1, median)
      MDI <- unname(geiger:::.area.between.curves(ltt, apply(dtt.sims, 
                                                             1, median), dtt.data, sort(mdi.range)))
      names(MDI) = disp
      colnames(dtt.sims) = NULL
      yy = range(dtt.sims)
      ylim = range(c(ylim, yy))
    }
  }
  if (is.numeric(nsim) & is.character(td$data)) {
    if (nsim > 0) {
      Q <- geiger:::.Qmatrix.from.gfit(fit_discrete_traits_all[[i]][[1]]) # use BM estimates of transition rates for null model simulations
      sims <- sim.char(td$phy, Q, nsim, model = "discrete")
      dtt.sims <- modified.dtt(td$phy, sims)
      mean.sims <- apply(dtt.sims, 1, mean)
      median.sims <- apply(dtt.sims, 1, median)
      MDI <- unname(geiger:::.area.between.curves(ltt, apply(dtt.sims, 
                                                             1, median), dtt.data, sort(mdi.range)))
      names(MDI) = disp
      colnames(dtt.sims) = NULL
      yy = quantile(dtt.sims, probs = c(0.001, 0.999))
      ylim = range(c(ylim, yy))
    }
  }
  if (plot) {
    plot(ltt, dtt.data, xlab = "relative time", ylab = "disparity", 
         ylim = ylim, bty = "n", type = "n")
    if (!is.null(dtt.sims)) {
      poly = geiger:::.dtt.polygon(dtt.sims, ltt, alpha = 1 - CI)
      polygon(poly[, "x"], poly[, "y"], col = geiger:::.transparency("lightgray", 
                                                                     0.5), border = NA)
      lines(ltt, median.sims, lty = 2)
    }
    lines(ltt, dtt.data, type = "l", lwd = 2)
  }
  res = list(dtt = dtt.data, times = ltt, sim = dtt.sims, MDI = MDI)
  drp = sapply(res, function(x) is.null(x))
  if (any(drp)) 
    res = res[-which(drp)]
  if (calculateMDIp) {
    pVal <- modified_getMDIp(res)
    res <- c(res, MDIpVal = pVal)
  }
  return(res)
}

modified_getMDIp<-function(dttRes) { # this is the version on github which fixes a big in estimating the P values
  foo<-function(x) {
    return(geiger:::.area.between.curves(x= dttRes$times, f1=x, f2=dttRes$dtt))
  }
  mdis<-apply(dttRes$sim,2,foo)
  
  #Two sided test
  p1<-length(which(mdis>=0))/length(mdis)
  p2<-length(which(mdis<=0))/length(mdis)
  
  pVal<-min(p1,p2)
  return(pVal)
}

####################################################
#
#	Generic function to compute rank envelope test
#
#	two tailed test: test="two.sided"
#	one sided tests: test="less" OR test="greater"
#
####################################################

rank_env_dtt<-function(x, Plot=F, test="two.sided"){
  
  spp_num<-length(x$times)		
  sims<-x$sim
  sims<-as.matrix(sims)
  
  s1<-sims[-c(1),]
  
  r<-x$times[-c(1)]
  
  r<-as.vector(r)
  
  obs<-as.vector(x$dtt)
  obs<-obs[-c(1)]
  
  c1<-list(r,obs, s1)
  names(c1)=c("r","obs","sim_m") 
  c2<-create_curve_set(c1)
  
  res<-rank_envelope(c2, alternative=test, savedevs = TRUE)
  
  if(Plot==T)
    plot(res, xlab="Relative time", ylab="Disparity", main="")	
  return(res)	
}

