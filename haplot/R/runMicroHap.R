
# Testing set:
haplo.sum<- readRDS("data/satro_sample/example_1.rds")  %>% mutate(id = as.character(id)) %>% tbl_df()
colnames(haplo.sum) <- c("group","id", "locus", "haplo", "depth", "logP.call", "logP.miscall", "pos", "allele.balance","rank")
collect.data <- RunGibbs(haplo.sum, "tag_id_1377", 200)

haplo.tbl <- haplo.sum
locus <- "tag_id_1377"
prior.model<-"uniform"
n.sam <- 20

## This function runs Gibbs sampling of observed haplotypes
#' @param haplo.tbl data frame. A haplotype data frame generated from running `runHaplo`.  The table contains "group" label, individual "id" label,
#'  "locus" label, observed "haplo"type, haplotype read "depth", "logP.call", "logP.miscall", "pos", "allele.balance","rank", etc.
#' @param locus string. the locus name. Required
#' @param n.sample integer. number of iterations. required
#' @param n.burn integer. number of burn-in cycle. n.burn < n.sample. Default sets as 0. Optional.
#' @param random.seed integer. Set the random seed number. Default sets as 43454. Optional
#' @param prior.model String. Choose two different prior models: "uniform" - prior set all prior haplotype weight to 1, or "empirical": the prior
#' alpha values are defined by the number of observed cases under refinement
RunGibbs <- function(haplo.tbl, locus, n.sam, n.burn=0, random.seed = 43454, prior.model="uniform")
{
  set.seed(random.seed)

  # reformat haplot.tbl
  haplo.tbl <- tidyHaplo(haplo.tbl, locus)

  #initialize parameters
  param <- setParam(haplo.tbl, n.sam, prior.model)

  #precompute all
  match.matrix <- PreComputeMatching(param, haplo.tbl)

  ##call updating for N.sam iterations
  for(i in 1:n.sam)
  {
    param$f <- UpdateF(param)
    param$pf <- UpdatePf(param, haplo.tbl)
    param$h <- UpdateH(param, match.matrix)

    param$save.freq[,i] <- param$f
    param$save.pfreq[,,i] <- param$pf
    param$save.hap[,,i] <- param$h
  }

  ### clean out the burn-in steps
  # param$save.freq <- param$save.freq[,-(1:n.burn)]
  # param$save.pfreq <- param$save.pfreq[,-(1:n.burn)]
  # param$save.hap <- param$save.hap[-(1:n.burn)]
  return(param)
}

tidyHaplo <- function(haplo.tbl, locus.select) {

  haplo.tbl <- haplo.tbl %>% dplyr::filter(locus == locus.select)
  n.sites <- length(unlist(strsplit(haplo.tbl$haplo[1],"")))
  haplo.tbl %>%
    #tidyr::separate(haplo, paste0("haplo.",1:n.sites), sep="(?!^)", extra="drop", remove=F) %>%
    #tidyr::separate(logP.call, paste0("logC.",1:n.sites), sep=",") %>%
    #tidyr::separate(logP.miscall, paste0("logW.",1:n.sites), sep=",") %>%
    dplyr::mutate(uniq.id = as.numeric(factor(id, levels=unique(haplo.tbl$id))))
}

setParam <- function(haplo.tbl, n.sam, prior.model){

  param <- NULL

  param$group <- unique(haplo.tbl$group)
  param$n.group <- length(param$group)
  param$n.indiv <- length(unique(haplo.tbl$id))
  param$n.sites <- length(strsplit(haplo.tbl$haplo[1],"")[[1]])

  param$grp.assoc.indiv <- haplo.tbl %>%
    dplyr::group_by(uniq.id) %>%
    dplyr::summarise(grp.indx = which(group[1]==param$group)) %>%
    dplyr::arrange(uniq.id) %>% dplyr::select(grp.indx) %>% unlist

  # define the sample space of all possible true haplotypes
  all.haplotype <- haplo.tbl %>%
    dplyr::filter(rank < 3 , allele.balance >0.3, depth > 10 ) %>%
    dplyr::group_by(haplo) %>%
    dplyr::summarise(count=n())

  param$haplo <- all.haplotype$haplo
  param$haplo.ct <- all.haplotype$count
  param$n.haplo <- length(all.haplotype$haplo)
  param$haplo.pair <- rbind(t(combn(1:param$n.haplo, 2)), matrix(rep(1:param$n.haplo,2), ncol=2))
  param$n.haplo.pair <- dim(param$haplo.pair)[1]

  ###cache the draw
  param$save.freq <- array(0, dim=c(param$n.haplo, n.sam))
  param$save.pfreq <- array(0, dim=c(param$n.group, param$n.haplo, n.sam))
  param$save.hap <- array(0, dim=c(param$n.indiv, param$n.haplo, n.sam))

  ## initialize all current parameters
  param$alpha <- param$haplo.ct
  if(prior.model == "uniform") param$alpha <- rep(1,param$n.haplo)
  param$f <- gtools::rdirichlet(1,param$alpha)

  param$pf <- t(sapply(1:param$n.group, function(i) {
    haplo.ct <- haplo.tbl %>%
      dplyr::filter(group==param$group[i]) %>%
      dplyr::group_by(haplo) %>%
      dplyr::summarise(count=n()) %>%
      dplyr::right_join(., data.frame("haplo"=param$haplo, stringsAsFactors = F),by="haplo")

    haplo.ct[is.na(haplo.ct)] <- 0
    gtools::rdirichlet(1, haplo.ct$count+1)
  }))


  param$h <- array(0, dim=c(param$n.indiv, param$n.haplo))

  haplo.select <- haplo.tbl %>%
    dplyr::group_by(uniq.id) %>%
    dplyr::summarise(haplo.1.indx = sample(1:param$n.haplo,
                                    1,
                                    prob=tabulate(match(haplo, param$haplo),
                                                     nbins=param$n.haplo)+
                                             param$pf[which(group[1]==param$group),]),
              haplo.2.indx = sample(1:param$n.haplo,
                                    1,
                                    prob=tabulate(match(haplo, param$haplo),
                                                  nbins=param$n.haplo)+
                                      param$pf[which(group[1]==param$group),]))

  param$h[cbind(haplo.select$uniq.id,
                haplo.select$haplo.1.indx)] <- 1
  param$h[cbind(haplo.select$uniq.id,
                haplo.select$haplo.2.indx)] <- param$h[cbind(haplo.select$uniq.id,
                                                             haplo.select$haplo.2.indx)]+1

  param$indic.combn <- matrix(0, nrow=param$n.haplo , ncol=param$n.haplo.pair)

  param$indic.combn[cbind(param$haplo.pair[,1],
                    1:param$n.haplo.pair)] <- 1
  param$indic.combn[cbind(param$haplo.pair[,2],
                    1:param$n.haplo.pair)] <- param$indic.combn[cbind(param$haplo.pair[,2],
                                                                1:param$n.haplo.pair)] + 1


  return(param)
}

PreComputeMatching <- function(param, haplo.tbl){

  sites.matrix <- strsplit(haplo.tbl$haplo, "") %>% unlist %>% matrix(ncol=param$n.sites, byrow=T)
  logC.matrix <- strsplit(haplo.tbl$logP.call, ",") %>% unlist %>% as.numeric() %>% matrix(ncol=param$n.sites, byrow=T)
  logI.matrix <- strsplit(haplo.tbl$logP.miscall, ",") %>% unlist %>% as.numeric() %>% matrix(ncol=param$n.sites, byrow=T)

  n.reads <- dim(haplo.tbl)[1]
  logP.read.match.ref <- sapply(1:param$n.haplo, function(i) {
    ref.matrix <- strsplit(param$haplo[i],"")%>% unlist %>% rep(., n.reads) %>% matrix(ncol=param$n.sites, byrow=T)
    rowSums(logC.matrix * (sites.matrix==ref.matrix) +
      logI.matrix * (sites.matrix!=ref.matrix))
  })

  exact.match <- sapply(1:param$n.haplo, function(i){
    haplo.tbl$haplo==param$haplo[i]
  })

  # impose that any observed haplotype that completely matches with one of the reference haplotype MUST
  # derived from that reference haplotype, thus ignore the possibility that observed haplotype
  # could generated from genotype error

  exact.match.pair.haplo <- exact.match %*% param$indic.combn > 0

  logP.match.by.indiv <- (exact.match.pair.haplo *
                            rowSums(logP.read.match.ref *exact.match) *
                            2 ) + (
                              (!exact.match.pair.haplo) *
                                (logP.read.match.ref %*% param$indic.combn))


  # Summing all log P by individual
  index.read.to.indiv <- haplo.tbl %>% dplyr::mutate(indx = row_number()) %>% dplyr::select(uniq.id, indx)%>% as.matrix(ncol=2,byrow=T)
  indic.matrix.read.by.indiv <- matrix(0, nrow=param$n.indiv, ncol=n.reads)
  indic.matrix.read.by.indiv[index.read.to.indiv] <- 1

  indic.matrix.read.by.indiv %*% logP.match.by.indiv
}



UpdateF <- function(param){
  gtools::rdirichlet(1,param$alpha)
}

UpdatePf <- function(param, haplo.tbl){
  t(sapply(1:param$n.group, function(i) {
    gtools::rdirichlet(1, colSums(param$h[param$grp.assoc.indiv==i,])+param$f)
  }))
}

UpdateH <- function(param, match.matrix){

  logP.contrib.pf <- param$pf %*% param$indic.combn %>%
    .[param$grp.assoc.indiv,] %>%
    log

  logP.h <- match.matrix + logP.contrib.pf

  # need to normalize to prevent numeric overflow
  prob.distrib <- t(apply(logP.h, 1, function(i) {
    exp(i-max(i))
  }))

  pair.indx <- sapply(1:param$n.indiv, function(i){
    sample(param$n.haplo.pair,1, prob=prob.distrib[i,])
  })

  h <- array(0, dim=c(param$n.indiv, param$n.haplo))
  h[cbind(1:param$n.indiv,
          param$haplo.pair[pair.indx,1]) ] <- 1

  h[cbind(1:param$n.indiv,
          param$haplo.pair[pair.indx,2]) ] <-
    h[cbind(1:param$n.indiv,
            param$haplo.pair[pair.indx,2]) ] + 1

  h

}


# prior alpha parameter: set as all 1 (weak prior), or set # based on what's observed


