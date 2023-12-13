#### Published code ####

## Requires data from NPS figshare
## NB Data from Bed B cannot be shared.
## NB Data shared only the first 100 replicates of each community due to data storage limitations

#### Dependencies ####
library(vegan)
library(ggplot2)
library(forams)
library(tidyverse)

#### Data ####

# Abundance data for composition analyses
sad <- read.csv("./Abun_Data/Surface_Abundance_Data.csv")

# Surface data for succession analyses
bish <- read.csv("./JK_Data/Bishops_Cove.csv")
br5 <- read.csv("./JK_Data/Brasier.csv")
bristy <- read.csv("./JK_Data/Bristy_Cove.csv")
fac <- read.csv("./JK_Data/H38.csv")
gm <- read.csv("./JK_Data/Goldmine.csv")
gp <- read.csv("./JK_Data/Green_Point.csv")
h14 <- read.csv("./JK_Data/H14.csv")
h5 <- read.csv("./JK_Data/H5.csv")
lmp <- read.csv("./JK_Data/LMP.csv")
mel <- read.csv("./JK_Data/Melrose.csv")
mpd <- read.csv("./JK_Data/D_Surface.csv")
mpe <- read.csv("./JK_Data/E_Surface.csv")
mpg <- read.csv("./JK_Data/G_Surface.csv")
piz <- read.csv("./JK_Data/Pizzeria.csv")
primo <- read.csv("./JK_Data/H26.csv")
sh <- read.csv("./JK_Data/Shingle_Head.csv")
sts <- read.csv("./JK_Data/St._Shotts.csv")

sur_df_list <- c(split(bish, f = bish$JK_ID), split(br5, f = br5$JK_ID), split(fac, f = fac$JK_ID),
                 split(gm , f = gm$JK_ID), split(gp, f = gp$JK_ID), split(h14, f = h14$JK_ID),
                 split(h5 , f = h5$JK_ID), split(lmp, f = lmp$JK_ID), split(mel, f = mel$JK_ID),
                 split(mpd, f = mpd$JK_ID), split(mpe, f = mpe$JK_ID), split(mpg, f = mpg$JK_ID),
                 split(piz, f = piz$JK_ID), split(primo, f = primo$JK_ID), split(sts, f = sts$JK_ID))

# Surface data for tiering analyses

bish_hd <- read.csv("./Height_Distributions/bish_hd.csv")
br5_hd <- read.csv("./Height_Distributions/br5_hd.csv")
bristy_hd <- read.csv("./Height_Distributions/bristy_hd.csv")
fac_hd <- read.csv("./Height_Distributions/fac_hd.csv")
gm_hd <- read.csv("./Height_Distributions/gm_hd.csv")
gp_hd <- read.csv("./Height_Distributions/gp_hd.csv")
h14_hd <- read.csv("./Height_Distributions/h14_hd.csv")
h5_hd <- read.csv("./Height_Distributions/h5_hd.csv")
lmp_hd <- read.csv("./Height_Distributions/lmp_hd.csv")
mel_hd <- read.csv("./Height_Distributions/mel_hd.csv")
mpd_hd <- read.csv("./Height_Distributions/mpd_hd.csv")
mpe_hd <- read.csv("./Height_Distributions/mpe_hd.csv")
mpg_hd <- read.csv("./Height_Distributions/mpg_hd.csv")
piz_hd <- read.csv("./Height_Distributions/piz_hd.csv")
primo_hd <- read.csv("./Height_Distributions/primo_hd.csv")
sh_hd <- read.csv("./Height_Distributions/sh_hd.csv")
sts_hd <- read.csv("./Height_Distributions/sts_hd.csv")

hd_list <- list(bish_hd, br5_hd, bristy_hd, fac_hd, gm_hd, gp_hd, h14_hd, h5_hd, lmp_hd,
             mel_hd, mpd_hd, mpe_hd, mpg_hd, piz_hd, primo_hd, sh_hd, sts_hd)


bishabuns <- read.csv("./Height_Distributions/bishabuns.csv")
br5abuns <- read.csv("./Height_Distributions/br5abuns.csv")
bristyabuns <- read.csv("./Height_Distributions/bristyabuns.csv")
facabuns <- read.csv("./Height_Distributions/facabuns.csv")
gmabuns <- read.csv("./Height_Distributions/gmabuns.csv")
gpabuns <- read.csv("./Height_Distributions/gpabuns.csv")
h14abuns <- read.csv("./Height_Distributions/h14abuns.csv")
h5abuns <- read.csv("./Height_Distributions/h5abuns.csv")
lmpabuns <- read.csv("./Height_Distributions/lmpabuns.csv")
melabuns <- read.csv("./Height_Distributions/melabuns.csv")
mpdabuns <- read.csv("./Height_Distributions/mpdabuns.csv")
mpeabuns <- read.csv("./Height_Distributions/mpeabuns.csv")
mpgabuns <- read.csv("./Height_Distributions/mpgabuns.csv")
pizabuns <- read.csv("./Height_Distributions/pizabuns.csv")
primoabuns <- read.csv("./Height_Distributions/primoabuns.csv")
shabuns <- read.csv("./Height_Distributions/shabuns.csv")
stsabuns <- read.csv("./Height_Distributions/stsabuns.csv")

abunslist <- list(bishabuns, br5abuns, bristyabuns, facabuns, gmabuns, gpabuns, h14abuns, h5abuns, lmpabuns,
                melabuns, mpdabuns, mpeabuns, mpgabuns, pizabuns, primoabuns, shabuns, stsabuns)

for (i in 1:length(abunslist)) {
  rownames(abunslist[[i]]) <- abunslist[[i]][1]
  abunslist[[i]] <- abunslist[[i]][-c(1)]
}

for (i in 1:length(abunslist)) {
  abunslist[[i]] <- abunslist[[i]][,order(colnames(abunslist[[i]]))]
}


#### Community composition analyses ####

mds_sad <- metaMDS(sad[c(1:3, 5:18), 2:31], try = 200, distance = "bray", k = 2, autotransform = T)
mds_sad_s <- as.data.frame(scores(mds_sad)$sites)
mds_sad_s$surfaces <- sad$Surface[c(1:3, 5:18)]
ggplot(mds_sad_s, aes(x = NMDS1, y = NMDS2, col = surfaces)) + geom_point(shape=21, size = 4) + theme_bw()

#### Succession analyses ####

abunsurfacelist <- list()

for (i in 1:length(sur_df_list)) {
  abunsurfacelist[[i]] <- getabun(sur_df_list[[i]])
}

relabunsurfacelist <- list()

for (i in 1:length(abunsurfacelist)) {
  relabunsurfacelist[[i]] <- getrelabun(abunsurfacelist[[i]])
}

W_metrics <- list()

for (i in 1:length(relabunsurfacelist)) {
  W_metrics[i] <- ABCplot(sur_df_list[[i]], relabunsurfacelist[[i]])[[2]]
}

sur_names <- list()

for(i in 1:length(sur_df_list)) {
  sur_names[i] <- sur_df_list[[i]]$Surface[1]
}

w_box <- as.data.frame(cbind(unlist(sur_names), unlist(W_metrics)))

ggplot(w_box, aes(x = V1, y = as.numeric(V2))) + geom_boxplot()

#### Tiering analyses ####

tiermetrics <- data.frame(dvsh = rep(1, length(hd_list)), dvsu = rep(1,length(hd_list)))

for (i in 1:length(hd_list)) {
  tiermetrics$dvsh[i] <- mean(count.sp.vs.comm2(get.height.table(hd_list[[i]][-c(1)]), abunslist[[i]]))
  tiermetrics$dvsu[i] <- mean(count.at.vs.comm2(hd_list[[i]][-c(1)], abunslist[[i]]))
}

#### Custom functions ####

getabun <- function(a) {
  test <- as.data.frame(matrix(ncol = nrow(as.data.frame(levels(as.factor(a$Sp_ID)))),
                               nrow = nrow(as.data.frame(levels((as.factor(a$Surface)))))))
  
  row.names(test) <- as.vector(levels((as.factor(a$Surface))))
  colnames(test) <- as.vector(levels(as.factor(a$Sp_ID)))
  
  tecy <- function(x, y, z){
    fishy <- table(a$Sp_ID[x$Surface == y])
    z[row.names(z) %in% y, (intersect((names(fishy)), names(z)))] <- as.data.frame(fishy)[,2]
    test <<- z
  }
  
  
  listy <- levels(as.factor(a$Surface))
  
  
  for(i in 1:length(listy)) {
    tecy(a, listy[1], test)
  }
  test[is.na(test)] <- 0
  return(test)
}


getrelabun <- function(x) {
  x$totab <- rowSums(x)
  test2 <- x
  for (i in 1:nrow(test2)) {
    test2[rownames(test2[i]),] <- (test2[rownames(test2[i]),]/x[rownames(test2[i]),ncol(x)])*100
  }
  test2 <- test2[1:(ncol(test2)-1)]
  test2 <- test2[complete.cases(test2),]
  return(test2)
}

x <- sur_df_list[[1]]
y <- relabunsurfacelist[[1]]

ABCdf <- function(x, y){
  b <- ddply(x,"Sp_ID",numcolwise(sum))
  b <- b[,c("Sp_ID", "Area")] 
  sar <- sum(b$Area)
  for (i in 1:nrow(b)){
    b$Area[i] <- (b$Area[i]/sar)*100
  }
  
  a <- as.data.frame(t(y))
  a <- cbind(a, rownames(a))
  colnames(a) <- c("relabun", "sp2")
  
  b <- b[order(-b$Area),]
  a <- a[order(-a$relabun),]
  
  return(as.data.frame(c(b, a)))
}

ABCplot <- function(x, y) {
  df1 <- ABCdf(x, y)
  df1$rank <- 1:nrow(df1)
  df1$cum_b <- cumsum(df1$Area)
  df1$cum_a <- cumsum(df1$relabun)
  gplabc <- ggplot() + geom_line(data = df1, aes(x = log(rank), y = log((1+cum_b)/(101-cum_b))), color = "blue") +
    geom_line(data = df1, aes(x = log(rank), y = log((1+cum_a)/(101-cum_a))), color = "red") +
    xlab("Species Rank") + ylab("Abundace (red) vs Biomass (blue)") + theme_bw()
  dub <- (abc(df1[c(2, 3)])@W.Stat[2])
  newlist <- list(gplabc, dub, df1)
  return(newlist)
}

## Tiering functions adapted from Mitchell and Kenchington, 2018, Nature E&E

count.sp.vs.comm2 <- function(c1, c2)
{
  c1<-t(c1)
  res1<-c(1:nrow(c1)) 
  for(j in 1:nrow(c1)) {
    r1<-seq(1,nrow(c1),1)
    r1<-r1[-j]
    res2<-c(1:ncol(c1))
    for(i in 1:ncol(c1)) {
      res2[i]<-ifelse(c1[j,i]>sum(c1[r1,i]),c1[j,i]-sum(c1[r1,i]),0)	
    }
    res1[j]<-sum(res2)/sum(c1[j,])
  }
  res1
  res3 <- rep(res1, c2)
  mean(res3)
  return(res3)
}

count.at.vs.comm2<-function(mat3, ab, breaks=10)#this tells us how much of each species has no overlap with anyother within the community.
{
  c1<-sp.to.zone(get.sp.active.table(mat3,breaks))
  #for(i in 1:length(c2)) {
  #  m <- matrix(nrow = nrow(c1[[i]]), ncol = ncol(c1[[i]]))
  #  m[1:nrow(c1[[i]]), 1:ncol(c1[[i]])] <- c2[[i]]
  #  
  #  c1[[i]] <- c1[[i]] / m
  #}
  res1<-c(1:length(c1)) #where nrow should correspond to the number of abudant species.res
  for(k in 1:length(c1)) {
    c2<-c1
    main.sp<-c1[[k]]
    c2[[k]]<-NULL #removes that species from the list of everything
    comm.mat<-Reduce('+',c2) #all the things added together
    sp.mat<-Reduce('-',list(main.sp,comm.mat)) # subtracts community from species
    res1[k]<-sum(sp.mat[sp.mat>0])/sum(main.sp) #gets the remaining, non overlapped taxa
    
  }
  names(res1)<-names(c1)
  res3 <- rep(res1, ab)
  mean(res3)
  return(res3)# return a matrix col = no of species, which gives the percentage of each taxa specimens which don't overlap with anything else in the comm.
}
