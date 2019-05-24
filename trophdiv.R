###########################################################################################################
# function to compute trophic diversity indices                                                           #
# following Villéger et al. 2008 (Marine Ecology Progress Series, 364: 135-146)                           #
#                                                                                                         #
#  inputs :                                                                                               #
#  - ab: a matrix (C x S) of the abundances of S species in the C communities of interest                 #
#  - tl: a vector (length S) with the trophic levels of the S species of interest                         #
#     -> names of species in tl and ab must be the same                                                   #
#     -> NA are not allowed in 'tl', and are ignored in 'ab'                                              #
#                                                                                                         #
# outputs: a dataframe of 10 columns with trophic diversity indices values for each community             #
#          - total abundance (abtot)                                                                      #
#          - number of species (nbsp)                                                                     #
#          - number of trophic levels (nbtl), trophic richness                                            #
#          - minimum trophic levels present (mintl)                                                       #
#          - maximum trophic levels present (maxtl)                                                       #
#          - range of trophic levels present (rgetl), trophic richness                                    #
#          - mean trophic level (meantl), MTI of Pauly & Watson 2005                                      #
#          - abundance weighted standard deviation of trophic levels (sdtl)                               #                                             
#          - FDvar (Mason et al. 2003), trophic divergence                                                #                            
#          - FROm (modifed from Mouillot et al. 2005), trophic evenness                                   #
#                                                                                                         #
#      FROm computation requires at least 3 different trophic levels                                      #
#                                                                                                         #
###########################################################################################################


trophdiv<-function(ab,tl) {


# check if dimensions of objects correspond
C<-dim(ab)[1]
S<-dim(ab)[2]
if (length(tl)!=S) stop("error : number of species different in 'ab' and 'tl'")
if (length(which(is.na(tl)==T))>0) stop("error : NA are not allowed in 'tl'")

# species names
if (length(which(names(tl)!=colnames(ab)))>0) stop("error : species names are different in 'ab' and 'tl'")

# vector for indices values
nbsp<-rep(NA,C) ; abtot<-rep(NA,C)
nbtl<-rep(NA,C) ; mintl<-rep(NA,C) ; maxtl<-rep(NA,C) ; rgetl<-rep(NA,C)
meantl<-rep(NA,C)
sdtl<-rep(NA,C) ; FDvar<-rep(NA,C)
FROm<-rep(NA,C)

# loop to compute trophic diversity on each of the C communities
for (k in 1:C) {

  # filter to keep only abundances of species present and corresponding names and TL values
  presk<-which(ab[k,]>0)
  abkp<-ab[k,presk]
  tlkp<-tl[presk]

  # total abundance
  abtot[k]<-sum(abkp)

  # number of species
  nbsp[k]<-length(abkp)

  # number of trophic levels
  nbtl[k]<-length(unique(tlkp))

  # minimum trophic level
  mintl[k]<-min(tlkp)

  # maximum trophic level
  maxtl[k]<-max(tlkp)

  # range of trophic level
  rgetl[k]<-maxtl[k]-mintl[k]

  # computation of relative abundances
  abrel<-abkp/sum(abkp)

  # abundance weighted mean trophic level
  meantl[k]<-round( tlkp%*%abrel, 3)

  # abundance weighted standard deviation of trophic level
  sdtl[k]<-round( ( ((tlkp)^2%*%abrel)-meantl[k]^2   )^0.5 , 3)

  # FDvar
  logtl<-log(tlkp)
  V<- (logtl^2%*%abrel)  -  (logtl%*%abrel)^2
  FDvar[k]<- round( 2/pi*atan(5*V) , 3)

  # FROm  (if at leat 3 different TL)
      if (nbtl[k]>2)   {
      o<-order(tlkp,decreasing=F) ; to<-tlkp[o] ; bo<-abkp[o]/sum(abkp) ; s<-length(bo) ; os<-1/(s-1)
      EW<- abs(to[-1]-to[-s]) / (bo[-1]+bo[-s]) ; PEW<-EW/sum(EW)
      minPEW<-sapply(PEW, function(x) { min(x,os) } )
      FROm[k]<-round( ( sum(minPEW)-os ) / ( 1- os ) , 3)     } # end of if nbtl>2

} # end of k

res<-data.frame(abtot=abtot, nbsp=nbsp, nbtl=nbtl, mintl=mintl, maxtl=maxtl, rgetl=rgetl, meantl=meantl, sdtl=sdtl, FDvar=FDvar, FROm=FROm)
row.names(res)<-row.names(ab)
invisible(res)

                        } # end of function trophdiv

################################################################################
# example 

showexample<-F

if (showexample==T) {
abex<-matrix(sample(c(round(runif(15,0,100)),rep(0,8),NA)),4,6,dimnames=list(paste("com",1:4,sep=""),paste("sp",1:6,sep="")))
tlex<-round(runif(6,2,4.5),1) ; names(tlex)<-paste("sp",1:6,sep="")
ex<-trophdiv(abex,tlex)
abex ; tlex ; ex
} # end of if showexample