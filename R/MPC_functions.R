##########################################################
### Rapid Evaluation of Multispecies Connectivity (Reconnect) Workflow
### MPC functions
### Functions for calculating MetaPopulation Capacity (MPC) metrics (depend on CONN_functions "eucdist_pm" and "costdist_pm")
### Jacqueline Oehri, jacqueline.oehri@gmail.com
### 18.11.2022
##########################################################

#######
#' @name MPC
#' @title Metapopulation capacity (MPC) core function
#' @description Calculates MPC for a given set of habitat patches. The function relies on the 'sf' and 'stars'packages
#' @param pa A numeric vector of areas of habitat patches 1-n
#' @param mdist A square matrix of pairwise distances between habitat patches 1-n
#' @param alpha A species-specific dispersal distance (e.g. average gap-crossing distance), default=317m
#' @param dispfun A species-specific dispersal survival function, currently 3 options can be chosen by parameter dispfop: "negex", "linear" and "log-sech"
#' @param ex A stochasticity parameter: if x >1, it becomes very unlikely for populations to go extinct after a critical patch size has been reached. If x <1, there exists no critical patch size and populations can even go extinct if they are large (Hanski 1994)
#' @param self Logical: should self-colonization of patches be modeled (default, TRUE) or not (FALSE)
#' @param evec Logical: should the dominant eigenvalue associated eigenvector be returned (default, TRUE) or not (FALSE)
#' @param savememory Logical: should objects that are not necessary be removed from the environment to save memory?
#' @return An 'MPC' object: a list containing a numeric MPC value (mpc), metapopulation density (mpcdens), the MPC-associated squared dominant eigenvector indicating patch importance (pimport) and the patch IDs in case they are indicated (pid)
#' @export
MPC   = function(pa,mdist,alpha,dispfun,ex,self,evec,symmetric,lower,savememory) {

  #matrix of dispersal survival probabilities
  dispr  = as.matrix(dispfun(d=mdist,alpha=alpha))                                                 # Haung et al. 2020: Convert distances measurements to probability of dispersal

  #v2: get names of patches and remove temporary objects to save memory
  pid    = colnames(mdist)
  if(savememory == TRUE) {
    rm(list=c("mdist"))
  }

  #add self-colonization option (1=on, 0=off;Strimas-Mackey & Brodie 2018)
  diag(dispr) = ifelse(self, 1, 0)                                                        # Strimas-Mackey & Brodie 2018: the self is only if you explicitly want to exclude the self; because all current distance functions yield 1 for dist =0...

  #create a matrix for patch areas (Huang et al. 2020)
  fragArea = matrix(pa, ncol=1)                                                           # Huang et al. 2020: actually take areas in km^2
  j.area   = matrix(fragArea,nrow=length(fragArea),ncol = length(fragArea),byrow = TRUE)

  #create the colonization matrix (Huang et al. 2020)
  col.matrix = dispr*j.area

  #v2: remove temporary objects to save memory
  rm(list=c("dispr","j.area"))

  #create a matrix for patch areas (Huang et al. 2020)
  i.area   = matrix(fragArea,nrow=length(fragArea),ncol = length(fragArea),byrow = FALSE)

  #v2: remove temporary objects to save memory
  rm(list=c("fragArea"))

  #create the metapopulation matrix (Huang et al. 2020)
  M = col.matrix*i.area^ex                                                                       # Cf. Eqn 1 in Huang et al. 2020

  #v2: remove temporary objects to save memory
  rm(list=c("col.matrix","i.area"))

  #v2: if the symmetric argument is TRUE, matrix calculation is forced to be fast & based only on one triangle (lower = TRUE by default)
  #v2: bear in mind: the eigs_sym function only works for matrices with at least 3 rows/cols!
  if(ncol(M)>=3){
  if(symmetric==TRUE) {
  MPC0    = RSpectra::eigs_sym(M, k = 1, which = "LM", opts = list(retvec = evec), lower = lower) # Faster eigenvalue calculation
  } else {
  ## calculate metapopulation capacity
  MPC0    = RSpectra::eigs(M, k = 1, which = "LM", opts = list(retvec = evec))                    # Faster eigenvalue calculation (with the largest modulus [euclidean norm])
         }                                                                                        # if matrix is not symmetric, this can result
  } else if(ncol(M)<3) {
    if(symmetric==TRUE) {
      MPC0    = eigen(M, symmetric=TRUE, only.values = !evec)                                     # Strimas-Mackey & Brodie 2018: create associated eigenvector if evec=TRUE
    } else {
      ## calculate metapopulation capacity
      MPC0   = eigen(M, only.values = !evec)                                                      # Strimas-Mackey & Brodie 2018: create associated eigenvector if evec=TRUE
    }                                                                                             # if matrix is not symmetric, this can result
  }

  ## calculate "raw" metapopulation capacity only once
  mpc = Mod(MPC0$values[1])

  ## generate the MPC object
  if (evec) {                                                                                  # Stott et al. 2010: Mod(x), return the Modulus, the simple root of the characteristic determinant of M
    MPC = list(mpc = mpc , mpcdens = mpc/sum(pa,na.rm=TRUE),
                pimport = abs(Mod(MPC0$vectors[, 1])), pid = pid)
                                                                                               # (Mod(x) is in contrast to Strimas-Mackey & Brodie 2018 who used Re(x) to only return the real part of the number)
  } else {                                                                                     # Take the absolute value, in order to (similar to Strimas-Mackey & Brodie 2018: exponent 2) highlight the absolute weight of the eigenvector
    MPC = list(mpc = mpc , mpcdens = mpc/sum(pa,na.rm=TRUE), pimport = NA, pid=NA)
  }

  #v2: remove temporary objects to save memory
  rm(list=c("M","MPC0"))

  return(MPC)
}


#######
#' @name MPC_fun
#' @title Metapopulation capacity (MPC) wrapper function
#' @description Calculates MPC for a given set of habitat patches. The function relies on the 'sf' and 'stars'packages
#' @param x A raster object where habitat patches are set to 1 and matrix area is set to NA OR a shapefile object containing habitat patches
#' @param id A character string, sfp geometry id column describing habitat patch id's
#' @param resi A resistance raster file, should have the same crs and extent as sf
#' @param pa A numeric vector of areas of habitat patches 1-n
#' @param mdist A square matrix of pairwise distances between habitat patches 1-n
#' @param alpha A species-specific dispersal distance (e.g. average gap-crossing distance), default=317m
#' @param dispfop An option to choose the dispersal survival function (dispfun parameter), currently 3 options implemented: "negex", "linear" and "log-sech"
#' @param ex A stochasticity parameter: if x >1, it becomes very unlikely for populations to go extinct after a critical patch size has been reached. If x <1, there exists no critical patch size and populations can even go extinct if they are large (Hanski 1994)
#' @param self Logical: should self-colonization of patches be modeled (default, TRUE) or not (FALSE)
#' @param evec Logical: should the dominant eigenvalue associated eigenvector be returned (default, TRUE) or not (FALSE)
#' @param symmetric Logical: is the pairwise distance matrix symmetric (TRUE) or not (FALSE, default). In case symmetric = TRUE, eigenvalue calculations are faster because only the lower triangle will be considered. If symmetric is not specified, the eigenvalue function will use the whole matrix as input.
#' @param lower Logical: In case the symmetric argument is TRUE: which triangle of the matrix should be used for calculation, the lower left (TRUE, default) or upper right (FALSE)
#' @param areafun A function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @param dist_type A character, distance type: either "edge" or "centroid" (default: "edge").
#' @param cost_type A character, algorithm applied to calculate distances: either "euclidean" (default), "least-cost" or "commute-time"
#' @param trans_fun A function expression, argument for the gdistance::transition() function, applied when moving between resi raster cells. Default for conductance: function(x) 1/mean(x)
#' @param neighb A numeric value, directions argument for the gdistance::transition() function: directions in which cells are connected (4, 8, 16, ...)
#' @param savememory Logical: should objects that are not necessary be removed from the environment to save memory?
#' @return the 'MPC' object created by the function 'MPC'
#' @export
MPC_fun  = function(x=NULL,id=NULL,resi=NULL,pa=NULL,mdist=NULL,alpha=500,dispfop="negex",
                    ex=0.5,self=TRUE,evec=TRUE,symmetric=FALSE,lower=TRUE,
                    areafun=sf::st_area,dist_type="edge",cost_type="euclidean",
                    trans_fun=function(x) 1/mean(x),neighb=8,
                    savememory=FALSE) {

  ## in case x is provided: prepare data: make objects pa (numeric vector of patch areas) and mdist (square matrix of interpatch distances)
   if(class(x)[1]!="NULL"){

    if(cost_type  %in% c("least-cost","commute-time")) {
      ## get the patch area and cost distance matrix object
      pamd = costdist_pm(x=x,resi=resi,id=id,areafun=areafun,
                         trans_fun=trans_fun,neighb=neighb,
                         dist_type=dist_type,cost_type=cost_type)
      pa    = pamd$pa
      mdist = pamd$mdist

    } else if (cost_type == "euclidean") {
      ## get the patch area and euclidean distance matrix object
      pamd  = eucdist_pm(x=x,areafun=areafun,dist_type=dist_type)
      pa    = pamd$pa
      mdist = pamd$mdist
    }

    ## remove to save memory
    rm(list="pamd")

  } # if a habitat raster or shapefile is provided


    ## if x is not specified but mdist and pa are, we can directly go ahead:
    if(class(pa)[1]=="numeric" & class(mdist)[1]=="matrix"& length(pa)==dim(mdist)[1] & dim(mdist)[1]==dim(mdist)[2]) {

    ## define the dispersal function
    dispfun = mkdispfun(option=dispfop)

    ## actually calculate the metapopulation capacity
    MPC_obj = MPC(pa=pa,mdist=mdist,alpha=alpha,dispfun=dispfun,ex=ex,self=self,evec=evec,
                   symmetric=symmetric,lower=lower,savememory=savememory)

  } else {
    print("--- error --- pa or mdist do not conform to one of the following ---")
    print('class(pa)[1]=="numeric" & class(mdist)[1]=="matrix"& length(pa)==dim(mdist)[1] & dim(mdist)[1]==dim(mdist)[2]')
    MPC_obj = "error"
  } # else

  return(MPC_obj)
}


########
#' @name MPCser
#' @title calculate a series of Metapopulation Capacity (MPC) - based metrics
#' @description it uses the same arguments as the MPC_fun with two additional ones: roi and a baseline MPCvalue
#' @param roi shapefile of region of interest
#' @param MPCbasl numeric value of baseline MPCvalue
#' @return a list with MPC based metrics: MPCraw = raw MPC value, MPCev = eigenvector for each patch ID, MPCdens=MPC density, MPCrmax=MPC relative to maximum potential.
#' @export
MPCser = function(roi=NULL,MPCbasl=NULL,
                  x=NULL,id=NULL,resi=NULL,pa=NULL,mdist=NULL,alpha=500,dispfop="negex",
                  ex=0.5,self=TRUE,evec=TRUE,symmetric=FALSE,lower=TRUE,
                  areafun=sf::st_area,dist_type="edge",cost_type="euclidean",
                  trans_fun=function(x) 1/mean(x),neighb=8,
                  savememory=FALSE){

                   mpc     = MPC_fun(x=x,id=id,resi=resi,pa=pa,mdist=mdist,alpha=alpha,dispfop=dispfop,
                                     ex=ex,self=self,evec=evec,symmetric=symmetric,lower=lower,
                                     areafun=areafun,dist_type=dist_type,cost_type=cost_type,
                                     trans_fun=trans_fun,neighb=neighb,
                                     savememory=savememory);
                   MPCraw  = mpc$mpc
                   MPCdens = mpc$mpcdens
                   MPCev   = mpc$pimport
                   MPCevid = mpc$pid
                   if(length(roi)>0){
                   MPCmax  = MPC_fun(x=roi,alpha=alpha,evec = FALSE)[["mpc"]]
                   MPCrmax = MPCraw/MPCmax
                   } else {
                   MPCmax  = NA
                   MPCrmax = NA}
                   if(length(MPCbasl)>0)  {
                   MPCrbasl = MPCraw/MPCbasl} else {MPCrbasl = NA}
                   MPCser = list(MPCraw  =MPCraw,
                                 MPCdens =MPCdens,
                                 MPCev   =MPCev,
                                 MPCevid =MPCevid,
                                 MPCmax  =MPCmax,
                                 MPCrmax =MPCrmax,
                                 MPCrbasl=MPCrbasl)
                   return(MPCser)
                   }



##########################################################
##### project:       Multispecies Connectivity Modelling
##### author:        Jacqueline Oehri (JO), jacqueline.oehri@gmail.com
##########################################################
##### References:
#####                 1) Hanski, Ilkka, and Otso Ovaskainen. 2000. “The Metapopulation Capacity of a Fragmented Landscape.” Nature 404 (6779): 755–58. doi:10.1038/35008063.
#####                 2) Schnell, Jessica K., Grant M. Harris, Stuart L. Pimm, and Gareth J. Russell. 2013. “Estimating Extinction Risk with Metapopulation Models of Large-Scale Fragmentation.” Conservation Biology 27 (3): 520–30. doi:10.1111/cobi.12047.
#####                 3) Strimas‐Mackey, M., & Brodie, J. F. (2018). Reserve design to optimize the long‐term persistence of multiple species. Ecological Applications, 28(5), 1354-1361. https://doi.org/10.1002/eap.1739; https://github.com/mstrimas/metacapa;
#####                 4) Huang, R., Pimm, S. L., & Giri, C. (2020). Using metapopulation theory for practical conservation of mangrove endemic birds. Conservation Biology, 34(1), 266-275. https://doi.org/10.1111/cobi.13364
#####                 5) Hanski, I. (1994). A practical model of metapopulation dynamics. Journal of animal ecology, 151-162.
#####                 6) Stott, I., Townley, S., Carslake, D., & Hodgson, D. J. (2010). On reducibility and ergodicity of population projection matrix models. Methods in Ecology and Evolution, 1(3), 242-252.
#########################################################


