##########################################################
### Rapid Evaluation of Multispecies Connectivity (Reconnect) Workflow
### CONN functions
### Functions related to landscape patch statistics and connectivity
### Jacqueline Oehri, jacqueline.oehri@gmail.com
### 18.11.2022
##########################################################


#######
#' @name costdist_pm
#' @title Pairwise cost distance matrix of habitat patches in a resistance landscape
#' @description this is a wrapper function combining the nice suite of [gdistance](https://cran.r-project.org/web/packages/gdistance/index.html) functions
#' @param x A raster object where habitat patches are set to 1 and matrix area is set to NA OR a shapefile object containing habitat patches
#' @param resi resistance raster file, should have the same crs and extent as sf
#' @param id character string, sfp geometry id column describing habitat patch id's
#' @param areafun function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @param trans_fun function expression, argument for the gdistance::transition() function, applied when moving between resi raster cells. Default for conductance: function(x) 1/mean(x)
#' @param neighb numeric, directions argument for the gdistance::transition() function: directions in which cells are connected (4, 8, 16, ...)
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @param cost_type character, gdistance algorithm applied to calculate distances: either "least-cost" or "commute-time" (default: "least-cost", the randomised shortest paths, "gdistance::rSPDistance()" is not implemented)
#' @return A list with an sf object of the habitat patches (sfpa) and a corresponding square matrix (cosdi) where pairwise cost distances among habitat patches are described according to the above defined criteria
#' @export
costdist_pm = function(x=NULL,resi=NULL,id=NULL,areafun=sf::st_area,
                       trans_fun=function(x) 1/mean(x),neighb=8,
                       dist_type="edge",cost_type="least-cost") {

  ## change resitance raster values depending on the distance type!
  if(dist_type=="edge") {habval = 0} else if(dist_type=="centroid") { habval = 1 }

  ## change resistance layer where habitat patches are, use a habitat raster or shapefile - depending on what is provided!
  if(class(x)[1]=="RasterLayer") {
    ## set habitat cells zero in the resistance layer! (Assuming habitat cells have a value >0)
    resi[x>0] = habval
    ## generate an sf object
    sfp   = rstosf(rs=x)
    rm(list="x")
  } else if(class(x)[1]=="sf") {
    sfp = x
    rm(list="x")
    ## extract cells of habitat and set these cells zero in the resistance raster!
    resi2   = resi
    resi2[] = NA
    cells   = raster::extract(resi2, as_Spatial(sfp$geometry),method="simple",cellnumbers=TRUE) %>% unlist
    if(length(cells)>0){
    resi[cells[is.finite(cells)]] = habval
    }
    rm(list="resi2")
  }

  ## in case no id is provided by the user, just use the name of the first column
  if(class(id)[1]!="character") {
    if(length(names(sfp)[names(sfp)!="geometry"])>0) {id = names(sfp)[names(sfp)!="geometry"][1]
    } else {sfp = sf::st_sf(pid=1:length(sfp$geometry),geometry=sfp$geometry); id = "pid"}
  } # if class(id)

  ## only conduct the following if sfp has more than 1 geometry!
  if(length(sfp$geometry)>1) {
  ## apply gdistance functions
  ## create transition matrix
  tr  = gdistance::transition(x=resi,transitionFunction=trans_fun,directions=neighb)
  ## apply different geocorrection for projected or lat lon rasters
  if(!raster::isLonLat(resi)) {
    trC  = gdistance::geoCorrection(tr, type="c", multpl=FALSE, scl=TRUE)
  } else if(raster::isLonLat(resi)) {
    trC  = gdistance::geoCorrection(tr, type="r", multpl=FALSE, scl=TRUE)
  }

  ## gives a sparse matrix (same as distr object above...)
  ## the accuracty of the resistance distance depends on the number of neighbours, the further away the more accurate get measurements...
  if(cost_type=="least-cost"){
    cosdi =  tryCatch({
                      result = gdistance::costDistance(trC,as_Spatial(sf::st_centroid(sfp$geometry)))
                      },
                      error = function(e) {
                        message('trycatch-error cosdi')
                        print(e)
                      }
                      )
  } else if(cost_type=="commute-time") {
    cosdi =  tryCatch({
                    result = gdistance::commuteDistance(trC,as_Spatial(sf::st_centroid(sfp$geometry)))
                      },
                    error = function(e) {
                    message('trycatch-error cosdi')
                    print(e)
                    }
                    )
  }

  ## if gdistance does not work, return NA
  if(class(cosdi)[1]!="simpleError") {
  ## turn into a matrix and add correct row and column names
  cosdi = as.matrix(cosdi)
  } else { cosdi = matrix(NA,nrow=length(sfp[[id]]),ncol=length(sfp[[id]]))}
  # sf has more than one patch
  } else { cosdi = as.matrix(0) }
  ## give names to distance matrix
  colnames(cosdi) = rownames(cosdi) = sfp[[id]]
  ## extract patch area vector
  sfpa = as.numeric(areafun(sfp))
  return(list(pa=sfpa,mdist=cosdi))
} # function end


#######
#' @name costdist_lm
#' @title Resistance distance to habitat in landscape
#' @description this is a wrapper function for the terra::costDist function
#' @param x A raster object where habitat patches are set to 1 and matrix area is set to NA OR a shapefile object containing habitat patches
#' @param resi resistance raster file, should have the same crs and extent as sf
#' @param scale numeric value: scale factor in the terra::costDist function - the cost distance is divided by this number
#' @param maxiter numeric value: maximum number of iterations used by the terra::costDist function
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @param cost_type character, gdistance algorithm applied to calculate distances: at the moment, only "least-cost" (default) available.
#' @return A list with an sf object of the habitat patches (sfpa) and a corresponding square matrix (cosdi) where pairwise cost distances among habitat patches are described according to the above defined criteria
#' @export
costdist_lm = function(x=NULL,resi=NULL,
                       scale=1,maxiter=100,
                       dist_type="edge",cost_type="least-cost") {

  ## make the habitat value
  habval = 0

  ## change resistance layer where habitat patches are, use a habitat raster or shapefile - depending on what is provided!
  if(class(x)[1]=="RasterLayer") {
    ## change resitance raster values depending on the distance type!
    if(dist_type=="edge") {
    ## set habitat cells zero in the resistance layer! (Assuming habitat cells have a value >0)
    resi[x>0] = habval
    }
    ## generate an sf object
    sfp   = rstosf(rs=x)
    rm(list="x")
    ##
  } else if(class(x)[1]=="sf") {
    sfp = x
    rm(list="x")
    if(dist_type=="edge") {
    ## set habitat cells zero in the resistance layer! (Assuming habitat cells have a value >0)
    ## extract cells of habitat and set these cells zero in the resistance raster!
    resi2   = resi
    resi2[] = NA
    cells   = raster::extract(resi2, as_Spatial(sfp$geometry),method="simple",cellnumbers=TRUE) %>% unlist
    resi[cells[is.finite(cells)]] = habval
    rm(list="resi2")
    }
  }

  if(dist_type=="centroid") {
    ## set habitat cells zero in the resistance layer! (Assuming habitat cells have a value >0)
    ## extract cells of habitat and set these cells zero in the resistance raster!
    resi2   = resi
    resi2[] = NA
    cells   = raster::extract(resi2, as_Spatial(sf::st_centroid(sfp$geometry)),method="simple",cellnumbers=TRUE) %>% unlist
    resi[cells[is.finite(cells)]] = habval
    rm(list="resi2")
  }

  ## apply cost distance function terra
  if(cost_type=="least-cost"){
  ## interestingly, terra package ALSO has the cost distance
  tr  = terra::costDist(terra::rast(resi), target=habval, scale=scale, maxiter=maxiter)
  #rm(list="resi")
  #tr  = as(tr, "Raster")
  #plot(tr)
  #text(tr)
  ## turn into a matrix and add correct row and column names (JO: I checked that in the terra raster it is the same sequence as in raster package!!)
  cosdi = values(tr)
  colnames(cosdi) = "costdist"
  rownames(cosdi) = 1:length(cosdi)
  #rm(list="tr")
  return(cosdi)
  }

} # function end


########
#' @name eucdist_pm
#' @title Pairwise euclidean distance between habitat patches
#' @description convert a raster or shapefile into patch area and mdist matrix objects
#' @param x a raster or shapefile describing habitat patches
#' @param id character string, sfp geometry id column describing habitat patch id's
#' @param areafun function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @return A list with a vector of habitat patch areas (pa) and a distance matrix (mdist)
#' @export
eucdist_pm = function(x=NULL,id=NULL,areafun=sf::st_area,dist_type="edge") {

  ## prepare data: make objects pa (numeric vector of patch areas) and mdist (square matrix of interpatch distances)
  if(class(x)[1]!="NULL"){
    ## change resistance layer where habitat patches are, use a habitat raster or shapefile - depending on what is provided!
    if(class(x)[1]=="RasterLayer") {
      ## generate an sf object
      sfp   = rstosf(rs=x)
      rm(list="x")
    } else if(class(x)[1]=="sf") {
      sfp = x
      rm(list="x")
    }

    ## in case no id is provided by the user, just use the name of the first column
    if(class(id)[1]!="character") {
      if(length(names(sfp)[names(sfp)!="geometry"])>0) {id = names(sfp)[names(sfp)!="geometry"][1]
      } else {sfp = sf::st_sf(pid=1:length(sfp$geometry),geometry=sfp$geometry); id = "pid"}
    } # if class(id)

    ## calculate patch areas
    pa    = as.numeric(areafun(sfp))

    ## calculate distance matrix
    if(dist_type=="edge") {
      mdist = sf::st_distance(sfp) } else if(dist_type=="centroid") {
        mdist = sf::st_distance(sf::st_centroid(sfp))
      }
    mdist = matrix(mdist, dim(mdist)[1], dim(mdist)[2])
    colnames(mdist) = rownames(mdist) = sfp[[id]]

    ## return the list of euclidean edge or centroid distances
    return(list(pa=pa,mdist=mdist))

  } # !is.null
} # function end


#######
#' @name eucdist_lm
#' @title Euclidean distance to habitat in landscape
#' @description this is a much simplified wrapper function for the terra::distance function
#' @param x a raster object where habitat patches are set to 1 and matrix area is set to NA
#' @param exclude numeric value of the cells that should not be considered for computing distances, used by the terra::distance function.
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @return A list with an sf object of the habitat patches (sfpa) and a corresponding square matrix (cosdi) where pairwise cost distances among habitat patches are described according to the above defined criteria
#' @export
eucdist_lm = function(x=NULL,
                       exclude=NULL,
                       dist_type="edge") {

  ## change resistance layer where habitat patches are, use a habitat raster or shapefile - depending on what is provided!
  if(class(x)[1]=="RasterLayer") {

    ## change raster values depending on the distance type!

    if(dist_type=="edge") {
      ## set habitat cells NA (Assuming habitat cells have a value >0)
      x[x<1] = NA
    }

    if(dist_type=="centroid") {
      ## set habitat cells zero in the resistance layer! (Assuming habitat cells have a value >0)
      ## generate an sf object
      sfp     = rstosf(rs=x)
      cells   = raster::extract(x, as_Spatial(sf::st_centroid(sfp$geometry)),method="simple",cellnumbers=TRUE) %>% unlist
      x[] = NA
      x[cells[is.finite(cells)]] = 1
    }

  } # raster

    ## interestingly, terra package ALSO has the cost distance
    tr  = terra::distance(terra::rast(x), target=NA,exclude=exclude)

    ## turn into a matrix and add correct row and column names (JO: I checked that in the terra raster it is the same sequence as in raster package!!)
    cosdi = values(tr)
    colnames(cosdi) = "eucdist"
    rownames(cosdi) = 1:length(cosdi)
    #rm(list="tr")
    return(cosdi)
  } # function end


#######
#' @name eucdist_hab
#' @title Euclidean distance from within habitat to habitat edge
#' @description this is a much simplified wrapper function for the terra::distance function
#' @param x a raster object where habitat patches are set to 1 and matrix area is set to NA
#' @param exclude numeric value of the cells that should not be considered for computing distances, used by the terra::distance function.
#' @return A list with an sf object of the habitat patches (sfpa) and a corresponding square matrix (cosdi) where pairwise cost distances among habitat patches are described according to the above defined criteria
#' @export
eucdist_hab = function(x=NULL,exclude=NULL) {
## change resistance layer where habitat patches are, use a habitat raster or shapefile - depending on what is provided!
if(class(x)[1]=="RasterLayer") {
  ## set all habitat to one value
  x[x>=1] = 10
  ## set all non-habitat values to 1
  x[x<1|is.na(x)]  = 1
  ## set all habitat to NA
  x[x==10] = NA

  ## interestingly, terra package ALSO has the cost distance
  tr  = terra::distance(terra::rast(x),target=NA,exclude=exclude)

  ## turn into a matrix and add correct row and column names (JO: I checked that in the terra raster it is the same sequence as in raster package!!)
  cosdi = values(tr)
  colnames(cosdi) = "eucdist"
  rownames(cosdi) = 1:length(cosdi)
  #rm(list="tr")
  return(cosdi)

  } # raster
} # function end


########
#' @name rstosf
#' @title convert habitat raster to shapefile with habitat patches
#' @description convert a binary habitat distribution raster file into a shapefile
#' @param rs a raster file describing habitat patches
#' @return A shapefile with patch ids in column "pid"
#' @export
rstosf = function(rs=NULL) {
  ## prepare data: make objects pa (numeric vector of patch areas) and mdist (square matrix of interpatch distances)
  if(class(rs)[1]!="NULL"){
    ## if rs is a raster layer, make an sf object out of it..
    if(class(rs)[1]=="RasterLayer") {
      # either make clumps or be aware that it works only if the raster has 1 and NA..clumps depend on landscapemetrixs, but advantage of unique IDs...
      if(max(values(rs),na.rm=TRUE)==1){
        rs  = rs
        rs[rs<1|!is.finite(rs)] = NA
        rs   = landscapemetrics::get_patches(rs, directions=4, class=1)[[1]][[1]]  ## the landscapemetrics function works only when non-habitat cells are set to NA...
      } # making sure to get habitat raster clumped into different patches
      rs   = stars::st_as_stars(rs) %>% sf::st_as_sf(merge = TRUE)    # if it is a raster, make a shapefile out of it
      rs   = sf::st_make_valid(rs,reason=TRUE)                        # make sure polygons are correctly extracted
    } # if raster
    names(rs)[1] = "pid"
    return(rs)
  } # !is.null rs
}


#######
#' @name PCnum
#' @title Probability of Connectivity numerator (PCnum)
#' @description Calculates the numerator of the probability of connectivity index (PC) which is the probability that two points randomly placed in a landscape fall into habitat areas that are interconnected/reachable from each other (cf. [Saura & Pascual-Hortal(2007)](https://www.sciencedirect.com/science/article/pii/S0169204607000709), [Saura and Rubio (2010)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1600-0587.2009.05760.x) ).
#' @param pa A numeric vector of areas of habitat patches 1-n
#' @param mdist A square matrix of pairwise distances between habitat patches 1-n
#' @param alpha A species-specific dispersal distance (e.g. average gap-crossing distance, in the units of the mdist values, default = 500)
#' @param dispfop An option to choose the dispersal survival function (dispfun parameter), currently 3 options implemented: "negex", "linear" and "log-sech".
#' @return A numeric value of the probability of connectivity numerator value "PCnum" in Saura & Rubio (2010).
#' @export
PCnum = function(pa=NULL,mdist=NULL,alpha=500,dispfop="negex",savememory=FALSE) {

  ## define the dispersal function
  dispfun = mkdispfun(option=dispfop)

  #matrix of dispersal survival probabilities
  #dispr       = as.matrix(dispfun(mdist))                                                 # Haung et al. 2020: Convert distances measurements to probability of dispersal
  #matrix of dispersal survival probabilities
  dispr  = as.matrix(dispfun(d=mdist,alpha=alpha))                                                 # Haung et al. 2020: Convert distances measurements to probability of dispersal

  #(diagonal is automatically = 1 for the above functions)
  #diag(dispr) = 1

  #v2: get names of patches and remove temporary objects to save memory
  pid    = colnames(mdist)

  ## remove matrix if potentially large
  if(savememory == TRUE) {rm(list=c("mdist"))}

  ## This notation was largely copied from the metapopulation capacity function, this notation allows to save RAM by sequentially removing objects not needed anymore..
  #create a matrix for patch areas (Huang et al. 2020)
  fragArea = matrix(pa, ncol=1)                                                           # Huang et al. 2020: actually take areas in km^2
  j.area   = matrix(fragArea,nrow=length(fragArea),ncol = length(fragArea),byrow = TRUE)
  #create a matrix for patch areas (Huang et al. 2020)
  i.area   = matrix(fragArea,nrow=length(fragArea),ncol = length(fragArea),byrow = FALSE)

  #v2: remove temporary objects to save memory
  rm(list=c("fragArea"))

  #create the colonization matrix (Huang et al. 2020)
  col.matrix = dispr*j.area

  #v2: remove temporary objects to save memory
  rm(list=c("dispr","j.area"))

  #create the outer product of i and j areas multiplied by maximum dispersal probability (i.e. dispersal probability of )
  M = col.matrix*i.area                                                                    # Cf. Eqn 1 in Huang et al. 2020

  #v2: remove temporary objects to save memory
  rm(list=c("col.matrix","i.area"))

  ## calculate the PCnum metric (Saura & Rubio, 2010)
  PCnum = sum(M)

  ## return
  return(PCnum)

}


#######
#' @name ECA_pc
#' @title Equivalent Connected Area index (ECA) based on the probability of connectivity (PCnum)
#' @description The ECA, calculated as the  square root of the numerator of the probability of connectivity index (PCnum), can be interpreted as the size of a single, maximally connected habitat patch that would provide the same probability of connectivity value than the provided habitat configuration (cf. [Saura, Estreguil, Mouton & Rodríguez-Freire 2011](https://doi.org/10.1016/j.ecolind.2010.06.011)).
#' @param pa A numeric vector of areas of habitat patches 1-n
#' @param mdist A square matrix of pairwise distances between habitat patches 1-n
#' @param alpha A species-specific dispersal distance (e.g. average gap-crossing distance, in the units of the mdist values, default = 500)
#' @param dispfop An option to choose the dispersal survival function (dispfun parameter), currently 3 options implemented: "negex", "linear" and "log-sech"
#' @return A numeric value of the Equivalent Connected Area index (ECA) based on the probability of connectivity (PCnum)
#' @export
ECA_pc = function(pa=NULL,mdist=NULL,alpha=500,dispfop="negex",savememory=FALSE,pcnum=NULL) {

  ## derive ECA directly from the numerator of the probability of connectivity index (cf. Saura, Estreguil, Mouton & Rodríguez-Freire 2011)

  if(class(pa)[1]!="numeric" & class(mdist)[1]!="matrix" & class(pcnum)[1]=="numeric") {
    ECA = sqrt(pcnum)
  } else if (class(pa)[1]=="numeric" & class(mdist)[1]=="matrix" & class(pcnum)[1]!="numeric") {
    ECA = sqrt(PCnum(pa=pa,mdist=mdist,alpha=alpha,dispfop=dispfop,savememory=savememory))
  }

  ## return
  return(ECA)
}


#######
#' @name PCECA_fun
#' @title Probability of Connectivity wrapper function
#' @description This function works either with a habitat raster or shapefile (x), can include a resistance raster (resi) or a vector with patch area's (pa) and a distance matrix (mdist) see [Saura & Pascual-Hortal(2007)](https://www.sciencedirect.com/science/article/pii/S0169204607000709)
#' @param x  habitat raster or shape file, should be binary if raster: 1= habitat and 0= non-habitat
#' @param id character string, sfp geometry id column describing habitat patch id's
#' @param resi resistance raster file, should have the same crs and extent as sf
#' @param pa numeric vector of areas of habitat patches 1-n
#' @param mdist square matrix of pairwise distances between habitat patches 1-n
#' @param AL numeric value, necessary if resi=NULL; should describe landscape area in units^2 of pa
#' @param areafun function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @param cost_type character, algorithm applied to calculate distances: either "euclidean" (default), "least-cost" or "commute-time"
#' @param trans_fun function expression, argument for the gdistance::transition() function, applied when moving between resi raster cells. Default for conductance: function(x) 1/mean(x)
#' @param neighb numeric, directions argument for the gdistance::transition() function: directions in which cells are connected (4, 8, 16, ...)
#' @param alpha A species-specific dispersal distance (e.g. average gap-crossing distance, in the units of the mdist values, default = 500)
#' @param dispfop An option to choose the dispersal survival function (dispfun parameter), currently 3 options implemented: "negex", "linear" and "log-sech"
#' @return A list of indices derived from the probability of connectivity numerator "PCnum" in [Saura & Rubio, 2010](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1600-0587.2009.05760.x). PC: Probability of Connectivity, ECA: Equivalent Connected Area index [Saura, Estreguil, Mouton & Rodríguez-Freire 2011](https://doi.org/10.1016/j.ecolind.2010.06.011), ECA_AL: ECA relative to maximum possible (if ECA_AL is multiplied with 100, it is the ProtConn index [Saura et al. (2017)](https://doi.org/10.1016/j.ecolind.2016.12.047)), ECA_AP: ECA relative to habitat amount in landscape (fraction of habitat area that is effectively connected), AL: Landscape area, AP: total patch area
#' @export
PCECA_fun= function(x=NULL,id=NULL,resi=NULL,pa=NULL,mdist=NULL,AL=NULL,areafun=sf::st_area,dist_type="edge",
                    cost_type="euclidean", trans_fun=function(x) 1/mean(x),neighb=8,
                    alpha=500,dispfop="negex",savememory=FALSE) {

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
      pamd = eucdist_pm(x=x,areafun=areafun,dist_type=dist_type)
      pa    = pamd$pa
      mdist = pamd$mdist
    }

    ## remove to save memory
    rm(list="pamd")

  } # if a habitat raster or shapefile is provided

  ## if x is not specified but mdist and pa are, we can directly go ahead:
  if(class(pa)[1]=="numeric" & class(mdist)[1]=="matrix"& length(pa)==dim(mdist)[1] & dim(mdist)[1]==dim(mdist)[2]) {

    ## derive the PCnum
    pcnum = PCnum(pa=pa,mdist=mdist,alpha=alpha,dispfop=dispfop,savememory=savememory)

    ## save the sum of patch areas
    AP = sum(pa)

    ##get the landscape area from the extent of the habitat patches file (raster or sf)
    if(class(AL)[1]!="numeric"){
      if(class(x)[1]=="RasterLayer"|class(x)[1]=="sf"){
        ## make a polygon of extent and extract center coordinates
        LA = sf::st_sf(pid = 1, geometry=sf::st_as_sfc(sf::st_bbox(x)))
        sf::st_crs(LA) = sf::st_crs(x)
        ## get area of whole landscape
        AL = as.numeric(areafun(LA))
      } else if(class(x)[1]!="RasterLayer"&class(x)[1]!="sf"){
        ## if AL is not provided and not derivable from habitat patches file, set it to -1
        AL = -1
      }
    }

    ## calculate the PC
    PC = pcnum/(AL)^2

    ## calculate ECA (Equivalent Connected Area index [Saura, Estreguil, Mouton & Rodríguez-Freire 2011](https://doi.org/10.1016/j.ecolind.2010.06.011))
    ECA = ECA_pc(pcnum=pcnum)

    ## calculate ECA_AL (ECA relative to maximum possible, if ECA_AL is multiplied with 100, it is the ProtConn index [Saura et al. (2017)](https://doi.org/10.1016/j.ecolind.2016.12.047) ).
    ECA_AL = ECA/AL

    ## calculate ECA_AP (ECA relative to habitat area, fraction of habitat area that is effectively connected).
    ECA_AP = ECA/AP

    ## return PCnum, PC, and ECA, with total area in m^2..
    return(list(PCnum=pcnum,PC=PC,ECA=ECA,ECA_AL=ECA_AL,ECA_AP=ECA_AP,AL=AL,AP=AP)) # function end

  } else {
    print("--- error --- pa or mdist do not conform to one of the following ---")
    print('class(pa)[1]=="numeric" & class(mdist)[1]=="matrix"& length(pa)==dim(mdist)[1] & dim(mdist)[1]==dim(mdist)[2]')
  } # else
}


#######
#' @name dI_fun
#' @title Patch Importance
#' @description calculates patch importance for any connectivity metric provided, according to the formulas in e.g. [Saura & Pascual-Hortal(2007)](https://www.sciencedirect.com/science/article/pii/S0169204607000709) and [Saura and Rubio (2010)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1600-0587.2009.05760.x).
#' @param fun a character naming the connectivity metric function name; e.g. fun = "PCnum"
#' @param x a list containing all objects needed by fun, e.g. x  = list(pa=pa,mdist=mdist), whereby pa and mdist are objects in the R-environment
#' @param expr a character string describing the fun statement. Importantly, the user needs to specify "k", i.e. how list objects should  be subsetted; e.g. expression = "pa=x[[1]][-k],mdist=x[[2]][-k,-k] ,alpha=500,dispfop='negex',savememory=FALSE"
#' @param pid character, give an id name, e.g. describing habitat patch id's, the length of pid should match the length of the first element of x (i.e. length(x[[1]]))
#' @param resnr numeric, in case the fun has several result outputs, specify the number of the result output one is interested
#' @param aslist logical: in case the results should  be returned as list and not as matrix (default = FALSE)
#' @return A matrix with patch importance values in rows for each connectivity indicator in columns
#' @export
dI_fun= function(fun=NULL,expr=NULL,x=NULL,pid=NULL,resnr=1,aslist=FALSE) {

  ## make a function statement with the given fun and expr
  ## make sure there are no spaces in the statement for gsub below!
  statement = paste0(fun,"(",gsub(" ","",expr),")")

  ## calculate total metric not excluding any patch
  ## therefore make the "zero statement" which is the same but without k indices
  statement0 = gsub(pattern="\\[[^][]*],",replacement=",",statement,perl=TRUE)
  ## the zero measure, i.e. when all patches are included
  allm0      = eval(parse(text = paste(statement0)))
  ## unlist if list
  if(is.list(allm0)) { allm0 = unlist(allm0)}
  ## select resultnr corresponding to index
  allm0 = allm0[resnr]

  ## apply all at once for all patches, orient yourself on first element in x!
  allk = sapply(1:length(x[[1]]),function(k){eval(parse(text = paste(statement)))})

  ## select resultnr, if only 1 result is selected, unlist allk in case it is still one...
  if(!is.null(dim(allk)[1])) { allk = allk[resnr,]; if(is.null(dim(allk)[1])) { allk = unlist(allk) }}

  if(length(allm0)==1){
    ## calculate patch importance - if there is a single result
    dI = 100*((allm0 - allk)/allm0)
    ## make a  matrix even though its only 1 result
    dI = matrix(dI,ncol=1,nrow=length(dI))
    if(is.null(names(allm0))) {names(allm0)=fun}
    ## name it with the pid if one is given, otherwise just numbering along length of x[[1]]
    # if(length(pid)==length(dI)) { names(dI) = pid }
  } else if(length(allm0)>1){
    try = try(allk[1,],silent=TRUE)
    ## JO: sapply if there are several results
    if(class(try)!="try-error") {
      dI = sapply(1:length(allm0),function(i){ allki = unlist(allk[i,]); 100*((allm0[i] - allki)/allm0[i])})
    } else {
      ## JO: add calculation if there is only 1 patch...(then patch importance should be 100!)
      dI = sapply(1:length(allm0),function(i){ allki = unlist(allk[i]); 100*((allm0[i] - allki)/allm0[i])})
      ## make a  matrix even though its only 1 result
      dI = matrix(dI,ncol=length(dI))
    }
  }

  ## give  names if there are some
  colnames(dI) = names(allm0)
  ## give patch id (pid) names if there are some and if they match
  if(length(pid)==nrow(dI)) {rownames(dI) = pid}

  ## make a list output to be consistent with the Reconnect framework
  if(aslist==TRUE) {
  dIL  = list()
  for(jj in 1:ncol(dI)) {dIL[[jj]] = dI[,jj];if(length(pid)==nrow(dI)) { names(dIL[[jj]]) = pid} }
  names(dIL) = colnames(dI)
  dI = dIL
  rm(list=c("dIL"))
  }

  ## return object as (named) numeric vector or matrix
  return(dI)

} # function end

####### JO: there is a bug when only 1 patch present! Error in allk[i, ] : incorrect number of dimensions
#' @name dI_fun_alpha
#' @title Patch Importance, for metrics where an alpha dispersal distance is indicated
#' @description calculates patch importance for any connectivity metric provided, according to the formulas in e.g. [Saura & Pascual-Hortal(2007)](https://www.sciencedirect.com/science/article/pii/S0169204607000709) and [Saura and Rubio (2010)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1600-0587.2009.05760.x).
#' @param fun a character naming the connectivity metric function name; e.g. fun = "PCnum"
#' @param x a list containing all objects needed by fun, e.g. x  = list(pa=pa,mdist=mdist), whereby pa and mdist are objects in the R-environment
#' @param expr a character string describing the fun statement. Importantly, the user needs to specify "k", i.e. how list objects should  be subsetted; e.g. expression = "pa=x[[1]][-k],mdist=x[[2]][-k,-k] ,alpha=500,dispfop='negex',savememory=FALSE"
#' @param alpha a numeric value of the alpha dispersal capacity
#' @param pid character, give an id name, e.g. describing habitat patch id's, the length of pid should match the length of the first element of x (i.e. length(x[[1]]))
#' @param resnr numeric, in case the fun has several result outputs, specify the number of the result output one is interested
#' @param aslist logical: in case the results should  be returned as list and not as matrix (default = FALSE)
#' @return A matrix with patch importance values in rows for each connectivity indicator in columns
#' @export
dI_fun_alpha = function(fun=NULL,expr=NULL,x=NULL,pid=NULL,resnr=1,alpha=500,aslist=FALSE) {

  ## first replace the alpha placeholder in expression with the corresponding one
  #expr = gsub("alpha=.+,",paste("alpha=",alpha,",",sep=""),gsub(" ","",expr))
  ## make a function statement with the given fun and expr
  ## make sure there are no spaces in the statement for gsub below!
  #statement = paste0(fun,"(",gsub(" ","",expr),paste0(",alpha=",alpha),")")
  statement = paste0(fun,"(",gsub(" ","",expr),")")

  ## calculate total metric not excluding any patch
  ## therefore make the "zero statement" which is the same but without k indices
  statement0 = gsub(pattern="\\[[^][]*],",replacement=",",statement,perl=TRUE)
  ## the zero measure, i.e. when all patches are included
  allm0      = eval(parse(text = paste(statement0)))
  ## unlist if list
  if(is.list(allm0)) { allm0 = unlist(allm0)}
  ## select resultnr corresponding to index
  allm0 = allm0[resnr]

  ## apply all at once for all patches, orient yourself on first element in x!
  allk = sapply(1:length(x[[1]]),function(k){eval(parse(text = paste(statement)))})

  ## select resultnr, if only 1 result is selected, unlist allk in case it is still one...
  if(!is.null(dim(allk)[1])) { allk = allk[resnr,]; if(is.null(dim(allk)[1])) { allk = unlist(allk) }}

  if(length(allm0)==1){
    ## calculate patch importance - if there is a single result
    dI = 100*((allm0 - allk)/allm0)
    ## make a  matrix even though its only 1 result
    dI = matrix(dI,ncol=1,nrow=length(dI))
    if(is.null(names(allm0))) {names(allm0)=fun}
    ## name it with the pid if one is given, otherwise just numbering along length of x[[1]]
    # if(length(pid)==length(dI)) { names(dI) = pid }
  } else if(length(allm0)>1){
    try = try(allk[1,],silent=TRUE)
    ## JO: sapply if there are several results
    if(class(try)!="try-error") {
    dI = sapply(1:length(allm0),function(i){ allki = unlist(allk[i,]); 100*((allm0[i] - allki)/allm0[i])})
    } else {
    ## JO: add calculation if there is only 1 patch...(then patch importance should be 100!)
    dI = sapply(1:length(allm0),function(i){ allki = unlist(allk[i]); 100*((allm0[i] - allki)/allm0[i])})
    ## make a  matrix even though its only 1 result
    dI = matrix(dI,ncol=length(dI))
    }
  }

  ## give  names if there are some
  colnames(dI) = names(allm0)
  ## give patch id (pid) names if there are some and if they match
  if(length(pid)==nrow(dI)) {rownames(dI) = pid}

  ## make a list output to be consistent with the Reconnect framework
  if(aslist==TRUE) {
    dIL  = list()
    for(jj in 1:ncol(dI)) {dIL[[jj]] = dI[,jj];if(length(pid)==nrow(dI)) { names(dIL[[jj]]) = pid} }
    names(dIL) = colnames(dI)
    dI = dIL
    rm(list=c("dIL"))
  }

  ## return object as (named) numeric vector or matrix
  return(dI)

} # function end


#######
#' @name centr_igraph
#' @title Igraph centrality measures
#' @description calculates betweenness centrality (BC), closeness centrality (CC) and Node Degree (ND). This function works either with a habitat raster or shapefile (x), can include a resistance raster (resi) or a vector with patch area's (pa) and a distance matrix (mdist) see [Saura & Pascual-Hortal(2007)](https://www.sciencedirect.com/science/article/pii/S0169204607000709)
#' @param x  habitat raster or shape file, should be binary if raster: 1= habitat and 0= non-habitat
#' @param id character string, sfp geometry id column describing habitat patch id's
#' @param resi resistance raster file, should have the same crs and extent as sf
#' @param pa numeric vector of areas of habitat patches 1-n
#' @param mdist square matrix of pairwise distances between habitat patches 1-n
#' @param areafun function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @param cost_type character, algorithm applied to calculate distances: either "euclidean" (default), "least-cost" or "commute-time"
#' @param trans_fun function expression, argument for the gdistance::transition() function, applied when moving between resi raster cells. Default for conductance: function(x) 1/mean(x)
#' @param neighb numeric, directions argument for the gdistance::transition() function: directions in which cells are connected (4, 8, 16, ...)
#' @param alpha A species-specific dispersal distance (e.g. average gap-crossing distance, in the units of the mdist values, default = 500)
#' @param dispfop An option to choose the dispersal survival function (dispfun parameter), currently 3 options implemented: "negex", "linear" and "log-sech"
#' @return A list of [igraph](https://igraph.org/) R-package derived graph centrality measures: betweenness centrality (BC), closeness centrality (CC) and node degree (ND)
#' @export
centr_igraph = function(x=NULL,id=NULL,resi=NULL,pa=NULL,mdist=NULL,
                        weighted=TRUE,cutoffpr=0.01,MST=FALSE,normalized=TRUE,mode=c("undirected"),
                        areafun=sf::st_area,dist_type="edge",
                        cost_type="euclidean", trans_fun=function(x) 1/mean(x),neighb=8,
                        alpha=500,dispfop="negex",savememory=FALSE) {

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
      pamd = eucdist_pm(x=x,areafun=areafun,dist_type=dist_type)
      pa    = pamd$pa
      mdist = pamd$mdist
    }

    ## remove to save memory
    rm(list="pamd")

  } # if a habitat raster or shapefile is provided

  ## if x is not specified but mdist and pa are, we can directly go ahead:
  if(class(pa)[1]=="numeric" & class(mdist)[1]=="matrix"& length(pa)==dim(mdist)[1] & dim(mdist)[1]==dim(mdist)[2]) {

    ## Igraph: specify the graph
    graph = igraph::graph_from_adjacency_matrix(
      adjmatrix= mdist,
      mode = mode,
      weighted = weighted,
      diag = TRUE,
      add.colnames = NULL,
      add.rownames = NA)

    ## if the minimum spanning tree argument is on, derive the minimum spanning tree
    ## this can save computing time and highlight central nodes better
    if(MST) {
      ## define weights
      if(weighted){  graph = igraph::mst(graph=graph,weights=edge_attr(graph, 'weight'))
      } else { weights = igraph::mst(graph=graph) }
    }

    ## if the cutoff probability is given, remove corresponding edges
    if(!is.null(cutoffpr)) {
      ## define the dispersal function
      dispfun = mkdispfun(option=dispfop)
      ## define cutoff distance for given alpha and dispfun
      cutoffd = bufferdist(alpha=alpha,dispfun=dispfun,prob=cutoffpr)["dist"] %>% as.numeric
      ## remove "isolated" patches
      iso = which(edge_attr(graph, 'weight')>=cutoffd)
      if(length(iso)>0) {graph = delete.edges(graph, iso)}
    } else if(is.null(cutoffpr)){
      cutoffd = -1
    }

    ## define weights
    if(weighted){ weights = edge_attr(graph, 'weight') } else { weights = NULL }

    ## calculate metrics..
    ### normalized: Logical scalar, whether to calculate the normalized closeness, i.e. the inverse average distance to all reachable vertices.
    ## igraph does have specific functions for large graphs - estimate_betweenness and estimate_closeness, which the manual says are not quadratic in runtime. You define a cutoff, which is the largest path length that will be included in the calculation. Traditionally, betweenness considers paths of any length. Defining a cutoff substantially cuts down the runtime:
    ## https://stackoverflow.com/questions/41753929/how-long-does-it-take-for-igraph-r-to-estimate-network-centrality-measures-for-a

    BC = igraph::estimate_betweenness(
      graph,
      v = V(graph),
      weights  = weights,
      cutoff   = cutoffd
      #,normalized=normalized
    )

    CC = igraph::estimate_closeness(
      graph,
      v = V(graph),
      weights = weights,
      cutoff=cutoffd
      #,normalized=normalized
    )

    ###JO: multiply by 100 to get larger values..
    CC = CC*100

    ### degreee
    ND = igraph::degree(graph=graph,v = V(graph), mode="all", normalized = FALSE)

    ### return the graph as well for potential further analyses?

    ## return PCnum, PC, and ECA, with total area in m^2..,
    return(list(BC=BC,CC=CC,ND=ND,graph=graph)) # function end

  } else {
    print("--- error --- pa or mdist do not conform to one of the following ---")
    print('class(pa)[1]=="numeric" & class(mdist)[1]=="matrix"& length(pa)==dim(mdist)[1] & dim(mdist)[1]==dim(mdist)[2]')
  } # else
}


#######
#' @name invCR_p
#' @title Inverse of cumulative resistance in patch network
#' @description invCR is the inverse of the average cumulative resistance of the shortest path between p randomly selected pairs of habitat patches k (pairs are drawn p times; cf. [Albert et al. (2017)](https://doi.org/10.1111/cobi.12943) )
#' @param x  habitat raster or shape file, should be binary if raster: 1= habitat and 0= non-habitat
#' @param id character string, sfp geometry id column describing habitat patch id's
#' @param resi resistance raster file, should have the same crs and extent as sf
#' @param pa numeric vector of areas of habitat patches 1-n
#' @param mdist square matrix of pairwise distances between habitat patches 1-n
#' @param pnr numeric, number of desired samples taken from all patches to calculate average cumulative resistance (default=100). Indicate 1 if all patches should be selected.
#' @param areafun function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @param cost_type character, algorithm applied to calculate distances: either "euclidean" (default), "least-cost" or "commute-time"
#' @param trans_fun function expression, argument for the gdistance::transition() function, applied when moving between resi raster cells. Default for conductance: function(x) 1/mean(x)
#' @param neighb numeric, directions argument for the gdistance::transition() function: directions in which cells are connected (4, 8, 16, ...)
#' @return A numeric value of the inverse of cumulative resistance in patch network (invCR_p)
#' @export
invCR_p = function(x=NULL,id=NULL,resi=NULL,pa=NULL,mdist=NULL,
                 pnr = 1,
                 areafun=sf::st_area,dist_type="edge",
                 cost_type="least-cost", trans_fun=function(x) 1/mean(x),neighb=8,
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

     if(pnr>1){
      pnr = min(pnr,ncol(mdist))
      ## sample mdist, number of pnr times, reduce pnr if nr of patches are lower..
      sample = sample(colnames(mdist),size=pnr)
      idx    = which(colnames(mdist)%in%sample)
      ## select only the sampled ones
      mdist  = mdist[idx,idx]
      pa     = pa[idx]
      }
    ## since distances are already least cost, we can simply take the average..
    diag(mdist) = NA
    ## calculate invCR, scale by 100 to make small numbers better visible
    invCR_p = (1/mean(mdist,na.rm=TRUE))*100

    ## return PCnum, PC, and ECA, with total area in m^2..
    return(invCR_p) # function end

  } else {
    print("--- error --- pa or mdist do not conform to one of the following ---")
    print('class(pa)[1]=="numeric" & class(mdist)[1]=="matrix"& length(pa)==dim(mdist)[1] & dim(mdist)[1]==dim(mdist)[2]')
  } # else
}



#######
#' @name invCR_l
#' @title Inverse of cumulative resistance in landscape
#' @description invCR_l is the inverse of the average cumulative resistance distance to habitat patches in a landscape
#' @param x  habitat raster or shape file, should be binary if raster: 1= habitat and 0= non-habitat
#' @param resi resistance raster file, should have the same crs and extent as sf
#' @param scale numeric value: scale factor in the terra::costDist function - the cost distance is divided by this number
#' @param maxiter numeric value: maximum number of iterations used by the terra::costDist function
#' @param cost_type character, algorithm applied to calculate distances: either "euclidean" (default) or "least-cost"
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @param exclude numeric value of the cells that should not be considered for computing distances, used by the terra::distance function.
#' @return A numeric value of the inverse of cumulative resistance in the landscape (invCR_l)
#' @export
invCR_l = function(x=NULL,resi=NULL,
                   scale=1,maxiter=100,cost_type="euclidean",
                   dist_type="edge",exclude=NULL) {

  ## in case x is provided: prepare data: make objects pa (numeric vector of patch areas) and mdist (square matrix of interpatch distances)
  if(class(x)[1]!="NULL"){

    if(cost_type  %in% c("least-cost")) {
      ## get the patch area and cost distance matrix object
      invCR_l = 1/mean(costdist_lm(x=x,resi=resi,
                         scale=scale,maxiter=maxiter,
                         dist_type=dist_type,cost_type=cost_type),na.rm=TRUE)
    } else if (cost_type == "euclidean") {
      ## get the patch area and euclidean distance matrix object
      invCR_l = 1/mean(eucdist_lm(x=x,exclude=exclude,dist_type=dist_type),na.rm=TRUE)

    }

    return(invCR_l)
  } # if a habitat raster file is provided

} # function end

#######
#' @name perimf
#' @title extract patch perimeter
#' @description get perimeters of each polygon in a shapefile
#' @param shp a shapefile with polygons
#' @return A numeric vector of perimeter of each polygon in shapefile
#' @export
perimf = function(shp) { shp %>% sf::st_cast("MULTILINESTRING") %>% sf::st_length() %>% return }


#######
#' @name NNDf
#' @title extract nearest neighbour distance
#' @description get nearest neighbour distance (NND) for each row in a distance matrix
#' @param mdist a square distance matrix
#' @return A numeric vector of nearest neigbour distance for each row in a distance matrix
#' @export
NNDf = function(mdist) {diag(mdist)=NA;
NND = as.numeric(apply(mdist,MARGIN=1,FUN=function(x){unique(x[which.min(x)])}));
#diag(mdist) = 0 # dont forget to set the mdist again correctly
return(NND)
} # function end


#######
#' @name LandStats
#' @title Make a selection of basic landscape statistical metrics
#' @description Calculate basic landscape statistics for euclidean and cost distance landscapes
#' @param x  habitat raster or shape file, should be binary if raster: 1= habitat and 0= non-habitat
#' @param id character string, sfp geometry id column describing habitat patch id's
#' @param resi resistance raster file, should have the same crs and extent as sf
#' @param pa numeric vector of areas of habitat patches 1-n
#' @param mdist square matrix of pairwise distances between habitat patches 1-n
#' @param areafun function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @param cost_type character, algorithm applied to calculate distances: either "euclidean" (default), "least-cost" or "commute-time"
#' @param trans_fun function expression, argument for the gdistance::transition() function, applied when moving between resi raster cells. Default for conductance: function(x) 1/mean(x)
#' @param neighb numeric, directions argument for the gdistance::transition() function: directions in which cells are connected (4, 8, 16, ...)
#' @return A named vector of landscape statistical metrics: total habitat area (tpa), average and standard deviation of habitat area (mnpa and sdpa), total number of patches (tnp), total nearest neighbour distances (tnnd), average and standard deviation of nearest neighbour distances (mnnnd and sdnnd)
#' @export
LandStats = function(x=NULL,id=NULL,resi=NULL,pa=NULL,mdist=NULL,
                 areafun=sf::st_area,dist_type="edge",
                 cost_type="euclidean", trans_fun=function(x) 1/mean(x),neighb=8,
                 savememory=FALSE) {

  ## in case x is provided as sf: convert
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

    ## Make habitat total area,
    tpa   = sum(pa,na.rm=TRUE)
    mnpa  = mean(pa,na.rm=TRUE)
    sdpa  = sd(pa,na.rm=TRUE)
    ## number of patches,
    tnp   = length(pa)
    ## and nearest neighbour distance
    nnd   = NNDf(mdist)
    ##
    tnnd  = sum(nnd,na.rm=TRUE)
    mnnnd = mean(nnd,na.rm=TRUE)
    sdnnd = sd(nnd,na.rm=TRUE)

    ## return named vector
    landstats = c(tpa,mnpa,sdpa,tnp,tnnd,mnnnd,sdnnd)
    ## names
    names(landstats) =  c("tpa","mnpa","sdpa","tnp","tnnd","mnnnd","sdnnd")
    ##
    return(landstats)

  } else {
    print("--- error --- pa or mdist do not conform to one of the following ---")
    print('class(pa)[1]=="numeric" & class(mdist)[1]=="matrix"& length(pa)==dim(mdist)[1] & dim(mdist)[1]==dim(mdist)[2]')
  } # else
}



#######
#' @name FragStats
#' @title Make a selection of FRAGSTATS landscape statistical metrics
#' @description calculate basic FRAGSTATS metrics, such as patch cohesion using the landscapemetrics [R-package](https://r-spatialecology.github.io/landscapemetrics/)
#' @param x  habitat Raster*Layer file, should be binary if raster: 1= habitat and 0= non-habitat
#' @param metrics character, abbreviation of metrics, cf. landscapemetrics [R-package](https://r-spatialecology.github.io/landscapemetrics/)
#' @param directions numeric, directions argument for directions in which cells are connected (4 or 8)
#' @param count_boundary A species-specific dispersal distance (e.g. average gap-crossing distance, in the units of the mdist values, default = 500)
#' @param edgedepth numeric, distance (measured in nr of cells) a cell has to be away from the patch boundary to be considered as core cell.
#' @return A named vector of landscape statistical metrics: default: number of patches (lsm_c_np),patch density (lsm_c_pd), aggregation index (lsm_c_ai), clumpiness index (lsm_c_clumpy), patch cohesion index (lsm_c_cohesion), normalized landscape shape index (lsm_c_nlsi), effective mesh size (lsm_c_mesh), percentage of like adjacencies (lsm_c_pladj)
#' @export
FragStats = function(x=NULL,metrics=c("lsm_c_np","lsm_c_pd","lsm_c_ai","lsm_c_clumpy","lsm_c_cohesion","lsm_c_nlsi","lsm_c_mesh","lsm_c_pladj"),
                     directions=4,count_boundary=FALSE,edgedepth=2) {


  ##########################################################
  ## check out the landscapemetrics R package! There are many easy summary stats doable with this!
  ## (e.g. number of patches, patch area etc)
  ## list some summary landscape metrics at class level
  #landscapemetrics::list_lsm(level = c("class"), type = "aggregation metric")
  #landscapemetrics::list_lsm(level = c("landscape"), type = "aggregation metric")
  #landscapemetrics::list_lsm(level = c("patch"), type = "aggregation metric")

  ### make sure raster habitat cells are = 1 and non-habitat are NA
  x[x>=1] = 1
  x[x<1]  = NA

  ### calculate some general landscape metricsf patches (compare raster and shapefile version..)
  ### these things are useful when we have a raster but these can differ from polygon circles...
  lmetr   =  landscapemetrics::calculate_lsm(x,
                                             what=metrics,
                                             directions = directions, # make sure directions are the same as for lclump
                                             count_boundary = count_boundary,
                                             edge_depth = edgedepth, neighbourhood = directions) # many options...


    ## return named vector
    landstats = lmetr$value
    ## names
    names(landstats) =  lmetr$metric
    ##
    return(landstats)

} # function end


##########################################################
##### project:       Multispecies Connectivity Modelling
##### author:        Jacqueline Oehri (JO), jacqueline.oehri@gmail.com
##########################################################
### References
#Albert, C. H., Rayfield, B., Dumitru, M., & Gonzalez, A. (2017). Applying network theory to prioritize multispecies habitat networks that are robust to climate and land‐use change. Conservation Biology, 31(6), 1383-1396.
#McGarigal, K., Cushman, S.A., and Ene E. 2012. FRAGSTATS v4: Spatial Pattern Analysis Program for Categorical and Continuous Maps. Computer software program produced by the authors at the University of Massachusetts, Amherst. Available at the following website: https://www.umass.edu/landeco/
#Saura, S., & Pascual-Hortal, L. (2007). A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and urban planning, 83(2-3), 91-103.
#Saura, S., & Rubio, L. (2010). A common currency for the different ways in which patches and links can contribute to habitat availability and connectivity in the landscape. Ecography, 33(3), 523-537.
#Saura, S., Estreguil, C., Mouton, C., & Rodríguez-Freire, M. (2011). Network analysis to assess landscape connectivity trends: application to European forests (1990–2000). Ecological Indicators, 11(2), 407-416.
#Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the world’s ecoregions: How well connected are they?. Ecological indicators, 76, 144-158.
