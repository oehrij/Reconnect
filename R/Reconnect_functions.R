##########################################################
### Rapid Evaluation of Multispecies Connectivity (Reconnect) Workflow
### Reconnect functions
### Functions that represent the core body of the Reconnect framework
### Jacqueline Oehri, jacqueline.oehri@gmail.com
### 18.11.2022
##########################################################

#######
#' @name mkdispfun
#' @title define a dispersal probability function
#' @description create a function that calculates the dispersal probability between two points separated by distance "d", for a hypothetical species with dispersal capacity "alpha"
#' @param option four options are currently available: "negex" (negative exponential, default), "linear" (linear distance decay), "log-sech" (long-tail distance decay) and "quantile", (negative exponential with alpha being median [quantile=0.5] or maximum [quantile=0.1] dispersal distance, default quantile value = 0.5, should be >0 and <=1)
#' @return a dispersal probability function describing the dispersal kernel for given distances (d) and dispersal capacities of species (alpha).
#' @export
mkdispfun = function(option="negex"){
if(option=="negex")         {dispfun= function(d,alpha=500) { exp( -(1/alpha) * abs(d) ) }}                             # alpha = usually average dispersal distance
else if(option=="linear")   {dispfun= function(d,alpha=500) { d = abs(d); ifelse(d < alpha, 1 - d/alpha, 0) }}          # alpha = usually maximum dispersal distance
else if(option=="log-sech") {dispfun= function(d,alpha=500,beta=1.77) { (2*atan((alpha/d)^(1/(1/(beta-1)))))/pi}}       # beta  = thickness of tail in fat-tailed distribution # alpha=usually average dispersal distance
else if(option=="quantile") {dispfun= function(d,alpha=500,q=0.5) { exp( (log(q)/alpha) * abs(d) ) }}                   # alpha = usually median (50% quantile, q = 0.5, default) or maximum (10%, q= 0.1) dispersal distance, assuming that 50% or 10% of dispersal events are beyond the alpha value.
return(dispfun)
  }


#######
#' @name bufferdist
#' @title identify distance of low dispersal probability (<= 0.1 %) for a given dispersal capacity (alpha) and a distance probability function (dispfun)
#' @description function to find the optimal distance (or buffer radius) to consider where dispersal probability <= 0.1 % for a given alpha and dispersal probability function
#' @param alpha numeric vector of average dispersal distance(s) in m
#' @param dispfun function for dispersal probability with a given distance (d) and dispersal probability (alpha)
#' @return a numeric vector of buffer distances to use for given alphas
#' @export
bufferdist = function(alpha=500,dispfun=function(d,alpha){exp(-(1/alpha)*abs(d))},prob=0.001){
                      res= 1; d=1; while(res>=prob) {res = dispfun(d,alpha); d=d+1};
                      return(data.frame(dist=d,prob=res)) }

## from Huang et al. 2020
#dispfun(75000,alpha=317)
#[1] 0.009459466
## with the negative exponential:
#dispfun(75000,alpha=317)
#1.773943e-103

#######
#' @name normalizef
#' @title normalizing function
#' @description normalize values of a vector between 0-1
#' @param x numeric vector of values to be normalized
#' @return a numeric vector of normalized values
#' @export
normalizef = function(x,na.rm=TRUE) {
  if(length(unique(x)[is.finite(unique(x))])>1){
  100*((x-min(x,na.rm=na.rm))/(max(x,na.rm=na.rm)-min(x,na.rm=na.rm)))
  } else if(length(unique(x)[is.finite(unique(x))])<=1){x} }


#######
#' @name readini
#' @title read in the excel Reconnect inifile
#' @description read in the Reconnect inifile
#' @param path character specifying path to the ini file ("C:/user/mypath/myfile.xlsx")
#' @param sheet numeric value of the inifile excel sheet, default = 1
#' @return a dataframe with the specified connectivity functions
#' @export
readini  = function(path=NULL,sheet=1) {
  remc   = readxl::read_excel(path=path,sheet=sheet)
  header = as.character(remc[(which(remc[,1]=="Reconnect ini file start")+1),])
  remc   = data.frame(remc[(which(remc[,1]=="Reconnect ini file start")+2):nrow(remc),])
  colnames(remc) = header
  return(remc)
}

#######
#' @name minarm
#' @title remove small patches
#' @description remove all habitat patches smaller than a minimum area
#' @param shp a shapefile with polygons
#' @param minsize numeric value of minimum patch area in units of shapefile
#' @return a new shapefile where small patches have been removed
#' @export
minarm = function(shp=NULL,minsize=0) {
  ## derive patch areas
  paas = sf::st_area(shp)
if(minsize > 0) {
  aidx = which(as.numeric(paas)>=minsize)
  shp  = shp[aidx,]
}
  return(shp)
}

########
#' @name mkcenters
#' @title make centerpoints of potential windows of selection with given distance apart, within extent "ext" and where "incr" is subtracted
#' @description get an sf object of by incr regularly spaced centerpoints in a given extent (ext)
#' @param ext a raster::extent object
#' @param incr numeric value defining distance among centerpoints
#' @return a shapefile containing a regular grid of centerpoints within extent ext
#' @export
mkcenters = function(ext=NULL,incr=NA) {

  if(is.finite(incr)){
    print("making centerpoints along extent...")
    ## derive centerpoints,
    csx        = seq(from=ext[1]+incr,to=ext[2]-incr,by=incr)
    csy        = seq(from=ext[3]+incr,to=ext[4]-incr,by=incr)
  } else if(!is.finite(incr)) {
    print("make only one centerpoint...")
    ## make just 1 centerpoint for incr = NA
    csx        = ext[1]+((ext[2]-ext[1])/2)
    csy        = ext[3]+((ext[4]-ext[3])/2)
  }
  coords     = expand.grid(x=csx,y=csy)
  coords$cid = 1:nrow(coords)
  ## make an sf object (has the advantage of already having a crs)
  #centers = sf::st_as_sf(coords,coords=c("x", "y"),crs=st_crs(ext))
  centers = sf::st_as_sf(coords,coords=c("x", "y"))
  return(centers)
}


########
#' @name mkbuffer
#' @title make a square or circular buffer for a given coordinate
#' @description make a buffer sf object around a given sf spatial points center point
#' @param centers a sf object (including crs, a projection with unit = meters) with coordinates for center points of buffer locations
#' @param radius the radius of the circular buffer, or half the side length of the square buffer in meters
#' @param form the form of the buffer: currently only circular or square
#' @return a sf object with corresponding buffers
#' @export
mkbuffer = function(centers,radius,form=c("square")) {

  crs =st_crs(centers)
  plist = list()

  if(form=="circle"){
    for(ii in 1:nrow(centers)) {
      #print(ii)
      cords  = sf::st_coordinates(centers$geometry[ii])
      circ   = dismo::circles(cords,d=radius,n=360,lonlat=FALSE) # if I increase the n (default 360, the circles are more precise later..)
      circ   = circ@polygons@polygons
      circ   = circ[[1]]@Polygons[[1]]
      coords = slot(circ,"coords")
      polyA  = list(x=coords[,1],y=coords[,2])
      p      = st_polygon(list(as.matrix(cbind(polyA$x,polyA$y))))
      plist[[ii]] = p
    } # centsp
    # circle
  } else if(form=="square"){
    for(ii in 1:nrow(centers)) {
      cords  = sf::st_coordinates(centers$geometry[ii])
      polyA  = list(x=c(cords[1]-(radius),cords[1]+(radius),cords[1]+(radius),cords[1]-(radius)),
                    y=c(cords[2]-(radius),cords[2]-(radius),cords[2]+(radius),cords[2]+(radius)))
      ## add the closing loop entry
      polyA$x[length(polyA$x)+1] = polyA$x[1]
      polyA$y[length(polyA$y)+1] = polyA$y[1]
      p      = st_polygon(list(as.matrix(cbind(polyA$x,polyA$y))))
      plist[[ii]] = p
    } # centsp
    # square
  }
  sps = st_as_sfc(plist,crs=crs)
  sps = st_sf(ID = c(1:length(sps)), geometry = sps)
  return(sps)
}


########
#' @name alignst
#' @title align raster "orig" with raster "target"
#' @description change the coordinate reference system, extent and resolution of "orig" to match the one of "target". The function relies on the 'sf' and 'terra'packages
#' @param orig a raster layer that should be modified
#' @param target a raster layer to which the orig should be matched
#' @param method a character value ("near","bilinear","cubic","cubicspline") for the method used for estimating new cell values in the terra::project function
#' @return a raster layer that matches to the "target" layer in coordinate reference system, extent and resolution.
#' @export
alignst = function(orig=NULL,target=NULL,method="near") {
  ## extract names
  orign = names(orig)
  ## make terra objects out of them (faster)
  orig      = terra::rast(orig)
  target    = terra::rast(target)
  ## align crs if necessary, i.e. reproject if necessary
  if(!(terra::crs(orig) == terra::crs(target))) {
    ## reproject everything to the target raster, use default settings of bilinear and near method, see helpfile
    orig  = terra::project(x=orig,y=target,align=TRUE,mask=TRUE,method=method)
  }
  # ## if rasters still not align
  # if(terra::ext(orig)!=terra::ext(target)) {
  ## resample because rasters have a different extent!
  orig = terra::resample(x=orig, y=target, method=method)
  #}
  orig  = as(orig, "Raster")
  ## add name just in case it got lost on the way
  names(orig) = orign
  return(orig)
}


########
#' @name Reconnect_core
#' @title Function within Reconnect_wrap: apply connectivity functions in the "Reconnect_inifile.xlsx" for a given landscape
#' @description this function is applied in the Reconnect_wrap function, to calculate connectivity indices for a buffer window around a given sf spatial points center point
#' @param cid character of the buffer-center shapefile id-column
#' @param pidname character describing the patch-id column name in the habitat patch shapefile (hshp)
#' @param hr  raster layer or raster stack object, containing a habitat layer and a resistance layer.
#' @param foc shapefile, of focal region of interest
#' @param hshp shapefile of habitat patches
#' @param minsizev numeric vector of minimum habitat patch areas in m2
#' @param cutpatches logical, should habitat patches be cut at the boundary of the focal region of interest (foc)? Default = TRUE
#' @param cost_type character, algorithm applied to calculate distances: either "euclidean" (default), "least-cost" or "commute-time"
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @param trans_fun function expression, argument for the gdistance::transition() function, applied when moving between resi raster cells. Default for conductance: function(x) 1/mean(x)
#' @param neighb numeric, directions argument for the gdistance::transition() function: directions in which cells are connected (4, 8, 16, ...)
#' @param areafun function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @param cond dataframe, connectivity functions ini file (cf. 1-Data folder & readini function)
#' @param nplot numeric value defining the frequency with which plotting is done for the buffer id's (default: 10, i.e. every 10th id will be plotted)
#' @param cols color ramp palette to be used (default: terrain.colors(100))
#' @param maxnpatch numeric, maximum number of patches to be treated in a landscape (so that connectivity functions don't get overwhelmed)
#' @param radius numeric, radius of buffer window
#' @param resolution numeric, resolution of raster grid
#' @param incr numeric, incremental distance among buffer windows
#' @param respath character string indicating file path where results should be saved
#' @return csv files saved to respath location, where results from the application of connectivity indicators specified in the cond file are saved (dfla = landscape-level, dfpa = patch-level, dfpi = pixel-level).
#' @export
Reconnect_core  = function(cid      = NULL,
                      pidname  = NULL,
                      hr       = NULL,
                      foc      = NULL,
                      hshp     = NULL,
                      minsizev = NULL,
                      cutpatches = TRUE,
                      cost_type  = "euclidean",
                      dist_type  = "edge",
                      trans_fun  = function(x) 1/mean(x),
                      neighb     = 8,
                      areafun  = sf::st_area,
                      cond     = NULL,
                      nplot    = 10,
                      cols     = terrain.colors(100),
                      maxnpatch = 28500,
                      radius   = NULL,
                      resolution = NULL,
                      incr       = NULL,
                      respath    = NULL,
                      habi = NULL
                      ) {

  ########################################
  ## Create a dataframe where patch IDs and corresponding connectivity measurements are stored!!
  dfpi  = data.frame(pid=NA,minsize=NA,alpha=NA,metric=NA,value=NA)
  dfpa  = data.frame(pid=NA,minsize=NA,alpha=NA,metric=NA,value=NA)
  dfla  = data.frame(pid=NA,minsize=NA,alpha=NA,metric=NA,value=NA)

  ########################################
  ## Generate shapefile for region of interest, include patches that are touched by window of interest (for correct minarea removal!)
  hshp1 = hshp[foc,]

  ########################################
  ## go through minsize argument
  for (minsiz in sort(minsizev,decreasing=FALSE)){
    ## minsiz = sort(minsizev,decreasing=FALSE)[1]
    print(sprintf("apply minsize %0.02fm...",minsiz))
    shp1 = minarm(shp=hshp1,minsize=minsiz)
    ## apply cutpatches argument only here, after minimum patches have been removed..!
    if(cutpatches==TRUE) {shp1 = sf::st_intersection(shp1,foc$geometry)}

    ### sometimes there is no habitat in focal area
    if(length(shp1[[pidname]])>0) {
      ### create new foc1 as polygon
      foc1  = sf::st_union(foc,shp1) %>% sf::st_union(by_feature=FALSE)
      foc1  = st_sf(geometry=foc1)
    } else if(length(shp1[[pidname]])==0) {
      foc1  = st_sf(geometry=foc$geometry)
    }

    ## try if extents overlap, if not go to next cid...
    try = try(raster::crop(hr,raster::extent(foc1)),silent=TRUE)
    if(class(try)[1]=="try-error") {print(sprintf("no overlap between cid %d and habitat data...",cid)); next}

    ### crop raster according to extent of foc (or should it be according to extent of cutpatches, which might give assymetrical landscapes??)
    hr1   = raster::crop(hr,raster::extent(foc1))
    hr1   = raster::mask(hr1,foc1)

    ### determine hr and resistance rasters
    lcr2 = hr1[["habitat"]]
    if(nlayers(hr1)==2) {re1 = hr1[["resistance"]]} else if(nlayers(hr1)<2) {re1 = NULL}

    ### set all habitat patches in raster (determined by pixel values) not contained in shp1 to NA!
    if(length(shp1[[pidname]])>1) {
      lcr2[!values(lcr2)%in%shp1[[pidname]]] = NA
    } else if(length(shp1[[pidname]])==1) {
      lcr2[values(lcr2)!=shp1[[pidname]]] = NA
    } else if(length(shp1[[pidname]])==0) {
      lcr2[] = NA
    }

    ### make 1 habitat raster for connfunctions, where habitat cells = 1 and the rest = NA
    lcr1 = lcr2
    lcr1[!is.na(lcr1)] = 1

    ### derive resolution from input data
    resolution = mean(res(hr1),na.rm=TRUE)

    ### short plot
    if ((cid/nplot)-floor(cid/nplot)==0) {
      plot(foc$geometry,main=sprintf("cid %s focal area \n habitat and resistance \n minsize %0.02fm %0.0fm res",cid,minsiz,resolution))
      if(length(re1)>0){raster::plot(re1,add=TRUE)}
      raster::plot(lcr2,add=TRUE,col="darkgreen",legend=FALSE)
      plot(foc1$geometry,add=TRUE,border="darkorange",lwd=2)
    }

    ########################################
    ## 2) Start Reconnect calculations for area of interest, foc1a, pa & mdist are always calculated
    ########################################
    ## determine total area of foc1 in any case dfla$foc_ha = foc1a/10000
    foc1a = as.numeric(areafun(foc1))
    dfl   = data.frame(pid="all",minsize=minsiz,alpha=NA,metric="foc_ha",value=foc1a/10000)
    dfla  = rbind(dfla,dfl)
    ## determine number of patches in landscape
    dfl   = data.frame(pid="all",minsize=minsiz,alpha=NA,metric="npatch",value=length(shp1[[pidname]]))
    dfla  = rbind(dfla,dfl)
    ## determine total and mean habitat patch area in landscape
    pa_ha = as.numeric(sf::st_area(shp1))
    dfl   = data.frame(pid="all",minsize=minsiz,alpha=NA,metric=c("tpa_ha","mnpa_ha"),value=c(sum(pa_ha,na.rm=TRUE)/10000,mean(pa_ha,na.rm=TRUE)/10000))
    dfla  = rbind(dfla,dfl)

    ########################################
    ## derive patch areas in case there are some habitat patches...
    if(length(shp1[[pidname]])>0 & length(shp1[[pidname]])<=maxnpatch){

      ########################################
      ## determine patch id names
      pid = shp1[[pidname]]

      ########################################
      ## distance matrix: is quite slow, only make the distance matrix once
      ## make the euclidean distance matrix
      if(cost_type == "euclidean") {
        ## get the patch area and euclidean distance matrix object
        pamd  = eucdist_pm(x=shp1,id=pidname,areafun=sf::st_area,dist_type=dist_type)
        pa1    = pamd$pa
        mdist1 = pamd$mdist
        ## remove to save memory
        rm(list="pamd")
      } # euclidean

      ## make the resitance distance matrix
        if(length(re1)>0 & cost_type == "least-cost") {
            ## get the patch area and cost distance matrix object
            pamd = costdist_pm(x=shp1,resi=re1,id=pidname,areafun=areafun,
                               trans_fun=trans_fun,neighb=neighb,
                               dist_type=dist_type,cost_type=cost_type)
            pa1    = pamd$pa
            mdist1 = pamd$mdist
            ## remove to save memory
            rm(list="pamd")
        } # least-cost

      ## add area of patches  (is quite fast..)

      ########################################
      ## Always add the patch area (important for summary after)
      dfp  = expand.grid(pid=pid,minsize=minsiz)
      dfp[["alpha"]]  = NA
      dfp[["metric"]] = "pa_ha"
      dfp[["value"]]  = pa1/10000
      dfpa = rbind(dfpa,dfp)

    ########################################
    ## Apply connectivity functions for specified alphas
    ########################################

    ## go through connectivity metrics
    for(conn in unique(cond$conname)) {
      ## conn    = unique(cond$conname)[1]
      condd     = cond[cond$conname==conn,]
      print(condd)
      resname   = gsub(" ","",strsplit(condd[,"resname"],split=",")[[1]])
      confun    = eval(parse(text=condd[,"confun"]))
      statement = condd[,"statement"]
      alphav    = strsplit(condd[,"alphav"],split=",")[[1]]
      alphav    = alphav[alphav!="NA"] %>% as.numeric()
      minssz    = strsplit(condd[,"minsizev"],split=",")[[1]]
      minssz    = minssz[minssz!="NA"] %>% as.numeric()
      level     = unique(strsplit(condd[,"level"],split=",")[[1]])
      resultnr  = strsplit(condd[,"resultnr"],split=",")[[1]]
      resultnr  = resultnr[resultnr!="NA"]%>% as.numeric()
      aggfun    = eval(parse(text=condd[,"aggfun"]))

      ## check if alpha's are used
      alphav = alphav[is.finite(alphav)]
      ## check if minsize is specified
      minssz = minssz[is.finite(minssz)]
      if(length(minssz)>0) {
      ## assign corresponding alphav if minsize is given for a specific alpha!
      alphav = alphav[which(minssz == minsiz)]
      }
      ## apply functions
      if(length(alphav)>0) {
        alphav = alphav[order(alphav,decreasing=FALSE)]
        ## try function first
        try = try(sapply(alphav,function(alpha1){eval(parse(text = paste(statement)))}),silent=FALSE)
        if(class(try)[1]=="try-error"){ conres = -99
        } else if(class(try)[1]!="try-error"){
          conres = sapply(alphav,function(alpha1){eval(parse(text = paste(statement)))})
          colnames(conres) = alphav
          if(length(conres)==0)  {conres = -99 }
          if(length(resultnr)>0) {conres = conres[resultnr,]}
        }
        if("pixel" %in% level){
          if(is.list(conres)){pp = unlist(conres)} else {pp = as.vector(conres)}
          dfpp = expand.grid(pid=1:(length(pp)/length(alphav)),minsize=minsiz,metric=resname,alpha=alphav)
          dfpp[["value"]]  = pp
          dfpp = dfpp[,c("pid","minsize","alpha","metric","value")]
          dfpi = rbind(dfpi,dfpp)
        }
        if("patch" %in% level){
          dfp = expand.grid(pid=pid,minsize=minsiz,metric=resname,alpha=alphav)
          dfp[["value"]]  = unlist(conres)
          dfp  = dfp[,c("pid","minsize","alpha","metric","value")]
          dfpa = rbind(dfpa,dfp)
        }
        if("landscape"%in%level){
          dfl            = expand.grid(pid="all",minsize=minsiz,metric=resname,alpha=alphav,value=NA)
          dfl[["value"]] = unlist(lapply(conres,aggfun))
          dfl = dfl[,colnames(dfla)]
          dfla = rbind(dfla,dfl)
        }
      } else if(length(alphav)==0){
        ## try function first

        try = try(eval(parse(text = paste(statement))),silent=FALSE)
        if(class(try)[1]=="try-error"){ conres = -99
        } else if(class(try)[1]!="try-error"){
          conres = eval(parse(text = paste(statement)))
          if(length(conres)==0)  {conres = -99 }
          if(length(resultnr)>0){conres = conres[[resultnr]]}
        }
        if("pixel" %in% level){
            pp = unlist(conres)
          dfpp = expand.grid(pid=1:length(pp),minsize=minsiz)
          dfpp[["alpha"]]  = NA
          dfpp[["metric"]] = resname
          dfpp[["value"]]  = pp
          dfpi = rbind(dfpi,dfpp)
        }
        if("patch" %in% level){
          dfp = expand.grid(pid=pid,minsize=minsiz)
          dfp[["alpha"]]  = NA
          dfp[["metric"]] = resname
          dfp[["value"]]  = unlist(conres)
          dfpa = rbind(dfpa,dfp)
        }
        if("landscape"%in%level){
          dfl  = expand.grid(pid="all",minsize=minsiz,metric=resname,alpha=NA,value=NA)
          if(length(conres)>nrow(dfl)){conres = aggfun(unlist(conres))} else if(length(conres)==0) { conres = -99 }
          dfl[["value"]] = conres
          dfl = dfl[,colnames(dfla)]
          dfla = rbind(dfla,dfl)
        }
      }
    } # for all conn functions


    } # apply confun only if nr of patches are > 0 and < maxnrpatch

  } # minsize

  ########################################
  ## finish datasets for window of interest, write to file..
  dfpi$level = "pixel"
  dfpa$level = "patch"
  dfla$level = "landscape"

  ## write to file if at least 1 row!
  if(nrow(dfpi)>1){
    dfpi = dfpi[-1,]
    dfpi$cid    = cid
    dfpi$radius = radius
    dfpi$resolution  = resolution
    dfpi$incr   = incr
    write.csv(dfpi,paste(respath,sprintf("dfpi_rad%0.0fm_res_%0.0fm_incr_%0.0fm_cid%d_%s.csv",radius,resolution,incr,cid,habi),sep="/"),row.names=FALSE)
  }

  ## write to file if at least 1 row!
  if(nrow(dfpa)>1){
    dfpa = dfpa[-1,]
    dfpa$cid    = cid
    dfpa$radius = radius
    dfpa$resolution  = resolution
    dfpa$incr   = incr
    write.csv(dfpa,paste(respath,sprintf("dfpa_rad%0.0fm_res_%0.0fm_incr_%0.0fm_cid%d_%s.csv",radius,resolution,incr,cid,habi),sep="/"),row.names=FALSE)
  }

  ## write to file if at least 1 row!
  if(nrow(dfla)>1){
    dfla = dfla[-1,]
    dfla$cid    = cid
    dfla$radius = radius
    dfla$resolution  = resolution
    dfla$incr   = incr
    write.csv(dfla,paste(respath,sprintf("dfla_rad%0.0fm_res_%0.0fm_incr_%0.0fm_cid%d_%s.csv",radius,resolution,incr,cid,habi),sep="/"),row.names=FALSE)
  }


  }# function end


########
#' @name Reconnect_wrap
#' @title wrapper function for the Reconnect_core function
#' @description subsets your habitat raster file into bufferzones, prepares data such that connectivity indices can be calculated in a moving buffer using the Reconnect_core function and connectivity functions specified in the Reconnect_ini file.
#' @param lcr raster layer object, containing a habitat patches (binary, habitat = 1, other = NA)
#' @param lres raster layer object, containing a resistance surface
#' @param roi shapefile, region of interest
#' @param cond dataframe, connectivity functions ini file (cf. 1-Data folder & readini function)
#' @param respath character string indicating file path where results should be saved
#' @param resolution numeric, resolution of raster grid
#' @param radius numeric, radius of circular moving window (if moving window is square, radius = 0.5*window sidelength
#' @param incr, numeric, incremental distance between moving window centers
#' @param habi character, habitat name
#' @param bufform character, form of buffer window used by the mkbuffer function (default: "square", alternative "circle")
#' @param minsizev numeric vector of minimum habitat patch areas in m2
#' @param cutpatches logical, should habitat patches be cut at the boundary of the focal region of interest (foc)? Default = TRUE
#' @param ncores numeric, number of cores that should be used in the parallel processing of moving buffers
#' @param oviewplot logical, should an overview plot be generated or not? (default: TRUE)
#' @param cost_type character, algorithm applied to calculate distances: either "euclidean" (default), "least-cost" or "commute-time"
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @param trans_fun function expression, argument for the gdistance::transition() function, applied when moving between resi raster cells. Default for conductance: function(x) 1/mean(x)
#' @param neighb numeric, directions argument for the gdistance::transition() function: directions in which cells are connected (4, 8, 16, ...)
#' @param areafun function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @param maxnpatch numeric, maximum number of patches to be treated in a landscape (so that connectivity functions don't get overwhelmed)
#' @param cst numeric, start of moving window indices that should be computed (if NULL as in default, cst will be set to 1)
#' @param cen numeric, end of moving window indices that should be computed (if NULL as in default, cen will be set to the maximum nr of moving windows existent)
#' @param hcols character string of colors used in habitat layer plotting
#' @param rcols character string of colors used in resistance layer plotting
#' @return csv files saved to respath location, where results from the application of connectivity indicators specified in the cond file are saved (dfla = landscape-level, dfpa = patch-level, dfpi = pixel-level).
#' @export
Reconnect_wrap = function(lcr=NULL,lres=NULL,roi=NULL,cond=NULL,respath=NULL,resolution=NULL,radius=NULL,incr=NA,
                     habi="test",bufform="square",minsizev=c(0),cutpatches=TRUE,
                     ncores=1,oviewplot=TRUE,
                     cost_type  = "euclidean",
                     dist_type  = "edge",
                     trans_fun  = function(x) 1/mean(x),
                     neighb     = 8,
                     areafun    = sf::st_area,
                     maxnpatch  = 28500,
                     cst  = NULL,
                     cen  = NULL,
                     hcols = rev(c("yellow","orange","red","purple","blue","lightblue")),
                     rcols = rev(c("purple","lightblue","yellow","olivedrab","darkgreen"))
                     ) {

### set warnings off..
options(warn=-1)

######
### 1.0) make a color palette of the given colour schemes!
hcols = colorRampPalette(hcols)
rcols = colorRampPalette(rcols)

  ######
  ### If resolution and radius are not indicated, make a default
  if(length(resolution)==0) {resolution = mean(res(lcr),na.rm=TRUE)[1]}
  if(length(radius)==0)     {radius = max((extent(lcr)[2]-extent(lcr)[1]),(extent(lcr)[4]-extent(lcr)[3]))/2}
  print(sprintf("resolution=%0.0fm, radius=%0.0fm, increment=%0.0fm...",resolution,radius,incr))

  ######
  ### assign a new hr object
  hr = lcr
  ### make sure that non-habitat cells are NA and not 0
  hr[hr<1] = NA
  ### name habitat raster
  names(hr) = "habitat"

  ### 2.1) resample habitat raster in case the desired resolution does not match the one of the raster at hand
  if(resolution != mean(res(lcr),na.rm=TRUE)[1]) {
    ### determine factor of aggregation
    fact = resolution/mean(res(lcr),na.rm=TRUE)[1]
    hr[is.na(hr)] = 0
    hr = raster::aggregate(x=hr,fact=fact,fun=function(x,na.rm=na.rm){sum(x)},expand=FALSE,na.rm=TRUE)
    ### Conservative approach: set all the pixels with less than 75% habitat to 0!
    hr[hr<=0.75*(fact^2)] = 0
    hr[hr>0.75*(fact^2)]  = 1
    hr[hr==0] = NA
  }

  ######
  ### 3.0) add resistance to the habitat raster if provided (and otherwise create one of resitance 1...)
  if(length(lres)>0){
    re = lres
    ## aggregate to larger resolution if necessary (is done automatically in alignst..)
    re  = alignst(orig=re,target=hr,method="bilinear")
  } else if(length(lres)==0) {
    ## create a resistance layer where all resistance is set to 1 (more compatible with overall functionality)
    re   = hr
    re[] = 1
  }
    ## make sure it is named correctly
    names(re) = "resistance"
    ## make sure that infinite cells are just NA
    re[!is.finite(re)] = NA
    ## stack habitat and resistance raster
    hr = raster::stack(hr,re)

  ######
  ### 4.0) mask to region of interest if roi is provided, and set extent (used for centers below)
  #####
  ### make the roi with a buffer distance!! (like this, the global extent will anyway be kept when incr = NA)
    if(length(roi)>0){
      if(sf::st_crs(hr)!=sf::st_crs(roi)) {roi = sf::st_transform(roi,crs = st_crs(hr))}
      } else if(length(roi)==0){
    roi = sf::st_sf(ID= 1, geometry=sf::st_as_sfc(sf::st_bbox(hr)))
    sf::st_crs(roi) = sf::st_crs(hr)
    }
    ## make the roi1 with a buffer distance: JO: make bufferdist with increment instead of radius!
    #if(!is.na(incr)){roi1 = sf::st_buffer(roi,dist=incr)}else if(is.na(incr)){roi1 = sf::st_buffer(roi,dist=resolution)}
    roi1 = sf::st_buffer(roi,dist=radius)
    ext  = raster::extent(roi1)
    hr   = raster::crop(hr,ext)
    hr   = raster::mask(hr,roi1)

  ######
  ### 5.0) identify raster clumps, i.e. generate IDs!
  hr1   = landscapemetrics::get_patches(hr[["habitat"]], directions=4, class=1)[[1]][[1]]
  print(sprintf("found %0.0f patches..",max(values(hr1),na.rm=TRUE)))

  ######
  ### 5.1) if it is a raster, make a shapefile out of it, to be able to exclude patches below minsize..
  hshp = stars::st_as_stars(hr1) %>% sf::st_as_sf(merge = TRUE)
  hshp = sf::st_make_valid(hshp,reason=TRUE)
  names(hshp)[[1]] = "pid"

  ######
  ### 5.2) replace habitat raster by patches raster
  hr[["habitat"]] = hr1
  rm(list="hr1")

  ######
  ### 5.3) write the patches to shapefile
  if(file.exists(paste(respath,sprintf("hshp_res%0.0f_rad%0.0f_incr%0.0f_habi%s.shp",resolution,radius,incr,habi),sep='/'))){
     file.remove(paste(respath,sprintf("hshp_res%0.0f_rad%0.0f_incr%0.0f_habi%s.shp",resolution,radius,incr,habi),sep='/')) }
  st_write(hshp, paste(respath,sprintf("hshp_res%0.0f_rad%0.0f_incr%0.0f_habi%s.shp",resolution,radius,incr,habi),sep='/'), driver="ESRI Shapefile")  # create to a shapefile

  ######
  ### 6.0) set centerpoints of buffers (& make & show buffers!)
  centers = mkcenters(ext=ext,incr=incr)

  #### if a region of interest is given, only cut out centers for region of interest
  ### 6.1) mask to region of interest if roi is provided, and set extent (used for centers below)
  if(sf::st_crs(centers)!=sf::st_crs(roi1)) {
      if(!is.na(sf::st_crs(centers))) {
          centers = sf::st_transform(centers,crs = sf::st_crs(roi1))
          } else if(is.na(sf::st_crs(centers))) {
            sf::st_crs(centers) =  sf::st_crs(roi1)
           }
         }
  centers = sf::st_intersection(centers,roi1)

  ### write centers shapefile to file, it contains the IDs and corresponding moving window center locations!
  if(file.exists(paste(respath,sprintf("centers_res%0.0f_rad%0.0f_incr%0.0f_habi%s.shp",resolution,radius,incr,habi),
                       sep='/'))) { file.remove(paste(respath,sprintf("centers_res%0.0f_rad%0.0f_incr%0.0f_habi%s.shp",resolution,radius,incr,habi),sep='/'))}
  sf::st_write(centers, paste(respath,sprintf("centers_res%0.0f_rad%0.0f_incr%0.0f_habi%s.shp",resolution,radius,incr,habi),sep='/'), driver="ESRI Shapefile")  # create to a shapefile

  ####################################
  ## quick overview plot
  if(oviewplot==TRUE) {
  pdf(paste(respath,sprintf("Reconnect_ov_rad%0.0fm_res_%0.0fm_incr%0.0fm_%s.pdf",radius,resolution,incr,habi),sep="/"))
  ## plot all data
  raster::plot(hr,col=rev(rcols(100)))
  ## take a small extent
  ext2 = ext*0.2
  #plot(ext2)
  if(nlayers(hr)==2){
  raster::plot(hr[["resistance"]],col=adjustcolor(rev(rcols(100)),alpha.f=0.5),legend=FALSE,main=sprintf("nr of buffers = %d",length(centers$geometry)),sub="yellow = roi")
  raster::plot(hr[["habitat"]],col=hcols(100),add=TRUE,legend=FALSE)
  } else if(nlayers(hr)==1) {raster::plot(hr[["habitat"]],col=hcols(100),legend=FALSE,sub="yellow = roi")}
  plot(centers$geometry,add=TRUE,pch=3,cex=0.4,lwd=1.2,col="black")
  plot(roi$geometry,add=TRUE,border="yellow",col=NA,lwd=2)
  ## make a polygon of extent and extract center coordinates
  p         = st_as_sf(as(ext2, 'SpatialPolygons'))
  st_crs(p) = st_crs(centers)
  centers1  = centers[p,]
  ### make circular buffers
  sps = mkbuffer(centers=centers1,radius=radius,form=bufform)
  ## controlplot
  raster::plot(hr[["habitat"]],col="darkgreen",legend=FALSE,
               main=sprintf("nr of buffers = %d",length(centers$geometry)),sub="yellow = roi")
  plot(sps$geometry,add=TRUE,col=adjustcolor(hcols(100),alpha.f=0.3),border=hcols(10),lwd=0.8)
  plot(roi$geometry,add=TRUE,border="yellow",col=NA,lwd=2)

  dev.off()
  } # overviewplot

  ######
  ### 7.0) start calculation of MULTIPLEX connectivity for connectivity indices x,y,z
  print(sprintf("start parallel moving windows %s...",habi))
  a = Sys.time()
  ## add packages to the environment in case not contained
  rpacks = unique(unlist(strsplit(cond$rpack,split=",")))
  rpacks = rpacks[rpacks!="NA"]
  ## add gdistance package to the packages if cost type != euclidean!
  if(cost_type!="euclidean"){rpacks = unique(c(rpacks,"gdistance"))}
  ## for parallelized treatment of patch IDs
  ## only registerdoparallel if cluster is not registered
  if(!foreach::getDoParRegistered()){
  doParallel::registerDoParallel(cores=ncores)
  }
  ## define nr of plots here
  #if(resolution<=60){nplot = 500} else if(resolution>60&resolution<=120){nplot=50} else {nplot=1}
  nplot=ceiling(0.1*length(unique(centers$cid)))
  ## define the centerids for which calculations should be done
  if(length(cst)==0) {cst=1}
  if(length(cen)==0) {cen=length(unique(centers$cid))}
  #apply function for each patch ID
  foreach::foreach(cid=unique(centers$cid)[cst:cen],
                    .packages= unique(c("sf","raster","stars","terra","dismo","landscapemetrics",rpacks)),
                    .export  = ls(globalenv())
                    #.combine = rbind
   ) %dopar% {
#for(cid in unique(centers$cid)) {
    # cid = centers$cid[5000]
    print(paste("center cid =",cid,sep=""))

    ########################################
    ## 7.1) Determine focal area of interest for connectivity meausurements
    ## make square buffers ## for testing: radius = radius/5
    foc = mkbuffer(centers=centers[centers[["cid"]]==cid,],radius=radius,form=bufform)

    ### extract raster information for buffer region of interest
    if(sf::st_crs(hr)!=sf::st_crs(foc)) {
      foc = sf::st_transform(foc,crs = sf::st_crs(hr))
    }

    # make plots only for some fraction of landscapes treated
    if ((cid/nplot)-floor(cid/nplot)==0) {
      pdf(paste(respath,sprintf("Reconnect_detail_rad%0.0fm_res_%0.0fm_cid%d_%s.pdf",radius,resolution,cid,habi),sep="/"))
      #sink(paste(respath,sprintf("LOG_Reconnect_detail_rad%0.0fm_res_%0.0fm_cid%d_%s.txt",radius,resolution,cid,habi),sep="/"))
    }

    ########################################
    ## Start of Reconnect_core function
    ########################################
    pidname = names(hshp)[[1]]
    cols    = hcols(100)

    Reconnect_core(cid=cid,pidname=pidname,hr=hr,foc=foc,hshp=hshp,minsizev=minsizev,
              cutpatches=cutpatches,cond=cond,nplot=nplot,cols=cols,
              cost_type  = cost_type,
              dist_type  = dist_type,
              trans_fun  = trans_fun,
              neighb     = neighb,
              areafun    = areafun,
              maxnpatch  = maxnpatch,
              radius     = radius,
              resolution = resolution,
              incr       = incr,
              respath    = respath,
              habi       = habi)

    print("Reconnect core done...")
    ## plot finish closeup
    if ((cid/nplot)-floor(cid/nplot)==0) {dev.off()
                                          #sink()
                                          } # plot

  } # centers dopar

  print(sprintf("all cid's finished..!"))
  print(Sys.time()-a)

} # Reconnect_wrap end


########
#' @name Reconnect_summary
#' @title Summarize moving buffer results from the Reconnect_wrap function into manageable formats
#' @description synthesize all moving buffer windows into 1 file for the grand region of interest (landscape & pixel level connectivity will be summarized as raster files, patch-level connectivity will be summarized as a shapefile)
#' @param sourcepath character string indicating file path where Reconnect_wrap outputs are saved
#' @param respath character string indicating file path where results should be saved
#' @param habnames character string, indicating the habitat names of interest (same as used in Reconnect_wrap function)
#' @param mosafun R-function used for generating a raster mosaic for the pixel-level connectivity functions
#' @param lcr0 raster layer, original large raster layer that was used as input in the Reconnect_wrap function
#' @param resolution numeric, resolution of raster grid
#' @param radius numeric, radius of circular moving window (if moving window is square, radius = 0.5*window sidelength
#' @param incr, numeric, incremental distance between moving window centers
#' @param cutpatches logical, should habitat patches be cut at the boundary of the focal region of interest (foc)? Default = TRUE
#' @param cols  color pallette to be used when plotting (default = topo.colors)
#' @param varset character string of landscape-level connectivity variables that should be set to another value than 0 (default) in case of NA (e.g. variables where connectivity is highest at low values, such as nearest neighbour distance NND)
#' @param valset character or numeric string of values (current options if character: "mean","min","max","radius") to which landscape-level connectivity variables in varset should be set (order of valset must match varset).
#' @return raster files saved to respath (landscape & pixel -level connectivity), shapefiles saved to respath (patch-level connectivity), and plots depicting connectivity at pixel-, patch- and landscape-level.
#' @export
Reconnect_summary = function(sourcepath = NULL,respath=NULL,habnames="test",resolution=NULL,radius=NULL,incr=NA,
                        mosafun=base::min,lcr0=NULL,cutpatches=TRUE,cols=topo.colors,varset=c("NND"),valset=c("max")) {

  ## set warnings off..
  options(warn=-1)
  #######
  ## 1) go through habnames
  for(habi in habnames) {
  #print(habi)
    ######
    ### If resolution and radius are not indicated, make a default
    if(length(resolution)==0) {resolution = mean(res(lcr),na.rm=TRUE)[1]}
    if(length(radius)==0)     {radius = max((extent(lcr)[2]-extent(lcr)[1]),(extent(lcr)[4]-extent(lcr)[3]))/2}
    print(sprintf("resolution=%0.0fm, radius=%0.0fm, increment=%0.0fm...",resolution,radius,incr))

  ########################
  ############ landscape level
  ########################
  ### 2.1) Load data outputs for scale, resolution, radius and habitat of interest
  dfla = load_data(path=sourcepath,
                    pattern=sprintf("dfla_rad%s.+_res_%s.+_incr_%s.+_cid.+_%s.csv",radius,resolution,incr,habi))

  if(length(dfla)>0) {
  print("\n found results at landscape-level... \n")
  ### JO: set all errors in Reconnect_core to 0 because they are summed afterwards, and no na.rm argument can be set.
  #dfla[dfla == -99] = NA
  ## JO 2.3.2023: it seems fasterize automatically excludes NA values!
  dfla[dfla$value == -99 & is.finite(dfla$value),"value"] = NA
  ## extract metrics to later make corresponding stacks
  metrics = unique(dfla$metric)
  ############
  ############
  ### 2.2) Landscape level aggregation (and pixel level, later): just average across moving buffer windows...
  centers   = sf::st_read(paste(sourcepath,sprintf("centers_res%0.0f_rad%0.0f_incr%0.0f_habi%s.shp",resolution,radius,incr,habi),sep="/"))
  polys     = mkbuffer(centers=centers,radius=radius,form=c("square"))
  polys$cid = centers$cid
  dfla      = dfla[dfla$level=="landscape",]
  ### create a wide landscape dataset
  dfla1 = reshape2::dcast(dfla[,c("cid","minsize","alpha","metric","value")], cid~metric+alpha+minsize)

  ## JO: 2.8.2023: in order to create smooth maps, remove NAs! assign a value of 0 for all connectivity metrics except the ones in varset (such as nearest neighbour that needs to be set to a maximum or similar).
  #nnd = grep("^NND_",colnames(dfla1))
  for (vsvr in varset) {
   ## default replaceval is NA
   replaceval = NA
   ## desired value for variable
   vsvl   = valset[which(varset==vsvr)]
   ## the actual variable names in dataset (potentially several with gap crossing distance and minsize argument)
   nnd    = colnames(dfla1)[grep(paste(paste0("^",vsvr,"_"),collapse="|"),colnames(dfla1))]
  if(length(nnd)>0) {
   for(nd in nnd) {
    if(is.numeric(vsvl)) {
       replaceval = vsvl
    } else if(is.character(vsvl)) {
        if (vsvl == "max") {
          replaceval = max(dfla1[,nd],na.rm=TRUE)
        } else if(vsvl == "min") {
          replaceval = min(dfla1[,nd],na.rm=TRUE)
          } else if(vsvl == "mean") {
          replaceval = mean(dfla1[,nd],na.rm=TRUE)
          }
      } # character
  dfla1[is.na(dfla1[,nd]),nd] =  replaceval
  } # nd
  } # nnd > 0
  } # varset

  ## In case no varset is specified, replace all NA's with a 0!
  dfla1[is.na(dfla1)] = 0

  ### merge it to the polys dataframe
  polys = merge(polys, dfla1, by="cid",all.x=TRUE)
  ## define an example raster
  lcla   = lcr0
  lcla[] = 0
  ## select fields onf interest
  allfields = colnames(dfla1)[-1]

  ## JO: 2.3.2023: group metrics, currently only implemented at landscape level, might be necessary for patch and pixel level as well?
  for (metr in metrics) {
    #metr = metrics[5]
    print(metr)
    fields = allfields[grep(sprintf("^%s_[0-9].+|^%s_NA.+",metr,metr),allfields)]
    ## stack all landscape metrics
    stla  = stack()
    ## make overview plot
    pdf(paste(respath,sprintf("Reconnectland_%s_rad%0.0fm_res_%0.0fm_incr_%0.0fm_%s.pdf",metr,radius,resolution,incr,habi),sep="/"))

  for(field in fields) {
    #field = fields[1]
    print(field)
    ### calculate mean for each raster cell (JO: 2.3.2023: only summarize non-na fields! - potentially interpolate with fillNA of SUHI package later..)
    polyr1  = fasterize::fasterize(sf=polys, raster=lcla, field=field,fun="sum")
    polyc1  = fasterize::fasterize(sf=polys, raster=lcla, field=field,fun="count")
    polymn1 = polyr1/polyc1
    names(polymn1) = field
    raster::plot(polymn1,main=names(polymn1)[1])
    stla = stack(stla,polymn1)
  } # field

  dev.off()

  ### write raster stack of land (takes quite long and should maybe only be done for certain fields)
  ## JO 2.3.2023: write as gridfile, to keep names! JO 1.8.2023: write as terra objects tif files to keep names
  #writeRaster(stla, paste(respath,sprintf("Reconnectland_%s_rad%0.0fm_res_%0.0fm_incr_%0.0fm_%s.grd",metr,radius,resolution,incr,habi),sep="/"),overwrite=TRUE)
  writeRaster(terra::rast(stla), paste(respath,sprintf("Reconnectland_%s_rad%0.0fm_res_%0.0fm_incr_%0.0fm_%s.tif",metr,radius,resolution,incr,habi),sep="/"),overwrite=TRUE)

  ### clean environment: remove not used objects
  gc()

  } # metrics

  rm(list=c("stla","polyc1","polymn1","polyr1","dfla","dfla1","polys","centers","lcla"))

  } # dfla

  ########################
  ############ patch level
  ########################
  ### 2.1) Load data outputs for scale, resolution, radius and habitat of interest
  dfpa = load_data(path=sourcepath,
                   pattern=sprintf("dfpa_rad%s.+_res_%s.+_incr_%s.+_cid.+_%s.csv",radius,resolution,incr,habi))


  if(length(dfpa)>0){
  print("\n found results at patch-level... \n")
  ### JO: set all errors in Reconnect_core to NA:
  dfpa[dfpa$value == -99 & is.finite(dfpa$value),"value"] = NA

  ### 2.3) Patch level aggregation, load shapefiles, normalize and make weighted means for each patch
  hshp1   = sf::st_read(paste(sourcepath,sprintf("hshp_res%0.0f_rad%0.0f_incr%0.0f_habi%s.shp",resolution,radius,incr,habi),sep="/"))
  dfpa    = dfpa[dfpa$level=="patch" & dfpa$pid!="no",]

  ############
  ### create a wide landscape dataset
  dfpa1 = reshape2::dcast(dfpa[,c("cid","pid","minsize","alpha","metric","value")], cid+pid~metric+alpha+minsize)

  ############
  ### prepare data for aggregation
  ## get the total patch area for each pid:
  hdf         = data.frame(hshp1)
  hdf$pid     = as.character(hdf$pid)
  hdf$area_ha = as.numeric(sf::st_area(hshp1)/10000)

  ############
  ### merge hdf data with the dfpa1 data
  dfpa1 = merge(dfpa1,hdf[,c("pid","area_ha")],by="pid",all.x=TRUE)

  ############
  ### get the number of cids per pids
  dfpa11 = dfpa1 %>% dplyr::group_by(pid) %>% dplyr::summarize(ngid = dplyr::n()) %>% data.frame

  ############
  ### merge ngid info with the dfpa1 data
  dfpa1 = merge(dfpa1,dfpa11[,c("pid","ngid")],by="pid",all.x=TRUE)

  ############
  ### create a weights column
  ### JO: make column selection more sophisticated? Currently assumes the lowest minsize arg. is listed first)
  ### JO: the division by ngid is not necessary for weighted mean function..Because of multiple overlaps, sum of paa/area_ha does not add up to 1!
  dfpa1$weights = (1/dfpa1$ngid)*dfpa1[[grep("pa_ha_NA",colnames(dfpa1))[1]]]/dfpa1$area_ha

  ############
  ### select fields of interest
  fields = colnames(dfpa1)[-which(colnames(dfpa1)%in%c("pid","cid","area_ha","ngid","weights"))]

  ############
  ### aggregate data for fields and pid using a weighted mean
  df2  = dfpa1[,c("pid",fields,"weights")]
  #df2 = data.frame(pid=c(1,1,1,2,2,3,3,3),a=c(1,2,1,4,5,9,10,10),b=c(1,1,1,1,1,1,1,1),c=c(0,0,0,0,0,0,1,1))
  #df2$weights = c(rep(1,3),c(0.5,1),c(1,2,2))
  ## somehow this throws an error
  #dfpa2   = df2 %>% group_by(pid) %>% summarise_all(sum) %>% data.frame
  dfpa2    = df2 %>% dplyr::group_by(pid) %>% dplyr::summarise_all(dplyr::funs(weighted.mean(., w=weights,na.rm=TRUE))) %>% data.frame

  ############
  ### reassign fieldnames because dplyr changed plus symbol to point symbol!!
  fields = colnames(dfpa2)[-which(colnames(dfpa2)%in%c("pid","cid","area_ha","ngid","weights"))]

  ############
  ### assign patch importance for different resolutions to original shapefile
  hshp2 = merge(hshp1,dfpa2,by="pid",all.x=TRUE)

  ############
  ### make a plot for overview
  pdf(paste(respath,sprintf("Reconnectpatch_rad%0.0fm_res_%0.0fm_incr_%0.0fm_%s.pdf",radius,resolution,incr,habi),sep="/"))
  for(field in fields) {
    print(field)
    #field = fields[28]
    #hshp3  = hshp2[is.finite(hshp2[[field]])]
    pavals = hshp2[[field]]
    ## JO: 2.3.2023 changed plotting here..
    plot(hshp2$geometry,border=NA,col="grey",main=field)
    if(length(pavals[is.finite(pavals)])>0) {
    plot(hshp2$geometry,border=NA,col=cols(105)[1+normalizef(pavals)],add=TRUE)
    #plot(hshp2$geometry,border="black",lwd=0.02,col=cols(105)[1+normalizef(pavals)],main=field)
    legval = unique(round(seq(range(pavals,na.rm=TRUE)[1],range(pavals,na.rm=TRUE)[2],length.out=5),4))
    add_legend("right",legend=legval,fill=cols(105)[1+normalizef(legval)],bty="n")
    }
  } # field
  dev.off()

  ############
  ### save shapefile
  ### JO 2.3.2023 extract the colnames and save it separately
  colnms = data.frame(colnms = colnames(hshp2))
  write.csv(colnms,paste(respath,sprintf("Reconnectpa%0.0fm%0.0fm%0.0fm%s_cols.csv",radius,resolution,incr,habi),sep='/'),row.names=FALSE)
  ### colnames need to have less than 10 characters - JO: the st_write does it automatically...
  #colnames(hshp2) = gsub("pimp","",gsub("km|ha","",gsub("1000","1m",gsub("_","",colnames(hshp2)))))
  ## not sure if driver needs to be indicated in st_write function: #driver="ESRI Shapefile"
  if(file.exists(paste(respath,sprintf("Reconnectpa%0.0fm%0.0fm%0.0fm%s.shp",radius,resolution,incr,habi),
                       sep='/'))){file.remove(paste(respath,sprintf("Reconnectpa%0.0fm%0.0fm%0.0fm%s.shp",radius,resolution,incr,habi),sep='/'))}
  sf::st_write(hshp2, paste(respath,sprintf("Reconnectpa%0.0fm%0.0fm%0.0fm%s.shp",radius,resolution,incr,habi),sep='/'))  # create to a shapefile

  rm(list=c("hshp1","hshp2","df2","dfpa1","dfpa2","dfpa"))

  } # nrow dfpa

  ########################
  ############ pixel level
  ########################
  ############
  ### 2.1) Load data outputs for scale, resolution, radius and habitat of interest
  dfpi = load_data(path=sourcepath,
                   pattern=sprintf("dfpi_rad%s.+_res_%s.+_incr_%s.+_cid.+_%s.csv",radius,resolution,incr,habi))

  if(length(dfpi)>0){
  print("\n found results at pixel-level... \n")
  ### JO: set all errors in Reconnect_core to NA:
  dfpi[dfpi$value == -99 & is.finite(dfpi$value),"value"] = NA

  ### 2.4) Pixel level aggregation (use mosaic? to be added)
  centers   = sf::st_read(paste(sourcepath,sprintf("centers_res%0.0f_rad%0.0f_incr%0.0f_habi%s.shp",resolution,radius,incr,habi),sep="/"))
  polys     = mkbuffer(centers=centers,radius=radius,form=c("square"))
  polys$cid = centers$cid
  rm(list="centers")
  ### 2.3) Patch level aggregation, load shapefiles, normalize and make weighted means for each patch
  hshp   = sf::st_read(paste(sourcepath,sprintf("hshp_res%0.0f_rad%0.0f_incr%0.0f_habi%s.shp",resolution,radius,incr,habi),sep="/"))

  ############
  ### define r0
  rr0   = lcr0
  rr0[] = NA
  ########################################
  ### same code as in Reconnect_core...
  ########################################
  ### 2.1) resample habitat raster in case the desired resolution does not match the one of the raster at hand
  if(resolution != mean(res(rr0),na.rm=TRUE)[1]) {
    ### determine factor of aggregation
    fact = resolution/mean(res(rr0),na.rm=TRUE)[1]
    rr0[is.na(rr0)] = 0
    rr0 = raster::aggregate(x=rr0,fact=fact,fun=function(x,na.rm=na.rm){sum(x)},expand=FALSE,na.rm=TRUE)
    ### Conservative approach: set all the pixels with less than 75% habitat to 0!
    rr0[] = NA
  }


  ############
  ### select fields of interest
  dfpi       =  dfpi[dfpi$level=="pixel",]
  ### create a field column
  dfpi$field = paste(dfpi$metric,dfpi$alpha,dfpi$minsize,sep="_")
  fields     = unique(dfpi$field)

  pdf(paste(respath,sprintf("Reconnectpix_rad%0.0fm_res_%0.0fm_incr_%0.0fm_%s.pdf",radius,resolution,incr,habi),sep="/"))

  for(field in fields) {
    #field = fields[1]
    print(field)

    ### extract right fieldname
    dfpi0 = dfpi[dfpi$field==field, ]
    ## extract minsize
    minsize = as.numeric(gsub("^.+_.+_(.+)$","\\1",field))

    ## make a raster mosaic list
    rlist = sapply(polys$cid,function(cid){
      #for(cid in cids){}
      #### go through cids
      #cid = 1
      foc = polys[polys$cid==cid,]
      ## JO 2.8.2023 not necessary: plot the raster...
      #raster::plot(rr0,main=field)
      ##plot(foc$geometry,add=TRUE)
      ########################################
      ### same code as in Reconnect_core...
      ########################################
      ## Generate shapefile for region of interest, include patches that are touched by window of interest (for correct minarea removal!)
      hshp1 = hshp[foc,]
      #print(sprintf("apply minsize %0.02fm...",minsize))
      shp1 = minarm(shp=hshp1,minsize=minsize)
      ## apply cutpatches argument only here, after minimum patches have been removed..!
      if(cutpatches) {shp1 = sf::st_intersection(shp1,foc$geometry)}
        ### sometimes there is no habitat in focal area
        if(length(shp1[[1]])>0) {
          ### create new foc1 as polygon
          foc1  = sf::st_union(foc,shp1) %>% sf::st_union(by_feature=FALSE)
          foc1  = st_sf(geometry=foc1)
        } else if(length(shp1[[1]])==0) {
          foc1  = st_sf(geometry=foc$geometry)
        }
        ## try if extents overlap, if not go to next cid...
        try = try(raster::crop(rr0,raster::extent(foc1)),silent=TRUE)
        if(class(try)[1]=="try-error") {print(sprintf("no overlap between cid %d and habitat data...",cid)); next}
        ### crop raster according to extent of foc (or should it be according to extent of cutpatches, which might give assymetrical landscapes??)
        rr1   = raster::crop(rr0,raster::extent(foc1))
        rr1   = raster::mask(rr1,foc1)
        ########################################
        ### extract pixel values for center ID
        dfpi1 = dfpi0[dfpi0$cid==cid,]
        ### fill raster with field values
        rr1[as.numeric(dfpi1$pid)] = dfpi1$value
        ## JO 2.8.2023: not necessary:
        #raster::plot(rr1,add=TRUE)
        ### make a terra version out of it
        #rr1 = terra::rast(rr1)
        return(rr1)
      #rlist[[which(cids)==cid]]
    })
    if(length(rlist)==1) {
    mosa = rlist[[1]]
    } else if(length(rlist)>1){
    ### define mosaic function
    rlist$fun   = mosafun
    rlist$na.rm = TRUE
    mosa = do.call(raster::mosaic, rlist)
    }

    ## Create a SpatRasterCollection
    # rsrc   = terra::sprc(rlist)
    # mosa   = terra::mosaic(rsrc,fun="mean")

    ### checkplot
    raster::plot(mosa,main=field)

    ### write raster stack of land (takes quite long and should maybe only be done for certain fields) JO 2.3.2023: write as gridfile, to keep names!
    ## JO 1.8.2023 write as terra tif file to keep layer names
    #writeRaster(mosa, paste(respath,sprintf("Reconnectpix_rad%0.0fm_res_%0.0fm_incr_%s_%s.grd",radius,resolution,field,habi),sep="/"), overwrite=TRUE)
    writeRaster(terra::rast(mosa), paste(respath,sprintf("Reconnectpix_rad%0.0fm_res_%0.0fm_incr_%s_%s.tif",radius,resolution,field,habi),sep="/"), overwrite=TRUE)

    rm(list = c("rlist","rsrc","mosa"))
  } # field

  dev.off()

  }# dfpi

  } # habi

} # function end




########
#' @name Reconnect_single
#' @title Function independent of Reconnect_wrap: apply connectivity functions in the "Reconnect_inifile.xlsx" for a given landscape
#' @description simple Reconnect core function: apply connectivity functions in the "Reconnect_inifile.xlsx" for a given landscape, without moving windows
#' @param cid character name or ID string
#' @param pidname character describing the patch-id column name in the habitat patch shapefile (hshp)
#' @param lcr  raster habitat layer where habitat > 1 and non-habitat = 0 or NA.
#' @param lres raster layer of resistance surface
#' @param foc shapefile, of focal region of interest (optional)
#' @param hshp shapefile of habitat patches (optional)
#' @param minsizev numeric vector of minimum habitat patch areas in m2
#' @param cutpatches logical, should habitat patches be cut at the boundary of the focal region of interest (foc)? Default = TRUE
#' @param cost_type character, algorithm applied to calculate distances: either "euclidean" (default), "least-cost" or "commute-time"
#' @param dist_type character, distance type: either "edge" or "centroid" (default: "edge").
#' @param trans_fun function expression, argument for the gdistance::transition() function, applied when moving between resi raster cells. Default for conductance: function(x) 1/mean(x)
#' @param neighb numeric, directions argument for the gdistance::transition() function: directions in which cells are connected (4, 8, 16, ...)
#' @param areafun function defining how patch areas are calculated in case the input is a raster or shapefile (default: sf::st_area)
#' @param cond dataframe, connectivity functions ini file (cf. 1-Data folder & readini function)
#' @param cols color ramp palette to be used (default: terrain.colors(100))
#' @param maxnpatch numeric, maximum number of patches to be treated in a landscape (so that connectivity functions don't get overwhelmed)
#' @param resolution numeric, resolution of raster grid
#' @param respath character string indicating file path where results should be saved
#' @return csv files saved to respath location, where results from the application of connectivity indicators specified in the cond file are saved (dfla = landscape-level, dfpa = patch-level, dfpi = pixel-level).
#' @export
Reconnect_single  = function(cid      = NULL,
                        pidname  = NULL,
                        lcr      = NULL,
                        lres     = NULL,
                        foc      = NULL,
                        hshp     = NULL,
                        minsizev = 0,
                        cutpatches = TRUE,
                        cost_type  = "euclidean",
                        dist_type  = "edge",
                        trans_fun  = function(x) 1/mean(x),
                        neighb     = 8,
                        areafun    = sf::st_area,
                        cond       = NULL,
                        cols       = terrain.colors(100),
                        maxnpatch  = 28500,
                        resolution = NULL,
                        respath    = NULL) {

  ###############################################################################################
  ## Generate raster or shapefile for region of interest, include patches that are touched by window of interest (for correct minarea removal!)
  if(length(lcr)==0 & length(hshp)==0) { stop("please either provide a habitat landscape raster or shape file")}
  ## change resolution to raster resolution if given
  if(length(lcr)!=0 & length(resolution)==0) { resolution = mean(res(lcr),na.rm=TRUE)[1] }
  ## Make sure there are habitat network shape and raster files
  if(length(hshp)==0) {
    ## generate an sf object
    hshp    = rstosf(rs=lcr)
    pidname = names(hshp)[1]
  } else if(length(lcr)==0) {
    ## generate a raster object
    if(length(pidname)==0) {pidname = names(hshp)[1]}
    bbx    = sf::st_bbox(hshp)
    lcr    = raster::raster(xmn=bbx[1], xmx=bbx[3], ymn=bbx[2], ymx=bbx[4], crs=sf::st_crs(hshp), resolution=resolution, vals=NA)
  }
  ## identify raster clumps and make sure patch Ids are the same in shape and raster file
  cells        = raster::extract(lcr, as_Spatial(hshp$geometry),method="simple",cellnumbers=TRUE)
  names(cells) = hshp[[pidname]]
  ## add patch numbers
  for(x in 1:length(cells)) {lcr[cells[[x]][,1]] = as.numeric(names(cells[x]))}
  ######
  ### assign a new hr object
  hr = lcr
  ### make sure that non-habitat cells are NA and not 0
  hr[hr<1] = NA
  ### name habitat raster
  names(hr) = "habitat"
  ######
  ### 3.0) add resistance to the habitat raster if provided (and otherwise create one of resitance 1...)
  if(length(lres)>0){
    re = lres
    ## aggregate to larger resolution if necessary
    if( (mean(res(lres),na.rm=TRUE)[1] != mean(res(lcr),na.rm=TRUE)[1]) | (raster::extent(lres) != raster::extent(lcr)) ) {
      re = terra::resample(x=terra::rast(re), y=terra::rast(hr), method="bilinear") %>% as("Raster")
    }
  } else if(length(lres)==0) {
    ## create a resistance layer where all resistance is set to 1 (more compatible with overall functionality)
    re   = hr
    re[] = 1
  }
  ## make sure it is named correctly
  names(re) = "resistance"
  ## make sure that infinite cells are just NA
  re[!is.finite(re)] = NA
  ## stack habitat and resistance raster
  hr = raster::stack(hr,re)

  ###############################################################################################
  ### mask to region of interest if foc is provided, and set extent
  if(length(foc)>0){
    if(sf::st_crs(hr)!=sf::st_crs(foc)) {foc = sf::st_transform(foc,crs=st_crs(hr))}
  } else if(length(foc)==0){
    foc = sf::st_sf(ID= 1, geometry=sf::st_as_sfc(sf::st_bbox(hr)))
    sf::st_crs(foc) = sf::st_crs(foc)
  }

  ###############################################################################################
  ###  write to file
  if(file.exists(paste(respath,sprintf("Single_hshp_ID%s_res%0.0f.shp",as.character(cid),resolution),sep='/'))){
    file.remove(paste(respath,sprintf("Single_hshp_ID%s_res%0.0f.shp",as.character(cid),resolution),sep='/')) }
  st_write(hshp, paste(respath,sprintf("Single_hshp_ID%s_res%0.0f.shp",as.character(cid),resolution),sep='/'), driver="ESRI Shapefile")  # create to a shapefile
  ###
  writeRaster(terra::rast(hr), paste(respath,sprintf("Single_hr_ID%s_res_%0.0fm.tif",as.character(cid),resolution),sep="/"),overwrite=TRUE)

  ###############################################################################################
  ## Create a dataframe where patch IDs and corresponding connectivity measurements are stored!!
  dfpi  = data.frame(pid=NA,minsize=NA,alpha=NA,metric=NA,value=NA)
  dfpa  = data.frame(pid=NA,minsize=NA,alpha=NA,metric=NA,value=NA)
  dfla  = data.frame(pid=NA,minsize=NA,alpha=NA,metric=NA,value=NA)

  ### in case foc is given, crop the shapefile and rasterfile with it..
  hshp1 = hshp[foc,]

  ########################################
  ## go through minsize argument
  for (minsiz in sort(minsizev,decreasing=FALSE)){
    ## minsiz = sort(minsizev,decreasing=FALSE)[1]
    print(sprintf("apply minsize %0.02fm...",minsiz))
    shp1 = minarm(shp=hshp1,minsize=minsiz)

    ## apply cutpatches argument only here, after minimum patches have been removed..!
    if(cutpatches==TRUE) {shp1 = sf::st_intersection(shp1,foc$geometry)}

    ### sometimes there is no habitat in focal area
    if(length(shp1[[pidname]])>0) {
      ### create new foc1 as polygon
      foc1  = sf::st_union(foc,shp1) %>% sf::st_union(by_feature=FALSE)
      foc1  = st_sf(geometry=foc1)
    } else if(length(shp1[[pidname]])==0) {
      foc1  = st_sf(geometry=foc$geometry)
    }

    ## try if extents overlap, if not go to next cid...
    try = try(raster::crop(hr,raster::extent(foc1)),silent=TRUE)
    if(class(try)[1]=="try-error") {print(sprintf("no overlap between cid %s and habitat data...",as.character(cid))); next}

    ### crop raster according to extent of foc (or should it be according to extent of cutpatches, which might give assymetrical landscapes??)
    hr1   = raster::crop(hr,raster::extent(foc1))
    hr1   = raster::mask(hr1,foc1)

    ### determine hr and resistance rasters
    lcr2 = hr1[["habitat"]]
    if(nlayers(hr1)==2) {re1 = hr1[["resistance"]]} else if(nlayers(hr1)<2) {re1 = NULL}

    ### set all habitat patches in raster (determined by pixel values) not contained in shp1 to NA!
    if(length(shp1[[pidname]])>1) {
      lcr2[!values(lcr2)%in%shp1[[pidname]]] = NA
    } else if(length(shp1[[pidname]])==1) {
      lcr2[values(lcr2)!=shp1[[pidname]]] = NA
    } else if(length(shp1[[pidname]])==0) {
      lcr2[] = NA
    }

    ### lcr1: make 1 habitat raster for connfunctions, where habitat cells = 1 and the rest = NA
    lcr1 = lcr2
    lcr1[!is.na(lcr1)] = 1

    ### derive resolution from input data
    resolution = mean(res(hr1),na.rm=TRUE)

    ### short plot
      plot(foc$geometry,main=sprintf("cid %s focal area \n habitat and resistance \n minsize %0.02fm %0.0fm res",cid,minsiz,resolution))
      if(length(re1)>0){raster::plot(re1,add=TRUE)}
      raster::plot(lcr2,add=TRUE,col="darkgreen",legend=FALSE)
      plot(foc1$geometry,add=TRUE,border="darkorange",lwd=2)

    ########################################
    ## 2) Start Reconnect calculations for area of interest, foc1a, pa & mdist are always calculated
    ########################################
    ## determine total area of foc1 in any case dfla$foc_ha = foc1a/10000
    foc1a = as.numeric(areafun(foc1))
    dfl   = data.frame(pid="all",minsize=minsiz,alpha=NA,metric="foc_ha",value=foc1a/10000)
    dfla  = rbind(dfla,dfl)
    ## determine number of patches in landscape
    dfl   = data.frame(pid="all",minsize=minsiz,alpha=NA,metric="npatch",value=length(shp1[[pidname]]))
    dfla  = rbind(dfla,dfl)
    ## determine total and mean habitat patch area in landscape
    pa_ha = as.numeric(sf::st_area(shp1))
    dfl   = data.frame(pid="all",minsize=minsiz,alpha=NA,metric=c("tpa_ha","mnpa_ha"),value=c(sum(pa_ha,na.rm=TRUE)/10000,mean(pa_ha,na.rm=TRUE)/10000))
    dfla  = rbind(dfla,dfl)

    ########################################
    ## derive patch areas in case there are some habitat patches...
    if(length(shp1[[pidname]])>0 & length(shp1[[pidname]])<=maxnpatch){

      ########################################
      ## determine patch id names
      pid = shp1[[pidname]]

      ########################################
      ## distance matrix: is quite slow, only make the distance matrix once
      ## make the euclidean distance matrix
      if(cost_type == "euclidean") {
        ## get the patch area and euclidean distance matrix object
        pamd  = eucdist_pm(x=shp1,id=pidname,areafun=sf::st_area,dist_type=dist_type)
        pa1    = pamd$pa
        mdist1 = pamd$mdist
        ## remove to save memory
        rm(list="pamd")
      } # euclidean

      ## make the resitance distance matrix
      if(length(re1)>0 & cost_type == "least-cost") {
        ## get the patch area and cost distance matrix object
        pamd = costdist_pm(x=shp1,resi=re1,id=pidname,areafun=areafun,
                           trans_fun=trans_fun,neighb=neighb,
                           dist_type=dist_type,cost_type=cost_type)
        pa1    = pamd$pa
        mdist1 = pamd$mdist
        ## remove to save memory
        rm(list="pamd")
      } # least-cost

      ## add area of patches  (is quite fast..)

      ########################################
      ## Always add the patch area (important for summary after)
      dfp  = expand.grid(pid=pid,minsize=minsiz)
      dfp[["alpha"]]  = NA
      dfp[["metric"]] = "pa_ha"
      dfp[["value"]]  = pa1/10000
      dfpa = rbind(dfpa,dfp)

      ########################################
      ## Apply connectivity functions for specified alphas
      ########################################

      ## go through connectivity metrics
      for(conn in unique(cond$conname)) {
        ## conn    = unique(cond$conname)[1]
        condd     = cond[cond$conname==conn,]
        print(condd)
        resname   = gsub(" ","",strsplit(condd[,"resname"],split=",")[[1]])
        confun    = eval(parse(text=condd[,"confun"]))
        statement = condd[,"statement"]
        alphav    = strsplit(condd[,"alphav"],split=",")[[1]]
        alphav    = alphav[alphav!="NA"] %>% as.numeric()
        minssz    = strsplit(condd[,"minsizev"],split=",")[[1]]
        minssz    = minssz[minssz!="NA"] %>% as.numeric()
        level     = unique(strsplit(condd[,"level"],split=",")[[1]])
        resultnr  = strsplit(condd[,"resultnr"],split=",")[[1]]
        resultnr  = resultnr[resultnr!="NA"]%>% as.numeric()
        aggfun    = eval(parse(text=condd[,"aggfun"]))

        ## check if alpha's are used
        alphav = alphav[is.finite(alphav)]
        ## check if minsize is specified
        minssz = minssz[is.finite(minssz)]
        if(length(minssz)>0) {
          ## assign corresponding alphav if minsize is given for a specific alpha!
          alphav = alphav[which(minssz == minsiz)]
        }
        ## apply functions
        if(length(alphav)>0) {
          alphav = alphav[order(alphav,decreasing=FALSE)]
          ## try function first
          try = try(sapply(alphav,function(alpha1){eval(parse(text = paste(statement)))}),silent=FALSE)
          if(class(try)[1]=="try-error"){ conres = -99
          } else if(class(try)[1]!="try-error"){
            conres = sapply(alphav,function(alpha1){eval(parse(text = paste(statement)))})
            colnames(conres) = alphav
            if(length(conres)==0)  {conres = -99 }
            if(length(resultnr)>0) {conres = conres[resultnr,]}
          }
          if("pixel" %in% level){
            if(is.list(conres)){pp = unlist(conres)} else {pp = as.vector(conres)}
            dfpp = expand.grid(pid=1:(length(pp)/length(alphav)),minsize=minsiz,metric=resname,alpha=alphav)
            dfpp[["value"]]  = pp
            dfpp = dfpp[,c("pid","minsize","alpha","metric","value")]
            dfpi = rbind(dfpi,dfpp)
          }
          if("patch" %in% level){
            dfp = expand.grid(pid=pid,minsize=minsiz,metric=resname,alpha=alphav)
            dfp[["value"]]  = unlist(conres)
            dfp  = dfp[,c("pid","minsize","alpha","metric","value")]
            dfpa = rbind(dfpa,dfp)
          }
          if("landscape"%in%level){
            dfl            = expand.grid(pid="all",minsize=minsiz,metric=resname,alpha=alphav,value=NA)
            dfl[["value"]] = unlist(lapply(conres,aggfun))
            dfl = dfl[,colnames(dfla)]
            dfla = rbind(dfla,dfl)
          }
        } else if(length(alphav)==0){
          ## try function first

          try = try(eval(parse(text = paste(statement))),silent=FALSE)
          if(class(try)[1]=="try-error"){ conres = -99
          } else if(class(try)[1]!="try-error"){
            conres = eval(parse(text = paste(statement)))
            if(length(conres)==0)  {conres = -99 }
            if(length(resultnr)>0){conres = conres[[resultnr]]}
          }
          if("pixel" %in% level){
            pp = unlist(conres)
            dfpp = expand.grid(pid=1:length(pp),minsize=minsiz)
            dfpp[["alpha"]]  = NA
            dfpp[["metric"]] = resname
            dfpp[["value"]]  = pp
            dfpi = rbind(dfpi,dfpp)
          }
          if("patch" %in% level){
            dfp = expand.grid(pid=pid,minsize=minsiz)
            dfp[["alpha"]]  = NA
            dfp[["metric"]] = resname
            dfp[["value"]]  = unlist(conres)
            dfpa = rbind(dfpa,dfp)
          }
          if("landscape"%in%level){
            dfl  = expand.grid(pid="all",minsize=minsiz,metric=resname,alpha=NA,value=NA)
            if(length(conres)>nrow(dfl)){conres = aggfun(unlist(conres))} else if(length(conres)==0) { conres = -99 }
            dfl[["value"]] = conres
            dfl = dfl[,colnames(dfla)]
            dfla = rbind(dfla,dfl)
          }
        }
      } # for all conn functions


    } # apply confun only if nr of patches are > 0 and < maxnrpatch

  } # minsize

  ########################################
  ## finish datasets for window of interest, write to file..
  dfpi$level = "pixel"
  dfpa$level = "patch"
  dfla$level = "landscape"


  ## write to file if at least 1 row!
  if(nrow(dfpi)>1){
    dfpi = dfpi[-1,]
    dfpi$cid    = cid
    dfpi$resolution  = resolution
    write.csv(dfpi,paste(respath,sprintf("Single_dfpi_res_%0.0fm_cid_%s.csv",resolution,as.character(cid)),sep="/"),row.names=FALSE)
  }

  ## write to file if at least 1 row!
  if(nrow(dfpa)>1){
    dfpa = dfpa[-1,]
    dfpa$cid    = cid
    dfpa$resolution  = resolution
    write.csv(dfpa,paste(respath,sprintf("Single_dfpa_res_%0.0fm_cid_%s.csv",resolution,as.character(cid)),sep="/"),row.names=FALSE)
  }

  ## write to file if at least 1 row!
  if(nrow(dfla)>1){
    dfla = dfla[-1,]
    dfla$cid    = cid
    dfla$resolution  = resolution
    write.csv(dfla,paste(respath,sprintf("Single_dfla_res_%0.0fm_cid_%s.csv",resolution,as.character(cid)),sep="/"),row.names=FALSE)
  }


}# function end


############
### Adapted from Hadley on stack overflow (https://stackoverflow.com/users/16632/hadley)
### Stack overflow: https://stackoverflow.com/questions/23190280/whats-wrong-with-my-function-to-load-multiple-csv-files-into-single-dataframe
#' @name load_data
#' @title load csv data files from source
#' @description import separate csv files as one coherent dataset
#' @param path character string describing path where csv data are stored
#' @param pattern character string describing regex pattern of csv files that should be loaded
#' @return a numeric vector of buffer distances to use for given alphas
#' @export
load_data = function(path=NULL,pattern=NULL) {
  fls     = dir(path, pattern = pattern, full.names = TRUE)
  tbls    = lapply(fls, read.csv)
  do.call(rbind, tbls)
}

############
## Add legend outside plot margins
## Adapted from Jan van der Laan on Stack Overflow, https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
#' @name add_legend
#' @title add legend
#' @description add legend outside of plot margins
#' @return a function enabling plotting on outside plot margins
#' @export
add_legend = function(...) {
  opar = par(fig=c(0.02, 0.97, 0.02, 0.97), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

##########################
#' @name fillNAr
#' @title fill NA values in raster as simple mean of surrounding cells
#' @description simple linear interpolation of raster datagaps using the raster::focal() function
#' @param x raster layer with NAs
#' @param roi shapefile of region of interest (optional)
#' @param nmax maximum nr of iterations (in case window size is small and NA gaps are large)
#' @param ws numeric value, must be an odd number, describing maximum window size applied in "focal" function of the raster package
#' @return a raster layer with filled NA gaps, derived from linear interpolation for NA cells in windows with iteratively decreasing sizes, starting with ws
#' @export
fillNAr = function(x=NULL,roi=NULL,nmax=100,ws=5) {

  ## make a general naidx1
  naidx1 = 1:ncell(x)

  ## if roi is given, only account for cells within it!
  if(exists("roi")) {
    if(length(roi)>0){
      if(raster::compareCRS(x,roi)==FALSE) {
        roi = sf::st_transform(roi,crs = sf::st_crs(x))
      }
      x = raster::crop(x,raster::extent(roi))
      naidx1 = raster::extract(x,roi,cellnumbers=TRUE)[[1]][,"cell"]
    }}

  ## setting up loop: identify all NA raster cells
  vals    = values(x)
  naidx   = which(is.na(vals))
  naidx   = naidx[naidx%in%naidx1]
  rm(list="naidx1")

  ## loop through window sizes starting with largest to fill gaps..
  for(wss in rev(seq(1,ws,by=2))){
    print(sprintf("window size %d...",wss))
    ## determine navals and length of NAs
    navals = x[naidx]
    diff   = ll = length(navals[is.na(navals)])
    nround = 0
    while(nround <=nmax){
      print(sprintf("filling round %d..",nround))
      x = raster::focal(x, w=matrix(1,nrow=wss, ncol=wss), fun=mean, NAonly=TRUE, na.rm=TRUE)
      navals = x[naidx]
      diff   = ll - length(navals[is.na(navals)])
      ll     = length(navals[is.na(navals)])
      if(diff== 0) {
        print("no further improvement with current window size")
        break
      }
      nround = nround+1
    } # while
  } ## wss

  if(exists("roi")) {
    if(length(roi)>0){
      x = raster::mask(x,roi)
    }}

  return(x)
}

##########################################################
##### project:       Multispecies Connectivity Modelling
##### author:        Jacqueline Oehri (JO), jacqueline.oehri@gmail.com
##########################################################
