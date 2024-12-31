##########################################################
### Neutral Landscapes for Multispecies Connectivity
### NLMC functions
### Functions to create a set of neutral or simulated landscapes
### Jacqueline Oehri, jacqueline.oehri@gmail.com
### 13.01.2022
##########################################################


########
#' @name mkcentersi
#' @title make centerpoints of potential windows of selection with given distance apart, within extent where "incr" is subtracted
#' @description get an sf object of by incr regularly spaced centerpoints in a given extent (ext)
#' @param ext a raster::extent object
#' @param incr numeric value defining distance among centerpoints
#' @return a shapefile containing a regular grid of centerpoints within ext
#' @export
mkcentersi = function(ext=NULL,incr=NA) {

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
#' @name mkbuffsim
#' @title make a square or circular buffer for a given coordinate
#' @description make a buffer sf object around a given sf spatial points center point
#' @param centers a sf object (including crs, a projection with unit = meters) with coordinates for center points of buffer locations
#' @param radii numeric vector of the desired radii of the circular buffer, or half the side length of the square buffer in meters
#' @param form the form of the buffer: currently only circular or square
#' @return a sf object with corresponding buffer
#' @export
mkbuffsim = function(centers,radii,form="circle") {

  crs   = st_crs(centers)
  plist = list()
  ## make equal radii in case radii have not the same length as number of centers
  if(length(radii)==1|length(radii)!=length(centers$geometry)){ radii = rep(mean(radii,na.rm=TRUE),length(centers$geometry))}

  ## actual start
  if(form=="circle"){
    for(ii in 1:nrow(centers)) {
      #print(ii)
      cords  = sf::st_coordinates(centers$geometry[ii])
      radius = radii[ii]
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
      radius = radii[ii]
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


#' @name stratsamp
#' @title stratified random sampling in strata of interest, with minimum distance argument
#' @description sample an equal number of replicates for strata of interest, such as urban and non-urban land covers, with minimum distance criterion
#' @param rast a raster or raster stack containing the strata from which to sample
#' @param strata a character string describing the stratum raster layer name
#' @param nsamp numeric value indicating the desired number of replicates for each stratum type
#' @param mindist numeric value indicating the minimum distance among samples withing a stratum
#' @param remedge logical -  should samples have a minimum distance (mindist) to the raster edge (TRUE) or not (FALSE, default)
#' @param maxiter numeric value indicating the maximum number of iterations for searching the optimal sample set per stratum
#' @return a list of 3 lists: idx = raster cell index values of samples per stratum, fracmind=fraction of samples undercutting mindist threshold, length=number of samples
#' @export
stratsamp = function(rast=NULL,strata=NULL,nsamp=10,mindist=100,remedge=FALSE,maxiter=500) {
  #### create dataframe for sampling! (add some spatial stratification?!)
  #### sample a larger subset than the number of replicates, so that you ideally could select a subset from that!
  dimy = dim(rast)[1]
  dimx = dim(rast)[2]
  ## create sample index
  indidx   =  c(1:ncell(rast))
  ## in case no strata is defined:
  if(length(strata)==0) { strata = names(rast[[1]])}
  ## create strata names index
  stratnam = values(rast[[strata]])

  ## remove indidx cellnumbers closer than mindist from edge
  if(remedge==TRUE){
    rowid   = ceiling(indidx/dimx)
    colid   = ceiling(indidx-((rowid-1)*dimx))
    remrow  = which(rowid <= (mindist*0.5) | rowid >= (dimy-(mindist*0.5)) )
    remcol  = which(colid <= (mindist*0.5) | colid >= (dimx-(mindist*0.5)) )
    torem   = unique(c(remrow,remcol))
    ## check:
    # indrem  = indidx[torem]
    # r = rast
    # r[] = 0
    # r[indrem] = 999
    # plot(r)
    ## make updated indidx
    indidx = indidx[-torem]
    stratnam = stratnam[-torem]
  }

  ## create an empty indexlist
  idxlist   = list()
  fracmind  = list()
  idxlength = list()

  for(s in unique(stratnam[is.finite(stratnam)])) {
    ## loop through strata
    print(s)
    ## set the new indidx to sample from
    indidx1 = indidx[which(stratnam%in%s)]
    ## set strata specific sampling
    nsamp1 = nsamp
    ## preparation
    nogo     = c()
    ind      = c()
    finished = FALSE
    nround   = 1

    ## sample function with mindist argument
    while(nround<=maxiter & finished==FALSE & nsamp1>0) {
      ## make sure sample size is not larger than ground sample
      maxl = length(indidx1[!indidx1%in%c(ind,nogo)])
      if(nsamp1>maxl){
        #print(sprintf("need to reduce sample size strata %s to %d",s,maxl))
        nsamp1 = maxl}
      ## take a sample from the indidx
      if(nsamp1 > 0) {
        ## only take a sample if size1 is not below or equal to 0
        size1    = (nsamp1-length(ind))
        if(size1 > 0) {
          ind      = c(ind,sample(x=indidx1[!indidx1%in%c(ind,nogo)],size=size1,replace=FALSE))
        }
        ## create a vector of centers for all
        allc     = data.frame(xyFromCell(rast[[strata]],cell=as.numeric(as.character(ind))))
        ## make an sf object (has the advantage of already having a crs)
        allc     = sf::st_as_sf(allc,coords=c("x", "y"),crs=st_crs(rast[[strata]]))
        mdist    = st_distance(allc)
        diag(mdist) = NA
        rownames(mdist) = colnames(mdist) = ind
        #print(mdist)
        ## identify indices that are too close to each other
        idrep = apply(mdist,MARGIN=1,FUN=function(x){names(which(x<(mindist)))}) %>% lapply(length) %>% unlist
        #print(idrep)

        ## replace half of the points with too close a distance..
        if(length(idrep[idrep>0])!=0 & length(idrep[idrep>0])<nsamp1){
          ## replace half of the ids that have the maximum number of small distance connections & that are above 0!
          torep1 = names(idrep[which(idrep>0&idrep==max(idrep))])
          idrep1 = sample(torep1,size=min(length(torep1),floor(nsamp1*0.5)))
          ## add the too close points to the nogo set
          nogo  = c(nogo,idrep1)
          ## sample the number of replicates from the cells that satisfy the distance criterion
          ind  = rownames(mdist)[which(!rownames(mdist)%in%nogo)]
          ## the case of
        } else if(length(idrep[idrep>0])>=nsamp1) {
          ## replace ids that have the maximum number of small distance connections & that are above 0!
          idrep1 = rownames(mdist)
          ## add the too close points to the nogo set
          nogo   = c(nogo,idrep1)
          ## sample the number of replicates from the cells that satisfy the distance criterion
          ind  = rownames(mdist)[which(!rownames(mdist)%in%nogo)]
        } else if(length(idrep[idrep>0])==0) {
          finished = TRUE
          print(sprintf("%d samples finished after %d rounds",nsamp1,nround))
        }
        nround= nround+1
      } else if (nsamp1<=0) {
        # keep all the settings from last round
        print(sprintf("strata %s has less than 3 samples for mindist %0.01f m...",s,mindist))
      }
    } #while

    ### add actual pixel indices to indexlist
    idxlist[[as.character(s)]]  = ind
    ### add fraction of connections that are closer than mindist
    mdist1 = as.numeric(mdist)
    fracmind[[as.character(s)]] = (length(mdist1[mdist1<mindist])-nsamp1)/(length(mdist1)-nsamp1)
    ### add length of index
    idxlength[[as.character(s)]]= length(ind)

    print(sprintf("finished strata %s going to next...",s))

  } # strata

  return(list(idx=idxlist,fracmind=fracmind,length=idxlength))

} # function end


########
#' @name simpcirc_gen
#' @title Generate simple circle raster landscapes
#' @description wrapper function applying the simpcirc function: simulate landscapes with defined number of patches, fraction of habitat (note: because of the random generation process, landscapes can slightly deviate from desired settings in hfr, npatch, sdpa and cf)
#' @param simname character, name of simulation run
#' @param ncs numeric value defining number of cells in x and y dimension of the square raster (default = 500)
#' @param hfrs numeric vector defining desired fraction of habitat in within the whole landscape (note that the larger hfr, the longer it takes to randomly sample habitat patch locations, especially if the clumping factor cf < 1, i.e. if they should be spread out)
#' @param npatches numeric vector defining desired number of habitat patches in the simulated landscapes
#' @param sdpas numeric vector defining desired standard deviation of habitat patch areas - if sd = 0 (default), all habitat patches have the same size
#' @param cfs numeric vector defining desired clumping factor of habitat patch distribution - if cf = 1 (default), habitat patches are non-overlapping; if cf<1 patches are sampled to be further apart, if cf>1,patches are sampled to be closer and overlapping, likely leading to a reduction in total habitat area, lower desired npatch and different forms of habitat patches.
#' @param reps numeric value indicating number of replicates per run
#' @param respath character string indicating file path where results should be saved
#' @param remedge logical -  should samples have a minimum distance (mindist) to the raster edge (TRUE) or not (FALSE, default)
#' @return a raster and/or a shapefile containing a landscape with habitat patches with above defined characteristics
#' @export
simpcirc_gen = function(simname="test",ncs=500,hfrs=0.17,npatches=c(1:10),sdpas=0,cfs=1,reps=1,respath=NULL,remedge=TRUE) {

  ## set warnings off..
  options(warn=-1)

  ## create repository to save landscapes of a certain simulation run
  if(!dir.exists(paste(respath,simname,sep="/"))){dir.create(paste(respath,simname,sep="/"))}

  ## to avoid problems with same names: add a "1" to all the variable names...
  remedge1 = remedge

  ### start generation of rasters
  for(nc1 in ncs) {
    #nc1      = ncs[1]
    print(nc1)
    for(hfr1 in hfrs) {
      #hfr1=hfrs[3]
      print(hfr1)
      ## create an empty log dataframe
      print("setting up log...")
      log     = data.frame(simname=simname,nc=nc1,hfr=hfr1,npatch=NA,sdpa=NA,cf=NA,rep=NA,metric=NA,value=NA)
      for(npatch1 in npatches) {
        #npatch1=npatches[1]
        for(sdpa1 in sdpas) {
          #sdpa1 = sdpas[1]
          for(cf1 in cfs) {
            #cf1=cfs[1]
            for(re1 in c(1:reps)) {
            ### make landscape simulation: raster and shapefile result
            res = simpcirc(dimx=nc1,dimy=nc1,hfr=hfr1,npatch=npatch1,sdpa=sdpa1,cf=cf1,form="circle",return="all",remedge=remedge1)
            ### unify new shape
            sps1 = sf::st_union(res$sps,by_feature=FALSE) %>% st_cast("POLYGON") %>% st_sf # st_cast(sps1,"POLYGON")
            ### log statistics
            ## total habitat area (slightly different if calculated from shapefile or raster...)
            ## total intended area
            ta      = sum(as.numeric(sf::st_area(res$sps))) # old area
            ## total new area
            #tna     = sum(res$rast[res$rast==1])
            tna    = sum(as.numeric(sf::st_area(sps1))) # new area
            ## total nr of new patches
            tnpatch = length(sps1$geometry)
            ## total habitat perimeter
            tnperim = sum(perimf(sps1))
            ## total nearest neighbour distance
            tnndf   = sum(NNDf(sf::st_distance(sps1)))
            ## add calculations to log statistics
            log1    = data.frame(simname=simname,nc=nc1,hfr=hfr1,npatch=npatch1,sdpa=sdpa1,cf=cf1,rep=re1
                                 ,metric=c("ta","tna","tnpatch","tnperim","tnndf"),value=c(ta,tna,tnpatch,tnperim,tnndf))
            log     = rbind(log,log1)

            # ### plot result
            # raster::plot(res$rast,main=sprintf("tna=%0.0f hfr = %0.02f \n npatch = %0.0f cf = %0.02f",tna,hfr1,npatch1,cf1),cex.main=0.8)
            # plot(sps1$geometry,add=TRUE,border="darkgreen")

            ### write to file ## Actually it takes much less space to save the 0 as 0 and saving them as -99 or NA takes up 1.5 times the space
            writeRaster(terra::rast(res$rast),paste(respath,simname,
                                       sprintf("%s_nc_%0.02f_hfr_%0.0f_npatch_%0.02f_cf_%d_rep.tif",simname,hfr1,npatch1,cf1,re1),sep="/"),overwrite=TRUE)

            ### write to shapefile
            # if(file.exists(paste(respath,simname,
            #                      sprintf("%s_nc_%0.02f_hfr_%0.0f_npatch_%0.02f_cf.shp",simname,hfr,npatch,cf),sep="/"))){
            #   file.remove(paste(respath,simname,
            #                     sprintf("%s_nc_%0.02f_hfr_%0.0f_npatch_%0.02f_cf.shp",simname,hfr,npatch,cf),sep="/")) }
            #   st_write(res$sps, paste(respath,simname,
            #                      sprintf("%s_nc_%0.02f_hfr_%0.0f_npatch_%0.02f_cf.shp",simname,hfr,npatch,cf),sep="/"), driver="ESRI Shapefile")  # create to a shapefile
            } # reps
          } # cfs
        } # sds
      } # npatchs
      ### write log for each nc and hfr
      log = log[-1,]
      write.csv(log,paste(respath,simname,sprintf("%s_%d_nc_%0.02f_hfr.csv",simname,nc1,hfr1),sep="/"),row.names=FALSE)
    } # hfrs
  } # ncs

} # function end


########
#' @name simpcirc
#' @title Simple circle habitat patch simulation
#' @description simulate landscapes with defined number of patches, fraction of habitat (note: because of the random generation process, landscapes can slightly deviate from desired settings in hfr, npatch, sdpa and cf)
#' @param dimx number of raster cells in the x dimension
#' @param dimy number of raster cells in the y dimension
#' @param hfr desired fraction of habitat in within the whole landscape (note that the larger hfr, the longer it takes to randomly sample habitat patch locations, especially if the clumping factor cf < 1, i.e. if they should be spread out)
#' @param npatch desired number of habitat patches in the landscape
#' @param sdpa desired standard deviation of habitat patch areas - if sd = 0 (default), all habitat patches have the same size
#' @param cf desired clumping factor of habitat patch distribution - if cf = 1 (default), habitat patches are non-overlapping; if cf<1 patches are sampled to be further apart, if cf>1,patches are sampled to be closer and overlapping, likely leading to a reduction in total habitat area, lower desired npatch and different forms of habitat patches.
#' @param form desired form of the habitat patches, currently available: "circle" (default) and "square".
#' @param return return a raster file ("raster",default), a shapefile ("shapefile") or a raster and shapefile ("all")?
#' @param remedge logical -  should samples have a minimum distance (mindist) to the raster edge (TRUE) or not (FALSE, default)
#' @return a raster and/or a shapefile containing a landscape with habitat patches with above defined characteristics
#' @export
simpcirc = function(dimx=500,dimy=500,hfr=0.17,npatch=5,sdpa=0,cf=1,form="circle",return="raster",remedge=TRUE) {

  ## create a raster with the specified dimensions, and cell resolution of 1
  rast = raster(nrows=dimy,ncols=dimx,xmn=0,xmx=dimx,ymn=0,ymx=dimy,crs="")
  ## extract nr of cells
  ncr = ncell(rast)
  ## fill raster with zeros
  rast[] = rep(0,ncr)
  ## determine total raster area
  totA     = ncr*(res(rast)[1]*res(rast)[2])
  ## determine the fractional area of habitat in landscape
  rfrA     = hfr*totA
  ## determine average area of habitat patches for a given habitat fraction
  mnpa     = rfrA/npatch
  ## determine standard deviation
  sdpa     = sdpa*mnpa
  ## make a vector for patch areas with given sd
  print(sprintf("creating patch areas with mean pa=%0.1f and sdpa=%0.1f",mnpa,sdpa))
  mpa      = rnorm(n=npatch, mean=mnpa, sd=sdpa)
  # print(mpa)
  ## in case sd is so big that some patches have a negative value...
  mpa      = mpa[mpa>0]
  ## sample landscape for final number of patches
  npatch   = length(mpa)
  print(sprintf("sample landscape for final number of patches = %d",npatch))
  ### assign radii corresponding to the areas:
  radii=sqrt(mpa/pi)
  ### actually create a minimum distance argument based on the radius of interest! (you could also take mean radius..)
  ### multiply with clumping factor!!
  mindist = mean(2*radii,na.rm=TRUE)*(1/cf)
  print(sprintf("minimum distance = %0.1f units",mindist))
  ### stratsamp function with mindist only works well when hfr<=0.5
  if(hfr>0.5){ print("habitat fraction is above 50% - set to max 50%"); hfr = 0.5}
  ####
  ### apply stratsamp function: with mindist argument! strata (in this case only 1) does not need to be defined
  ind0 = stratsamp(rast=rast,strata=NULL,nsamp=npatch,mindist=mindist,remedge=remedge,maxiter=500)
  ### only extract index values, in case stratsamp did not yield a result,  sample randomly..
  if(length(ind0$idx[[1]])>=1){ind = ind0$idx[[1]]} else if (length(ind0$idx[[1]])<1){
    #ind = (((dimy/2)*dimx)-(dimx/2))
    ind = sample(1:ncell(rast),size=npatch,replace=FALSE)}
  ## create a vector of centers for all
  allc      = data.frame(xyFromCell(rast,cell=as.numeric(as.character(ind))))
  allc$name = ind
  ## make an sf object (has the advantage of already having a crs)
  allc     = sf::st_as_sf(allc,coords=c("x", "y"),crs=st_crs(rast))
  # mdist    = st_distance(allc)
  # mdist    = matrix(mdist, dim(mdist)[1], dim(mdist)[2])
  # rownames(mdist) = colnames(mdist) = ind
  ####
  ### use the mkbuffsim function to actually make polygons
  sps = mkbuffsim(centers=allc,radii=radii,form="circle")
  ### make a raster clump where circles are!
  vals    = unique(unlist(raster::extract(rast, as_Spatial(sps$geometry),method="simple",cellnumbers=TRUE)))
  vals    = vals[is.finite(vals) & vals>0]
  rast[vals] = 1
  # plot(rast,col=c(NA,"black"),main=sprintf("hfr = %0.02f-npatch = %0.0f-sdpa = %0.02f",hfr,npatch,sdpa),cex.main=0.8)
  # plot(sps$geometry,add=TRUE,border="red")

  if(return == "raster"){
  return(rast)
  } else if (return=="shapefile"){
  return(sps)
  } else if (return =="all"){
  return(list(rast=rast,sps=sps))
  }
  } # function end


########
#' @name randclust_gen
#' @title Generate random cluster landscapes
#' @description wrapper function applying the random cluster algorithm developped by developped by [Saura and Martinez-Millan 2000](https://link.springer.com/article/10.1023/A:1008107902848) and implemented in the [NLMR package](https://github.com/ropensci/NLMR/) by [Sciaini et al. 2018](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13076).
#' @param simname character string, name of simulation run
#' @param ncs numeric value defining number of cells in x and y dimension of the square raster (default = 500, algorithm runs faster the smaller ncs)
#' @param hfrs numeric vector defining desired fraction of habitat in the landscape
#' @param cfs numeric vector defining the desired "clumping factor" realized by the proportion of elements randomly selected to form clusters: small values give more fine scale heterogeneity
#' @param reps numeric value indicating number of replicates per run
#' @param respath character string indicating file path where results should be saved
#' @param ncores number of cores to be used when parallel processing
#' @return a raster file containing a landscape with habitat patches with above defined characteristics, saved to respath/simname
#' @export
randclust_gen = function(simname="test",ncs=500,hfrs=0.17,cfs=1,reps=1,respath=NULL,ncores=1) {

  ## create repository to save landscapes of a certain simulation run
  if(!dir.exists(paste(respath,simname,sep="/"))){dir.create(paste(respath,simname,sep="/"))}

  ## start generation of rasters
  for(nc1 in ncs) {
    #nc1      = ncs[1]
    print(nc1)
    print(sprintf("starting parallel ncores = %d..",ncores))
    ## for parallelized treatment of patch IDs
    doParallel::registerDoParallel(cores=ncores)
    #apply function for each patch ID
    foreach::foreach(hfr1=hfrs,
                     .packages= unique(c("sf","raster","stars","RandomFields","NLMR","landscapetools","landscapemetrics")),
                     .export  = ls(globalenv())
                     #.combine = rbind
    ) %dopar% {
    # for(hfr1 in hfrs) {
      #hfr1=hfrs[1]
      print(hfr1)
      ## create an empty log dataframe
      print("setting up log...")
      log     = data.frame(simname=simname,nc=nc1,hfr=hfr1,cf=NA,rep=NA,metric=NA,value=NA)
          for(cf1 in cfs) {
            #cf1=cfs[1]
            print(cf1)
            for(re1 in c(1:reps)) {
            ### make landscape simulation: raster and shapefile result
            a   = Sys.time()
            res =  NLMR::nlm_randomcluster(ncol=nc1,nrow=nc1,
                                     resolution = 1,
                                     p=cf1,                   # small values give more fine scale heterogeneity & large values give more larger scale heterogeneity..gives the proportion of elements randomly selected to form clusters -> if you have high proportion (e.g.0.9) then equal ai proportions are less likely to be realized
                                     ai = c((1-hfr1),hfr1),   # gives the fractional cover of each of the land cover types (if you list 3 then you will have 3 types)
                                     neighbourhood = 4, rescale = TRUE)
            Sys.time()-a

            ### identify raster clumps, i.e. generate IDs!
            res1   = landscapemetrics::get_patches(res, directions=4, class=1)[[1]][[1]]
            print(sprintf("found %0.0f patches..",max(values(res1),na.rm=TRUE)))

            ######
            ### 2.3) if it is a raster, make a shapefile out of it, to be able to exclude patches below minsize..
            sps1 = stars::st_as_stars(res1) %>% sf::st_as_sf(merge = TRUE)
            sps1 = sf::st_make_valid(sps1,reason=TRUE)
            names(sps1)[[1]] = "pid"

            ### log statistics
            ## total habitat area (slightly different if calculated from shapefile or raster...)
            ## total intended area
            ta = hfr1*(nc1*nc1)
            ## total new area
            #tna   = sum(res[res==1])
            tna    = sum(as.numeric(sf::st_area(sps1))) # new area
            ## total nr of new patches
            tnpatch = length(sps1$geometry)
            ## total habitat perimeter
            tnperim = sum(perimf(sps1))
            ## total nearest neighbour distance
            tnndf   = sum(NNDf(sf::st_distance(sps1)))
            ## add calculations to log statistics
            log1    = data.frame(simname=simname,nc=nc1,hfr=hfr1,cf=cf1,rep=re1
                                 ,metric=c("ta","tna","tnpatch","tnperim","tnndf"),value=c(ta,tna,tnpatch,tnperim,tnndf))
            log     = rbind(log,log1)

            ### plot result
            # raster::plot(res,main=sprintf("tna=%0.0f hfr = %0.02f \n cf = %0.02f",tna,hfr1,cf1),cex.main=0.8)
            # plot(sps1$geometry,add=TRUE,border="darkgreen")

            ### write to file ## Actually it takes much less space to save the 0 as 0 and saving them as -99 or NA takes up 1.5 times the space
            writeRaster(terra::rast(res),paste(respath,simname,
                                       sprintf("%s_nc_%0.02f_hfr_%0.02f_cf_%d_rep.tif",simname,hfr1,cf1,re1),sep="/"),overwrite=TRUE)

            ### write to shapefile
            # if(file.exists(paste(respath,simname,
            #                      sprintf("%s_nc_%0.02f_hfr_%0.0f_npatch_%0.02f_cf.shp",simname,hfr,npatch,cf),sep="/"))){
            #   file.remove(paste(respath,simname,
            #                     sprintf("%s_nc_%0.02f_hfr_%0.0f_npatch_%0.02f_cf.shp",simname,hfr,npatch,cf),sep="/")) }
            #   st_write(res$sps, paste(respath,simname,
            #                      sprintf("%s_nc_%0.02f_hfr_%0.0f_npatch_%0.02f_cf.shp",simname,hfr,npatch,cf),sep="/"), driver="ESRI Shapefile")  # create to a shapefile
            } # reps
            } # cfs
      ### write log for each nc and hfr
      log = log[-1,]
      write.csv(log,paste(respath,simname,sprintf("%s_%d_nc_%0.02f_hfr.csv",simname,nc1,hfr1),sep="/"),row.names=FALSE)
    } # hfrs
  } # ncs

} # function end


########
#' @name reallands_gen
#' @title Generate real-world landscapes with habitat distribution from an input raster file
#' @description wrapper function sampling landscapes of a certain size and for a given land cover value
#' @param simname character, name of simulation run
#' @param rast land cover raster file containing codes corresponding to the "Code" column in lcdf.
#' @param ncs numeric value defining number of cells in x and y dimension of the square raster to be extractedx
#' @param samplesize numeric value indicating the number of desired sample landscapes
#' @param overlap numeric value indicating the allowed overlap among sampled landscapes (1=no overlap,if <1 overlap will be bigger, if >1 distance among sampled landscapes will be larger)
#' @param respath character string indicating file path where results should be saved
#' @param habnames character string defining land cover names to be used to define habitat pixels in rast, habnames must be contained in the lccol of the lcdf dataframe.
#' @param lcdf data frame containing a "Code" column describing codes contained in rast; and a lccol with associated habitat names
#' @param lccol character string defining column name of lcdf where habitat names are contained
#' @return a raster and/or a shapefile containing a landscape with habitat patches with above defined characteristics
#' @export
reallands_gen = function(simname="test",rast=NULL,ncs=100,samplesize=50,overlap=0.5,respath=NULL,
                         lcdf=NULL,habnames=NULL,lccol="hab") {

  ## create repository to save landscapes of a certain simulation run
  if(!dir.exists(paste(respath,simname,sep="/"))){dir.create(paste(respath,simname,sep="/"))}

  if(length(lcdf)==0) {
    valls = unique(raster::values(rast))
    valls = valls[is.finite(valls)&valls>0]
    habnames = paste(names(rast),letters[1:length(valls)],sep="_")
    lcdf = data.frame(hab=habnames,Code=valls)
  }

  ## derive extent
  ext  = raster::extent(rast)
  ## define increment related to desired cellsize (half - overlapping windows in this case...)
  radius  = 0.5*ncs*mean(res(rast),na.rm=TRUE)
  incr    = overlap*2*radius

  ## make centers with increment
  centers = mkcentersi(ext=ext,incr=incr)
  ## add crs
  sf::st_crs(centers) = sf::st_crs(rast)

  ## plot for check
  raster::plot(rast,main=names(rast))
  #plot(centers$geometry,add=TRUE)

  ## remove edge centers - not necessary anymore with the mkcentersi
  ## remove centers in NA pixels
  vals   = raster::extract(rast,centers)

  ## JO 3.8.2023: dont know why I did this with 0 here?:
  # identify and remove centers in NA pixels
  #torem  = which(vals==0|is.na(vals))
  torem  = which(is.na(vals))
  if(length(torem)>0) {
  outcent= centers[torem,]
  centers= centers[-torem,]
  #plot(centers$geometry,add=TRUE,col="red")
  ## remove also the centers that are too close to the outcents
  disto   = st_distance(centers,outcent)
  ## remove centers too close to the outcents
  torem   = which(apply(disto,MARGIN=1,FUN=function(x){any(x<=2*incr)}))
  ## make final centers
  centers = centers[-torem,]
  #plot(centers$geometry,add=TRUE,col="yellow")
  } # length torem

  ## sample samplesize from centers (make sure its not larger than the number of centers)
  samplesize = min(c(samplesize,length(centers$cid)),na.rm=TRUE)
  print(sprintf("samplesize = %d",samplesize))

  ## sample from centers
  centers = centers[sample(1:length(centers$cid),size=samplesize),]
  plot(centers$geometry,add=TRUE,col="black",bg="darkorange",pch=21)

  ### extract windows
  for(cid in centers$cid){
    # cid = centers$cid[1]
    foc = mkbuffer(centers=centers[centers[["cid"]]==cid,],radius=radius,form="square")
    foc = sf::st_sf(geometry=foc$geometry)
    sf::st_crs(foc) = sf::st_crs(rast)
    plot(foc$geometry,add=TRUE,border="blue",col=adjustcolor("blue",alpha.f=0.3))
    ## crop raster accordingly (masking not needed here, since shape is always square)
    rast1   = raster::crop(rast,raster::extent(foc))
    ## make different habitat rasters and save to file!

    for(hab in habnames) {
      #hab = habnames[1]
      ### set all non-habitat pixels to NA!
      res     =  rast1
      lccodes = lcdf[lcdf[[lccol]]==hab,"Code"]
      if(length(lccodes)==1) {
      res[res!=lccodes]  = NA
      res[res==lccodes]  = 1
      } else if(length(lccodes)>1){
      res[!res%in%lccodes]= NA
      res[res%in%lccodes] = 1
      } else if(length(lccodes)==0){
      print(sprintf("%s not found..",hab))
      res[]= NA
      }
      #raster::plot(res,col=c(NA,"darkgreen"),main=sprintf("cid: %d lc: %s",cid,hab),legend=FALSE)

      evl = any(raster::values(res)==1)

      if(evl>0 & !is.na(evl)){
        ### run log calculations and extract raster
        ### identify raster clumps, i.e. generate IDs!
        res1   = landscapemetrics::get_patches(res, directions=4, class=1)[[1]][[1]]
        print(sprintf("%s %0.0f patches..",hab,max(values(res1),na.rm=TRUE)))

        ### if it is a raster, make a shapefile out of it, to be able to exclude patches below minsize..
        sps1 = stars::st_as_stars(res1) %>% sf::st_as_sf(merge = TRUE)
        sps1 = sf::st_make_valid(sps1) # JO 31.12.2024: ,reason=TRUE
        names(sps1)[[1]] = "pid"

        ### log statistics
        ## total habitat area (slightly different if calculated from shapefile or raster...)
        ta     = NA ## there is no initial ta
        ## total new area
        #tna   = sum(res[res==1])
        tna    = sum(as.numeric(sf::st_area(sps1))) # new area
        ## total nr of new patches
        tnpatch = length(sps1$geometry)
        ## total habitat perimeter
        tnperim = sum(perimf(sps1))
        ## total nearest neighbour distance
        tnndf   = sum(NNDf(sf::st_distance(sps1)))

        ### plot result
        plot(sps1$geometry,add=TRUE,col="yellow")
        #raster::plot(res,main=sprintf("cid = %d hab = %s \n tna=%0.0f",cid,hab,tna),cex.main=0.8)


      } else {
        print(sprintf("%s %0.0f patches..",hab,0))

        ta      = NA
        tna     = 0
        tnpatch = 0
        tnperim = 0
        tnndf   = 0
      }

      ## add calculations to log statistics
      log    = data.frame(simname=simname,ncs=mean(dim(res)[1:2]),cid=cid,hab=hab
                          ,metric=c("ta","tna","tnpatch","tnperim","tnndf"),value=c(ta,tna,tnpatch,tnperim,tnndf))
      ### write log
      write.csv(log,paste(respath,simname,sprintf("%s_ncs_%d_cid_%d_hab_%s.csv",simname,ncs,cid,hab),sep="/"),row.names=FALSE)

      ### write to file ## Actually it takes much less space to save the 0 as 0 and saving them as -99 or NA takes up 1.5 times the space
      writeRaster(terra::rast(res),paste(respath,simname,
                            sprintf("%s_ncs_%d_cid_%d_hab_%s.tif",simname,ncs,cid,hab),sep="/"),overwrite=TRUE)

      ### write to shapefile
      # if(file.exists(paste(respath,simname,
      #                      sprintf("%s_ncs_%0.02f_hfr_%0.0f_npatch_%0.02f_cf.shp",simname,hfr,npatch,cf),sep="/"))){
      #   file.remove(paste(respath,simname,
      #                     sprintf("%s_ncs_%0.02f_hfr_%0.0f_npatch_%0.02f_cf.shp",simname,hfr,npatch,cf),sep="/")) }
      #   st_write(res$sps, paste(respath,simname,
      #                      sprintf("%s_ncs_%0.02f_hfr_%0.0f_npatch_%0.02f_cf.shp",simname,hfr,npatch,cf),sep="/"), driver="ESRI Shapefile")  # create to a shapefile


    } # hab

  } # cid

} # function end


########
#' @name summplots
#' @title Generate overview plots of simulation log data
#' @description plot relationships among metrics in simulated landscapes, grouped by groups
#' @param simname character, name of simulation run
#' @param makepdf logical, should a pdf of the plots be created? (default=FALSE)
#' @param rdf a results data frame, containing the column names "metric", "value", "simname" as well as names specified in the "groups" variable
#' @param metrics a character string defining metrics in the "metric" column that should be summarized
#' @param groups a character string indicating column names in "rdf" that serve as grouping variables when summarizing the data
#' @param cols  color pallette to be used when plotting (default = topo.colors)
#' @param simpcirc logical - are data coming from the simpcirc simulation (simpcirc_gen function output)? Default = FALSE
#' @param respath character string indicating file path where results should be saved
#' @return Plots depicting relationships among metrics, grouped by "groups" and saved as a pdf file in respath (if makepdf=TRUE)
#' @export
summplots= function(simname="test",respath=NULL,makepdf=FALSE,rdf=NULL, metrics=NULL, groups=NULL,cols=topo.colors, simpcirc=FALSE) {


  ## make pdf or not
  if(makepdf == TRUE) {
    ## create repository to save landscapes of a certain simulation run
    if(!dir.exists(paste(respath,simname,sep="/"))){dir.create(paste(respath,simname,sep="/"))}
    pdf(paste(respath,simname,sprintf("%s_LogStats.pdf",simname),sep="/"),width = 5, height=5)
  }
  par(mar=c(5,5,3,3))

  if(simpcirc == TRUE){
    ## first check the closeness of simulation to desires
    ## desired area vs. final area
    rdf1 = rdf[rdf$metric=="ta",]
    rdf2 = rdf[rdf$metric=="tna",]
    plot(rdf2$value~rdf1$value,bg=cols(max(rdf1$npatch))[rdf1$npatch], pch=21,
         xlab="total desired \n patch area (m2)",ylab="total realized \n patch area (m2)")
    legend("topleft",title="npatch",legend=round(seq(1,max(rdf1$npatch),length.out=5),0),fill=cols(max(rdf1$npatch))[round(seq(1,max(rdf1$npatch),length.out=5),0)],bty="n")

    ## desired nr of patches vs. actual nr. of patches
    rdf1 = rdf[rdf$metric=="tnpatch",]
    plot(rdf1$value~rdf1$npatch,bg=cols(length(levels(rdf1$cf)))[rdf1$cf], pch=21,
         xlab="total desired \n number of patches",ylab="total realized \n number of patches")
    abline(0,1)
    legend("topleft",title="cf",legend=levels(rdf1$cf),fill=cols(length(levels(rdf1$cf)))[1:length(levels(rdf1$cf))],bty="n")
  } ## simpcirc specific summary

  ## first check desired area vs. final area
  for(m in metrics) {
    print(m)
    rdf1 = rdf[rdf$metric==m,]
    for(mm in metrics){
      print(mm)
      rdf2 = rdf[rdf$metric==mm,]
      ## merge dataframes to be sure
      rdf3 = merge(rdf1,rdf2,by=colnames(rdf1)[!colnames(rdf1)%in%c("metric","value")],all=TRUE)
      if(any(is.finite(unique(rdf3$value.x))) & any(is.finite(unique(rdf3$value.y)))) {
        for(g in groups){
          print(g)
          ## plot mm as a function of m
          plot(rdf3$value.y~rdf3$value.x,ylab=unique(rdf3$metric.y),xlab=unique(rdf3$metric.x),las=2,bg=NA,col=NA,pch=21,
               main=sprintf("group = %s",as.character(g)),cex.main=1)
          if(is.numeric(rdf3[[g]])){
            points(rdf3$value.y~rdf3$value.x,bg=adjustcolor(cols(max(rdf3[[g]],na.rm=TRUE))[rdf3[[g]]],alpha.f=0.5),pch=21,cex=1.5)
            legend("topleft",legend=round(seq(1,max(rdf3[[g]],na.rm=TRUE),length.out=5),0), bty="n",
                   fill=cols(max(rdf3[[g]],na.rm=TRUE))[round(seq(1,max(rdf3[[g]],na.rm=TRUE),length.out=5),0)],title=g)
          } else if (is.factor(rdf3[[g]])){
            points(rdf3$value.y~rdf3$value.x,bg=adjustcolor(cols(length(levels(rdf3[[g]])))[rdf3[[g]]],alpha.f=0.5),pch=21,cex=1.5)
            legend("topleft",legend=levels(rdf3[[g]]),fill=cols(length(levels(rdf3[[g]]))),title=g, bty="n")
          }

        } # groups
      } # finite values
    } # mm metrics
  } # m metrics

  if(makepdf == TRUE) {
    dev.off() }

  print("finished..!")

} # function end


##########################################################
##### project:       Multispecies Connectivity Modelling
##### author:        Jacqueline Oehri (JO), jacqueline.oehri@gmail.com
##########################################################
