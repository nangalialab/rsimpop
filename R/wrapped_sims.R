### This file contains wrapped calls to the simulator that serve as examples
#' @title 
#' Runs a single compartment/single driver scenario
#' @description 
#' run_selection_sim models the evolution of a population of cells with a selective advantage given to a randomly
#' chosen cell, once the population size reaches the specified equilibrium size.
#' @param initial_division_rate - Rate of symmetric cell division during development
#' @param final_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param target_pop_size - Size of target population
#' @param nyears_driver_acquisition - When driver is acquired.
#' @param nyears - Total number of years to run the simulation
#' @param fitness - relative fitness advantage.  Cells carrying this divide at a rate=(1+fitness)*baserate
#' @param minprop - Minimum aberrant cell fraction at nyears for driver to be regarded as have "taken"
#' @param maxtry  - Maximum number of attempts to introduce a driver
#' @param max_driver_count - If positive simulation will stop when number drivers matches this.  This evaluated at daily time snaps.
#' @return simpop object.
#' @export
#' @examples
#' selsim=run_selection_sim(0.05,1/365,target_pop_size = 5e4,nyears = 50,fitness=0.3)
#' @details 
#' The selective advantage is introduced as a driver event. The driver is subject to stochastic 
#' extinction and so is repeatedly dropped in until it "takes".The prevailing state when the driver is introduced is
#' saved and reinstated with each attempted introduction.
run_selection_sim=function (initial_division_rate = 0.1, final_division_rate = 1/365, 
                            target_pop_size = 1e+05, nyears_driver_acquisition = 15, 
                            nyears = 40, fitness = 0.2, minprop = 0.001, mindriver = 1, 
                            maxtry = 40, max_driver_count=-1) 
{
  cfg = getDefaultConfig(target_pop_size, rate = initial_division_rate, 
                         ndriver = 1, basefit = fitness)
  params = list(n_sim_days = nyears_driver_acquisition * 365, 
                b_stop_at_pop_size = 1, b_stop_if_empty = 0,max_driver_count=max_driver_count)
  growthphase = sim_pop(NULL, params = params, cfg)
  #browser()
  ndivkeep = 0
  mdivkeep = 0
  gdivkeep = 0
  tree0 = get_tree_from_simpop(growthphase)
  if (growthphase$status == 0) {
    tree1 = tree0
    params[["n_sim_days"]] = nyears * 365
    params[["b_stop_if_empty"]] = 1
    dc = 0
    tries = 0
    tree1_tmp = tree1
    while (dc < max(minprop * target_pop_size, mindriver)) {
      if (tries >= maxtry) {
        return(NULL)
      }
      cat("No driver found: tries=", tries, "\n")
      tries = tries + 1
      tree1_tmp = addDriverEvent(tree1, tree1$cfg, 1, fitness = fitness)
      print(tree1_tmp$cfg$info)
      params[["b_stop_at_pop_size"]] = 1
      adult2 = sim_pop(tree1_tmp, params = params, tree1_tmp$cfg)
      tree2 = get_tree_from_simpop(adult2)
      params[["b_stop_at_pop_size"]] = 0
      cfg = tree2$cfg
      cfg$compartment$rate[2] = final_division_rate
      cfg$compartment$popsize[2] = target_pop_size
      lastsim = sim_pop(tree2, params = params, cfg)
      dc = rsimpop:::getSingleDriverCount(lastsim)
    }
  }
  else {
    gdivkeep = mean(nodeHeights(tree0)[which(tree0$edge[, 
                                                        2] <= length(tree0$tip.label)), 2])
    cfg$compartment$rate[2] = final_division_rate
    cfg$compartment$popsize[2] = target_pop_size
    params[["n_sim_days"]] = nyears_driver_acquisition * 
      365
    params[["b_stop_at_pop_size"]] = 0
    adult1 = sim_pop(tree0, params = params, cfg)
    tree1 = get_tree_from_simpop(adult1)
    params[["n_sim_days"]] = nyears * 365
    params[["b_stop_if_empty"]] = 1
    dc = 0
    tries = 0
    tree1_tmp = tree1
    while (dc < max(minprop * target_pop_size, mindriver)) {
      if (tries >= maxtry) 
        return(NULL)
      cat("No driver found: tries=", tries, "\n")
      tries = tries + 1
      tree1_tmp = addDriverEvent(tree1, tree1$cfg, 1, fitness = fitness)
      if (TRUE) {
        ndivkeep = nodeHeights(tree1_tmp)[which(tree1_tmp$edge[, 
                                                               2] == tree1_tmp$events$node[3]), 2]
        mdivkeep = mean(nodeHeights(tree1_tmp)[which(tree1_tmp$edge[, 
                                                                    2] <= length(tree1_tmp$tip.label)), 2])
      }
      print(tree1_tmp$cfg$info)
      adult2 = sim_pop(tree1_tmp, params = params, tree1_tmp$cfg)
      dc = getSingleDriverCount(adult2)
    }
    lastsim = adult2
  }
  fulltree = get_tree_from_simpop(lastsim)
  fulltree$tries = tries
  fulltree$ndivkeep = ndivkeep
  fulltree$mdivkeep = mdivkeep
  fulltree$gdivkeep = gdivkeep
  return(fulltree)
}

#' Runs a multiple driver scenario
#'
#' The user specifies the rate at which drivers are introduced and the distribution of the selection coefficients are drawn from.
#' Note that unlike in \code{\link{run_selection_sim}} the drivers are allowed to stochastically die out.
#' The drivers arrive at the specified rate with the gap between successive introductions exponentially distributed.
#' @param initial_division_rate - Rate of symmetric cell division during development
#' @param final_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param target_pop_size - Size of target population
#' @param drivers_per_year - The expected number of drivers per year.
#' @param nyears - Total number of years to run the simulation
#' @param fitness_gen - Function that returns a single value for the selection coefficient or vector with 200,000 selection coefficients.
#' @param bForceReseed - Whether to reseed.
#' @param offset - time offset for choosing seed.
#' @param user_trajectory - data.frame specifying target population size and division rate change points.  The timestamp (ts) column is in days, division_rate is in divs per day: ts, target_pop_size,division_rate, compartment
#' @param simpop - existing rsimpop object.  Created if null - non-null not currently tested.
#' @param drivers_per_cell_per_day - Rate of driver acquisition, per cell, per day. This is the preferred specification for the frequency of driver acquisition.
#' @return simpop object.
#' @export
#' @examples
#' fitnessGen=function(){
#'  trials=rexp(100,rate=30)
#'  idx=which(trials>0.05)
#'  if(length(idx)==0){
#'   fitnessGen()
#'  }else{
#'   trials[idx[1]]
#'  }
#' }
#' dps=run_driver_process_sim(0.1,1/(2*190),target_pop_size = 1e5,nyears = 10,fitness=fitnessGen,drivers_per_year = 5)
#' trajectory=data.frame(ts=365*(1:80),target_pop_size=5e4+100*(1:80),division_rate=1/(2*190),compartment=1)
#' trajectory$target_pop_size[15:20]=2*trajectory$target_pop_size[15:20]
#' trajectory$target_pop_size[21:25]=0.2*trajectory$target_pop_size[21:25]
#' dps=run_driver_process_sim(initial_division_rate=0.1,user_trajectory=trajectory, target_pop_size = 1e5,nyears = 80,fitnessGen=fitnessGen,drivers_per_cell_per_day = 1/(365*1e5))
run_driver_process_sim=function( 
                                 initial_division_rate=0.1,
                                 final_division_rate=1/365,
                                 target_pop_size=1e5,
                                 drivers_per_year=0.1,
                                 nyears=40,
                                 fitnessGen=fitness_vector,
                                 bForceReseed=FALSE,
                                 offset=0,
                                 user_trajectory=NULL,
                                 simpop=NULL,
                                 drivers_per_cell_per_day=NA
                                 )
{
  
  if(bForceReseed){
    delay=(offset/10 +1)
    Sys.sleep(delay)
    cat("delay=",delay,"\n")
    initSimPop(-1,bForce = TRUE)
  }
  ## setup trajectory
  if(!is.null(user_trajectory) & (class(user_trajectory) %in% c("data.frame", "data.table"))){
    message("cell population will be modelled based on provided trajectory")
    if(ncol(user_trajectory) != 4){
      stop("trajectory data frame should have 4 columns named as 'ts', 'target_pop_size', 'division_rate', 'compartment'")
    }
    if(any(user_trajectory$division_rate>1)){
      stop("Supplied division rate too high; should be less than 1.0")
    }
    if(any(user_trajectory$target_pop_size>20e6)){
      stop("Supplied population too high; should be less than 20e6")
    }
    user_trajectory=user_trajectory[order(user_trajectory$ts,user_trajectory$compartment),]
    if(length(unique(user_trajectory$compartment))>1){
      stop("Multiple compartments in trajectory not currently supported.")
    }
  }else{
    nsd=nyears*365
    user_trajectory = data.frame(ts=c(nsd, nsd), target_pop_size=c(target_pop_size, target_pop_size),
                                 division_rate=c(final_division_rate, final_division_rate), compartment=c(1,1))
  }
  nsd=user_trajectory$ts[nrow(user_trajectory)]
  
  ## setup drivers per cell per day
  if(is.na(drivers_per_cell_per_day)){
    dpcpd=(drivers_per_year/365)/user_trajectory$target_pop_size[1]
    warning("Specification of driver frequency is preferably done using drivers_per_cell_per_day rather than drivers_per_year")
  }else{
    dpcpd=drivers_per_cell_per_day
  }
  maxdrivers=dpcpd*get_max_popsize_from_trajectory(simpop,user_trajectory)*nsd
  ##Add 10 sd as insurance  (assume poisson variance)
  maxdrivers=max(1000,ceiling(maxdrivers+10*sqrt(maxdrivers)))
  
  
  ## setup fitness coefficients
  fitness_gen=fitnessGen ## Keeping backward compatibility with parameter neame
  drivers_fitness=fitness_gen
  if(class(drivers_fitness)=="function"){
    ##backward compatibility
    #warning("Deprecation warning:drivers_fitness should now be a vector of fitness coefficients")
    test=drivers_fitness()
    if(length(test)==1){
      fitness_gen=sapply(1:maxdrivers,function(dummy) drivers_fitness())
    }else{
      stop("bad drivers_fitness supplied")
    }
  }else{
    stop("Please supply function for generating a single selection coefficient")
    #if(length(drivers_fitness)<2e5){
    #driversFitness=c(driversFitness,rep(0,1e5-length(driversFitness)))
    #stop("Please provide at least 200,000 driver fitness values")
    #}
  }
  
  
  
  if(is.null(simpop)){
    
    params=list(n_sim_days=user_trajectory$ts[1],  ##This is the longest that the simulation will continue
                b_stop_at_pop_size=0,
                b_stop_if_empty=0,
                driver_rate_per_cell_per_day=dpcpd)
    cfg=getDefaultConfig(user_trajectory$target_pop_size[1],rate=user_trajectory$division_rate[1],
                         ndriver=1,basefit = 0)
    cfg2=getDefaultConfig(user_trajectory$target_pop_size[1],rate=initial_division_rate,
                          ndriver=1,basefit = 0)
    params2=params
    params2$b_stop_at_pop_size=1
    tree0=sim_pop(NULL,params=params2,cfg2,driversFitness = fitness_gen)
    #if(tree0$status==1){
    maxd=max(1,tree0$maxDriverID)
    params$n_sim_days=max(tree0$timestamp)
    tree0=sim_pop(tree0,params=params,cfg,b_verbose = FALSE, trajectory=user_trajectory,#[2:nrow(user_trajectory),],
                  driversFitness=fitness_gen[-(1:maxd)])
    #}
    
  }else{
    #message("simpop object alredy exist! use continue_driver_process_sim() function instead")
    # I should call continue_driver_process_sim() function instead
    #continue_driver_process_sim(simpop,nyears,fitness_gen, user_trajectory,drivers_per_year)
    params=list(n_sim_days=max(simpop$timestamp),  ##This is the longest that the simulation will continue
                b_stop_at_pop_size=0,
                b_stop_if_empty=0,
                driver_rate_per_cell_per_day=dpcpd)
    cfg=simpop$cfg
    user_trajectory=user_trajectory %>% filter(ts>max(simpop$timestamp) & compartment %in% cfg$compartment$val)
    if(dim(user_trajectory)[1]==0){
      warning("Supplied user trajectory finishes before simulation")
      return(simpop)
    }
    
    utf=user_trajectory %>% group_by(compartment) %>% summarise(rate=division_rate[1],popsize=target_pop_size[1])
    idx=match(utf$compartment,cfg$compartment$val)
    if(length(idx)>0){
      cfg$compartment$rate[idx]=utf$rate
      cfg$compartment$popsize[idx]=utf$popsize
    }
    #browser()
    tree0=sim_pop(simpop,params=params,cfg,b_verbose = FALSE, trajectory=user_trajectory,#[2:nrow(user_trajectory),],
                  driversFitness=fitness_gen)
  }
  #  }else{
  #    
  
  #    nsd=nyears*365
  #    user_trajectory = data.frame(ts=c(nsd, nsd), target_pop_size=c(target_pop_size, target_pop_size),
  #                                 division_rate=c(initial_division_rate, final_division_rate), compartment=c(1,1))
  #    
  #    cfg=getDefaultConfig(user_trajectory$target_pop_size[1],rate=user_trajectory$division_rate[1],
  #                         ndriver=1,basefit = 0)
  #    params=list(n_sim_days=user_trajectory$ts[1],
  #                b_stop_at_pop_size=1,
  #                b_stop_if_empty=0,
  #                driver_rate_per_cell_per_day=dpcpd)
  #    drivers_fitness=fitness_gen
  #    # fitnessD=fitness_gen
  #    tree0=sim_pop(NULL,params=params,cfg,driversFitness = fitness_gen[1:100000])
  #    if(tree0$status==1){
  #      ##  Ended at initial pop size..  
  #      params$b_stop_at_pop_size=0
  #      cfg$compartment$rate[2]=final_division_rate
  #      tree0=sim_pop(tree0,params=params,cfg, trajectory=user_trajectory[2:nrow(user_trajectory),],
  #                  driversFitness = fitness_gen[-(1:100000)])
  #    }else{
  #      browser()
  #    }
  #  }
  
  tree0
  return(tree0)
}


#' Restarts a multiple driver scenario
#'
#' The user specifies the rate at which drivers are introduced and the distribution the selection coefficients are drawn from.
#' Note that unlike in \code{\link{run_selection_sim}} the drivers are allowed to stochastically die out.
#' The drivers at the specified rate with the gap between successive introductions exponentially distributed.
#' @param insim - Result of previous sim
#' @param final_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param target_pop_size - Size of target population
#' @param drivers_per_year - The expected number of drivers per year
#' @param nyears - Total number of years to run the simulation
#' @param bForceReseed - Whether to reseed.
#' @param offset - time offset for choosing seed.
#' @return simpop object.
#' @export
#' @examples
#' fitnessGen=function(){
#'  trials=rexp(100,rate=30)
#'  idx=which(trials>0.05)
#'  if(length(idx)==0){
#'   fitnessGen()
#'  }else{
#'   trials[idx[1]]
#'  }
#' }
#' dps=run_driver_process_sim(0.1,1/(2*190),target_pop_size = 1e5,nyears = 30,fitness=fitnessGen,drivers_per_year = 5)
#' ## Do stuff - subsample, plot trees etc and the restart to simulate the next 10 years:
#' dps2=continue_driver_process_sim(dps,40,fitnessGen = fitnessGen)
#' ## Do stuff and simulation the next 10 years
#' dps2=continue_driver_process_sim(dps2,50,fitnessGen = fitnessGen)
continue_driver_process_sim=function(insim,nyears,fitness_gen, user_trajectory=NULL,drivers_per_year){
  
  if(!is.null(user_trajectory)){
    maxt=max(insim$timestamp)
    idx=which(user_trajectory$ts>maxt)
    if(length(idx)>0){
      idx=idx[1]
      user_trajectory = user_trajectory[idx:nrow(user_trajectory),]
      dpcpd=(drivers_per_year/365)/user_trajectory$target_pop_size[1]
      params=list(n_sim_days=user_trajectory$ts[1],  ##This is the longest that the simulation will continue
                  b_stop_at_pop_size=1,
                  b_stop_if_empty=0,
                  driver_rate_per_cell_per_day=dpcpd)
      cfg=getDefaultConfig(user_trajectory$target_pop_size[1],rate=user_trajectory$division_rate[1],
                           ndriver=1,basefit = 0)
      
    }else{
      idx=1
    }
    if( insim$maxDriverID>0){
      fitness_gen=drivers_fitness[-(1:insim$maxDriverID)]
    }
    sp=insim
    tree0 = sim_pop(get_tree_from_simpop(sp), params=params,cfg,b_verbose = FALSE, driversFitness=fitness_gen,
                    trajectory_year=user_trajectory$ts[2:nrow(user_trajectory)], trajectory_pop = user_trajectory$target_pop_size[2:nrow(user_trajectory)],
                    trajectory_div_rate=user_trajectory$division_rate[2:nrow(user_trajectory)])
  }else{
    tree0=insim
    params=insim$params
    params$maxt=NULL
    params$n_sim_days=nyears*365
    if(length(which(tree0$cfg$info$fitness==0 & tree0$cfg$info$population>0))>2){
      browser()
    }
    
    if(max(tree0$timestamp)>=nyears*365){
      tree0
      return(tree0)
    }
    #driversFitness=fitnessGen
    if(class(fitness_gen)=="function"){
      ##backward compatibility
      warning("Deprecation warning:driversFitness should now be a vector of fitness coefficients")
      test=driversFitness()
      if(length(test)==1){
        fitnessGen=sapply(1:200000,function(dummy) driversFitness())
      }else{
        stop("bad driversFitness supplied")
      }
    }else{
      if(length(fitness_gen)<2e5){
        #driversFitness=c(driversFitness,rep(0,1e5-length(driversFitness)))
        stop("Please provide at least 200,000 driver fitness values")
      }
    }
    fitnessD=fitness_gen
    
    #while(max(tree0$timestamp)<nyears*365)
    #{
      adult=sim_pop(tree0,params=params,tree0$cfg,driversFitness = fitnessD)
    #  if( adult$maxDriverID>0){
    #    fitnessD=fitnessD[-(1:adult$maxDriverID)] #fitnessGen=driversFitness[-(1:adult$maxDriverID)]
    #  }
    #  if(adult$maxDriverID>100000){  
    #    stop("too many iterations!")
    #  }
      tree0=adult
      # adult=sim_pop(tree0,params=params,tree0$cfg,b_verbose = FALSE)
      # tree0=adult
    #}
  }
  
  tree0
  return(tree0)
}



#' Runs a single compartment/single driver scenario with transient selection
#'
#' @description 
#' run_transient_selection models the evolution of a population of cells with a transient selective advantage given to a randomly
#' chosen cell, once the population size reaches the specified equilibrium size.
#' @param initial_division_rate - Rate of symmetric cell division during development
#' @param final_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param target_pop_size - Size of target population
#' @param nyears_driver_acquisition - When driver is acquired.
#' @param nyears_transient_end - When selective advantage is set to 0.
#' @param nyears - Total number of years to run the simulation
#' @param fitness - relative fitness advantage.  Cells carrying this divide at a rate=(1+fitness)*baserate
#' @param minprop - Minimum aberrant cell fraction at nyears for driver to be regarded as have "taken"
#' @param maxtry  - Maximum number of attempts to introduce a driver
#' @details 
#' This function is similar to the run_selection_sim() function. However, the driver is switch off after the time specified in the nyears_transient_end argument.  
#' @return simpop object.
#' @export
#' @examples
#' tselsim=run_transient_selection(0.05,1/365,target_pop_size = 5e4,nyears_driver_acquisition=15,
#' nyears_transient_end=30,nyears=50,fitness=0.5)
run_transient_selection=function( initial_division_rate,
                                  final_division_rate,
                                  target_pop_size=1e5,
                                  nyears_driver_acquisition=15,
                                  nyears_transient_end=30,
                                  nyears=40,
                                  fitness=0.2, ## relative fitness
                                  minprop=0.05)
{
  selsim=run_selection_sim(initial_division_rate,
                           final_division_rate,
                           target_pop_size,
                           nyears_driver_acquisition,
                           nyears_transient_end,
                           fitness,minprop)
  ##switch off
  cfg=selsim$cfg
  cfg$info$fitness[3]=0
  params=selsim$params
  params[["n_sim_days"]]=nyears*365
  params[["maxt"]]=NULL
  final=sim_pop(selsim,params=params,cfg)
  # final=combine_simpops(selsim,final)
  return(final)
}

#' Runs a simple neutral simulation
#'
#' @param initial_division_rate - Rate of symmetric cell division during development
#' @param final_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param target_pop_size - Size of target population
#' @param nyears - Total number of years to run the simulation
#' @return simpop object.
#' @export
#' @examples
#' neutsim=run_neutral_sim(0.05,1/365,target_pop_size = 5e4,nyears=80)
run_neutral_sim=function( initial_division_rate,
                          final_division_rate,
                          target_pop_size=1e5,
                          nyears=40)
{
  cfg=getDefaultConfig(target_pop_size,rate=initial_division_rate,ndriver=1,basefit = 0)
  params=list(n_sim_days=nyears*365,
              b_stop_at_pop_size=1,
              b_stop_if_empty=0
  )
  growthphase=sim_pop(NULL,params=params,cfg)
  cfg$compartment$rate[2]=final_division_rate
  cfg$compartment$popsize[2]=target_pop_size
  params[["b_stop_at_pop_size"]]=0
  adult1=sim_pop(growthphase,params=params,cfg)
  return(adult1)  #return(combine_simpops(growthphase,adult1))
}


add_driver=function(simpop,params,acceptance_threshold=0.05){
  tree1=get_tree_from_simpop(simpop)
  params[["b_stop_if_empty"]]=1
  tree1=simpop
  dc=0
  tries=0
  tree1_tmp=tree1
  while(dc/target_pop_size<0.05){
    cat("No driver found: tries=",tries,"\n")
    idx=sample(length(tree1$tip.label)-1,1)+1
    celltype=rep(NA,length(tree1$tip.label))
    celltype[idx]=-1
    tree1_tmp=assign_celltype(tree1,celltype,tree1$cfg)
    simpop2=sim_pop(tree1_tmp,params=params,tree1_tmp$cfg)
    
    dc=getSingleDriverCount(simpop2)
    #print(adult2$cfg$info)
    tries=tries+1
    if(tries>max_tries){
      stop("Unable to add driver.. Too many attempts")
    }
  }
  simpop=combine_simpops(simpop,simpop2)
  simpop
}

#' Runs a simple single compartment neutral simulation with a specified trajectory
#'
#' @param simpop - Rate of symmetric cell division during development
#' @param initial_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param trajectory - data.frame - with fields ts(timestamp in days),target_pop_size,division_rate
#' @param nyears - Total number of years to run the simulation
#' @return simpop object.
#' @export
#' @examples
#' trajectory=data.frame(ts=365*(1:80),target_pop_size=5e4+100*(1:80),division_rate=1/(2*190))
#' trajectory$target_pop_size[5:10]=2*trajectory$target_pop_size[5:10]
#' trajectory$target_pop_size[11:15]=0.2*trajectory$target_pop_size[11:15]
#' sp=run_neutral_trajectory(NULL,0.5,trajectory)
#' plot(sp)
run_neutral_trajectory=function(simpop,initial_division_rate,trajectory){
  if(any(trajectory$division_rate>1)){
    stop("Supplied division rate too high; should be less than 1.0")
  }
  if(any(trajectory$target_pop_size>20e6)){
    stop("Supplied population too high; should be less than 20e6")
  }
  if(any(trajectory$target_pop_size<2)){
    stop("Supplied population too low; should be more than 1")
  }
  if(any(is.na(trajectory$target_pop_size)) || any(is.na(trajectory$ts))|| any(is.na(trajectory$division_rate))){
    stop("trajectory should contain no missing data")
  }
  params=list(n_sim_days=trajectory$ts[1],##This is the longest that the simulation will continue
              b_stop_at_pop_size=1,
              b_stop_if_empty=0
  )
  #browser()
  if(is.null(simpop)){
    cfg=getDefaultConfig(target_pop_size = trajectory$target_pop_size[1],rate = initial_division_rate,basefit = 0,ndriver = 1)
    
    
    sp=sim_pop(NULL,params=params,cfg,b_verbose = FALSE)
    idx=1
  }else{
    maxt=max(simpop$timestamp)
    idx=which(trajectory$ts>maxt)
    if(length(idx)>0){
      idx=idx[1]
    }else{
      idx=1
    }
    
    sp=simpop
  }
  
  params[["b_stop_at_pop_size"]]=0
  for(i in idx:(length(trajectory$ts)-1)){
    cfg=list(compartment=data.frame(val=c(0,1),rate=c(-1,trajectory$division_rate[i]),popsize=c(1,trajectory$target_pop_size[i])),
             info=sp$cfg$info)
    cfg=getDefaultConfig(target_pop_size = trajectory$target_pop_size[i],rate = trajectory$division_rate[i],basefit = 0,ndriver = 1)
    params[["n_sim_days"]]=max(sp$timestamp+1,trajectory$ts[i+1]) #occasionally the current simulation finishes after the next target date - skip here.
    spx=sim_pop(get_tree_from_simpop(sp),params=params,cfg,b_verbose = FALSE)
    sp=spx  #sp=combine_simpops(sp,spx)
  }
  sp
}


run_2compartment_model=function(){
  ##Initialise LT-HSC
  #
  
}

simple_continuous_selection=function(){
  
}



#' Converts annual growth selection coefficient to per division selection coefficient
#' @param S input annual selection coefficient
#' @param divrate division rate parameter specified in divisions per day.
#' @return s : per division selection coefficient to be supplied to rsimpop
#' @export
#'
convert_annualS_to_s=function(S,divrate){
  log(1 + S)/(365 * divrate)
}


getSingleDriverCount=function(sim){
  if(length(which(sim$events$driverid>0))>0){
    sim$cfg$info$population[3]
  }else{
    0
  }
}
