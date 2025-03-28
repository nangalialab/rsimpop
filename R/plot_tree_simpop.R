#' Plots simpop events on a the tree
#'
#' @param tree phylo/simpop
#' @param legpos legend position
#' @param fmode  0=no selection coeff annotation;1=per driver selection coeff;2=composite driver selection coeff
#' @param events.keep  Output from previous plot_tree_events call where it is desired to keep the same colour for the nodes.
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> . Parameters that are passed to plot_tree
#'
#' @return simpop object
#' @export
#' @examples
#' cfg=getDefaultConfig(target_pop_size =5e4,rate=0.1)
#' params=list(n_sim_days=2*365,b_stop_at_pop_size=1,b_stop_if_empty=0)
#' growthphase=sim_pop(NULL,params=params,cfg)
#' stree=get_subsampled_tree(growthphase,50)
#' plot_tree_events(stree,fmode=1)
plot_tree_events=function (tree, legpos = "topleft", fmode=0,events.keep=NULL,use.desc=FALSE,...)
{
  if (is.null(tree$events)) {
    stop("tree has no events")
  }
  events = tree$events
  events = events[order(events$ts), ]
  events = events %>% left_join(tree$cfg$compartment,by=c("value"="val"))
  if(use.desc){
    events=events %>% mutate(value=desc)
  }
  tree = plot_tree(tree, ...)
  idx.child = match(tree$edge[, 2], tree$edge[, 1])
  N = dim(tree$edge)[1]
  TT = max(tree$timestamp)
  duration = ifelse(is.na(idx.child), TT - tree$tBirth, tree$tBirth[idx.child] -
                      tree$tBirth)
  events$fitness=tree$cfg$drivers$fitness[match(events$driver,tree$cfg$drivers$driver)]
  events$fitness[is.na(events$fitness)]=0
  if(fmode==2){
    warning("WARNING: fmode=2 under development gives wrong results for multiple drivers on single branch")
    s=sapply(1:length(events$driverid),function(i){
      if(events$driverid[i]>0){

        parents=intersect(rsimpop:::get_parents(events$node[i],tree$edge),events$node)
        idx=which(events$node %in% parents)
        #sum(events$fitness[idx])
        sum(events$fitness[idx])
      }else{
        0
      }
    })
    events$fitness=s
  }
  if(fmode>0){
    events$key=sprintf("%s:%s:%3.2f", events$value,
                     events$driverid,events$fitness)
  }else{
    events$key=sprintf("%s:%s", events$value,
                       events$driverid)
  }
  keepkey=character(0)
  if(!is.null(events.keep)){
    if(fmode>0){
      keepkey=with(events.keep,sprintf("%s:%s:%3.2f", value,driverid,fitness))
    }else{
      keepkey=with(events.keep,sprintf("%s:%s", value,driverid,fitness))
    }
    keepkey=sort(unique(keepkey))
  }
  otherkeys=setdiff(events$key,keepkey)
  df = data.frame(uval = c(keepkey,unique(sort(otherkeys))), stringsAsFactors = FALSE)
  cols=c(RColorBrewer::brewer.pal(9,"Set1"))
  df$col=rep(cols,26)[1:length(df$uval)]
  df$pch=rep(c(19,17,15,25,setdiff(0:25,c(19,17,15,25))),each=length(cols))[1:length(df$uval)]
  df=df %>% filter(uval %in% events$key)
  
  events = events %>% left_join(df, by = c(key = "uval"))

  events$idx = match(events$node, tree$edge[, 2])
  fracs = lapply(1:N, function(x) c())
  cols = lapply(1:N, function(x) c())

  for (i in 1:length(events$node)) {
    thisnode = events$node[i]
    idx = events$idx[i]
    info = rsimpop:::get_edge_info(NULL, tree, thisnode)
    frac = (events$ts[i] - tree$tBirth[events$idx[i]])/(duration[events$idx[i]]+1e-6)##Fix zero duration
    points(x = info$x, y = info$yt - frac * (info$yt - info$yb),
           pch = events$pch[i], col = events$col[i], cex = 1)
    fracs[[idx]] = c(fracs[[idx]], frac)
    cols[[idx]] = c(cols[[idx]], events$col[i])
    kids = get_all_node_children(thisnode, tree)
    for (k in match(kids, tree$edge[, 2])) {
      fracs[[k]] = 0
      cols[[k]] = events$col[i]
    }
  }
  for (j in 1:N) {
    info = rsimpop:::get_edge_info(NULL, tree, tree$edge[j, 2])
    ff = c(fracs[[j]][-1], 1)
    segments(x0 = info$x, y1 = info$yt - fracs[[j]] * (info$yt -
                                                         info$yb), y0 = info$yt - ff * (info$yt - info$yb),
             col = cols[[j]], lwd = 1)
  }
  if (!is.null(legpos)) {
    legend(legpos, legend = df$uval, col = df$col, pch = df$pch,
           cex = 1)
  }
  tree$events=events
  tree
}

get_tree_overlap=function(tree,selsim){
  st=tree
  parents=match(selsim$edge[,1],selsim$edge[,2])
  child1=match(selsim$edge[,2],selsim$edge[,1])
  schild1=match(st$edge[,2],st$edge[,1])
  #st=get_subsampled_tree(selsim,40);
  st2=get_elapsed_time_tree(st)
  overlaps=lapply(1:dim(st$edge)[1],function(i){
    if(is.na(schild1[i])){
      ##terminal branch.. nothing we can do with this!
      ##TODO Fix up get_subsampled_tree to store parent tree IDs.
      NULL
    }else{
      thisBirth=st$tBirth[i]
      parent=which(selsim$tBirth==st$tBirth[schild1[i]])[1]
      if(length(parent)==0){
        stop("Cant find in parent tree!")
      }
      ts=selsim$tBirth[parent]
      tslist=data.frame(ts=rep(NA,1000),dcount=rep(-1,1000),node=st2$edge[i,2])#rep(NA,100)
      k=0

      while(!is.na(parent) && ts>=thisBirth){
        if(k>0){
          tslist$ts[k]=ts
          tslist$dcount[k]=selsim$dcount[parent]
        }
        k=k+1
        parent=parents[parent]
        ts=selsim$tBirth[parent]

      }
      tslist[which(!is.na(tslist$ts)),]
    }
  })
  st2$overlaps=overlaps
  st2$countdf=do.call("rbind",overlaps)
  st2$countdf=with(st2,countdf[order(countdf$node,countdf$ts),])
  st2$totalpop=sum(selsim$cfg$info$population)-1
  st2
}


##TODO:Redo using plot_tree_annots framework
plot_tree_overlaps=function(tree,legpos="topleft",scale=1,b.add.count=TRUE,...){
  events=tree$countdf
  if(is.null(tree$events)){
    stop("tree has no events")
  }
  tree=plot_tree(tree,...)
  idx.child=match(tree$edge[,2],tree$edge[,1])
  N=dim(tree$edge)[1]
  TT=max(tree$timestamp)
  duration=ifelse(is.na(idx.child),TT-tree$tBirth,tree$tBirth[idx.child]-tree$tBirth)
  df=data.frame(uval=unique(sort(sprintf("%s:%s",events$node,events$ts))),stringsAsFactors = FALSE)
  #df$col=c(RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(12,"Set3"))[1:length(df$uval)]
  allcol=c(RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(12,"Set3"))
  df$col="black"
  events$idx=match(events$node,tree$edge[,2])
  events$col="red"
  fracs=lapply(1:N,function(x) c())
  cols =lapply(1:N,function(x) c())
  acfs =lapply(1:N,function(x) c())
  dcounts=lapply(1:N,function(x) c())


  for(i in 1:N){
    thisBirth=tree$tBirth[i]
    overlap=tree$overlaps[[i]]

    if(!is.null(overlap)){
      #if(dim(overlap)[1]==1){
      #  fracs[[i]]=0
      #  acfs[[i]]=overlap$dcount/tree$totalpop
      #}else{
      overlap=overlap[order(overlap$ts),]
      fracs[[i]]=(overlap$ts-thisBirth)/duration[i]
      #if(length(fracs[i])>0){
      #  fracs[[i]]=c(0,diff(fracs[[i]]))
      #}
      acfs[[i]]=overlap$dcount/tree$totalpop
      dcounts[[i]]=overlap$dcount
      cols[[i]]=allcol[1:dim(overlap)[1]]
      #}
    }
  }


  ##browser()
  for(j in 1:N){
    #if(j==36){
    #  browser()
    #}
    info=get_edge_info(NULL,tree,tree$edge[j,2])
    ff=fracs[[j]]
    ff2=c(ff[-1],1)
    ##cat(info$yt-fracs[[j]]*(info$yt-info$yb),info$yt-ff*(info$yt-info$yb))
    ##if(length(ff)>1){

    if(!is.null(ff)){
      yt=info$yt-ff*(info$yt-info$yb)
      yb=info$yt-ff2*(info$yt-info$yb)
    rect(xleft=info$x-scale*acfs[[j]],xright=info$x+scale*acfs[[j]],ytop=yt,
             ybottom=yb,col=cols[[j]],border=NA)
      if(b.add.count){
        text(info$x,0.5*(yt+yb),sprintf("%d",dcounts[[j]]),col="black",cex=0.4)
      }
    }
    ##}

  }
  ##Add the Events..
  for(i in 1:dim(tree$events)[1]){

    if(tree$events$driverid[i]>0){
      #browser()
      j=which(tree$edge[,2]==tree$events$node[i])
      info=get_edge_info(NULL,tree,tree$edge[j,2])
      frac=(tree$events$ts[i]-tree$tBirth[j])/duration[j]
      points(x=info$x,y=info$yt-frac*(info$yt-info$yb),pch=4,col="black",cex=2)
    }
  }


  if(!is.null(legpos)){
    ##legend(legpos,legend=df$uval,col=df$col,pch=19,cex=1)
  }
  tree
}

#' Adds a vector "dcount" to a phylo object giving the number of descendent tips for each child node
#'
#' @param tree id of specified node
#' @param tree phylo
#' @return phylo with dcount vector in parallel to edge.length
#' @export
add_descendent_counts_to_tree=function(tree){
  N=length(tree$tip.label)
  idx=which(tree$edge[,2]<=N)
  dcount=rep(0,length(tree$edge.length))
  dcount[idx]=1
  ##First we index the edges...
  parents=match(tree$edge[,1],tree$edge[,2])
  for(i in idx){
    parent=parents[i]
    while(!is.na(parent)){
      dcount[parent]=dcount[parent]+1
      parent=parents[parent]
    }
  }
  tree$dcount=dcount
  tree
}


annotate_tree_with_pop_coalescences=function(subtree,selsim){
  st=subtree
  parents=match(selsim$edge[,1],selsim$edge[,2])
  #st=get_subsampled_tree(selsim,40);
  st2=get_elapsed_time_tree(st);tree=plot_tree_events(st2)
  parent=which(selsim$tBirth==st$tBirth[which(st$edge[,1]==st$events$node[3])[1]])[1]
  while(!is.na(parent)){
    ts=selsim$tBirth[parent]
    cat(selsim$edge[parent,],selsim$dcount[parent],ts,"\n")
    abline(h=tree$top-ts)
    parent=parents[parent]
  }
}

