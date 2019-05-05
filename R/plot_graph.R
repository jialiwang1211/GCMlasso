#' Plot graphical model
#'
#' Plot the conditional dependence between variables on the latent scale.
#' The assocations between variables are computed from the precision matrices
#' from the MCMC, and the percentage of edges to be included in the graph
#' is determined by the \code{edge_perc} value.
#'
#' @param GCMlasso_obj \code{GCMlasso} object.
#'
#' @param var variables in the graph.
#'
#' @param edge_perc cut-off quantile value of edge inclusion.
#'
#' @param seed a random integer.
#'
#' @examples plot_graph(GCMlasso_obj,var=1:15,edge_perc=0.65)
#'
#' @export
#' @import igraph
plot_graph<-function(GCMlasso_obj=CMlasso_obj,var,edge_perc,seed=1){
  set.seed(seed)
  data<-GCMlasso_obj$data_ordered
  Omega<-GCMlasso_obj$Omega.st[var,var,]
  numsamp=dim(Omega)[3]
  Prec.mcmc<-Omega[,,1:numsamp]
  Prec<-apply(Prec.mcmc,c(1,2),mean)
  quantile_val<-quantile(abs(Prec),edge_perc)
  Prec[abs(Prec)<quantile_val]=0

  colnames(Prec)<-rownames(Prec)<-names(data)[var]
  network=graph_from_adjacency_matrix(Prec, weighted=T, mode="lower", diag=F)

  V(network)$size <- 5
  E(network)[E(network)$weight>0]$color <- "blue"
  E(network)[E(network)$weight< -0]$color <- "red"

  plot.igraph(network, layout=layout_with_gem)
}
