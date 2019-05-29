returnDirectlyConnectedNodes <- function(node, linkList){
  tmp_link <- linkList[linkList$source %in% node,]
  tmp_link2 <- linkList[linkList$target %in% node,]
  
  tmp_nodes <- union(tmp_link$target, tmp_link2$source)
  tmp_nodes <- union(tmp_nodes, node)
  return(tmp_nodes)
}



returnAllConnectedNodes <- function(node, linkList){
  for (i in 1:5) {
    node <- returnDirectlyConnectedNodes(node, linkList)
  }
  return(node)
}


drawNetworkGraph <- function(linkList, Nodes, overlay){
  
  rownames(Nodes) <- 1:nrow(Nodes)
  
  ref <- 1:nrow(Nodes)
  names(ref) <- Nodes$nodes
  linkList$source1 <- ref[as.character(linkList$source)] - 1
  linkList$target1 <- ref[as.character(linkList$target)] - 1
  
  linkList$Value <- 1
  #linkList$Colour <- c("#CD6155", "#566573")[as.numeric(linkList[,3] > 0) + 1]
  
  a <- forceNetwork(Links = linkList, Nodes = Nodes, zoom = T,opacityNoHover = 0.5, 
                    Source = "source1", Target = "target1", arrows = T,
                    NodeID = "nodes", Value ="Value" , #linkColour = linkList$Colour,
                    Group = overlay, opacity = 0.9)
  
  return(a)
  
}

#' Draw the Network
#'
#' This function generate a force-directed network based on the link list and nodes information. 
#' @param celltag.obj A CellTag object with link list and nodes filled
#' @param tag Which tags would you like to plot?
#' @param overlay What information would you like to overlay with the network? This should be one of the column names of the node information.
#' @return A CellTag object with the attribute (network) modified.
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' drawSubnet(bam.test.obj, "CellTagV1_2", "Cluster")
#' 
drawSubnet <- function(celltag.obj, tag, overlay){
  # e.g. tag; "celltag2.1_698"
  # e.g. color: "cluster" or "tag" or "SuperClone"
  Nodes <- celltag.obj@nodes
  linkList <- celltag.obj@network.link.list
  
  no <- returnAllConnectedNodes(tag, linkList)
  sub_link <- linkList[(linkList$source %in% no) | (linkList$target %in% no),]
  sub_Nodes <- Nodes[Nodes$nodes %in% no ,]
  
  a <- drawNetworkGraph(sub_link, sub_Nodes, overlay)
  
  celltag.obj@network <- a
  return(celltag.obj)
}

