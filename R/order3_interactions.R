#' @title Look at shared taxa between samples of different types
#' @author Kara J.M. Taylor
#' @description Basically a specific format for a distance matrix between two input mattrices
#' @param X Input community matrix for first sample type
#' @param Y Input community matrix for second sample type
#' @param rare Logical or numeric; Rarefy the data first? If numeric, the number to rarefy to.
#' @param commList Output from `as.bipartite.net`
#' @param method Distance metric; See `vegan::vegdist` for options.
#' @param binary Logical; Whether the metric should use the binary version or not.
#' @returns A distance matrix formatted for downstream analysis.
#' @export


as.bipartite.net <- function(X, Y, rare=F){

  rows <- deparse(substitute(X))
  cols <- deparse(substitute(Y))

  if (is.numeric(rare)){
    X <- X[rowSums(X)>rare,]
    Y <- Y[rowSums(Y)>rare,]
  }

  Y <- Y[which(row.names(Y)%in%row.names(X)),]
  X <- X[which(row.names(X)%in%row.names(Y)),]

  Y <- Y[order(row.names(Y)),]
  X <- X[order(row.names(X)),]

  if (!identical(row.names(X), row.names(Y))){
    stop("Samples are not ordered the same way. Check input `row.names`.")}

  if (is.numeric(rare)){
    X <- vegan::rrarefy(X, rare)
    Y <- vegan::rrarefy(Y, rare)
  }
  bind <- rbind(X,Y)
  bind[bind>0] <- 1

  if (identical(X, Y)){
    noZeros <- which(colSums(bind)<=1)
  } else {noZeros <- which(colSums(bind)<=1)}
  X <- X[,-noZeros]
  Y <- Y[,-noZeros]

  row.names(X) <- paste0(row.names(X),rows)
  row.names(Y) <- paste0(row.names(Y),cols)

  return(list(rows=X, cols=Y))
}

#' @export
cross.comparison <- function(commList, method=NULL, binary=NULL){

  commX <- commList$rows
  commY <- commList$cols

  if (!identical(dim(commX),dim(commY))){
    stop("Input matrices should have the same dimensions")}

  if (!is.null(method)&is.null(binary)){
    stop("If using a distance metric, should the data be binary?")}

  if (is.null(method)){

    commX[commX>0] <- 1
    commY[commY>0] <- 1

    pull_overlap <- function(rowX, matY){

      each_rowY <- function(rowY, rowX){
        Y <- rowX+rowY
        Z <- which(Y==2)
        return(length(Z))
      }

      counts <- apply(matY, MARGIN=1, FUN=each_rowY, rowX=rowX)
    }

    count_overlap <- apply(commX, MARGIN=1, FUN=pull_overlap, matY=commY)

    reorder <- match(names(sort(rowSums(count_overlap), decreasing=T)), row.names(count_overlap))
    count_overlap <- count_overlap[reorder,reorder]

    return(t(count_overlap))

  } else {

    bind <- rbind(commX, commY)
    bind <- bind[,-which(apply(bind,2,function(x){length(which(x==0))})>=(nrow(commX)+nrow(commY)-1))]

    veg <- vegan::vegdist(rbind(commX, commY), method=method, diag = T, upper=T, binary=binary)
    veg <- as.matrix(veg)

    #veg <- 1-veg # similarity matrix

    overlap <- veg[row.names(commX),row.names(commY)]

    return(overlap)
  }

}
