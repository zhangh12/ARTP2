makeVector <- function (x) 
{
    d <- dim(x)
    if (is.null(d)) 
        return(x)
    nn <- NULL
    if (d[1] == 1) {
        nn <- colnames(x)
    }
    else if (d[2] == 1) {
        nn <- rownames(x)
    }
    dim(x) <- NULL
    if (!is.null(nn)) 
        names(x) <- nn
    if ((!is.vector(x)) && (is.list(x))) 
        x <- unlist(x)
    x
}
