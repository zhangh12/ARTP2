removeWhiteSpace <- function (str, leading = 1, trailing = 1) {

    if ((leading) && (trailing)) {
        ret <- gsub("^\\s+|\\s+$", "", str, perl = TRUE)
    }
    else if (leading) {
        ret <- gsub("^\\s+", "", str, perl = TRUE)
    }
    else if (trailing) {
        ret <- gsub("\\s+$", "", str, perl = TRUE)
    }
    else {
        ret <- str
    }
    ret

}
