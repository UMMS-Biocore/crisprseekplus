#' fileInputFunc
#'
#' If input file is empty, use sample file 
#'
#' @note \code{fileInputFunc}
#' @return If no file is uploaded, use sample file
#' @param fileInputFunc, enter correct files
#'
#' @examples
#'     x<- fileInputFunc()
#'
#'
#' @export
#

fileInputFunc <- function(input, sampleFile) {
  if(is.null(input)) {
    sampleFile
  }
  else {
    input
  }
}