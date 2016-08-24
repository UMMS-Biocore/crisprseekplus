#' fileInputFunc
#'
#' If input file is empty, use sample file 
#'
#' @note \code{fileInputFunc}
#' @return If no file is uploaded, use sample file
#' @param input, enter correct files
#' @param sampleFile, sampleFile
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