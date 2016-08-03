#' disableDownload
#'
#' Enable or disable to download button depending on if analysis is complete
#'
#' @note \code{disableDownload}
#' @return the download button either enabled or disabled;
#' @param disableDownload, disable the download button
#'
#' @examples
#'     x<- disableDownload()
#'
#' @importFrom shiny toggleState
#'
#' @export
#' 
disableDownload <- function(input){
  toggleState(
    id = "downloadData",
    condition =  input > 0
  )
}