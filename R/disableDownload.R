#' disableDownload
#'
#' Enable or disable to download button depending on if analysis is complete
#'
#' @note \code{disableDownload}
#' @return the download button either enabled or disabled;
#' @param input, disable the download button
#'
#' @examples
#'     x<- disableDownload()
#'
#' @importFrom shinyjs toggleState
#' @export
#' 
disableDownload <- function(input = NULL){ 
  if (is.null(input)) return(NULL)
  toggleState(
    id = "downloadData",
    condition =  input > 0
  )
}