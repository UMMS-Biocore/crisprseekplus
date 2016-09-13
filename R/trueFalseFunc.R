#' trueFalseFunc
#'
#' If radio button input == 1, the function returns true
#' and if the radio button value is 2, the function returns false
#'
#' @note \code{trueFalseFunc}
#' @return true or false depending on input
#' @param input, user inputs
#'
#' @examples
#'     x<- trueFalseFunc()
#'
#'
#' @export
#
trueFalseFunc <- function(input = NULL) {
if (is.null(input)) return(NULL)
if(input == 1) {
        TRUE
    }
else if(input == 2) {
        FALSE
    }
}
