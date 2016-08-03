
#' startcrisprseekplus
#'
#' Starts the crisprseekplus to be able to run.
#'
#' @note \code{startcrisprseekplus}
#' @return the app
#'
#' @examples
#'     startcrisprseekplus()
#'
#' @export
#'

startcrisprseekplus <- function(){
if (interactive()) {
    #the upload file size limit is 1GB
    options( shiny.maxRequestSize = 1000 * 1024 ^ 2)
    
    addResourcePath(prefix = "www", directoryPath =
                      system.file("extdata", "www", 
                                  package = "crisprseekplus"))
    
    environment(cspServer) <- environment()
    app <- shinyApp( ui = shinyUI(cspUI),
                     server = shinyServer(cspServer))
    runApp(app)
  }
}


