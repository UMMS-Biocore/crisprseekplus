#' getLoadingMsg
#'
#' @note \code{getLoadingMsg}
#' @return loading msg
#' @examples
#'     x <- getLoadingMsg()
#' @export
#'
getLoadingMsg <- function() {
  imgsrc <- "www/images/loading.gif"
  a <- list(
  tags$head(tags$style(type = "text/css", "
                         #loadmessage {
                         position: fixed;
                         top: 0px;
                         left: 200px;
                         width: 70%;
                         height: 100;
                         padding: 5px 0px 5px 0px;
                         text-align: center;
                         font-weight: bold;
                         font-size: 100%;
                         color: #000000;
                         opacity: 0.8;
                         background-color: #FFFFFFF;
                         z-index: 100;
}")),
    conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
    tags$div("Please wait! Loading...", id = "loadmessage",
                tags$img(src = imgsrc))),
    
    conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",
                tags$div("Analysis complete! Click 'Download Output'", 
                id = "loadmessage" )))
    }

#' getLogo
#'
#'
#' @note \code{getLogo}
#' @return return logo
#' @examples
#'     x <- getLogo()
#' @export
#'
getLogo <- function(){
    imgsrc <- "www/images/logo.png"
    a<-list(img(src=imgsrc, align = "right"))
}

