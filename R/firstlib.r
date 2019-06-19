.onAttach <- function (lib, pkg) {
  ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"),"Version")
  ver <- as.character(ver)
  packageStartupMessage("\nESS ",ver," loaded.\n", 
                      # "Please cite as:\n   ",format(citation("ESS"), style = "text") ,
                         domain = NULL,  appendLF = TRUE)
}

.onLoad <- function(...) {
}
