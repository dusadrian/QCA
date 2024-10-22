`runGUI` <- function(
    x
) {
    
    if (missing(x)) {
        x <- system.file("gui", package = "QCA")
    }
    
    ## Future feature: local settings
    
    # if (file.exists(file.path(getwd(), "noedit"))) {
        
        # settings <- lapply(strsplit(readLines(file.path(x, "noedit")), split="="), admisc::trimstr)
        # opts <- unlist(lapply(settings, "[[", 1))
        # vals <- unlist(lapply(settings, "[[", 2))
        
        ## see which other options might be added
        
        # sink(file.path(getwd(), "noedit"))
        # for (i in seq(length(opts))) {
        #     cat(opts[i], " = ", vals[i], "\n", sep="")
        # }
        # sink()
    # }
    # else {
        # sink(file.path(getwd(), "noedit"))
        
        ## add other options if needed
        
        # sink()
    # }
    
    Sys.setenv(userwd = getwd())
    
    runApp(x, launch.browser = TRUE)
    
}
