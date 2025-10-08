### Configuration
################################################################################
### R Options
nodeInfo = system(paste0("scontrol show node ", Sys.info()["nodename"]), intern = TRUE)
cpus = as.integer(unlist(strsplit(unlist(strsplit(nodeInfo[grep("CfgTRES", nodeInfo)], "cpu="))[2], ","))[1])
memInKB = as.integer(unlist(strsplit(unlist(strsplit(nodeInfo[grep("CfgTRES", nodeInfo)], "mem="))[2], "M,"))[1]) * 1024^2

options(stringsAsFactors=FALSE,
        dplyr.summarise.inform=FALSE, 
        future.globals.maxSize=min(memInKB, 20 * 1024^3),
        mc.cores=min(cpus,1),
        future.fork.enable=TRUE, future.plan="multicore",
        future.rng.onMisuse="ignore")

### Fuctions
# Git directory and files to source must be done first, then all helper functions can be sourced
git_files_to_source = c("functions_analysis.R",
                        "functions_annotation.R",
                        "functions_biomart.R", 
                        "functions_degs.R", 
                        "functions_io.R", 
                        "functions_plotting.R", 
                        "functions_util.R")
git_files_to_source = file.path(param$path_to_git, "R", git_files_to_source)
file_exists = purrr::map_lgl(git_files_to_source, file.exists)
if (any(!file_exists)) stop(paste("The following files could not be found:", paste(git_files_to_source[!file_exists], collapse=", "), ". Please check the git directory at '", param$path_to_git, "'.!"))
invisible(purrr::map(git_files_to_source, source))

### Debugging mode
# 'default_debugging' for default, 'terminal_debugger' for debugging without X11, 'print_traceback' for non-interactive sessions 
param$debugging_mode = "default_debugging"
switch (param$debugging_mode, 
        default_debugging=on_error_default_debugging(), 
        terminal_debugger=on_error_start_terminal_debugger(),
        print_traceback=on_error_just_print_traceback(),
        on_error_default_debugging())

