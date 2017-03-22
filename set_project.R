## ExtractR for ECAD climate data

# ExtractR
#		- R/
# 		- data/
# 		- doc/
# 		- figs/
# 		- output/
#		- analysis/

#======================================

# R 

project_name <- "WessexLULCchange"

if (!dir.exists(normalizePath(file.path(Sys.getenv("HOME"), "OneDrive - Natural Environment Research Council",project_name)))) {
	dir.create(normalizePath(file.path(Sys.getenv("HOME"), "OneDrive - Natural Environment Research Council",project_name)))
}

setwd(normalizePath(file.path(Sys.getenv("HOME"), "OneDrive - Natural Environment Research Council",project_name)))

dir.create("R")
dir.create("data")
dir.create("doc")
dir.create("figs")
dir.create("output")	
dir.create("analysis")