### Rmarkdown configuration
################################################################################

### Create output directories
if (!file.exists(file.path(param$path_out, "figures"))) dir.create(file.path(param$path_out, "figures"), recursive=TRUE, showWarnings=FALSE)
if (!file.exists(file.path(param$path_out, "data"))) dir.create(file.path(param$path_out, "data"), recursive=TRUE, showWarnings=FALSE)


### R Options
options(citation_format="pandoc", 
        knitr.table.format="html",
        kableExtra_view_html=TRUE)


### Knitr default options
knitr::opts_chunk$set(echo=TRUE,                     # output code
                      cache=FALSE,                   # do not cache results
                      message=TRUE,                  # show messages
                      warning=TRUE,                  # show warnings
                      tidy=FALSE,                    # do not auto-tidy-up code
                      fig.width=10,                  # default fig width in inches
                      class.source='fold-hide',      # by default collapse code blocks
                      dev=c('png', 'svg', 'tiff'),          # create figures in png, tiff, and svg; the first device (png) will be used for HTML output
                      dev.args=list(png = list(type="cairo"),     # png: use cairo - works on cluster
                                    tiff = list(type="cairo", compression = 'zip')),  # tiff: use cairo - works on cluster, transparent background, compression type "zip" is ‘deflate (Adobe-style)’ 
                      dpi=300,                       # figure resolution
                      fig.retina=2,                 # retina multiplier
                      fig.path=paste0(file.path(param$path_out, "figures"), "/")  # Path for figures in png and pdf format (trailing "/" is needed)
)

### Set output hooks
knitr::knit_hooks$set(message=format_message, warning=format_warning)

### Required libraries
library(magrittr)

### Set the bibliographic environment
# Clear previous bibliography
knitcitations::cleanbib()

# If the DOI publication servers cannot be reached, there will be no citations, knitcitations will not write a references.bib file and pandoc will stop. This makes sure that there is always at least one citation.
rcug_ref = knitcitations::read.bibtex(file.path(param$path_to_git,"assets","rcug_references.bib"))
invisible(knitcitations::citep(rcug_ref))

### Figure heights
# Single figure, landscape
fig_standard_height = 4
# Two plots alongside (e.g. umaps)
fig_standard2_height = 5
# Three plots alongside (e.g. umaps)
fig_standard3_height = 4
# Four plots alongside (e.g. umaps)
fig_standard4_height = 2.5
# Four plots 2x2 (e.g. umaps)
fig_patchwork4_height = fig_standard2_height * 2
# Six plots 2x3 (e.g. umaps)
fig_patchwork6_height = fig_standard3_height * 2
# Eight plots 4x2 (e.g. umaps)
fig_patchwork8_height = fig_standard4_height * 2
# Twelve plots 4x3 (e.g. umaps)
fig_patchwork12_height = fig_standard4_height * 3
# Sixteen plots 4x4 (e.g. umaps)
fig_patchwork16_height = fig_standard4_height * 4




