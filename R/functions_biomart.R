#' Wrapper around the biomaRt::useEnsembl function to cope with unavailable Ensembl mirrors. Tries different Ensembl mirrors and returns a mart object with the mirror that works.
#' 
#' @param biomart A biomaRt database name. Possible database names can be retrieved with the function listEnsembl().
#' @param dataset Dataset you want to use. Possible dataset names can be retrieved with the function listDatasets(mart_obj).
#' @param mirror Specify an Ensembl mirror to connect to. The valid options here are 'www', 'uswest', 'useast', 'asia'. If no mirror is specified, then the first mirror that works will be used. Will be ignored if a host is specified.
#' @param version Ensembl version to connect to when wanting to connect to an archived Ensembl version.
#' @param host Host to connect to. Only needs to be specified if different from www.ensembl.org. 
#' @return A biomaRt object.
GetBiomaRt = function(biomart, dataset, mirror=NULL, version=NULL, host=NULL) {
  
  # Which mirrors to test
  if (is.null(mirror)) {
    mirrors_to_test = c("www", "uswest", "useast", "asia")
  } else {
    mirrors_to_test = c(mirror)
  }
  
  mart_obj = NULL
  if(is.null(host)) {
    # Test and if a mirror is not available, check the next one
    for(m in mirrors_to_test) {
      mart_obj = tryCatch({
        biomaRt::useEnsembl(biomart=biomart, dataset=dataset, mirror=m, version=version)
      },
      error=function(cond) {
        return(NULL)
      })
      
      if(!is.null(mart_obj)) break
    }
    
  } else {
    # Use specific host
    mart_obj = tryCatch({
      biomaRt::useEnsembl(biomart=biomart, dataset=dataset, host=host, version=version)
    },
    error=function(cond) {
      return(NULL)
    })
  }
  
  return(mart_obj)
}


#' Returns the mirror of a biomaRt object.
#' 
#' @param mart_obj A biomaRt object obtained by GetBiomaRt or useEnsembl name.
#' @return The mirror of the biomaRt object. Can be 'www', 'uswest', 'useast' or 'asia'.
GetBiomaRtMirror = function(mart_obj) {
  mirrors_to_test = c("uswest", "useast", "asia")
  mirror = "www"
  
  for(m in mirrors_to_test){
    if(grepl(pattern=m, x=mart_obj@host)){
      mirror = m
      break
    }
  }
  
  return(mirror)
}


#' Translate between Ensembl id, gene symbol, seurat-compatible unique rowname, and Entrez
#' 
#' @param annot_ensembl A gene annotation table with Ensembl ids.
#' @param annot_main A vector with database names ("ensembl_gene_id", "external_gene_name", "entrezgene_accession").
#' @return List and tables with diverse ID translations (ensembl_to_seurat_rowname, seurat_rowname_to_ensembl, seurat_rowname_to_entrez, annot_ensembl).

TranslateIDs = function(annot_ensembl, annot_main=NULL) {
  
  # Which mirrors to test
  if (is.null(annot_main)) {
    param_annot_main = c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
  } else {
    param_annot_main = c(annot_main)
  }
  
  # Create translation tables
  ensembl = param_annot_main["ensembl"]
  symbol = param_annot_main["symbol"]
  entrez = param_annot_main["entrez"]
  
  # Ensembl id to gene symbol (new)
  ensembl_to_symbol = unique(annot_ensembl[, c(ensembl, symbol)])
  ensembl_to_symbol[, symbol] = ifelse(nchar(ensembl_to_symbol[, symbol]) == 0, NA, ensembl_to_symbol[, symbol])
  ensembl_to_symbol = setNames(ensembl_to_symbol[, symbol], ensembl_to_symbol[, ensembl])
  
  # Ensembl id to seurat-compatible unique rowname (new)
  ensembl_to_seurat_rowname = unique(annot_ensembl[, c(ensembl, symbol)])
  ensembl_to_seurat_rowname[, symbol] = ifelse(nchar(ensembl_to_seurat_rowname[, symbol]) == 0, ensembl_to_seurat_rowname[, ensembl], ensembl_to_seurat_rowname[, symbol])
  ensembl_to_seurat_rowname[, symbol] = make.unique(gsub(pattern="_", replacement="-", x=ensembl_to_seurat_rowname[, symbol], fixed=TRUE))
  ensembl_to_seurat_rowname = setNames(ensembl_to_seurat_rowname[, symbol], ensembl_to_seurat_rowname[, ensembl])
  
  # Seurat-compatible unique rowname to ensembl id
  seurat_rowname_to_ensembl = setNames(names(ensembl_to_seurat_rowname), ensembl_to_seurat_rowname)
  
  # Ensembl to Entrez
  ensembl_to_entrez = unique(annot_ensembl[, c(ensembl, entrez)])
  ensembl_to_entrez[, entrez] = ifelse(nchar(ensembl_to_entrez[, entrez]) == 0, NA, ensembl_to_entrez[, entrez])
  ensembl_to_entrez = split(ensembl_to_entrez[, entrez], ensembl_to_entrez[, ensembl])
  
  # Seurat-compatible unique rowname to Entrez
  seurat_rowname_to_ensembl_match = match(seurat_rowname_to_ensembl, names(ensembl_to_entrez))
  names(seurat_rowname_to_ensembl_match) = names(seurat_rowname_to_ensembl)
  seurat_rowname_to_entrez = purrr::map(seurat_rowname_to_ensembl_match, function(i) {unname(ensembl_to_entrez[[i]])})
  
  # Entrez IDs is duplicating Ensembl IDs in annot_ensembl
  # Therefore, we remove Entrez IDs from the annotation table, after generating all required translation tables
  # Set rownames of annotation table to Ensembl identifiers
  annot_ensembl = as.data.frame(unique(annot_ensembl[, -match(entrez, colnames(annot_ensembl))]))
  rownames(annot_ensembl) = annot_ensembl[, ensembl]

  IDs_out <- list(ensembl_to_seurat_rowname, seurat_rowname_to_ensembl, seurat_rowname_to_entrez, annot_ensembl) 
  return(IDs_out)
}


