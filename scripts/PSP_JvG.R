###Background
#Here are functions related to phosphosite plus e.g. loading in data


directory_test <- function(){
  
  setwd(getSrcDirectory(function(x) {x}))
  print(getwd())
}



###Function: General function that loads psp data
#directory: directory where psp data is e.g. "data/psp". If NULL, directory is figured out based on location of this script (not recommended if you are not using this on JvG's computer)
load_psp_data <- function(filename,
                          directory = NULL){
  current_wd <- getwd()
  
  #Load to default location if directory is null
  if (is.null(directory)){
    #Set wd to script
    setwd(getSrcDirectory(function(x) {x}))
    psp_data <- read.delim(paste("..\\..\\data/biol_databases/phosphosite_plus", filename, sep = "/"))
  } else {
    psp_data <- read.delim(paste(directory, filename, sep = "/"))
  }
  setwd(current_wd)
  return(psp_data)
}



###Function: Uniprot site to site group ID
#uniprot_sites: vector of sites to convert of the form e.g. P81122_S604
#psp_phos_data: Dataset of phosphosites to use. If NULL, will load in itself
#restrict_to_phos: Only use phosphorylation modifications (probably redundant)
uniprotsite_to_sitegroupid <- function(uniprot_sites,
                                       psp_phos_data = NULL,
                                       restrict_to_phos = TRUE){
  if (is.null(psp_phos_data)){
    psp_phos_data <- load_psp_data("PSP__phos_dataset_20211007")
  }
  if (restrict_to_phos){
    psp_phos_data <- psp_phos_data[grep("p$", psp_phos_data$MOD_RSD), ]
  }
  #Make uniprot_site in psp_phos_data
  psp_phos_data$uniprot_site <- paste(psp_phos_data$ACC_ID, 
                                      gsub("-.+", "", psp_phos_data$MOD_RSD), 
                                      sep = "_")
  site_group_ids <- lapply(uniprot_sites, function(x) psp_phos_data$SITE_GRP_ID[which(psp_phos_data$uniprot_site == x)])
  names(site_group_ids) <- uniprot_sites
  site_group_ids[vapply(site_group_ids, function(x) {length(x) == 0}, logical(1))] <- NA
  return(site_group_ids)
}
  
  
  
  












