# These functions are designed to expand GO annotations obtained from e.g Uniprot in order to obtain *all* ancestor terms as well.
# The rationale is that e.g tRNA binding proteins may not be annotated as RNA binding protiens which is obviously unintuitive!

suppressMessages(library(GO.db))

##################################################################################
# Input:
# ------
# term: a single GO.ID
# ontologies: named list where names=GO.IDs and values=Ontologies, e.g CC, MF, BP,
#
# Return:
# --------
# The function required to get all the ancestors for the GO.ID
##################################################################################

determineAncesterFunction <- function(term, ontologies){
  ontology <- ontologies[term]
  if (is.na(ontology)){
    return(NA)
  }
  if (ontology == "MF"){
    return(GOMFANCESTOR)
  }
  else if (ontology == "CC"){
    return(GOCCANCESTOR)
  }
  else if (ontology == "BP"){
    return(GOBPANCESTOR)
  }
  else{return(NA)}
}


determineOffspringFunction <- function(term, ontologies){
  ontology <- ontologies[term]
  if (is.na(ontology)){
    return(NA)
  }
  if (ontology == "MF"){
    return(GOMFOFFSPRING)
  }
  else if (ontology == "CC"){
    return(GOCCOFFSPRING)
  }
  else if (ontology == "BP"){
    return(GOBPOFFSPRING)
  }
  else{return(NA)}
}

##################################################################################
# Input:
# ------
# go_ids: list of GO.IDs
# ontologies: named list where names=GO.IDs and values=Ontologies, e.g CC, MF, BP,
#
# Return:
# --------
# Named list where names=GO.IDs and values=all ancestor GO IDs
##################################################################################
getAllMappings <- function(go_ids, ontologies, verbose=TRUE, direction="ancester"){
  go2relations <- NULL
  
  if(direction=="ancester"){
    determineFunction <- determineAncesterFunction
  }
  
  else if(direction=="offspring"){
    determineFunction <- determineOffspringFunction
  }
  
  else{
    stop("direction must be `ancester` or `offspring`")
  }
  
  if(verbose){
    print(sprintf("Getting all %s GO terms for %i observed terms. This may take a while!", direction, length(go_ids)))
    pb <- txtProgressBar(min = 0, max = length(go_ids), style=3)
  }
  for (i in 1:length(go_ids)){
    query <- go_ids[i]
    
    go_relations <- tryCatch(get(query, determineFunction(query, ontologies)), error = function(e) e)
    
    if(class(go_relations)=="character"){
      go2relations[[go_ids[i]]] <- go_relations
    }
    else{
      if(verbose){
        print(sprintf("Could not get %ss for GO.ID=%s", direction, query))
      }
      go2relations[[go_ids[i]]] <- c(NA)
      
    }
    if(verbose){
      setTxtProgressBar(pb, i)
    }
  }
  
  return(go2relations)
}


##################################################################################
# Input:
# ------
# go_df: a data.frame for a single protein, with a column=="GO.ID"
# go2Ancester: Named list where names=GO.IDs and values=all ancestor GO IDs,
#              as returned by getAllMappings
#
# Return:
# --------
# Data.frame containing all ancestors for a protein
##################################################################################
ExpandTerms <- function(go_df, go_col, go2Ancester){
  
  observed_go_ids <- go_df %>% pull(go_col) %>% unique() # unique() shouldn't be req. but does not harm
  unprocessed_ids <- observed_go_ids  
  all_ancestors <- observed_go_ids
  
  while (length(unprocessed_ids) > 0){
    query_ancestors <- unique(go2Ancester[unprocessed_ids[1]])
    all_ancestors <- c(all_ancestors, unlist(query_ancestors))
    unprocessed_ids <- unprocessed_ids[-1] # removed processed go id
    unprocessed_ids <- base::setdiff(unprocessed_ids, query_ancestors) # and any which were found in ancestors (we don't want to re-query these)
  }
  all_ancestors <- setdiff(unique(all_ancestors), "all")
  
  return(data.frame("GOID"=all_ancestors))
}


##################################################################################
# Input:
# ------
# go_df: a data.frame with all initial GO terms for all proteins,
#        with columns==[($$option$$:feature_col), "GO.ID"]
# feature_col: the name of the column with the features, e.g UNIPROTKB
# go_col: the name of the column with the GO ids, e.g "GO.ID"
#
# Return:
# --------
# Data.frame: containing all ancestors for all proteins
#             plus term decriptions and ontologies
##################################################################################
getAllGOTerms <- function(go_df, feature_col="UNIPROTKB", go_col="GO.ID", return_early_debug=FALSE, verbose=TRUE){

  print(head(go_df, 2))
  print(head(go_df[[go_col]], 2))
  
  all_observed_go <- unique(go_df[[go_col]])
  all_observed_go <- all_observed_go[!is.na(all_observed_go)]
  
  ontologies <- AnnotationDbi::select(GO.db, all_observed_go, columns = c('ONTOLOGY'), keytype='GOID')
  ontologies <- setNames(ontologies$ONTOLOGY, ontologies$GOID)
  
  go2Ancester <- getAllMappings(all_observed_go, ontologies, direction="ancester", verbose=verbose)
  
  
  if(return_early_debug){
    return(list("go_df"=go_df, "go2Ancester"=go2Ancester))
  }
  
  print("Expanding GO terms to include all ancestors for all entries")
  full_go_df <- go_df %>% filter_(paste0("!is.na(", go_col, ")")) %>%
    group_by_(feature_col) %>% do(ExpandTerms(., go_col, go2Ancester))
  
  full_go_details <- AnnotationDbi::select(GO.db, as.character(unique(full_go_df$GOID)),
                            columns = c('TERM', 'ONTOLOGY'), keytype='GOID')
  
  full_go_df <- merge(full_go_df, full_go_details, by="GOID", all.x=TRUE)
  
  full_go_df <- full_go_df[,c(feature_col, 'GOID', 'TERM', 'ONTOLOGY')]
  colnames(full_go_df)[colnames(full_go_df)=="GOID"] <- go_col

  return(full_go_df)
}



