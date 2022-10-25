#' Check gene names against NCBI database and fix nomenclature
#'
#' @param gene_vector Required. Character vector of input gene names. These can be slightly incorrect in their capitalization or be gene aliases.
#' @param gene_db Required. Tibble data frame of the GENE database downloaded from NCBI. See below for details.
#'
#' How to download gene_db:
#' url <- "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
#' url <- "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz"
#' download.file(url, basename(url), method = "wget")
#' gene_db <- read_tsv(gzfile(basename(url)))
#' @return None
#'
#' @examples
#' \dontrun{
#' genes_clean <- check_gene_db(genes_dirty, gene_db)
#' }
#'
#' @export
check_gene_db <- function(gene_vector, gene_db) {
    .check_gene_db <- function(input_symbol, gene_db) {
        input_symbol_camel <- str_replace_all(snakecase::to_sentence_case(input_symbol), " ", "")
    
        if (input_symbol %in% gene_db$Symbol) {
            final_name <- input_symbol
            note <- "original_input_symbol"
            possible_df <- NULL
            possible_symbol <- gene_db$Symbol[which(gene_db$Symbol %in% input_symbol)[1]]
            description <- gene_db$description[which(gene_db$Symbol %in% input_symbol)[1]]
            synonyms <- gene_db$Synonyms[which(gene_db$Symbol %in% input_symbol)[1]]
        } else if (input_symbol_camel %in% gene_db$Symbol) {
            final_name <- input_symbol_camel
            note <- "converted_to_camel_case"
            possible_df <- gene_db[which(gene_db$Symbol %in% input_symbol_camel), c("Symbol", "Synonyms", "description")] %>% 
                    dplyr::mutate(input_symbol = input_symbol) %>%
                    dplyr::select(input_symbol, everything())
            if (dim(possible_df)[1] == 1) {
                possible_symbol <- possible_df$Symbol
                description <- possible_df$description
                synonyms <- possible_df$Synonyms
            } else {
                possible_symbol <- NA
                description <- NA
                synonyms <- NA
            }
        } else {
            if (length(which(grepl(input_symbol, gene_db$Synonyms))) >= 1) {
                final_name <- NA
                note <- "fuzzy_match_original_input_symbol"
                possible_df <- gene_db[which(grepl(input_symbol, gene_db$Synonyms)), c("Symbol", "Synonyms", "description")] %>% 
                    dplyr::mutate(input_symbol = input_symbol) %>%
                    dplyr::select(input_symbol, everything())
                if (dim(possible_df)[1] == 1) {
                    possible_symbol <- possible_df$Symbol
                    description <- possible_df$description
                    synonyms <- possible_df$Synonyms
                } else {
                    possible_symbol <- NA
                    description <- NA
                    synonyms <- NA
                }
            } else {
                final_name <- NA
                note <- "unknown"
                possible_df <- gene_db[agrep(input_symbol, gene_db$Synonyms), c("Symbol", "Synonyms", "description")] %>% 
                    dplyr::mutate(input_symbol = input_symbol) %>%
                    dplyr::select(input_symbol, everything())
                if (dim(possible_df)[1] == 1) {
                    possible_symbol <- possible_df$Symbol
                    description <- possible_df$description
                    synonyms <- possible_df$Synonyms
                } else {
                    possible_symbol <- NA
                    description <- NA
                    synonyms <- NA
                }
            } 
        }
        out <- list(input_symbol = input_symbol, final_name = final_name, possible_symbol = possible_symbol, note = note, description = description, synonyms = synonyms, possible_df = possible_df)
        return(out)
    }
    # Get results and combine
    y <- lapply(gene_vector, .check_gene_db, gene_db)

    a <- unlist(lapply(y, '[[', 1), recursive = FALSE)
    b <- unlist(lapply(y, '[[', 2), recursive = FALSE)
    c <- unlist(lapply(y, '[[', 3), recursive = FALSE)
    d <- unlist(lapply(y, '[[', 4), recursive = FALSE)
    e <- unlist(lapply(y, '[[', 5), recursive = FALSE)
    f <- unlist(lapply(y, '[[', 6), recursive = FALSE)
    g <- lapply(y, '[[', 7)
    results <- tibble(input_symbol = a, final_symbol = b, possible_symbol = c, note = d, description = e, synonyms = f, possible_df = g)
    return(results)
}





