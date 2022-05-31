#' Clonotype by CDR3 amino acid or nucleotide seq
#'
#' Use the cellranger-determined clonotypes (based on full nt seq of VDJ region) and the CDR3 amino acid sequences to add additional clonotype classifications to each original clonotype. Two or more cellranger (raw) clonotypes might have the same exact CDR3 amino acid sequence -- this function would give these different raw clonotypes the same "cdr3aa" clonotype id. The same process is done for CDR3 nt sequences. 
#'
#' @param clonotypes_and_contigs_list [Required] A list of length == 2. The first element is a tibble data frame of clonotypes (derived from the cellranger clonotypes.csv file). The second element is a tibble data frame of contig annotations (derived from the cellranger filtered_contig_annotations.csv file).
#'
#' 
#'
#' @return A list of length 2: list[[1]] are the clonotypes; list[[2]] are the contigs.
#'
#'
#'
#'
#' @examples
#' \dontrun{
#' clonotypes <- read_csv("outs/clonotypes.csv", col_types = cols(.default = "c"))
#' filtered_contig_annotations <- read_csv("outs/filtered_contig_annotations.csv", col_types = cols(.default = "c"))
#' out <- tk_add_nt_and_aa_clonotypes(list(clonotypes, filtered_contig_annotations))
#' }
#'
#' @export
add_nt_and_aa_clonotypes <- function(clonotypes_and_contigs_list = list(clonotypes, filtered_contig_annotations)) {

    # The "clonotypes.csv" tibble
    cl <- clonotypes_and_contigs_list[[1]]
    # The "filtered_contig_annotations.csv" tibble
    vdj <- clonotypes_and_contigs_list[[2]]


    # Make sure cdr3_aa is sorted properly
    cdr3s_aa_ord <- vector()
    cdr3s_nt_ord <- vector()
    for (i in seq_len(dim(cl)[1])) {
        curr_aa <- cl %>%
            dplyr::slice(i) %>%
            pull(cdr3s_aa) %>%
            str_split(., ";")

        cdr3s_aa_ord[i] <- curr_aa[[1]][order(curr_aa[[1]])] %>%
            paste(., collapse = ";")
        
        curr_nt <- cl %>%
            dplyr::slice(i) %>%
            pull(cdr3s_nt) %>%
            str_split(., ";")

        cdr3s_nt_ord[i] <- curr_nt[[1]][order(curr_nt[[1]])] %>%
            paste(., collapse = ";")
    }

    cl2 <- cl %>%
        add_column(cdr3s_aa_ord, cdr3s_nt_ord)


    # Find the distinct CDR3 aa and nt clonotypes
    cl3 <- cl2 %>%
        dplyr::select(cdr3s_aa_ord) %>%
        distinct() %>%
        dplyr::mutate(cdr3s_aa_ord_random_id = dplyr::row_number())
    
    cl4 <- cl2 %>%
        dplyr::select(cdr3s_nt_ord) %>%
        distinct() %>%
        dplyr::mutate(cdr3s_nt_ord_random_id = dplyr::row_number())
    


    cl2_aa_ids <- cl2 %>%
        left_join(cl3, by = c("cdr3s_aa_ord"))

    cl2_aa_ids_nt_ids <- cl2_aa_ids %>%
        left_join(cl4, by = c("cdr3s_nt_ord"))



    # Add random CDR3 aa and nt ids back to contigs table
    cl2_aa_ids_nt_ids_for_join <- cl2_aa_ids_nt_ids %>%
        dplyr::select(clonotype_id, cdr3s_aa_ord, cdr3s_nt_ord, cdr3s_aa_ord_random_id, cdr3s_nt_ord_random_id)

    vdj2 <- vdj %>%
        dplyr::left_join(cl2_aa_ids_nt_ids_for_join, by = c("raw_clonotype_id" = "clonotype_id"))


    # Get a list of each unique CDR3aa clonotype
    cdr3aa_clonotype_unique <- cl2_aa_ids_nt_ids %>%
        dplyr::select(cdr3s_aa_ord_random_id, cdr3s_aa_ord) %>%
        distinct()

    # Give a list of just each unique CDR3nt clonotype
    cdr3nt_clonotype_unique <- cl2_aa_ids_nt_ids %>%
        dplyr::select(cdr3s_nt_ord_random_id, cdr3s_nt_ord) %>%
        distinct()

 
    # This tally mimics the web_summary.html -- but at the CDR3aa clonotype level
    gem_tally_by_cdr3aa_clonotype <- vdj2 %>%
        distinct(barcode, .keep_all = TRUE) %>%
        group_by(cdr3s_aa_ord_random_id) %>%
        tally() %>%
        arrange(desc(n)) %>%
        dplyr::rename(number_of_gems_with_same_cdr3aa_clonotype = "n") %>%
        ungroup() %>%
        dplyr::left_join(cdr3aa_clonotype_unique, by = "cdr3s_aa_ord_random_id")


    # This tally mimics the web_summary.html -- but at the CDR3nt clonotype level
    gem_tally_by_cdr3nt_clonotype <- vdj2 %>%
        distinct(barcode, .keep_all = TRUE) %>%
        group_by(cdr3s_nt_ord_random_id) %>%
        tally() %>%
        arrange(desc(n)) %>%
        dplyr::rename(number_of_gems_with_same_cdr3nt_clonotype = "n") %>%
        ungroup() %>%
        dplyr::left_join(cdr3nt_clonotype_unique, by = "cdr3s_nt_ord_random_id")


    # Rename CDR3aa and CDR3nt clonotypes ids based on the contig frequency observed 
    gem_tally_by_cdr3aa_clonotype$clonotype_id_cdr3aa <- paste0("clonotype_cdr3aa", as.numeric(factor(gem_tally_by_cdr3aa_clonotype$cdr3s_aa_ord_random_id, levels = unique(gem_tally_by_cdr3aa_clonotype$cdr3s_aa_ord_random_id))))
    gem_tally_by_cdr3nt_clonotype$clonotype_id_cdr3nt <- paste0("clonotype_cdr3nt", as.numeric(factor(gem_tally_by_cdr3nt_clonotype$cdr3s_nt_ord_random_id, levels = unique(gem_tally_by_cdr3nt_clonotype$cdr3s_nt_ord_random_id))))


    # Add new names from clonotypes file to contigs file
    aa_for_join <- gem_tally_by_cdr3aa_clonotype %>%
        dplyr::select(cdr3s_aa_ord_random_id, clonotype_id_cdr3aa)

    nt_for_join <- gem_tally_by_cdr3nt_clonotype %>%
        dplyr::select(cdr3s_nt_ord_random_id, clonotype_id_cdr3nt)

    vdj3 <- vdj2 %>%
        dplyr::left_join(aa_for_join, by = "cdr3s_aa_ord_random_id") %>%
        dplyr::left_join(nt_for_join, by = "cdr3s_nt_ord_random_id") %>%
        dplyr::select(-cdr3s_aa_ord_random_id, -cdr3s_nt_ord_random_id) %>%
        dplyr::rename(cdr3s_aa = "cdr3s_aa_ord", cdr3s_nt = "cdr3s_nt_ord")


    cl5 <- cl2_aa_ids_nt_ids %>%
        dplyr::left_join(aa_for_join, by = "cdr3s_aa_ord_random_id") %>%
        dplyr::left_join(nt_for_join, by = "cdr3s_nt_ord_random_id") %>%
        dplyr::select(clonotype_id, frequency, proportion, cdr3s_aa_ord, cdr3s_nt_ord, inkt_evidence, mait_evidence, clonotype_id_cdr3aa, clonotype_id_cdr3nt) %>%
        dplyr::rename(cdr3s_aa = "cdr3s_aa_ord", cdr3s_nt = "cdr3s_nt_ord")
     

    return(list(clonotypes = cl5, contigs = vdj3))
}



