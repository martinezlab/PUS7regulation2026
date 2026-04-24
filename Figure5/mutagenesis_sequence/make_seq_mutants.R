#!/usr/bin/env Rscript
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

find_delpos <- function(sequence, psipos) {
  # 1. Get the part of the sequence before psipos
  prefix <- str_sub(sequence, 1, psipos - 1)
  # 2. Find a run of 'T's at the very end of this prefix
  t_run <- str_extract(prefix, "T+$")
  # 3. Calculate the length of the run. Replace NAs with 0.
  run_length <- nchar(t_run)
  run_length <- ifelse(is.na(run_length), 0, run_length)
  # 4. Calculate delpos
  delpos <- psipos - run_length
  # 5. Handle NAs
  delpos <- ifelse(is.na(sequence) | is.na(psipos), NA_integer_, delpos)
  
  return(delpos)
}

TGTAG_mutation <- function(sequence, fivemer, psipos) {
  if (is.na(fivemer) || is.na(psipos)) return(rep(NA_character_, 8))
  start_pos <- psipos - 2
  # Check if the 5-mer matches
  if (str_sub(sequence, start_pos, start_pos + 4) != fivemer) {
    warning("The 5-mer at the specified position does not match the provided fivemer.")
    return(list()) 
  }
  mutations <- c("TGTAG", "TGTAA", "TATAG", "TCTAG", "TTTAG", "TATAA", "TCTAA", "TTTAA")
  mutated_sequences <- rep(sequence, length(mutations))
  str_sub(mutated_sequences, start_pos, start_pos + 4) <- mutations
  
  return(setNames(as.list(mutated_sequences), paste0(mutations, "_seq")))
}
generate_mutated_rows <- function(df, mutation_fn, new_seq_type) {
  if (nrow(df) == 0) return(tibble())

  # 1 & 2. Generate mutations and unnest
  mutations_df <- df %>%
    mutate(mutations = pmap(list(seq, psipos_5mer, psipos), mutation_fn)) %>%
    filter(lengths(mutations) > 0) %>% 
    unnest_longer(mutations, values_to = "new_seq", indices_to = "mutation_id") %>% 
    filter(!is.na(new_seq))
  
  if (nrow(mutations_df) == 0) return(tibble())
  
  # 3. Update columns for the new rows
  mutations_df %>%
    mutate(
      seq = new_seq,
      sequence_type = new_seq_type,
      name = if (new_seq_type %in% c("UNUAR", "NNUNN")) {
        paste(parent_name, str_remove(mutation_id, "_seq"), sep = "_")
      } else {
        paste(parent_name, new_seq_type, sep = "_")
      },
      psipos_5mer = str_sub(seq, psipos - 2, psipos + 2),
      delpos = find_delpos(seq, psipos) 
    ) %>%
    dplyr::select(any_of((colnames(df))))
}

targets_df <- read.csv("target_sites_seq_mutations.csv", header = TRUE)

unuar_mutations <- generate_mutated_rows(
  df = targets_df,
  mutation_fn = TGTAG_mutation,
  new_seq_type = "UNUAR"
)

write.csv(unuar_mutations, "mutated_seqs.csv", row.names = FALSE)
