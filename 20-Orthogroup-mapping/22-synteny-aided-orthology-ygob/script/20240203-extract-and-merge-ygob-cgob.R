# title: extract and merge YGOB and CGOB orthology relationships for four yeasts
# author: Bin He
# date: 2024-02-03

# --------------
# Load libraries
# --------------
require(tidyverse)

# --------------
# Read in data
# --------------

# Read in YGOB
tmp <- read.table("../input/YGOB-Pillars.tab", na.string = "---")
# according to the YGOB_README:
# Col8 : C. glabrata 1
# Col12: S. cereevisiae 1
# Col16: K. lactis
# Col22: S. cerevisiae 2
# Col26: C. glabrata 2
ygob.select <- c(
  "Scer1" = 12,
  "Scer2" = 22,
  "Cgla1" = 8,
  "Cgla2" = 26,
  "Klac"  = 16
)
colnames(tmp)[ygob.select] = names(ygob.select)
ygob <- as_tibble(tmp[, ygob.select]) %>% 
  filter(
    # drop rows with all missing data
    !if_all(everything(), ~ is.na(.x))
    # note: there are still CEN, tRNA and other non protein coding
    #       gene elements in this table
  ) %>% 
  # if Scer1 is NA and Scer2 is not, move Scer2 to Scer1. order doesn't matter
  # in the Pillars.tab file (so Scer1 and Cgla1 are not necessarily paired)
  mutate(
    Scer1 = ifelse(is.na(Scer1) & !is.na(Scer2), Scer2, Scer1),
    Scer2 = ifelse(Scer1 == Scer2, NA, Scer2)
  ) %>% 
  # compute the number of genes each row
  mutate(
    across(everything(), is.na, .names = "na_{.col}"),
    n.NAs = rowSums(pick(starts_with('na')))
  ) %>% 
  select(-starts_with("na"))


# Read in CGOB
tmp <- read.table("../input/CGOB-Pillars-edited.tab", na.strings = "---")
# according to the CGOB_README:
# Col1 : C. albicans SC5314
# Col18: S. cerevisiae

cgob.select <- c(
  "Scer" = 18,
  "Calb" = 1
)
colnames(tmp)[cgob.select] = names(cgob.select)
cgob <- as_tibble(tmp[, cgob.select]) %>% 
  filter(
    # drop rows with all missing data
    !if_all(everything(), ~ is.na(.x)),
    # remove tRNA entries etc. in S. cerevisiae, most have no matches
    !str_sub(Scer, 1, 3) %in% c("it_", "t_s", "yit"),
    # turn of tRNA and predicted ORFs in C. albicans, most have no matches
    is.na(Calb) | str_sub(Calb, 1, 3) == "orf" 
  )

# test merging the two
merge0 <- ygob %>% 
  # remove rows where Scer1 is NA and there is only one column with non NA
  filter(!(is.na(Scer1) & n.NAs == 4)) %>% 
  # treating NA as non-matches (keep them in ygob but ignore in cgob)
  left_join(cgob, by = c("Scer1" = "Scer"), na_matches = "never") %>% 
  mutate(n.NAs = n.NAs + is.na(Calb)) %>% 
  relocate(n.NAs, .after = last_col()) %>% 
  arrange(desc(n.NAs))

# export the result
write_tsv(merge0, file = "../output/20240206-4sps-ygob-cgob-based-orthology-map.tsv")

# filter for the rows containing exactly one ortholog per species
single.gene <- merge0 %>% 
  filter(
    # only one Scer ortholog
    xor(is.na(Scer1), is.na(Scer2)),
    # only one Cgla ortholog
    xor(is.na(Cgla1), is.na(Cgla2)),
    # Klactis ortholog
    !is.na(Klac),
    # Calbicans ortholog
    !is.na(Calb)
  )

# export the result
write_tsv(merge0, file = "../output/20240206-4sps-ygob-cgob-one-gene-per-sps.tsv")
