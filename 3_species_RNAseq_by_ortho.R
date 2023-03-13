library(tidyverse)

# don't load DESeq2 yet! (will throw an error; load later)

# import RNAseq data and reorganize columns in temporal order
Lv_rna_counts <- read_tsv("Lv_counts_frac_clean.txt")
Lv_rna_counts <- relocate(Lv_rna_counts, starts_with("LVEGG"), .after = L_var_ID)
Lv_rna_counts <- relocate(Lv_rna_counts, starts_with("LV4"), .after = LVEGGC)
Lv_rna_counts <- relocate(Lv_rna_counts, starts_with("LVMB"), .after = LV32cellC)

Ht_rna_counts <- read_tsv("Ht_counts_frac_clean.txt")
Ht_rna_counts <-relocate(Ht_rna_counts, starts_with("HTEGG"), .after = H_tub_ID)
Ht_rna_counts <-relocate(Ht_rna_counts, starts_with("HT4"), .after = HTEGGC)
Ht_rna_counts <-relocate(Ht_rna_counts, starts_with("HTMB"), .after = HT32cellC)

He_rna_counts <- read_tsv("He_counts_frac_clean.txt")
He_rna_counts <-relocate(He_rna_counts, starts_with("HEEGG"), .after = H_ery_ID)
He_rna_counts <-relocate(He_rna_counts, starts_with("HE4"), .after = HEEGGC)
He_rna_counts <-relocate(He_rna_counts, starts_with("HEMB"), .after = HE32cellC)

# rename columns 
Lv_rna_counts <- rename(Lv_rna_counts, LVeggA = LVEGGA, LVeggB = LVEGGB, LVeggC = LVEGGC)
Lv_rna_counts <- rename(Lv_rna_counts, LV4cA = LV4cellA, LV4cB = LV4cellB, LV4cC = LV4cellC)
Lv_rna_counts <- rename(Lv_rna_counts, LV16cA = LV16cellA, LV16cB = LV16cellB, LV16cC = LV16cellC)
Lv_rna_counts <- rename(Lv_rna_counts, LV32cA = LV32cellA, LV32cB = LV32cellB, LV32cC = LV32cellC)
Lv_rna_counts <- rename(Lv_rna_counts, LVblA = LVMBA, LVblB = LVMBB, LVblC = LVMBC)
Lv_rna_counts <- rename(Lv_rna_counts, LVgaA = LVGA, LVgaB = LVGB, LVgaC = LVGC)
Lv_rna_counts <- rename(Lv_rna_counts, LVprA = LVPRA, LVprB = LVPRB, LVprC = LVPRC)

Ht_rna_counts <- rename(Ht_rna_counts, HTeggA = HTEGGA,HTeggB = HTEGGB,HTeggC = HTEGGC)
Ht_rna_counts <- rename(Ht_rna_counts, HT4cA = HT4cellA,HT4cB = HT4cellB,HT4cC = HT4cellC)
Ht_rna_counts <- rename(Ht_rna_counts, HT16cA = HT16cellA,HT16cB = HT16cellB,HT16cC = HT16cellC)
Ht_rna_counts <- rename(Ht_rna_counts, HT32cA = HT32cellA,HT32cB = HT32cellB,HT32cC = HT32cellC)
Ht_rna_counts <- rename(Ht_rna_counts, HTblA = HTMBA,HTblB = HTMBB,HTblC = HTMBC)
Ht_rna_counts <- rename(Ht_rna_counts, HTgaA = HTGA,HTgaB = HTGB,HTgaC = HTGC)
Ht_rna_counts <- rename(Ht_rna_counts, HTprA = HTPRA,HTprB = HTPRB,HTprC = HTPRC)

He_rna_counts <- rename(He_rna_counts, HEeggA = HEEGGA,HEeggB = HEEGGB,HEeggC = HEEGGC)
He_rna_counts <- rename(He_rna_counts, HE4cA = HE4cellA,HE4cB = HE4cellB,HE4cC = HE4cellC)
He_rna_counts <- rename(He_rna_counts, HE16cA = HE16cellA,HE16cB = HE16cellB,HE16cC = HE16cellC)
He_rna_counts <- rename(He_rna_counts, HE32cA = HE32cellA,HE32cB = HE32cellB,HE32cC = HE32cellC)
He_rna_counts <- rename(He_rna_counts, HEblA = HEMBA,HEblB = HEMBB,HEblC = HEMBC)
He_rna_counts <- rename(He_rna_counts, HEgaA = HEGA,HEgaB = HEGB,HEgaC = HEGC)
He_rna_counts <- rename(He_rna_counts, HEprA = HEPRA,HEprB = HEPRB,HEprC = HEPRC)

# round counts to nearest integer
Lv_rna_counts <- Lv_rna_counts |> mutate_if(is.numeric, round)
Ht_rna_counts <- Ht_rna_counts |> mutate_if(is.numeric, round)
He_rna_counts <- He_rna_counts |> mutate_if(is.numeric, round)

# collapse counts from alternate transcripts
#Ht_rna_counts$H_tub_ID<- gsub("\\..*", "", Ht_rna_counts$H_tub_ID)
Lv_rna_counts$L_var_ID <- str_replace(Lv_rna_counts$L_var_ID, "\\..*", "")
Lv_rna_counts<-Lv_rna_counts %>% group_by(L_var_ID) %>% summarise_each(funs(sum))

# save intermediate tables 
write_tsv(Lv_rna_counts, "./intermediate_files/Lv_rna_raw_counts.txt") 
write_tsv(Ht_rna_counts, "./intermediate_files/Ht_rna_raw_counts.txt") 
write_tsv(He_rna_counts, "./intermediate_files/He_rna_raw_counts.txt") 

# import orthogroup tables
Lv_ortho <- read_tsv("ortho_Lv_ready.tsv")
Ht_ortho <- read_tsv("ortho_Ht_ready.tsv")
He_ortho <- read_tsv("ortho_He_ready.tsv")

# rename column and delete header of orthogroup table
Lv_ortho <- rename(Lv_ortho, L_var_ID = L_var.longestprot)
Lv_ortho <- Lv_ortho[-1, ]

# join on gene ID column to append orthogroups
Lv_orth_rna_counts <- left_join(Lv_rna_counts, Lv_ortho)


library(DESeq2)


