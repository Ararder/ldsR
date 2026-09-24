# Builds inst/extdata/rsid_map.parquet. Two steps: (1) locally, write the RSIDs of the
# bundled LD score + weight files; (2) on dardel, look them up in the tidyGWAS dbSNP155 reference.

# step 1 (local, from package root) ----------------------------------------
library(tidyverse)
bind_rows(
  arrow::read_parquet("inst/extdata/eur_w_ld.parquet", col_select = "SNP"),
  arrow::read_parquet("inst/extdata/1000G_Phase3_weights_hm3_no_MHC.parquet", col_select = "SNP")
) |>
  distinct() |>
  arrow::write_parquet("ldsc_snps.parquet")
# rsync ldsc_snps.parquet to `out` on dardel

# step 2 (dardel) -----------------------------------------------------------
library(tidyverse)
ref <- "~/ki-pgi-storage/Data/downstreamGWAS/reference/dbSNP155"
out <- "/cfs/klemming/scratch/a/arvhar/ldsR_rsid_map"
snps <- arrow::read_parquet(file.path(out, "ldsc_snps.parquet")) |>
  mutate(id = as.integer(str_remove(SNP, "^rs")))

merged <- arrow::read_parquet(file.path(ref, "refsnp-merged/part-0.parquet")) |>
  filter(old_RSID %in% snps$SNP) |>
  transmute(old_id = as.integer(str_remove(old_RSID, "^rs")), new_id = as.integer(str_remove(RSID, "^rs")))
cat("merged rsids among ldsc snps:", nrow(merged), "\n")

# lookup = current dbSNP id; SNP keeps the rsid used by the LD score reference
key <- snps |>
  left_join(merged, by = c("id" = "old_id")) |>
  transmute(SNP, lookup = coalesce(new_id, id))

db <- arrow::open_dataset(file.path(ref, "v155"))
res <- map(as.character(1:22), \(chr) {
  hit <- db |>
    filter(CHR == chr, RSID %in% key$lookup) |>
    select(RSID, POS37 = POS_37, POS38 = POS_38, REF = REF_38, ALT = ALT_38) |>
    collect()
  cat("chr", chr, ":", nrow(hit), "\n")
  mutate(hit, CHR = as.integer(chr))
}) |> list_rbind()

map_tbl <- key |>
  inner_join(res, by = c("lookup" = "RSID"), relationship = "many-to-many") |>
  separate_longer_delim(ALT, ",") |>
  select(CHR, POS37, POS38, SNP, REF, ALT) |>
  arrange(CHR, POS38)

cat("rows:", nrow(map_tbl), " snps:", n_distinct(map_tbl$SNP), "/", nrow(snps), "\n")
print(summarise(map_tbl, na37 = sum(is.na(POS37)), na38 = sum(is.na(POS38)), multi_rs = sum(duplicated(paste(SNP, ALT)))))
arrow::write_parquet(map_tbl, file.path(out, "rsid_map.parquet"), compression = "zstd", compression_level = 19)
