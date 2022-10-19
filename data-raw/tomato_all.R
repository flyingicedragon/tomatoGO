## code to prepare `tomato_all` dataset goes here
go_basic <- read.table(
  "data-raw/go-basic.tb",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)[, 1:2]
names(go_basic) <- c("gsid", "name")
tomato_allgene <- read.table("data-raw/tomato_go.tsv", sep = "\t")
names(tomato_allgene) <- c("gene", "gsid", "type")
tomato_all <- subset(tomato_allgene, type != "")
tomato_bp <- subset(tomato_allgene, type == "biological_process")
tomato_cc <- subset(tomato_allgene, type == "cellular_component")
tomato_mf <- subset(tomato_allgene, type == "molecular_function")
tomato_all <- gson::gson(
  gsid2gene = tomato_all[, c(2, 1)],
  gsid2name = go_basic,
  species = "tomato",
  gsname = "GO_all",
  version = "4.0",
  accessed_date = "2022/10/14",
  info = "All GO data in one"
)
tomato_bp <- gson::gson(
  gsid2gene = tomato_bp[, c(2, 1)],
  gsid2name = go_basic,
  species = "tomato",
  gsname = "GO_all",
  version = "4.0",
  accessed_date = "2022/10/14",
  info = "Biological process GO data in one"
)
tomato_cc <- gson::gson(
  gsid2gene = tomato_cc[, c(2, 1)],
  gsid2name = go_basic,
  species = "tomato",
  gsname = "GO_all",
  version = "4.0",
  accessed_date = "2022/10/14",
  info = "Cellular component GO data in one"
)
tomato_mf <- gson::gson(
  gsid2gene = tomato_mf[, c(2, 1)],
  gsid2name = go_basic,
  species = "tomato",
  gsname = "GO_all",
  version = "4.0",
  accessed_date = "2022/10/14",
  info = "Molecular function GO data in one"
)

usethis::use_data(tomato_all, overwrite = TRUE)
usethis::use_data(tomato_bp, overwrite = TRUE)
usethis::use_data(tomato_cc, overwrite = TRUE)
usethis::use_data(tomato_mf, overwrite = TRUE)