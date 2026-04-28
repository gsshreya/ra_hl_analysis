library("GEOquery")

gse_ra <- getGEO("GSE55235", GSEMatrix = TRUE) #downloading data for RA
ra <- gse_ra[[1]]
meta_ra <- pData(ra)

group_raw <- meta_ra$`disease state:ch1`
keep_ra <- group_raw %in% c("healthy control", "rheumatoid arthritis") #filtering conditions

expr_ra_filt <- exprs(ra)[, keep_ra]
meta_ra_filt <- meta_ra[keep_ra, ]

group_ra <- ifelse(meta_ra_filt$`disease state:ch1` == "healthy control", "Control", "RA") #defining groups

gse_hl <- getGEO("GSE12453", GSEMatrix =TRUE) #downloading data for HL
hl <- gse_hl[[1]]

meta_hl <- pData(hl)

group_hl_raw <- meta_hl$characteristics_ch1

keep_hl <- grepl("hodgkin", group_hl_raw, ignore.case = TRUE) &
           grepl("classical", group_hl_raw, ignore.case = TRUE) |
           grepl("centroblasts|centrocytes", group_hl_raw, ignore.case = TRUE)

expr_hl_filt <- exprs(hl)[, keep_hl]
meta_hl_filt <- meta_hl[keep_hl, ]

group_hl <- ifelse(
  grepl("hodgkin", meta_hl_filt$characteristics_ch1, ignore.case = TRUE),
  "HL",
  "Control"
) #defining groups

write.csv(expr_ra_filt, "data/processed/ra_expression.csv")
write.csv(expr_hl_filt, "data/processed/hl_expression.csv")
write.csv(meta_ra_filt, "data/processed/ra_metadata.csv")
write.csv(meta_hl_filt, "data/processed/hl_metadata.csv")
write.csv(group_ra, "data/processed/ra_groups.csv", row.names = FALSE)
write.csv(group_hl, "data/processed/hl_groups.csv", row.names = FALSE) #saving 3 types of csvs for RA and HL each