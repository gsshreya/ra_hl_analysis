library(limma)

expr_ra <- read.csv("data/processed/ra_expression.csv", row.names=1) #to get expression data for ra
group_ra <- read.csv ("data/processed/ra_groups.csv")

group_ra <- group_ra[,1]

expr_ra <- as.matrix(expr_ra)
expr_ra <- log2(expr_ra +1)

group_ra <- factor(group_ra)

design_ra <- model.matrix (~ group_ra) #settin up ra cases vs controls
colnames(design_ra) <- c("Intercept", "RA_vs_Control")

fit_ra <- lmFit(expr_ra, design_ra)
fit_ra <- eBayes(fit_ra)

res_ra <- topTable(fit_ra, coef = "RA_vs_Control", number = Inf)

write.csv(res_ra, "data/processed/ra_deg_results.csv")

expr_hl <- read.csv("data/processed/hl_expression.csv", row.names=1)
group_hl <- read.csv("data/processed/hl_groups.csv")

group_hl <- group_hl[,1]

expr_hl <- as.matrix(expr_hl)
expr_hl <- log2(expr_hl + 1)

group_hl <- factor(group_hl)

design_hl <- model.matrix(~ group_hl)
colnames(design_hl) <- c("Intercept", "HL_vs_Control")

fit_hl <- lmFit(expr_hl, design_hl)
fit_hl <- eBayes(fit_hl)

res_hl <- topTable(fit_hl, coef = "HL_vs_Control", number = Inf)

write.csv(res_hl, "data/processed/hl_deg_results.csv")