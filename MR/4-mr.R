#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(ggplot2)

library(future)
library(future.apply)

library(TwoSampleMR) # mr
# remotes::install_github("MRCIEU/TwoSampleMR")
library(ieugwasr)
library(plinkbinr)
#library(LDlinkR)


root <- "/vol/projects/CIIM/meta_cQTL/out"
cohort <- "2000HIV-EU-discovery"
cov <- "main"
bfile <- paste0(root, "/", cohort, "/genotype/allchr")

#n_samples <- 267       # validation 2000HIV-EU
n_samples <- 1064       # discovery 2000HIV-EU

pheno2 <- "cmv"   # outcome
pheno2 <- "cytokines"  # comment out one or the other

e <- "FCRL6_Inflammation"
e <- "ENSG00000181036"  # comment out one or the other
root_in <- "/home/nvanunen/nvanunen/extra_projects/hua_gwas_mr/out_trans_maf_0.01"
root_out <- "/home/nvanunen/nvanunen/extra_projects/hua_gwas_mr/out_tier1/"

file1 <- paste0(root_in, "/1-clump/", e, ".tsv")
file2 <- paste0(root_in, "/2-outcome/", e, "_", pheno2, ".tsv")

testing <- FALSE

outpath <- paste0(root_out, "/", e, "_", pheno2)
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)


# clumping is done ahead with genomewide_clumping.R
trait1 <- file1 %>% basename() %>% str_remove("\\.tsv")
message(trait1)

# distinct cause for some eqtl it gives duplicate rows, same snp, beta, everything
exposure_dat <- fread(file1) %>% distinct()
outcome_all <- fread(file2) %>% distinct()


#selected <- c("chr2:206202754:C:T;rs17223074", "chr8:6342145:G:A;rs73661710", "chr16:64703819:G:T;rs74958584", "chr18:41148643:G:T;rs192206089", "chr20:60824146:G:A;rs34788939")
exposure_dat <- exposure_dat %>% dplyr::filter(ID %in% tier1$SNP)

full_sub <- data.frame()

trait2 <- unique(outcome_all$gene)[1]
for (trait2 in unique(outcome_all$gene)) {
    outcome <- outcome_all %>% filter(gene == trait2)
    nrow(outcome) # 1

    outcome_dat <- outcome %>%
        dplyr::rename(pval.outcome = "pval", beta.outcome = "beta",
        chr.outcome = "chr", pos.outcome = "location", effect_allele.outcome = "effect_allele",
        other_allele.outcome = "other_allele", se.outcome = "se", eaf.outcome = "eaf.exposure", 
        outcome = "gene", id.outcome = "pheno")
    
    #outcome_dat <- outcome_dat %>% dplyr::filter(SNP %in% (tier1 %>% separate(SNP, c("extra", "SNP"), sep = ";")))
        
    # sometimes this version works instead:
    #outcome_dat <- outcome %>% rename(pval = "pval.outcome", beta = "beta.outcome", chr = "chr.outcome", 
    #    location = "pos.outcome", effect_allele = "effect_allele.outcome", other_allele = "other_allele.outcome",
    #    se = "se.outcome", eaf.exposure = "eaf.outcome", gene = "outcome", pheno = "id.outcome")

    dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
    dat$samplesize.outcome <- n_samples
    dat$samplesize.exposure <- n_samples

    if (nrow(dat) == 0) {
        #failed <- rbind(failed, data.frame(trait1 = trait1, trait2 = trait2, 
        #    reason = "No snps after harmonising"))
        message("No snps after harmonising, skipping ...")
        next
    }
    # filter out mr_keep = FALSE
    dat <- dat %>% filter(mr_keep == TRUE)
    #if (nrow(dat) < 3) {
        #failed <- rbind(failed, data.frame(trait1 = trait1, trait2 = trait2, 
        #    reason = "Less than 3 snps after harmonising"))
    #    message("Less than 3 snps after harmonising, skipping ...")
    #    next
    #}

    outname <- paste0(outpath, "/", trait2) #ivw_pval %>% format(scientific = TRUE, digits = 3), "_",
    dir.create(outname, showWarnings = FALSE, recursive = TRUE)

    # save data
    fwrite(dat, paste0(outname, "/data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    if (nrow(dat) == 1) {
        # perform wald ratio
        res <- mr_singlesnp(dat)[1,]
    } else {
        # perform ivw MR
        all_res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_simple_median", "mr_weighted_median"))
        res <- all_res %>% filter(method == "Inverse variance weighted") %>%
            mutate(pleio = NA, het = NA, lou = NA)
        ivw_pval <- res %>% pull(pval)

        # if significant, perform sensitivity checks and plots
        if (ivw_pval < 0.05) {
            # single snp  & forest
            s_res <- mr_singlesnp(dat)
            tryCatch({
                frplot <- mr_forest_plot(s_res)
                ggsave(plot = frplot[[1]], filename = paste0(outname, "/forest.pdf"), width = 4, height = 5)
            }, error = function(e) {
                # e.g. Error in `levels<-`(`*tmp*`, value = as.character(levels)): factor level [5] is duplicated
                # happens when directionality also cannot be calculated
                # and there are warnings on perfect fits... 
                # e.g. In summary.lm(mod) : essentially perfect fit: summary may be unreliable
                message("Error in forest plot ...\n", e)
            })
            
            # Plot MR
            out_plot <- paste0(outname, "/scatter.pdf")
            mrplot <- mr_scatter_plot(all_res, dat)[[1]] + 
                theme_bw() + 
                labs(title = paste0("MR ", trait1, " vs ", trait2), 
                    subtitle = paste0("IVW p-value: ", ivw_pval %>% format(scientific = TRUE, digits = 3))
                                    #"\nDirectionality p-value: ", dir$steiger_pval %>% format(scientific = TRUE, digits = 3),
                                    #"\nHeterogeneity p-value: ", het$Q_pval %>% format(scientific = TRUE, digits = 3), 
                                    #"\nPleiotropy p-value: ", pleio$pval %>% format(scientific = TRUE, digits = 3),
                                    #"\nLeave-one-out max p-value: ", max(lou$p) %>% format(scientific = TRUE, digits = 3)
                                    ) +
                    geom_vline(aes(xintercept = 0), linetype = "dashed") + 
                    geom_hline(aes(yintercept = 0), linetype = "dashed")
            ggsave(plot = mrplot, filename = out_plot, width = 10, height = 5)
            message("Saved plot to ", out_plot)

            # sensitivity checks
            pleio <- mr_pleiotropy_test(dat)
            res$pleio <- pleio$pval
            het <- mr_heterogeneity(dat, method_list = "mr_ivw")
            res$het <- het$Q_pval
            lou <- mr_leaveoneout(dat)
            if (!is.null(lou)) {
                if (any(is.na(lou$p))) {
                    res$lou <- NA
                } else {
                    res$lou <- max(lou$p)
                }
            }
        }
    }
    full_sub <- rbind(full_sub, res)
}

outfile <- paste0(outpath, "/all_mr.tsv")
fwrite(full_sub, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved results to ", outfile)

full_sub %>% dplyr::filter(pval < 0.05)

