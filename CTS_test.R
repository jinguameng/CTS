## ============================================================================================================================
## Script: CTS.R
## Description: This script run the longitudinal GWAS with Conditional Two-Step approach
##              Reference paper: https://www.nature.com/articles/ejhg20151#Abs1
##
## Usage:
##		Rscript CTS.R <pheno_file> <geno_file> <frq_file> <id_col> <time_col> <y_col> <fixedEffects> <outpath> <outprefix>
##
##
## Arguments:
##		pheno_file      - Full path to the phenotype file, longitudinal dataset with long-format
##                      -- RDS format
##                      -- covariates should be included for LMM analysis
##                      -- FID and IID columns should also be included
##                      -- Example pheno_file: /data/h_vmac/zhanm32/SPAREAD/DATA/ADNI/CTS_ADNI.rds
##		geno_file       - Full path to the genotype file, file prefix included (PLINK .bim format)
##    frq_file        - Full path and full name of the SNP frequency file (e.g. /data/h_vmac/zhanm32/SPAREAD/DATA/ADNI/ADNI_EUR_snpweights_baseline_frq.frq)
## 		id_col          - Column name in the phenotype file representing subject IDs
##		time_col        - Column name in the phenotype file representing time variable
##		y_col           - Outcome variable column name
##    fixedEffects    - Comma-separated list of fixed effects for LMM model after CTS (e.g., "Age,Sex,Education")
##    outpath         - Full path to save the CTS outputs
##    outprefix       - Prefix of output
##
##
## Example Usage:
##		Rscript /data/h_vmac/zhanm32/SPAREAD/CODE/scripts/CTS_test.R \
##		"/data/h_vmac/zhanm32/SPAREAD/DATA/ADNI/CTS_ADNI.rds" \
##		"/data/h_vmac/zhanm32/SPAREAD/DATA/ADNI/ADNI_EUR_snpweights_baseline" \
##    "/data/h_vmac/zhanm32/SPAREAD/DATA/ADNI/ADNI_EUR_snpweights_baseline_frq.frq" \
##		"ID" \
##    "interval_years" \
##    "SPARE_AD" \
##    "Sex,PHC_Age_T1_MUSE,Education" \
##    "/data/h_vmac/zhanm32/SPAREAD/OUTPUT/longitudinalGWAS/" \
##    "EUR_ADNI"
##
## Important Note:
##    - Required packages: data.table, lme4, lmerTest, tidyverse, qqman
##		- Ensure that the phenotype file is cleaned before running this script:
##      1. handle missingness
##      2. make sure at least two time points per individuals
##      3. make sure the pheno_file include FID and IID columns that matching the geno file
## ============================================================================================================================




###############################################################################################################
#### Load packages
###############################################################################################################
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(lme4))
suppressPackageStartupMessages(require(lmerTest))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(qqman))

`%!in%` <- Negate(`%in%`)

###############################################################################################################
#### Check arguments and verify file existence
###############################################################################################################
args <- commandArgs(TRUE)

## Check if at least 5 arguments (without covariates) are provided
if (length(args) < 9) {
  stop("Error: Missing required arguments.\nUsage: Rscript CTS.R <pheno_file> <geno_file> <frq_file> <id_col> <time_col> <y_col> <fixedEffects> <outpath> <outprefix>\n")
}

## Required arguments
pheno_file    <- args[1] # Full path and full name of the phenotype file.
geno_file     <- args[2] # Full path and prefix of the genotype file.
frq_file      <- args[3] # Full and full name of the SNP frequency file.
id_col        <- args[4] # ID column name in the phenotype file. You have to make sure the IDs are consistent with the FID column in the genotype file.
time_col      <- args[5] # Time variable column in the phenotype file
y_col         <- args[6] # Outcome variable column in the phenotype file
fixedEffects  <- args[7] # Fixed effects
outpath       <- args[8] # Output path
outprefix     <- args[9] # Output prefix

## Check if the files exist before proceeding
if (!file.exists(pheno_file)) {
  stop(paste("Error: The phenotype file", pheno_file, "does not exist. Please provide a valid file."))
}

bimfile <- paste0(geno_file,".bim")
bedfile <- paste0(geno_file,".bed")
famfile <- paste0(geno_file,".fam")

if (!(file.exists(bimfile) && file.exists(bedfile) && file.exists(famfile))) {
  stop(paste("Error: The genotype files", geno_file, "do not exist. Please provide a valid file path and prefix."))
}


###############################################################################################################
#### Load phenotype
###############################################################################################################
## Load phenotype data
pheno_data <- readRDS(pheno_file)

## load fixed effects
fixedEffects <- strsplit(fixedEffects, ",")[[1]]  # Split covariates by comma without space, Ex. "Age,Sex,Education"

## Validate required columns in the phenotype file
required_cols <- c(id_col, time_col, y_col,fixedEffects)

missing_cols <- required_cols[!required_cols %in% colnames(pheno_data)]

if (length(missing_cols) > 0) {
  stop(paste("Error: The following required columns are missing from the phenotype file:", paste(missing_cols, collapse = ", ")))
}

message("=======================================================================")
message("Phenotype file is validated successfully!")
message("Pheno file: ", pheno_file)
message("ID Column: ", id_col)
message("Time Column: ", time_col)
message("Outcome Variable: ", y_col)
message("Fixed Effects for LMM: ", fixedEffects)
message("Number of observations in pheno data: ", nrow(pheno_data))
message("=======================================================================")
message(" ")
message(" ")

###############################################################################################################
#### Fit Conditional LMM without SNP
###############################################################################################################

##========================
## Step 0: prepare data
##========================

## subset phenotype file to remain only necessary columns
mydat <- pheno_data[,c(id_col,time_col,y_col)]
names(mydat) <- c("id","Time","y")
mydat$id <- as.character(mydat$id)

## sort based on individual and time
mydat <- mydat %>% arrange(id,Time)

## extract unique subject IDs
ids = unique(mydat$id)

##=================================================================================
## Step 1: data transformation: map the time-stationary part of the model to zero
##=================================================================================

transdata = NULL

for (i in ids) {
  xi = mydat[mydat$id == i, c("Time", "y")]  # Extract relevant columns
  if (nrow(xi) < 2) next # Skip subjects with fewer than 2 time points
  xi = as.matrix(xi)
  A = cumsum(rep(1, nrow(xi)))  # Sequence 1,2,3,...
  A1 = poly(A, degree = length(A) - 1)  # Orthogonal polynomials
  transxi = t(A1) %*% xi  # Apply transformation
  transxi = cbind(i, transxi)  # Add subject ID
  transdata = rbind(transdata, transxi)
}

transdata = as.data.frame(transdata)  # Convert matrix back to data frame
names(transdata) <- c("id","Time","y")
transdata$Time <- as.numeric(transdata$Time)
transdata$y <- as.numeric(transdata$y)

##=================================================================================
## Step 2: Fit the reduced mixed model and extract random slopes
##=================================================================================

mod2 = lmer(y ~ + Time - 1 + (Time-1|id), data = transdata, 
            control = lmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 100000)))

blups = ranef(mod2)$id # Extract random slopes (BLUPs)

##=================================================================================
## Step 3: prepare file for PLINK2
##=================================================================================

ids <- pheno_data[,c(id_col,"FID","IID")]

blups[[id_col]] <- rownames(blups)

blups <- merge(blups,ids,by = id_col)

blups <- blups[,c("FID","IID","Time")]

fwrite(blups,paste0(outpath,outprefix,"_CTSrandomSlope.pheno"),col.names = T, row.names = F, quote = F, sep = "\t")

###############################################################################################################
#### Association Test with PLINK2
###############################################################################################################

phenofile <- paste0(outpath,outprefix,"_CTSrandomSlope.pheno")

system(paste0("plink2 --bfile ",geno_file, " --linear allow-no-covars --pheno ",phenofile,
              "  --pheno-name Time --out ",outpath,outprefix,"_CTS"))

message("=======================================================================")
message(" ")
message(" ")

###############################################################################################################
#### Produce QQ plot and Manhattan plot
###############################################################################################################

message("Now generating the QQ plots and Manhatton Plots")

dat <- fread(paste0(outpath,outprefix,"_CTS.Time.glm.linear"))
dat <- as.data.frame(dat)

## calculate Lambda
z=qnorm(dat$P/2)
lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
message("The Lambda is ",lambda)
message(" ")

## QQ plot
png(filename = paste0(outpath,outprefix,"_CTS_QQ.png"),
    width = 8, height = 6, units = 'in', res = 300)
qqman::qq(dat$P)
dev.off()

message("QQ plot is located at: ",paste0(outpath,outprefix,"_CTS_QQ.png"))
message(" ")

## Manhattan Plot with APOE
png(filename = paste0(outpath,outprefix,"_CTS_Manhattan_withAPOE.png"),
    width = 8, height = 6, units = 'in', res = 300)
qqman::manhattan(dat, chr="#CHROM",bp="POS",p="P",snp="ID", annotateTop = T, annotatePval = 5e-8, col = c("darkgreen","darkblue"))
dev.off()

message("Manhatton plot with APOE region is located at: ",paste0(outpath,outprefix,"_CTS_Manhattan_withAPOE.png"))
message(" ")

## Manhattan Plot without APOE
APOE.snps <- dat$ID[dat$`#CHROM` == 19 & dat$POS>=44656625 & dat$POS<=45159250]
dat_noPAOE <- subset(dat,dat$ID %!in% APOE.snps)
png(filename = paste0(outpath,outprefix,"_CTS_Manhattan_noAPOE.png"),
    width = 8, height = 6, units = 'in', res = 300)
qqman::manhattan(dat_noPAOE, chr="#CHROM",bp="POS",p="P",snp="ID", annotateTop = T, annotatePval = 1e-6, col = c("darkgreen","darkblue"))
dev.off()
message("Manhatton plot without APOE region is located at: ",paste0(outpath,outprefix,"_CTS_Manhattan_noAPOE.png"))
message(" ")


message("=======================================================================")
message(" ")
message(" ")

###############################################################################################################
#### Extract SNPs with P-value < 1e-5 and perform LMM on those SNPs
###############################################################################################################
snplist <- data.frame(SNP=dat$ID[dat$P< 1e-5])

message("There are ", length(snplist), "SNPs has P value < 1e-5 from the CTS results")
message(" ")

snplistpath <- paste0(outpath,outprefix,"_CTS_P1e5.snps")
fwrite(snplist,snplistpath,col.names = F, row.names = F, quote = F, sep = "\t")

## Use PLINK2 subset the SNPs and convert to 0/1/2 format
system(paste0("plink2 --bfile ",geno_file," --extract ",snplistpath," --recode A --out ",outpath,outprefix,"_CTS_P1e5_SNPs"))

message("After subset the SNPs and convert to 0/1/2 format, the .raw file has been save to: ", paste0(outpath,outprefix,"_CTS_P1e5_SNPs.raw"))
message(" ")

## read subset file and merge with pheno file
genosub <- fread(paste0(outpath,outprefix,"_CTS_P1e5_SNPs.raw"))
genosub <- as.data.frame(genosub)
selected_cols <- names(genosub)[c(2, 7:ncol(genosub))]
genosub <- genosub[,selected_cols]

testsnps <- names(genosub)[c(2:ncol(genosub))]

## merge pheno_covs and genosub file
alldat <- merge(pheno_data,genosub,by = "IID")

LMM <- function(snp){
  # Construct the formula correctly
  if(length(fixedEffects) == 0){
    formula_str <- paste0(y_col,' ~ ',time_col,' * `', snp, 
                          '` + (1+',time_col,'|',id_col,')')
  }else{
    formula_str <- paste0(y_col,' ~ ',paste(fixedEffects,collapse = "+"),'+',time_col,' * `', snp, 
                          '` + (1+',time_col,'|',id_col,')')
  }
  
  model_formula <- as.formula(formula_str)
  
  # Fit the mixed-effects model
  mod <- lmer(model_formula, data = alldat, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  
  # Extract coefficients
  mod_coef <- data.frame(summary(mod)$coefficients)
  
  out <- mod_coef[nrow(mod_coef),c(1,2,4,5)]
  names(out) <- c("BETA","SE","T_STAT","P")
  
  return(out)
}

message("Now running the classic LMM for those selected SNPs....")
message(" ")

## run LMM for the test SNPs
testsnpsout <- do.call(rbind, lapply(testsnps, function(snp) {
  res <- LMM(snp)
  res$ID <- snp  # Add SNP name to the result
  return(res)
}))

## update the SNPs names 
testsnpsout$ID <- sapply(strsplit(testsnpsout$ID, "_"), `[`, 1)

## update the P-values for the test SNPs in the dat
idx <- match(dat$ID, testsnpsout$ID) # First, match the ID values from testsnpsout to dat
matched <- !is.na(idx) # Now, update only the rows that match (non-NA)

dat$BETA[matched]    <- testsnpsout$BETA[idx[matched]]
dat$SE[matched]      <- testsnpsout$SE[idx[matched]]
dat$T_STAT[matched]  <- testsnpsout$T_STAT[idx[matched]]
dat$P[matched]       <- testsnpsout$P[idx[matched]]


message("The update for the CTS results file for those SNPs is complete.")
message("=======================================================================")
message(" ")
message(" ")

###############################################################################################################
#### Prepare file for GWAMA
###############################################################################################################

## prepare data for GWAMA
## we need columns:  CHR SNP BP A1 TEST NMISS BETA SE L95 U95 STAT P
options(scipen = 999)

frqdat <- fread(frq_file,stringsAsFactors=F)
frqdat <- as.data.frame(frqdat)
frqdat <- frqdat[,c("SNP","A1","A2","MAF")]
names(frqdat) <- c("MARKERNAME","EA","NEA","EAF")
frqdat$STRAND <- "+"

dat <- dat[,c("ID","BETA","SE","OBS_CT")]
names(dat) <- c("MARKERNAME","BETA","SE","N")
dat$ORIGINAL_ORDER <- seq_len(nrow(dat))

dat <- merge(dat, frqdat, by = "MARKERNAME", all.x = TRUE)
dat <- dat[order(dat$ORIGINAL_ORDER), ] # Reorder based on the original order
dat <- dat[,c("MARKERNAME","EA","NEA","BETA","SE","N","EAF","STRAND")]

## save the updated association results
datpath <- paste0(outpath,outprefix,"EUR_ADNI_CTS_LMM_clean_GWAMA")
fwrite(dat,datpath,col.names = T, row.names = F, quote = F, sep = "\t")

message("The association results file was cleaned for GWAMA and file was saved to: ",datpath)
message(" ")

message("=======================================================================")
message(" ")
message(" ")
