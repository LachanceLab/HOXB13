# ---- Install & load ----
install.packages("wtest")      # run once
library(wtest)

# ---- Load your data ----
# Expecting a CSV with one row per individual and columns:
#   case      : 0=control, 1=case
#   HOXB13_go : 0=no gain-of-function variant, 1=carrier of ANY GoF variant
#   RFX6_Tdos : number of T alleles at the SNP (0, 1, or 2)
# Replace "your_file.csv" or build the data.frame directly if you already have it in R.
dat <- read.csv("HOXB13RFX6Epistasis_CaP.csv", stringsAsFactors = FALSE)

# Basic checks / recoding
stopifnot(all(dat$case %in% c(0,1)))
stopifnot(all(dat$HOXB13_go %in% c(0,1)))
stopifnot(all(dat$RFX6_Tdos %in% c(0,1,2)))

# ---- Build the genotype matrix in the format wtest expects ----
# wtest wants a matrix/data.frame with genotype columns (0/1/2 or 0/1), subjects in rows.
# We'll name the columns with the locus IDs (these appear in the output).
geno <- data.frame(
  HOXB13_go = as.integer(dat$HOXB13_go),   # binary “any rare GoF variant present”
  RFX6_Tdos = as.integer(dat$RFX6_Tdos)    # 0,1,2 T-allele copies
)

# Phenotype vector (0/1)
y <- as.integer(dat$case)

# ---- Step 1: Estimate h,f parameters for interactions (w.order = 2) ----
# B is the bootstrap count for estimating the null distribution parameters.
# Authors recommend B >= 400; increase if you want more stable estimates.
set.seed(1)
hf2 <- hf(data = geno, w.order = 2, B = 400)

# (Optional) diagnostic check that the estimated distribution behaves well
# In big studies set n.rep higher (e.g., 200); this plots histograms vs. theoretical curves.
# w.diagnosis(geno, w.order = 2, n.rep = 100, hf2 = hf2)

# ---- Step 2: Run the pairwise epistasis test on your specific pair ----
# which.marker are the column indices of the two variables to test (here 1 and 2).
w_pair <- wtest(data = geno, y = y, w.order = 2, which.marker = c(1, 2), hf2 = hf2)

# Look at results
print(w_pair$results)

# Columns returned (for w.order = 2, single pair):
# 1: SNP1 name, 2: SNP2 name, 3: W-value (interaction), 4: k, 5: p-value (interaction)
# 6-8: main-effect stats for SNP1 (W, k, p)
# 9-11: main-effect stats for SNP2 (W, k, p)

# ---- (Optional) Odds ratios for interpretability ----
# Single-SNP ORs (main effects)
OR_hoxb13 <- odds.ratio(geno, y, w.order = 1, which.marker = 1)
OR_rfx6   <- odds.ratio(geno, y, w.order = 1, which.marker = 2)

# Logistic OR for the interaction pair (note: this is from a logistic model inside wtest helper)
OR_pair   <- odds.ratio(geno, y, w.order = 2, which.marker = c(1, 2))

OR_hoxb13; OR_rfx6; OR_pair

# ---- (Optional) Q-Q plot for interaction p-values (useful if you test many pairs) ----
# With only one pair this isn’t necessary, but for completeness:
# w.qqplot(geno, y, w.order = 2, hf2 = hf2)