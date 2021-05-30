# In this code, we are going to practive analyzing files containing copy number variations.

# The dataset has 6 columns: (i) sample - represents the patient ID; (ii) chromosome represents the chromosome number; (iii) start -
# represents the start position of the segmented window; (iv) end - represents the end position of the segmented window;
# (v) num_probes - represents the number of probes in the segmented window; (vi) segment_mean - represents the mean copy number
# estimate of this particular segment.

# When segment_mean is greater than 0, there is amplification, and when it is less than 0, there is deletion.


# Loading the data
cnvLogs <- read.table(
    'data/BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt',
    header = TRUE,
    fill = TRUE
)
dim(cnvLogs)
head(cnvLogs)


# Data understanding
summary(cnvLogs)


# Copy number alterations of a single patient
# For a single patient, we are going to:
# 1. find how many potential copy number variations he/she has (just counting how many rows he/she has in this dataset)
# 2. find how many amplifications (segment_mean > 0) he/she has
# 3. find how many deletions (segment_mean < 0) he/she has
# 4. find the average segment_mean he/she has to get a global picture
patient <- "TCGA-A8-A06R-01A-11D-A011-01"
nrow(subset(cnvLogs, Sample == patient))
nrow(subset(cnvLogs, Sample == patient & Segment_Mean > 0))
nrow(subset(cnvLogs, Sample == patient & Segment_Mean < 0))
mean(subset(cnvLogs, Sample == patient)[["Segment_Mean"]])