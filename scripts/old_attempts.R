######################old attempts
#step 1: load the count matrix
countdata = read.csv("C:/R_projects/bulk RNA seq/data/counts.csv", row.names = 1, check.names = FALSE)
head(countdata)

countdata = read.delim("C:/R_Projects/Bulk RNA seq/data/counts.csv", row.names = 1, sep = "\t", check.names = FALSE)

head(countdata)
colnames(countdata)[1:5]  # Show first 5 sample names
#it showed that the column names don't just contain sample names, but other info as ful identifiers, so we will need to match and extract the core sample IDs from the meta data (series matrix)

#step 2: load the seriesmatrix
# Read the raw metadata 
meta_raw <- readLines("C:/R_Projects/Bulk RNA seq/data/GSE303918_series_matrix.txt")

# Step 2: Extract all lines that start with "!Sample_characteristics_ch1"
characteristics <- grep("!Sample_characteristics_ch1", meta_raw, value = TRUE)

# Step 3: Remove the label and split each line into fields
characteristics_clean <- lapply(characteristics, function(line) {
  fields <- strsplit(sub("!Sample_characteristics_ch1\t", "", line), "\t")[[1]]
  return(fields)
})

# Step 4: Transpose the list into a data frame
meta_matrix <- as.data.frame(t(do.call(rbind, characteristics_clean)), stringsAsFactors = FALSE)

# Step 5: Add sample IDs
sample_ids <- grep("!Sample_geo_accession", meta_raw, value = TRUE)
sample_ids <- strsplit(sub("!Sample_geo_accession\t", "", sample_ids), "\t")[[1]]
meta_matrix$sample_id <- sample_ids

dim(meta_matrix)
head(meta_matrix)

colnames(meta_matrix)

#filter for SOX1
sox1_ko_rows <- grepl("SOX1", meta_matrix$V4) & grepl("day 14", meta_matrix$V6)
sox1_ko <- meta_matrix[sox1_ko_rows, ]

sox1_ko$sample_id
nrow(sox1_ko)
#result showed that this dataset has no SOX1 sampes at day 14
#to check which TF appear at day 14
# Get all samples at day 14
day14_rows <- grepl("day 14", meta_matrix$V6)
day14_samples <- meta_matrix[day14_rows, ]

# See which TFs were knocked out
table(day14_samples$V4)
