# set working directory
setwd <- ("C:/Users/rofia/OneDrive/Desktop/Alignment")

-----------------------------------------------
# STEP 1: Install and load the required packages
-----------------------------------------------
install.packages("rentrez")
install.packages("msa")
install.packages("ggplot2")
install.packages("ggmsa")
install.packages("pdflatex")
install.packages("Biostrings")

library(rentrez)
library(msa)
library(ggplot2)
library(ggmmsa)
library(pdflatex)
library(Biostrings)
library(pdflatex)


-------------------------------------------
# STEP 2a: create a list of the identifiers
-------------------------------------------
  # Multiple identifiers separated by commas
ids <- c("ACD88516.1", "AJW32121.1", "ADC45736.1", "AAG01758.1", "ACZ92133.1",
         "ACS92895.1", "BCZ08606.1", "BCZ08605.1", "BCZ08604.1", "AAG01785.2",
         "AFH00351.1", "ABG80194.1", "ACD88520.1", "AAG01767.1", "AAK63824.1",
         "AAK63822.1", "BCZ08611.1", "BCZ08610.1", "BCZ08609.1", "AAA16879.1",
         "BCZ08615.1", "BCZ08614.1", "BCZ08613.1", "BCZ08612.1", "AAG01776.1",
         "BCZ08608.1", "BCZ08607.1")

-------------------------------------------------------------
# STEP 2b: download the protein sequences from Entrez database
--------------------------------------------------------------
sequences <- entrez_fetch(db = "protein", id = ids, rettype = "fasta")

# Save the sequences to a file
writeLines(sequences, "C:/Users/rofia/OneDrive/Desktop/Alignment/multipleSequences.fasta")

----------
# STEP 2c: assign the fasta sequences to an object
---------  
mySequenceFile <- "C:/Users/rofia/OneDrive/Desktop/Alignment/multipleSequences.fasta"
mySequences <- readAAStringSet(mySequenceFile)
mySequences


-------------------------------------
# STEP 3a: perform Alignment with msa
-------------------------------------
myFirstAlignment <- msa(mySequences)
myFirstAlignment

---------------------------------------
# STEP 3b: identify the masked aligment
---------------------------------------
myMaskedAlignment <- myFirstAlignment
colM <- IRanges(start=1, end=100)
colmask(myMaskedAlignment) <- colM
myMaskedAlignment


-------------------------------------
# STEP 3c: identify unmasked alignment
-------------------------------------
unmasked(myMaskedAlignment)


------------------------------------
# STEP 3d: calculate consensus matrix
------------------------------------
  conMat <- consensusMatrix(myFirstAlignment)
dim(conMat)
conMat[, 101:110]


---------------------------
# STEP 3e: export aligmnent
--------------------------
# Check if pdflatex is accessible
pdflatex_path <- "C:\\PROGRA~1\\MiKTeX\\miktex\\bin\\x64\\pdflatex.exe"
if (file.exists(pdflatex_path)) {
  print("pdflatex is accessible.")
} else {
  stop("pdflatex is not accessible. Please check the path and TeX installation.")
}

# Specify the file path for the output PDF directly
tmpFile <- "C:/Users/rofia/OneDrive/Desktop/alignment/MSA_Alignment.pdf"

# Print the temporary file path
print(tmpFile)

# Run the msaPrettyPrint function with the specified file path and verbose set to TRUE
msaPrettyPrint(myFirstAlignment, file = tmpFile, output = "pdf",
               showNames = "left", showNumbering = "none", showLogo = "top",
               showConsensus = "bottom", logoColors = "rasmol",
               verbose = TRUE, askForOverwrite = FALSE)


---------------------------------------------------
# STEP 4: calculate conservation score with BLOSUM62
---------------------------------------------------
conservationScore<- msaConservationScore(myFirstAlignment, BLOSUM62, gapVsGap=0,
                                         type="upperlower", thresh=c(50, 20))

# create a file path to save the conservation scores
file_path <- "C:/Users/rofia/OneDrive/Desktop/Alignment/msaConservation_BLOSUM62.csv"

# Extract conservation scores and corresponding regions for the entire alignment
positions <- 1:length(conservationScore)
conservation_scores <- conservationScore
conserved_regions <- which(conservationScore >= 50)

# Combine the data into a data frame
# We will use the length of the shortest vector to create the data frame
length_shortest <- min(length(positions), length(conservation_scores), length(conserved_regions))
conservation_data <- data.frame(Position = positions[1:length_shortest], 
                                Conservation_Score = conservation_scores[1:length_shortest], 
                                Conserved_Region = conserved_regions[1:length_shortest])

# Write the data frame to an Excel file 
write.csv(conservation_data, file = "conservation_result.csv", row.names = FALSE)



---------------------------------------
# STEP 5a: plot the conservation metric
---------------------------------------
# Threshold for highly conserved residues
  threshold <- 50

# Plot the conservation metric as a line graph with a horizontal line for the threshold
ggplot(conservation_data, aes(x = Position, y = Conservation_Score)) +
  geom_line(color = "steelblue", size = 0.5) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  labs(x = "Position", y = "Conservation Score", title = "Conservation Metric of Protein Alignment") +
  theme_minimal()
ggsave("conservation_plot.pdf", width = 10, height = 8)


-------------------------------------------------------------------
# STEP 5b: plot a graph of the alignment with the conservation score
-------------------------------------------------------------------
protein_sequences <- "C:/Users/rofia/AppData/Local/Temp/RtmpYn1tsK/seq1df46e1a564a.fasta"
alignment_plot <- ggmsa(protein_sequences, start = 1, end = 600, char_width = 1.5, seq_name = T) + geom_seqlogo() + geom_msaBar()

# Save the plot to a PDF file
ggsave(output_file, plot = alignment_plot, device = "pdf", width = 10, height = 6)

# Print a message to confirm that the plot is saved
cat("The plot has been saved to:", output_file, "\n")


------------------------------------
# STEP 5c: plot the conservation ratio
------------------------------------
# Calculate the total number of positions in the alignment
total_positions <- length(conservation_data$Position)
# Calculate the number of highly conserved positions (score >= 50)
highly_conserved_positions <- sum(conservation_data$Conservation_Score >= 50)
# Calculate the number of less conserved positions (score < 50)
less_conserved_positions <- sum(conservation_data$Conservation_Score < 50)
# Calculate the percentage of highly conserved positions
percentage_highly_conserved <- (highly_conserved_positions / total_positions) * 100
# Calculate the percentage of less conserved positions
percentage_less_conserved <- (less_conserved_positions / total_positions) * 100
# Create a data frame for the histogram
hist_data <- data.frame(Category = c("Highly Conserved", "Less Conserved"),
                        Percentage = c(percentage_highly_conserved, percentage_less_conserved))

# Plot the combined histogram
ggplot(hist_data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", color = "black", width = 0.5) +
  labs(x = "Conservation Category", y = "Percentage", title = "Percentage of Conservation Categories") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



-----------------------------------------------------------
# STEP 5d: calculate conservation score with shannon entropy
------------------------------------------------------------
# Convert the msa object to a matrix
alignment_matrix <- as.matrix(myFirstAlignment)

# Calculate Shannon Entropy for each position in the alignment
shannon_entropy <- apply(alignment_matrix, 2, function(column) {
  column_freq <- table(column)
  column_prob <- column_freq / sum(column_freq)
  -sum(column_prob * log2(column_prob))
})
  

--------------------------------------
# STEP 5e: plot the conservation score
--------------------------------------
# Create a data frame for the Shannon entropy values and positions
entropy_data <- data.frame(Position = 1:length(shannon_entropy),
                           Shannon_Entropy = shannon_entropy)

# Plot the Shannon entropy as a line graph
ggplot(entropy_data, aes(x = Position, y = Shannon_Entropy)) +
  geom_line(color = "black", size = 0.5) +
  labs(x = "Position", y = "Shannon Entropy", title = "Shannon Entropy of Protein Alignment") +
  theme_minimal()

  