# GeneDifferencialExpression
This code is an R script that loads RNA-Seq data and performs data manipulation, analysis, and visualization. Here is a step-by-step breakdown of the code:
Sets the working directory.
Loads several libraries including dplyr, tidyverse, GEOquery, ggplot2, and tibble.
Reads in the RNA-Seq data from a CSV file "DATA-1.csv".
Retrieves metadata from a Gene Expression Omnibus (GEO) database using the GEOquery package and extracts relevant information.
Renames and selects columns in the metadata.
Reshapes the RNA-Seq data into a long format.
Joins the metadata with the RNA-Seq data by a common column.
Writes the joined data to a CSV file "DATA-Long.csv".
Subsets the data into two separate data frames based on tissue type.
Calculates the mean of same FPKM (fragments per kilobase of transcript per million mapped reads) for each gene in the two separate data frames.
Separates the genes into four categories based on their FPKM values.
Eliminates the FPKM column from the four gene data frames.
Joins the columns that are repeated in the two separate data frames.
Calculates the overexpression and underexpression of genes in the two categories.
Writes the final data frames to CSV files "NUTO.csv" and "NOTU.csv".
Plots a bar graph using the ggplot2 package to visualize the overexpression and underexpression of genes in the two categories.
