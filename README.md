# PhalFNP
Analyses of the _Panicum hallii_ fast neutron mutagenesis population

Access the interactive depth images at the following link:
https://nolanbentley.github.io/PhalFNP/index.html

# Documenting the process
## Accessing the input files
These sequence mapping files and more originate from mapping reads to the _Panicum hallii_ V3.2 houased on the Phytozome (https://phytozome-next.jgi.doe.gov/info/Phallii_v3_2) database.

These files may be accessed via the JGI genome portal (https://genome.jgi.doe.gov/portal/ResaPanihallii/ResaPanihallii.info.html). 

This repo originally focused on the analysis of the bam files from this project (see data/bamFilesList.txt) as well as downstream variant files produced through the analysis of GATK and DELLY of collaborators. These downstream files are available upon request.

## Downloading the bam files
bam files were downloaded via JGI's curl-based API. This process was controlled by a modified version of rScriptToControlJamo.R. 

## Calculating depths
bam files were then analyzed across a series of genomic intervals to calculate coverage statistics using samtools depth with commands controlled by bamToDepth_FNP_v1.1.R.

Input: A full listing of the locally stored files can be found in data/bamFilesInDepth.csv

## Summarize the depth files to create the "preliminary" aggregate dataset
This script (depthToPrelimAggregate.R) looped across each interval and sample combination in order to generate a summary of each interval.This was valuable as it massively reduced the computational burden by  simplifying the genome into 1 record per sample and interval location. 

Input: A full listing of the locally stored files analyzed can be found in the examples directory.

## Convert the coverage to copy number and find trends across adjacent intervals
This script (prelimToWindowedCopyNumber.R) was used to filter and analyze the intervals in order to correlate coverage with copy number. The calculation used here is that the mean coverage per interval * 2 / peak coverage per interval should estimate the number of copies across that region. 

The principal assumption of this analysis is that, for the majority of the genome, these samples are diploid. This should be generally true for our samples, although technically, if any samples had whole genome duplication, they would not necessarily be detected through this method as the most frequent coverage would no longer correspond to a diploid state. 

After transformation, the median CN across 5 interval windows was calculated to generate the primary statistic used in downstream analyses.

For intervals where quality control dictated the removal of that locus and/or record, the next closest interval was used. The first and last two intervals were removed. 

Input: An example preliminary interval analysis file (prelimAggregateDepths_TIMESTAMP.csv) of this script can be found in the examples directory.

## Produce plots visualizing copy-number variants
This script (windowedCNToPlots.R) produced a series of ggplot2 -> plotly -> htmlwidget derived interactive html files that visualize the windowed copy number analysis for easy exploration. The current URL for this tool is  https://nolanbentley.github.io/PhalFNP/index.html 

Inputs: 
- An example windowed copy number file (windowedNs5x.csv) of this script can be found in the examples directory.
- The geneDf.csv file used to draw the gene's on single-chromosome plots can be found at data/geneDf.csv

## Produce html files to facilitate access to interative plots
This script produced a series of html files that link together the various plots into an easy to explore format. 

Input: An example windowed-CN plotted analysis data.frame (plottedAggDf_n.csv) can be found in the examples directory.






