# EditRBatch

This is a fork of [EditR](https://github.com/MoriarityLab/EditR) that can analyze multiple Sanger sequencing samples at a time.

## System Requirements

### Software Dependencies and Operating Systems
- R (latest version recommended)
- Dependencies installed via `dependencies.R`. The specific versions can be found in `manifest.json`. 

### Versions Tested On
- R version 4.2.3 on macOS

### Required Non-Standard Hardware
- No special hardware required beyond a standard desktop computer.

## Installation Guide

### Instructions
1. Ensure R is installed on your system.
2. Clone the repository:
   ```sh
   git clone https://github.com/your-repo/EditRBatch.git
   cd EditRBatch
   ```
2. Install dependencies by running:
   ```r
   source("dependencies.R")
   ```
3. Run Shiny app by executing:
   ```r
   library(shiny)
   runApp()
   ```

### Typical Install Time
- Less than 5 minutes on a standard desktop computer.

### Using the Shiny App

1. Create a metadata .xlsx file which has the following columns:
    - `Sample ID`: This should match the name of your .ab1 file.
    - `gRNA guide sequence`: The gRNA guide sequence used in your experiment.
    - `Reverse Y/N`: Y or N. Whether the sample was reverse sequenced. If so, the reverse complement of the sequence will be used for analysis.
2. In the Shiny app, select the .ab1 files you wish to analyze and the metadata file you created.
3. After the files are uploaded, the app will automatically analyze the samples and output the results in a table.
