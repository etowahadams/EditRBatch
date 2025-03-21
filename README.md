# EditRBatch
[![DOI](https://zenodo.org/badge/594460101.svg)](https://doi.org/10.5281/zenodo.15060671)

This is a fork of [EditR](https://github.com/MoriarityLab/EditR) that can analyze multiple Sanger sequencing samples at a time.

## System Requirements

### Software Dependencies and Operating Systems
- R (latest version recommended)
- Dependencies installed via `dependencies.R`. The specific versions can be found in `manifest.json`. 

### Versions Tested On
- R version 4.2.3, using RStudio (Version 2023.03.0+386) on OSX Sonoma. 

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
2. Open the folder in your R environment, and install dependencies by running:
   ```r
   source("dependencies.R")
   ```
3. Run Shiny app by executing:
   ```r
   library(shiny)
   runApp()
   ```

### Typical Install Time
- Less than 10 minutes on a standard desktop computer.

### Using the Shiny App

1. Create a metadata .xlsx file which has the following columns:
    - `Sample ID`: This should match the name of your .ab1 file.
    - `gRNA guide sequence`: The gRNA guide sequence used in your experiment.
    - `Reverse Y/N`: Y or N. Whether the sample was reverse sequenced. If so, the reverse complement of the sequence will be used for analysis.
2. In the Shiny app, select the .ab1 files you wish to analyze and the metadata file you created.
3. After the files are uploaded, the app will automatically analyze the samples and output the results in a table.

### Demo

In the `demo_data` folder, we have provided a sample metadata file and .ab1 files for you to test the app with.

If you select the ab1 files in the `demo_data/ab1_files` folder and the `demo_data/230105_ExampleData_Key.xlsx` file, you should see the following output, which shows the results of the EditR analysis for each sample:

![Demo Output](image.png)

It should take less than 5 minutes generate the output. 

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.
We retain the original license from the EditR project.