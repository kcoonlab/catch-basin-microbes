## Overview 
**Raw data and scripts for:**
Microbiota composition associates with mosquito productivity outcomes in belowground larval habitats

## Authors 
* Serena Y. Zhao
* Kerri L. Coon - kerri.coon@wisc.edu

## Analysis overview 
Scripts for each analysis are written in R. Each directory contains necessary files and code to recreate each figure, table, and associated statistical analyses reported in the manuscript. To repeat the analysis, clone the repository, and then run each script. Do not `cd` into the cloned repository. 

**Example**
Create a new project in RStudio. To run the script to recreate Fig. 1 in the manuscript: 
* Navigate to the terminal window and `git clone https://github.com/kcoonlab/catch-basin-microbes`.
* Open the script `Fig1.R` from the files panel window.
* Install required packages. 
* `cmd enter` from line `1`.

**Before getting started**
* Run the script `phyloseq-object.R` (code from line `1` to `39`) to generate a phyloseq object from the appropriate qiime artifacts. This object will be referenced in some downstream analyses.
* Note: this script also contains the code used to estimate alpha diversity (ASV richness and Shannon's H index) in each sample (see lines `41` to `47`)

## Recreate the manuscript figures, tables, and associated statistical analyses
Once the repository has been cloned (above), recreate each figure/table/data analysis as follows: 

**Fig. 1: Locations of collection sites**
* Fig 1A - Script: `Fig1.R`: run code from line `1` to `15`
* Fig 1B - Script: `Fig1.R`: run code from line `1` to `27`

**Fig. 2: Bacterial diversity in water sampled from study catch basins**
* Fig 2A - Script: `Fig2.R`: run code from line `1` to `29`
* Fig 2B - Script: `Fig2.R`: run code from line `1` to `50`
* Fig 2C - Script: `Fig2.R`: run code from line `1` to `69`
* Fig 2D - Script: `Fig2.R`: run code from line `1` to `88`

**Fig. 3: Catch basin microbiota biotypes identified by PAM clustering**
* Fig 3A - Script: `Fig3.R`: run code from line `1` to `33`
* Fig 3B - Script: `Fig3.R`: run code from line `1` to `54`
* Fig 3C - Script: `Fig3.R`: run code from line `1` to `80`
* Fig 3D - Script: `Fig3.R`: run code from line `1` to `103`
* Fig 3E - Script: `Fig3.R`: run code from line `1` to `135`
* Fig 3F - Script: `Fig3.R`: run code from line `1` to `166`

**Fig. 4: Bacterial taxa significantly associated with different catch basin variables**
* Fig 4 - Script: `Fig4.R`: run code from line `1` to `193`

**Fig. 5: Bacterial community differences by pupal occurrence**
* Fig 5A - Script: `Fig5.R`: run code from line `1` to `38`
* Fig 5B - Script: `Fig5.R`: run code from line `1` to `60`
* Fig 5C - Script: `Fig5.R`: run code from line `1` to `82`
* Fig 5D - Script: `Fig5.R`: run code from line `1` to `104`
* Fig 5E - Script: `Fig5.R`: run code from line `1` to `109`
* Fig 5F - Script: `Fig5.R`: run code from line `1` to `166`

**Table 1: Effects of water quality on mosquito productivity**
* Table 1A - Script: `Table1.R`: run code from line `1` to `53`
* Table 1B - Script: `Table1.R`: run code from line `1` to `102`
* Table 1C - Script: `Table1.R`: run code from line `1` to `152`

**Table 2: Effects of sampling date, water quality, and mosquito productivity on microbiota diversity**
* Table 2A - Script: `Table2.R`: run code from line `1` to `67`
* Table 2B - Script: `Table2.R`: run code from line `1` to `85`
  
**Table 3: Effects of microbiota diversity on mosquito productivity**
* Table 3 - Script: `Table3.R`: run code from line `1` to `26`

#### Miscellaneous data analyses
* Script: `misc-analyses.R`: run code from line `1` to `227`

#### Supplementary information 
* Fig. S1 - Script: `supp-info.R`: run code from line `1` to `31`
* Fig. S2A - Script: `supp-info.R`: run code from line `1` to `61`
* Fig. S2B - Script: `supp-info.R`: run code from line `1` to `80`
* Fig. S2C - Script: `supp-info.R`: run code from line `1` to `98`
* Fig. S3 - Script: `supp-info.R`: run code from line `1` to `109`
* Fig. S4A - Script: `supp-info.R`: run code from line `1` to `128`
* Fig. S4B - Script: `supp-info.R`: run code from line `1` to `143`
* Fig. S5 - Script: `supp-info.R`: run code from line `1` to `185`
* Fig. S6 - Script: `supp-info.R`: run code from line `1` to `260`
* Table S1A - Script: `supp-info.R`: run code from line `1` to `308`
* Table S1B - Script: `supp-info.R`: run code from line `1` to `347`
* Table S2 - File: `input-files/metadata.txt`
* Table S3A - Script: `supp-info.R`: run code from line `1` to `424`
* Table S3B - Script: `supp-info.R`: run code from line `1` to `445`
* Table S4A - Script: `supp-info.R`: run code from line `1` to `460`
* Table S4B - Script: `supp-info.R`: run code from line `1` to `473`
* Table S4C - Script: `supp-info.R`: run code from line `1` to `489`
* Table S5A - Script: `supp-info.R`: run code from line `1` to `510`
* Table S5B - Script: `supp-info.R`: run code from line `1` to `529`
* Table S6 - Script: `supp-info.R`: run code from line `1` to `576`
* Table S7A - Script: `Table2.R`: run code from line `1` to `67`
* Table S7B - Script: `Table2.R`: run code from line `1` to `85`

## Citation 
Zhao SY, Sommer AJ, Bartlett D, Harbison JE, Irwin P, Coon KL. Microbiota composition associates with mosquito productivity outcomes in belowground larval habitats.

