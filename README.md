## Overview 
**Raw data and scripts for:**
Microbiota composition associates with mosquito productivity outcomes in belowground larval habitats

## Authors 
* Serena Y. Zhao
* Kerri L. Coon - kerri.coon@wisc.edu

## Analysis Overview 
Scripts for each analysis are written in R. Each directory contains necessary files and code to recreate each figure of the manuscript. To repeat the analysis, clone the repository, and the run each script. Do not `cd` into the cloned repository. 

**Example**
Create a new project in RStudio. To run the script to recreate Figure 1 in the manuscript: 
* Navigate to the terminal window and `git clone https://github.com/kcoonlab/catch-basin-microbes`.
* Open the script `Fig1.R` from the files panel window.
* Install required packages. 
* `cmd enter` from line `1`.

## Recreate the Manuscript Figures
Once the repository has been cloned (above), recreate each figure as follows: 

**Fig 1: Locations of collection sites**
* Fig 1 - Script: `04_Sankey_Diagram/Metadata_Snakey_Diagram.R`: run code from line `1` to `127`

**Fig 2: Bacterial diversity in water sampled from study catch basins**
* Fig 2A - Script: `03_MosAIC_Phylogeny/plot_tree_metadata.R`: run code from line `1` to `257`
* Fig 2B - Script: `01_Genome_QC/checkM_analysis.R`: run code from line `1 to `60`
* Fig 2C - Script: `01_Genome_QC/checkM_analysis.R`: run code from line `1` to `239`
* Fig 2D - Script: `01_Genome_QC/checkM_analysis.R`: run code from line `1` to `136`
* Fig 2E - Script: `01_Genome_QC/checkM_analysis.R`: run code from line `1` to `159`
* Fig 2F - Script: `02_GTDB_Drep_Summary/gtdbtk_drep_stat.R`: run code from line `1` to `112`

**Fig 3: Catch basin microbiota biotypes identified by PAM clustering**
* Fig 3 - Script: `05_Virulence_Factor_Analysis.R`: run code from line `1` to `192`

**Fig 4: Bacterial taxa significantly associated with different catch basin variables**
* Fig 4A-C - Script: `06b_EnterobacterPopulationStructure/Enterobacter_Pop_Struc.R`: run code from line `1` to `426`
* Fig 4D-F - Script: `06a_SerratiaPopulationStructure/Serratia_Genus_Pop_Structure.R`: run code from line `1` to `370`
* Fig 4G-H - Script: `06c_ElizabethkingiaPopulationStructure/ElizabethkingiaPopStruc.R`: run code from line `1` to `300`

**Fig 5: Bacterial community differences by pupal occurrence**
* Fig 5A - Script: `07b_EnterobacterPangenome/EnterobacterPangenomeTree.R`: run code from line `1` to `249` 
* Fig 5B - Script: `07a_SerratiaPangenome/SerratiaMPangenome.R`: run code from line `1` to `210`
* Fig 5C - Script: `07c_ElizabethkingiaPangenome/Elizabethkingia_Anophelis_Pangenome.R`: run code from line `1` to `177`

#### Supplementary Figures 
* Fig S1 - Script: `01_Genome_QC/plot_QUAST.R`: run code from line `1` to `55`

## Citation 


