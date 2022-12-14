# Dengue_human_movement_model
 
Data and `R` code to support Lee S.A. *et al.*, (in prep.). The contribution of human movement to dengue expansion differs between regions in Brazil.


--------------------------------------------------------------------------------

Dengue transmission has been expanding across Brazil since its re-introduction into areas that were previously thought to be protected. Previous studies showed that geographical barriers to dengue transmission have been eroded in South Brazil and the western Amazon. In this study, we quantify the relative contribution of human movement to the geographical expansion of dengue outbreaks between 2001 and 2020 in Brazil using a Bayesian hierarchical model. 

--------------------------------------------------------------------------------

A description of each file and folder is provided below:

  **R/00_human_movement_collate.R:** `R` script to clean and explore commuting data from the 2010 Brazilian census. Note that this code requires the Brazilian municipality shapefile which can be downloaded from the IBGE website: https://www.ibge.gov.br/geociencias/organizacao-do-territorio/malhas-territoriais/15774-malhas.html?edicao=27415&t=acesso-ao-produto

  **01_run_model.R:** `R` script to fit spatial models to the data.

  **02_model_outputt.R:** `R` script to explore and visualise model output.
  
  **Data:** a folder containing databases necessary to create the visualisations provided in the manuscript and fit the spatial models.
  
  **Output:** a folder to save the output generated by the `R` scripts.

The analysis was performed using R version 4.2.0 (2022-04-22).

For any issues with the code please contact [Sophie Lee](https://www.lshtm.ac.uk/aboutus/people/lee.sophie).

