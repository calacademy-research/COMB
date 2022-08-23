
# 13 April 2022, Sarah Jacobs
# 
# This README describes steps taken to summarize, transform, and write raw data to usable format for downstream visualization and analysis.
# This script begins with six .xlsx files exported from an access database, managed by USFS, that holds annual plot visitation data for the Caples Watershed.
# Various summary statistics and transformations are performed and unused variables dropped. Data are then combined into a single dataframe and written to a .csv file
# The .csv output file serves as the input for the Shiny app (COMB/viz/app.R) to explore and visualize the plot based variables. 
# NOTE: if using this script to generate a new .csv file -- for example, if new data are incorporated into the dataset -- plese refer to the Caples_Access_readme to ensure data are exported correctly AND make sure any newly exported datasheets from the Access database are appropriately named. Previous protocol has appended dates of export to file names (YYYYMMDD)
