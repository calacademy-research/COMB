# COMB
**C** ombined **O** ccupancy **M** odels for **B** irds

## Repository Structure (draft): 

**_Considerations on Mixed Computation for Collaboration_**

**Goal**: Collaboration with groups used to different computational tools and approaches.

**Issues for development and production**:

- Building our work in discrete '**chunks**' in our favorite programming platform is going to be fastest.
- Outputs can be shared as **text** or compressed text for interoperability. 
- This will allow different programmers to collaborate on distinct stages of analysis.
- Althought development can occur in IDEs or in environments supporting multiple tools (e.g. https://jupyter.org/) ideally code could be executable by anyone at the command line.
- Data can be read in via code (the locations for data can be shared and documented _in the code_).
- We propose placing data read in via code in your local /input/ directory
- And digested data in your local /output/ directory
- Both directories are listed in .gitignore

Solution:  Adopt GitHub for _code_ and Google Drive / Sheets for _data_ (Academy has infinite Google space).

Proposed structure

File or Directory     | Description
------------- | -------------
README.md     | This markdown document: start here for orientation to the project space
\/**chunk_x**\/ | coding each discrete '**chunk or step**' in the \/**chunk_x**\/ directories as follows:
..\/input         | all raw data retrieved from Google Drive (large files) or Google Sheets ('by hand' metadata).
..\/src           | all code: R, Python, Awk, etc. scripts for reading in, cleaning, exploring, and analyzing data
..\/output        | all digested data (for next steps) as well as figures, tables, etc. for reports, manuscripts, presentations
..\/note          | all notebook analyses (jupyter, R-notebooks, Markdown docs, metadata summaries)
..\/hand          | all by-hand step descriptions how to reproduce  (readme_by_hand.txt, with links to tools)
..\/models        | copies of .txt files of occupancy models (used in JAGS)  not this can be moved to \\manuscript\ eventually
\/**chunk_y**\/  | ...

This proposal borrows shamelessly from the Human Rights Data Analysis Group's Patrick Ball. If you have time, you can see his explanation of 'Principled Data Processing' here:

https://www.youtube.com/watch?v=ZSunU9GQdcI

Where he explains his rationale on how to organize computational work for ‘self-documenting’ reproducible science.

Whatever we decide on, it will be great if our pipeline process self-documents with 'stepwise' inputs, code and outputs in the directory hierarchy and intermediate files are saved and automatically archived and work can be restarted from any step.  

Eventually we may combine these in a shared computing environment that everyone understands, perhaps using JupyterLab for notebooks that can collate our results.  If we do that, we might want to create a stable set of tools for a mixed R/python environment using Conda so we are all on the same page.
