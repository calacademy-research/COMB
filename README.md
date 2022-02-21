# COMB
**C** ombined **O** ccupancy **M** odels for **B** irds

To sift through sound to identify birds and predict their occupancy!

## Repository Structure (draft):

**_Considerations on Mixed Computation for Collaboration_**

**Goal**: Collaboration with groups used to different computational tools and approaches.

**Issues for development and production**:

- Building our work in discrete '**chunks**' in our favorite programming platform is going to be fastest.
- Outputs can be shared as **text** or compressed text for interoperability.
- This will allow different programmers to collaborate on distinct stages of analysis.
- Althought development can occur in IDEs or in environments supporting multiple tools (e.g. https://jupyter.org/) ideally code could be executable by anyone at the command line.
- Data can be read in via code (the locations for data can be shared and documented _in the code_).
- We propose that each discrete computational **chunk** has its own directory.
- Underneath this directory, all code would exist in the \/src/ subdirectory.
- Then we would use code to move key data into a /input/ subdirectory for each **chunk**.
- Similarly, digested data that is necessary for subsequent **chunks** will be written to a /output/ subdirectory
- For R both /input/ AND /output/ directories are created by code and listed in .gitignore

This way we can adopt GitHub for _code_ and use Google Drive / Sheets for _data_ (Academy has infinite Google space).

Proposed structure

File or Directory     | Description
------------- | -------------
README.md     | This markdown document: start here for orientation to the project space
\/**chunk_x**\/ | coding each discrete '**chunk or step**' in the \/**chunk_x**\/ directories as follows:
..\/src           | all code: R, Python, Awk, etc. scripts for reading in, cleaning, exploring, and analyzing data
..\/input         | all raw data retrieved from Google Drive (large files) or Google Sheets ('by hand' metadata).
..\/output        | all digested data (for next steps) as well as figures, tables, etc. for reports, manuscripts, presentations
..\/note          | all notebook analyses (jupyter, R-notebooks, Markdown docs, metadata summaries)
..\/hand          | all by-hand step descriptions how to reproduce  (readme\_by\_hand.txt, with links to tools)
..\/models        | copies of .txt files of occupancy models (used in JAGS)  not this can be moved to \\manuscript\ eventually
\/**chunk_y**\/  | ...

This proposal borrows shamelessly from the Human Rights Data Analysis Group's Patrick Ball. If you have time, you can watch his [YouTube video](https://www.youtube.com/watch?v=ZSunU9GQdcI) on 'Principled Data Processing', where he explains his rationale on how to organize computational work for ‘self-documenting’ reproducible science.

Whatever we decide on, it will be great if our pipeline process self-documents with 'stepwise' inputs, code and outputs in the directory hierarchy and intermediate files are saved and automatically archived and work can be restarted from any step.

Eventually we may combine these in a shared computing environment that everyone understands, perhaps using JupyterLab for notebooks that can collate our results.  If we do that, we might want to create a stable set of tools for a mixed R/python environment using Conda so we are all on the same page.

Thoughts welcome, meanwhile we will start to populate the git with some of our code that _points to data_ and _processes it_.
