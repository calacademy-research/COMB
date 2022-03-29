# Developing in COMB

## Guidelines

### Formatting

[styler](https://www.tidyverse.org/blog/2017/12/styler-1.0.0/) is an R code
formatter. For usage instructions follow the link.

Before making a pull request containing .R source files, it is strongly
recommended to style all that have changed under styler's default settings, as
this can meaningfully reduce the rate of merge conflicts due to formatting.
(Differences due to editor handling of tab expansion and trailing whitespace are
the most common.)

## Development Tools

### Jupyter Notebook

(Optional) If you would like to be able to run R .ipynb files from other
contributors or to author your own, you can run a local Jupyter notebook. The
IRkernel package is already installed in the project R environment by renv. For
the rest:

1. [Install Jupyter](https://jupyter.org/install).
2. Run the R script ```infrastructure/install_jupyter_r_kernel.R``` from the
   project R environment.
3. Run ```jupyter notebook``` or ```jupyter lab``` according to your preference
   and which are installed.
4. In your browser, find or create the .ipynb you are interested in.
