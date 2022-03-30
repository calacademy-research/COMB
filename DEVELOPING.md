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

### renv

#### First-time setup

[renv](https://rstudio.github.io/renv/articles/renv.html) is a dependency
management system for R. Using it is optional but highly recommended as a way of
keeping all contributors' package versions in sync. There is not much to do to
use it. When you run R from the project directory, you will notice a prompt
suggesting

```renv::restore()```

to sync with the lockfile. Running this command will install the
lockfile-specified package versions into renv/library. This can take about half
an hour the first time, but after that R will start up as fast as before, and
syncs as new packages are added will be quicker.

#### Adding packages

If you add a library() line for a new package during development and find that
you needed to install a package, please also run ```renv::snapshot()``` and
include the changes to the lockfile in the commit that adds the library()
statement.

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