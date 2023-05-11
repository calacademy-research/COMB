# convert_ipynb_2_rmarkdown.R

library(rmarkdown)
library(here)

nb_file = here("spatial","note","SALO_DIRECT.ipynb")

# convert to R Markdown
nb_rmd <- rmarkdown:::convert_ipynb(nb_file)
xfun::file_string(nb_rmd)
