# install_jupyter_r_kernel.R
#
# This simple script installs the R kernel for this project into a Jupyter
# environment that is assumed to have already been installed on the system.
#
# (It exists as a script rather than just documentation to trick renv into
# counting IRkernel as a dependency.)

library(IRkernel)

IRkernel::installspec()
