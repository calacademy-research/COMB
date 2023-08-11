library(rjags)
# install.packages("MCMCvis")
library(MCMCvis)


MCMCsummary(jagsResult)

MCMCsummary(jagsResult, 
            params = 'mu', 
            Rhat = TRUE, 
            n.eff = TRUE, 
            round = 2)

MCMCpstr(jagsResult, 
         params = 'mu', 
         func = mean,
         type = 'summary')
