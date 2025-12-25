Code to go from classifier outputs to structured data that we can feed into an occupancy model. 

There are three main inputs into the combined model: ARU, PC, and spatial. For each data source, we provide a polars dataframe containing the data along with parameters that specify the columns to use. 

The output is a dictionary that maps the name of the data (burn, y_pc, y_aru, score, ...) to the data. The values in the output are either an integer (e.g. number of sites) or a numpy array (e.g. y_aru).

The shapes of the output data are all (species, year, site, visit), except in the case of covariates where we remove irrelevant dimensions but keep the order. So covariates have shape (year, site). 

The script `classifier_to_aru.py` is used to take the outputs from the `perch-agile-modeling` and put them into a format suitable for the combining class. 

The only classes users should use are `CombinedParams` and `CombinedData`, as this class handles each data type (ARU, PC, spatial) behind the scenes. 