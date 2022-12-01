#quick table of vegetation and bird points by burn severity

#depends upon veg_plot_salo_comp.R

#table the variables

Caples_PlotData_20220225 %>%
  dplyr::filter(SamplingTime == 0) %>% 
  dplyr::select(Caples_Severity_Class) %>% 
  count(Caples_Severity_Class)
