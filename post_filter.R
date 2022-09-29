cols <- colnames(dataML)[16:104]

filter_values <- paste(cols, ">-2", "| ") %>% 
  paste(collapse = "") %>% 
  str_sub(end = -4)

small_dataML <- dataML %>% filter(rlang::eval_tidy(rlang::parse_expr(filter_values)))

dataML %>%
  filter(dusfly > -.5) %>% 
  ggplot(aes(x=1, y= dusfly))+
  geom_violin()
