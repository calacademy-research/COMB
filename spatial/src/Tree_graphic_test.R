library(ggplot2)
library(data.table)

dt.triangle <- data.table(group = c(1,1,1), polygon.x = c(1,2,1.75), polygon.y = c(1,1,2))
p <- ggplot()
p <- p + geom_polygon(
  data = dt.triangle
  ,aes(
    x=polygon.x
    ,y=polygon.y
    ,group=group
  )
)
p


library(ggplot2)

col_col <- c("#000000",'#000000','#000000')
col_fill <- c("#5cb85c","#f9f9f9","#f9f9f9")
d=data.frame(x=c(1,2,2, 1.5,1.5,2,2), y=c(1,1,2, 1.375,1.25,1.25,1.75), t=c('a', 'a', 'a',  'b', 'b', 'b', 'b'), r=c('x','z','y', 4,5,6,7))
p <- ggplot(data = d, aes(x = x, y = y, col = factor(t), fill = factor(t))) + geom_polygon(data = d, alpha = .75) + # geom_point() +
  scale_color_manual(values = col_col) + scale_fill_manual(values = col_fill) 
p + geom_point(data = d[1:3,]) + geom_text(data = d[1:3,], aes(x=x, y=y, label=r), hjust=0, vjust=1, size=4) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none')


library(ggplot2)

col_col <- c("#000000",'#000000','#000000')
col_fill <- c("#964B00","#808080","#5cb85c")

d=data.frame(x=c(0,0,2,0,0,1,0,1,1,0), y=c(1,5,5,0,1,1,0,0,-.25,-.25), cat=c('t', 't', 't','f','f','f','s','s','s','s'), r=c('bh','ch','cc','gnd','bh','lf','gnd','sf','',''))
p <- ggplot(data = d, aes(x = x, y = y, col = factor(cat), fill = factor(cat))) + geom_polygon(data = d, alpha = .75) + # geom_point() +
  scale_color_manual(values = col_col) + scale_fill_manual(values = col_fill) 
p + geom_point(data = d[1:9,]) + geom_text(data = d[1:9,], aes(x=x, y=y, label=r), hjust=0, vjust=1, size=4) +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none')

