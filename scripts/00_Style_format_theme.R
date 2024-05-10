# See below for theme
## Theme_format ----         
Style_format_theme <- theme(
  
  ## Title of the axis
  axis.title.x       = element_text(color="black", size=28),
  axis.title.y       = element_text(color="black", size=28),
  ## Text in the axis
  axis.text          = element_text(color="black", size=22),
  axis.text.x        = element_text(color="black", size=22),
  axis.text.y        = element_text(color="black", size=22),
  ## Line of the axis
  axis.line          = element_line(color="black", size= 1),
  axis.line.x        = element_line(color="black", size= 1),
  axis.line.y        = element_line(color="black", size= 1),
  ## Ticks of axis
  axis.ticks         = element_line(color="black", size= 1),
  axis.ticks.x       = element_line(color="black", size= 1),
  axis.ticks.y       = element_line(color="black", size= 1),
  axis.ticks.length  = unit(0.2,"cm"),
  ## Legend
  legend.title       = element_text(size=20),
  legend.text        = element_text(size=20),
  #legend.background  = element_rect( colour ="black", size=0.8, linetype="solid"),
  #legend.position    ="none",
  legend.position="right",
  #legend.position = c(0.6, 0.9),
  ## Panel
  panel.grid.major   = element_line(color="white"),
  panel.grid.minor   = element_line(color="white"),
  panel.background   = element_blank() ,
  ## Margins
  plot.margin       = unit(c(0.5,0.5,0.5,0.5), "cm"))
####