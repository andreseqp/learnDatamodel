# ----------- Aesthetic paramters ------------------#

# Define Colours -------------------------------------------------------

library("wesanderson")

colours <- c(rgb(red = .0, green = 0, blue = 0.8, alpha = 0.5),
             rgb(red = .8, green = 0.8, blue = 0, alpha = 0.5))

colboxes<- c('#d7191c','#fdae61','#2b83ba',"black")

palette('default')

colorValues <- c("#1b9e77","#d95f02","#7570b3",
                 "#e7298a","#66a61e","#e6ab02")[6:1]

paletteMeans <- colorRampPalette(c('#d73027','#fc8d59','#fee090',
                                   '#e0f3f8','#91bfdb','#4575b4')[6:1],
                                 alpha=TRUE)

paletteVar <- colorRampPalette(c('#d8b365','#f5f5f5','#5ab4ac'),alpha=TRUE)

multDiscrPallet <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

