## build optimized biome packages of data
## dhemerson.costa@ipam.org.br || wallace.silva@ipam.org.br

## set csv directory
root  <- './table_biomes/' #

## get filenames
files <- list.files(root)

for (i in 1:length(files)) {
  print(files[i])
  
  ## read data
  data <- read.csv(paste0(root, files[i]))
  
  ## get L8
  L8 <- subset(data, sat == 8)
  
  ## adjust divide by 10000
  L8$blue <- L8$blue / 10000
  L8$green <- L8$green / 10000
  L8$red <- L8$red / 10000
  L8$nir <- L8$nir / 10000
  L8$swir1 <- L8$swir1 / 10000
  L8$swir2 <- L8$swir2 / 10000
  
  ## get others
  others <- subset(data, sat != 8)
  
  ## merge
  recipe <- rbind(L8, others)
  
  ## export
  write.csv(recipe, paste0('./table_biomes2/', files[i]))

}
