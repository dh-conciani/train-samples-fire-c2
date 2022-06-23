## build optimized biome packages of data
## dhemerson.costa@ipam.org.br || wallace.silva@ipam.org.br

## set csv directory
root <- './table'

## get filenames
files <- list.files(root)

## get biomes
### change "-" by "_"
basenames <- gsub('-', '_', files)
### get biome names
biomes <- unique(sapply(strsplit(basenames, split='_', fixed=TRUE), function(x) (x[6])))

## for each biome
for (i in 1:length(biomes)) {
  print(paste0('processing biome: ', biomes[i]))
  ## get files
  files_i <- list.files(root, pattern= biomes[i])
  print(paste0(length(files_i), ' files detected'))
  ## build recipe
  recipe <- as.data.frame(NULL)
  for (j in 1:length(files_i)) {
    print(paste0('file ', j, ' of ', length(files_i)))
    ## read file
    file_ij <- read.csv(paste0(root, '/', files_i[j]), encoding= 'UTF-8')[-1][-2][-4][-13]
    
    ## compute indexes
    file_ij$ndvi <- (file_ij$nir - file_ij$red) / (file_ij$nir + file_ij$red)
    file_ij$baim <- 1 / ((0.05 - file_ij$nir)^2 + (0.2 - file_ij$swir1)^2)
    file_ij$mirbi = 10 * file_ij$swir1 - 9.8 * file_ij$nir + 2
    
    ## bind into recipe
    recipe <- rbind(recipe, file_ij)
  }
  ## export
  write.csv(recipe, paste0('./table_biomes/', 'samples_', biomes[i], '.csv'))
}

