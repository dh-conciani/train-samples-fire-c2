## build optimized biome packages of data
## dhemerson.costa@ipam.org.br || wallace.silva@ipam.org.br

## set csv directory
root  <- './table' #

## get filenames
files <- list.files(root)

## get biomes
### change "-" by "_"
basenames <- gsub('-', '_', files)
### get biome names
biomes <- unique(sapply(strsplit(basenames, split='_', fixed=TRUE), function(x) (x[7])))[1:6]

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
    file_ij <- read.csv(paste0(root, '/', files_i[j]), encoding= 'UTF-8')
    
    ## remove unnecessary columns
    file_ij <- file_ij[ , -which(names(file_ij) %in% c("system.index",
                                                       "cena",
                                                       "region",
                                                       "fogo",
                                                       "fire_str",
                                                       "id",
                                                       "path",
                                                       "name",
                                                       "type",
                                                       ".geo"))]
    
    ## skip if colums dont contains spectral signatures
    if (sum(colnames(file_ij) == 'nir') == 0) next
    ## standardize names
    colnames(file_ij) <- gsub('fogo', 'fire', colnames(file_ij))
    
    ## use regex to remove underscore when it is the first character in shortname string
    file_ij$shortname <- gsub("^\\_","",file_ij$shortname)
    
    ## parse region name from string 
    file_ij$region <- sapply(strsplit(file_ij$shortname, split='_', fixed=TRUE), function(x) (x[2]))
    
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
