# ** @description get spectral signatures from fire and non-fire classes 
# Mapbiomas Fogo - GT
# Wallace Silva, Dhemerson Conciani e Soltan Galeno

## read libraries
library(rgee)
library(rgeeExtra)
library(reticulate)
library(rlist)

## intialize
ee_Initialize()

## get block list module
#blockList_landsat <- module("users/geomapeamentoipam/GT_Fogo_MapBiomas:00_Tools/module-blockList")$landsat()

## set function to compute minNBR
addBand_NBR <- function (image) {
  # define expression
  exp <- '( b("nir") - b("swir2") ) / ( b("nir") + b("swir2") )'
  ## get mosaic
  minimoNBR <- image$
    expression(exp)$
    multiply(-1)$
    add(1)$
    multiply(1000)$
    int16()$
    rename('nbr')
  
  return(image$addBands(minimoNBR))
}

## define function to clip scenes
clipBoard_Landsat <- function(image) {
  return(image$updateMask(ee$Image()$
                            paint(image$geometry()$buffer(-3000))$eq(0))
         )
}

## function to standardize landsat [TM asnd ETM+] collection 2
corrections_LS57_col2 <- function(image) {
  opticalBands <- image$select('SR_B.')$multiply(0.0000275)$add(-0.2)
  thermalBands <- image$select('ST_B.*')$multiply(0.00341802)$add(149.0)
  
  ## stack image
  image <- image$addBands(opticalBands, NULL, TRUE)$
    addBands(thermalBands, null, true)
  
  ## cloud mask
  cloudShadowBitMask <- bitwShiftL(1, 3)
  cloudsBitMask <- bitwShiftL(1, 5)
  
  ## aplly QA mask
  qa <- image$select('QA_PIXEL');
  mask <- qa$bitwiseAnd(cloudShadowBitMask)$eq(0)$
    and(qa$bitwiseAnd(cloudsBitMask)$eq(0))
  
  ## set radiometric noise mask
  bitwiseExtract <- function(value, fromBit, toBit) {
    if (length(toBit) == 0) {
      toBit = fromBit
    }
    
    maskSize <- ee$Number(1)$add(toBit)$subtract(fromBit)
    mask <- ee$Number(1)$leftShift(maskSize)$subtract(1)
    
    return(value$rightShift(fromBit)$bitwiseAnd(mask))
  }
  
  clear <- bitwiseExtract(qa, 6) ## 1 if clear
  water <- bitwiseExtract(qa, 7) ## 1 if water
  
  radsatQA <- image$select('QA_RADSAT');
  band5Saturated <- bitwiseExtract(radsatQA, 4)  ##  0 if band 5 is not saturated
  anySaturated <- bitwiseExtract(radsatQA, 0, 6) ##  0 if no bands are saturated
  
  mask_saturation <- clear$
    or(water)$
    and(anySaturated$not())
  
  ## is visible bands with negative reflectance? 
   negative_mask <- image$select('SR_B1')$gt(0)$and(
     image$select('SR_B2')$gt(0))$and(
       image$select('SR_B3')$gt(0))$and(
         image$select('SR_B4')$gt(0))$and(
           image$select('SR_B5')$gt(0))$and(
             image$select('SR_B7')$gt(0))
   
   image = image$
     updateMask(mask)$
     updateMask(mask_saturation)$
     updateMask(negative_mask)
  
   ## rename bands
   oldBands <- c('SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7')
   newBands <- c('blue', 'green','red',  'nir',  'swir1','swir2')
   
   image <- image$select(oldBands, newBands)
   
   return (addBand_NBR(image))
}

## function to standardize landsat collection 2
corrections_LS8_col2 <- function(image) {
  opticalBands <- image$select('SR_B.')$multiply(0.0000275)$add(-0.2)
  ## rectfy to dark corpse reflectance == -0.0000000001
  opticalBands <- opticalBands$multiply(10000)$subtract(0.0000275 * 0.2 * 1e5 * 100)$round()
  thermalBands <- image$select('ST_B.*')$multiply(0.00341802)$add(149.0)
  
  ## stack image
  image <- image$addBands(opticalBands, NULL, TRUE)$
    addBands(thermalBands, NULL, TRUE)
  
  #  If the cloud bit (3) is set and the cloud confidence (9) is high
  #   or the cloud shadow bit is set (3), then it's a bad pixel.
  qa <- image$select('QA_PIXEL')
  cloud <- qa$bitwiseAnd(bitwShiftL(1, 3))$
    and(qa$bitwiseAnd(bitwShiftL(1, 9)))$
    or(qa$bitwiseAnd(bitwShiftL(1, 4)))
  
  # If the clear bit (6) is set 
  # or water bit is set (7), then it's a good pixel 
  good_pixel <- qa$bitwiseAnd(bitwShiftL(1, 6))$
    or(qa$bitwiseAnd(bitwShiftL(1, 7)))
  
  # read radsat 
  radsatQA = image$select('QA_RADSAT')
  
  #  is any band saturated? 
  saturated = radsatQA$bitwiseAnd(bitwShiftL(1, 0))$
    or(radsatQA$bitwiseAnd(bitwShiftL(1, 1)))$
    or(radsatQA$bitwiseAnd(bitwShiftL(1, 2)))$
    or(radsatQA$bitwiseAnd(bitwShiftL(1, 3)))$
    or(radsatQA$bitwiseAnd(bitwShiftL(1, 4)))$
    or(radsatQA$bitwiseAnd(bitwShiftL(1, 5)))$
    or(radsatQA$bitwiseAnd(bitwShiftL(1, 6)))
    
  ## is visible bands with negative reflectance? 
  negative_mask <- image$select('SR_B1')$gt(0)$and(
    image$select('SR_B2')$gt(0))$and(
      image$select('SR_B3')$gt(0))$and(
        image$select('SR_B4')$gt(0))$and(
          image$select('SR_B5')$gt(0))$and(
            image$select('SR_B7')$gt(0))

  image <- image$
    updateMask(cloud$not())$
    updateMask(good_pixel)$
    updateMask(saturated$not())$
    updateMask(negative_mask)
  
  ## rename bands
  oldBands <- c('SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7')
  newBands <- c('blue', 'green','red',  'nir',  'swir1','swir2')
  
  image <- image$select(oldBands, newBands)
  
  return (addBand_NBR(image))
}

satellites <- ee$List(list(
  list(
    'type'='ImageCollection',
    'name'= 'LT05_C02_T1_L2',
    'id'='LANDSAT/LT05/C02/T1_L2',
    'years'= c(seq(1985, 2011)),
    'allProcess' = function(col) {
      #col <- ee$ImageCollection(col)$filter(ee$Filter$inList('system:index', blockList_landsat)$not())
      return(ee$ImageCollection(col)$map(function(image) {
        image <- clipBoard_Landsat(image)
        image <- corrections_LS57_col2(image)
        return (image)
        })
      )
    },
    'reduceProcess' = function (col) {
      return(col$qualityMosaic('nbr'))
    }
  ),
  
  list(
      'type'='ImageCollection',
      'name'='LE07_C02_T1_L2',
      'id'='LANDSAT/LE07/C02/T1_L2',
      'years'= c(seq(1999, 2011)),
      'allProcess' = function(col) {
        #col <- ee$ImageCollection(col)$filter(ee$Filter$inList('system:index', blockList_landsat)$not())
        return(ee$ImageCollection(col)$map(function(image) {
          image <- clipBoard_Landsat(image)
          image <- corrections_LS57_col2(image)
          return (image)
        })
        )
      },
      'reduceProcess' = function (col) {
        return(col$qualityMosaic('nbr'))
      }
  ),
  
  list(
    'type'='ImageCollection',
    'name'='LC08_C02_T1_L2',
    'id'='LANDSAT/LC08/C02/T1_L2',
    'years'= c(seq(2013, 2022)),
    'allProcess' = function(col) {
      #col <- ee$ImageCollection(col)$filter(ee$Filter$inList('system:index', blockList_landsat)$not())
      return(ee$ImageCollection(col)$map(function(image) {
        image <- clipBoard_Landsat(image)
        image <- corrections_LS57_col2(image)
        return (image)
      })
      )
    },
    'reduceProcess' = function (col) {
      return(col$qualityMosaic('nbr'))
    }
  ),
  
  list(
    'type'='ImageCollection',
    'name'='LC09_C02_T1_L2',
    'id'='LANDSAT/LC09/C02/T1_L2',
    'years'= 2022,
    'allProcess' = function(col) {
      #col <- ee$ImageCollection(col)$filter(ee$Filter$inList('system:index', blockList_landsat)$not())
      return(ee$ImageCollection(col)$map(function(image) {
        image <- clipBoard_Landsat(image)
        image <- corrections_LS57_col2(image)
        return (image)
      })
      )
    },
    'reduceProcess' = function (col) {
      return(col$qualityMosaic('nbr'))
    }
  )
))


## set sample destination
samples <- ee$List(list(
 list(
    'type'='FOLDER',
    'id'='projects/mapbiomas-workspace/FOGO/AMOSTRAS_COLECAO1',
    'name'='AMOSTRAS_COLECAO1'
  ),
  
  list(
    'type'='FOLDER',
    'id'='projects/mapbiomas-workspace/FOGO/AMOSTRAS_COLECAO2',
    'name'='AMOSTRAS_COLECAO2'
  ),
  
  list(
    'type'='FOLDER',
    'id'='projects/mapbiomas-workspace/FOGO/AMOSTRAS_SENTINEL',
    'name'='AMOSTRAS_SENTINEL'
  ),
  
  list(
    'type'='FOLDER',
    'id'='projects/mapbiomas-workspace/FOGO/AVALIACOES_FOGO',
    'name'='AVALIACOES_FOGO'
  )
))



samples <- samples$filter(function(obj) {
  return(obj$name == 'AMOSTRAS_COLECAO2')
  })$map(function(obj) {
    obj$assets = ee$data$listAssets(obj$id)$assets$
      map(function(o) {
        o$year = o$id$slice(-4)
        return(ee$List(
          list(
          'type' = o$type,
          'id' = o$id,
          'name'= o$name$split('/')[7],
          'path'= o$path$split('/')[6],
          'year'= o$id$slice(-4),
          'sat'= o$id$split('_l')[1]$split('_')[0]
          ))
        )
      })
    
    return(obj)
    })$map(
      function(obj) {
        obj$assets = obj$assets$
          map(function(o){
            featureCollection <- ee$FeatureCollection(o$id)
            recipe <- ee.ImageCollection()
            
            satellites$map(function(satellite) {
              satellite$years = satellite$years$map(function(year) {
                return(paste('', year))
              })
              return(satellite)
            })$filter(function(satellite) {
              return(satellite$years$indexOf(o$year) != -1)
            })$map(function(satellite){
              start = ee$Date(paste0('', year, "-01-01"))$millis()
              end = ee$Date('1971-01-01')$millis()
              end = start$add(end)
              
              collection = satellite$allProcess(satellite$id)$
                filterBounds(featureCollection)$
                filterDate(start, end)
              
              recipe = recipe$merge(collection)
            })
            
            mosaic = recipe$qualityMosaic('nbr')
            
            sampleRegions = mosaic.sampleRegions(
                         collection=featureCollection,
                         scale=30, 
                         geometries=FALSE)
            
            o$sampleRegions = sampleRegions
             return (o)
          })$map(
            function(o) {
              description = paste0('assign_spectral_data-v0_20220604-',o$name$replace('train_test_fire_nbr_',''))
              
              task <- ee$batch$Export.table.toDrive(
                            collection=o$sampleRegions,
                            description=description,
                            folder='mapbiomas-fogo',
                            fileNamePrefix=description,
                            fileFormat='csv',
                            )
              
              task.start()

            }
          )
        return(obj)
      }
    )

