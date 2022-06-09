/** @description coletando comportamento espectral das amostras de fogo e não fogo do Mapbiomas-Fogo
  Grupo de trabalho de mapeamento de fogo no Brasil - MapBiomas Fogo
  
  2022-06-04 - Wallace Silva, Dhemerson Conciani e Soltan Galeno
 * 
*/ 
//--- --- --- funçoes e variaveis auxiliares para o dataset
var blockList_landsat = require('users/geomapeamentoipam/GT_Fogo_MapBiomas:00_Tools/module-blockList').landsat();

function addBand_NBR (image){
  var exp = '( b("nir") - b("swir2") ) / ( b("nir") + b("swir2") )';
  var minimoNBR = image
    .expression(exp)
    // -> na formula da USGS as cicatrizes ocupam os menores valores e multiplicamos o resultado por -1 para que 
    // as cicatrizes ocupem os valores maximos utilizados como referencia no processamento dos mosaicos de qualidade
    .multiply(-1)
    // -> adequações legadas
    .add(1)
    .multiply(1000)
    .int16()
    .rename("nbr");
  return image
    .addBands(minimoNBR);
}
// recortando bordas de cenas landsat
function clipBoard_Landsat(image){
  return image
    .updateMask(ee.Image().paint(image.geometry().buffer(-3000)).eq(0));
}

// correção radiometrica das imagens landsat coleção 2
function corrections_LS57_col2 (image){
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  // - return 
  
  image = image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
              
  // mascara de nuvem
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);


  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(qa.bitwiseAnd(cloudsBitMask).eq(0));

  // mascara de ruidos, saturação radiométrica
  function bitwiseExtract(value, fromBit, toBit) {
    if (toBit === undefined)
      toBit = fromBit;
    var maskSize = ee.Number(1).add(toBit).subtract(fromBit);
    var mask = ee.Number(1).leftShift(maskSize).subtract(1);
    return value.rightShift(fromBit).bitwiseAnd(mask);
  }

  var clear = bitwiseExtract(qa, 6); // 1 if clear
  var water = bitwiseExtract(qa, 7); // 1 if water

  var radsatQA = image.select('QA_RADSAT');
  var band5Saturated = bitwiseExtract(radsatQA, 4); // 0 if band 5 is not saturated
  var anySaturated = bitwiseExtract(radsatQA, 0, 6); // 0 if no bands are saturated

  var mask_saturation = clear
    .or(water)
    .and(anySaturated.not());
  
  // is visible bands with negative reflectance? 
  var negative_mask = image.select(['SR_B1']).gt(0).and(
    image.select(['SR_B2']).gt(0)).and(
      image.select(['SR_B3']).gt(0)).and(
        image.select(['SR_B4']).gt(0)).and(
          image.select(['SR_B5']).gt(0)).and(
            image.select(['SR_B7']).gt(0));
  
  // - return
  image = image
    .updateMask(mask)
    .updateMask(mask_saturation)
    .updateMask(negative_mask);
        

  var oldBands = ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7',];
  var newBands = ['blue', 'green','red',  'nir',  'swir1','swir2'];
  // - return
  image = image.select(oldBands,newBands);
  
  // -
  return addBand_NBR(image);
}

function corrections_LS8_col2 (image){

  // - radiometric correction
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  // rectfy to dark corpse reflectance == -0.0000000001
  opticalBands = opticalBands.multiply(10000).subtract(0.0000275 * 0.2 * 1e5 * 100).round();
  
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  
  // - return 
  image = image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
    
  // - masks
  // If the cloud bit (3) is set and the cloud confidence (9) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var qa = image.select('QA_PIXEL');
      var cloud = qa.bitwiseAnd(1 << 3)
      .and(qa.bitwiseAnd(1 << 9))
      .or(qa.bitwiseAnd(1 << 4));
  
  // If the clear bit (6) is set 
  // or water bit is set (7), then it's a good pixel 
  var good_pixel  = qa.bitwiseAnd(1 << 6)
      .or(qa.bitwiseAnd(1 << 7));

  // read radsat 
  var radsatQA = image.select('QA_RADSAT');
  // Is any band saturated? 
  var saturated = radsatQA.bitwiseAnd(1 << 0)
    .or(radsatQA.bitwiseAnd(1 << 1))
      .or(radsatQA.bitwiseAnd(1 << 2))
        .or(radsatQA.bitwiseAnd(1 << 3))
          .or(radsatQA.bitwiseAnd(1 << 4))
            .or(radsatQA.bitwiseAnd(1 << 5))
              .or(radsatQA.bitwiseAnd(1 << 6));

  // is any band with negative reflectance? 
  var negative_mask = image.select(['SR_B1']).gt(0).and(
    image.select(['SR_B2']).gt(0)).and(
      image.select(['SR_B3']).gt(0)).and(
        image.select(['SR_B4']).gt(0)).and(
          image.select(['SR_B5']).gt(0)).and(
            image.select(['SR_B7']).gt(0));

  
  // -return 
  image = image
  .updateMask(cloud.not())
  .updateMask(good_pixel)
  .updateMask(saturated.not())
  .updateMask(negative_mask);
  
  
  // correction bandnames to default
  var oldBands = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7',];
  var newBands = ['blue', 'green','red',  'nir',  'swir1','swir2'];
  // -
  // - return
  image = image.select(oldBands,newBands);
  
  // -
  return addBand_NBR(image);
}

var satellites = [
  {
    'type':'ImageCollection',
    'id':'LANDSAT/LT05/C02/T1_L2',
    'years':[
      1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011
    ],
    allProcess:function(col){
      col = ee.ImageCollection(col).filter(ee.Filter.inList('system:index', blockList_landsat).not());
      return ee.ImageCollection(col).map(function(image){
        image = clipBoard_Landsat(image);
        image = corrections_LS57_col2(image);
        
        return image;
      });
    },
    reduceProcess:function(col){
      return col.qualityMosaic('nbr');
    }
  },
  {
    'type':'ImageCollection',
    'name':'LE07_C02_T1_L2',
    'id':'LANDSAT/LE07/C02/T1_L2',
    'years':[1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,/*2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022*/],
     allProcess:function(col){
      col = ee.ImageCollection(col).filter(ee.Filter.inList('system:index', blockList_landsat).not());
     return ee.ImageCollection(col).map(function(image){
        image = clipBoard_Landsat(image);
        image = corrections_LS57_col2(image);
        
        return image;
      });
    },
    reduceProcess:function(col){
      return col.qualityMosaic('nbr');
    }
  },
  {
    'type':'ImageCollection',
    'name':'LC08_C02_T1_L2',
    'id':'LANDSAT/LC08/C02/T1_L2',
    'years':[
      2013,2014,2015,2016,2017,2018,2019,2020,2021,2022
    ],
    allProcess:function(col){
      col = ee.ImageCollection(col).filter(ee.Filter.inList('system:index', blockList_landsat).not());
      return ee.ImageCollection(col).map(function(image){
        image = clipBoard_Landsat(image);
        image = corrections_LS8_col2(image);
        
        return image;
      });
    },
    reduceProcess:function(col){
      return col.qualityMosaic('nbr');
    }
  },
  {
    'type':'ImageCollection',
    'name':'LC09_C02_T1_L2',
    'id':'LANDSAT/LC09/C02/T1_L2',
    'years':[
      2022
    ],
    allProcess:function(col){
      col = ee.ImageCollection(col).filter(ee.Filter.inList('system:index', blockList_landsat).not());
      return ee.ImageCollection(col).map(function(image){
        image = clipBoard_Landsat(image);
        image = corrections_LS8_col2(image);
        
        return image;
      });
    },
    reduceProcess:function(col){
      return col.qualityMosaic('nbr');
    }
  },
];

var samples = [
  {
    'type':'FOLDER',
    'id':'projects/mapbiomas-workspace/FOGO/AMOSTRAS_COLECAO1',
    'name':'AMOSTRAS_COLECAO1',
  },
  {
    'type':'FOLDER',
    'id':'projects/mapbiomas-workspace/FOGO/AMOSTRAS_COLECAO2',
    'name':'AMOSTRAS_COLECAO2',
  },
  {
    'type':'FOLDER',
    'id':'projects/mapbiomas-workspace/FOGO/AMOSTRAS_SENTINEL',
    'name':'AMOSTRAS_SENTINEL',
  },
  {
    'type':'FOLDER',
    'id':'projects/mapbiomas-workspace/FOGO/AVALIACOES_FOGO',
    'name':'AVALIACOES_FOGO',
  },
  
];

var samples = samples
  .filter(function(obj){
    return obj.name === 'AMOSTRAS_COLECAO2';
  })
  .map(function(obj){
    
    obj.assets = ee.data.listAssets(obj.id)
      .assets
      // .slice(0,2)
      .map(function(o){
        o.year = o.id.slice(-4);
        return {
          'type':o.type,
          'id':o.id,
          'name':o.name.split('/')[7],
          'path':o.name.split('/')[6],
          'year':o.id.slice(-4),
          'sat':o.id.split('_l')[1].split('_')[0]
        };
      });

    return obj;})
  .map(function(obj){
    obj.assets = obj.assets
      .map(function(o){
        var featureCollection = ee.FeatureCollection(o.id);
    
        var recipe = ee.ImageCollection([]);
  
        satellites
        .map(function(satellite){
          satellite.years = satellite.years.map(function(year){return ''+year});
          return satellite;
        }).filter(function(satellite){return satellite.years.indexOf(o.year) !== -1}) // filtro segundo a lista de anos nos objetos dos satelites
        .forEach(function(satellite){
            var start = ee.Date(''+o.year+'-01-01').millis();
            var end = ee.Date('1971-01-01').millis(); // os sistema Unix, armazenam o tempo utilizando o ano de 1970 como referencia de 0 na contagem, '
            end = start.add(end);                     // ao solicitar os millisegundos de 1971 ele retorna os milissegundos de um ano inteiro
            
            var collection = satellite.allProcess(satellite.id)
              .filterBounds(featureCollection)
              .filterDate(start,end);
            
            recipe = recipe.merge(collection);
            // var mosaic = satellite.reduceProcess(collection);
            // print(collection,start,end);
            // Map.addLayer(mosaic);
          });
          
          var mosaic = recipe.qualityMosaic('nbr');
          
          var landcover = ee.Image('projects/mapbiomas-workspace/public/collection6/mapbiomas_collection60_integration_v1')
            .select(['classification_'+o.year]['classification']);
          
          mosaic = mosaic.addBands(landcover);
          
          var sampleRegions = mosaic.sampleRegions({
            collection:featureCollection,
            // properties:{},
            scale:30, 
            // projection:,
            // tileScale:,
            geometries:false
          });  

          o.sampleRegions = sampleRegions;//.limit(10);
          
          return o;
          // Map.addLayer(mosaic);
    })
      .forEach(function(o){
        var description = 'assign_spectral_data-v0_20220604-'+o.name.replace('train_test_fire_nbr_','');
        // print('o.sampleRegions',o.sampleRegions);
        Export.table.toDrive({
          collection:o.sampleRegions,
          description:description,
          folder:'mapbiomas-fogo',
          fileNamePrefix:description,
          fileFormat:'csv',
          // selectors:,
          // maxVertices:
        });

      });
    
    // obj.years = obj.assets.map(function(o){return o.year})
    //   .sort();
    // obj.years = obj.years.filter(function(este,i){return obj.years.indexOf(este) === i;});
    
    return obj;
  });


