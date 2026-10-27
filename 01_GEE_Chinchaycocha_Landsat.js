// Lago Chinchaycocha (Junín) visto desde el espacio - Análisis 10 años  (2016–2026)
// PARÁMETROS INICIALES 
var ANIO_INICIO = 2016;
var ANIO_FIN    = 2026;
var MES_SECO_INI = 6;     // junio
var MES_SECO_FIN = 9;     // setiembre
var UMBRAL_MNDWI = 0;    // MNDWI  0  = agua

// definimos área de estudio 
var caja = ee.Geometry.Rectangle([-76.30, -11.25, -75.90, -10.80]);
// Punto en el centro del lago (para ubicar la cuenca)
var centroLago = ee.Geometry.Point([-76.10, -11.02]);

Map.centerObject(caja, 10);
Map.addLayer(caja, {color 'yellow'}, 'Zona de estudio', false);


// 1. MÁSCARA DEL LAGO
// Usamos el máximo histórico de agua del JRC (1984–2021) para que solo
// contemos píxeles que alguna vez fueron agua (evita confundir sombras
// de cerros o suelos oscuros con agua).

var jrc = ee.Image('JRCGSW1_4GlobalSurfaceWater');
var zonaLago = jrc.select('max_extent').eq(1).clip(caja);
Map.addLayer(zonaLago.selfMask(), {palette ['cyan']}, 'Máx. extensión histórica (JRC)', false);

//  2.  IMÁGENES LANDSAT 8 + 9 (Colección 2, Nivel 2 = reflectancia)

function prepararLandsat(img) {
   // 2A. Quitar nubes, sombras y nieve con la banda de calidad QA_PIXEL
  var qa = img.select('QA_PIXEL');
  var nubeDilatada = qa.bitwiseAnd(1  1).neq(0);
  var nube         = qa.bitwiseAnd(1  3).neq(0);
  var sombra       = qa.bitwiseAnd(1  4).neq(0);
  var nieve        = qa.bitwiseAnd(1  5).neq(0);
  var limpio = nubeDilatada.or(nube).or(sombra).or(nieve).not();

  // 2B. Pasar números digitales a reflectancia (0–1)
  var refl = img.select('SR_B.').multiply(0.0000275).add(-0.2);

  // 2C. MNDWI = (Verde − SWIR1)  (Verde + SWIR1)  =  B3 y B6 en Landsat 89
  var mndwi = refl.normalizedDifference(['SR_B3', 'SR_B6']).rename('MNDWI');

  return refl.addBands(mndwi).updateMask(limpio)
             .copyProperties(img, ['systemtime_start']);
}

var landsat = ee.ImageCollection('LANDSATLC08C02T1_L2')
  .merge(ee.ImageCollection('LANDSATLC09C02T1_L2'))
  .filterBounds(caja)
  .filter(ee.Filter.calendarRange(ANIO_INICIO, ANIO_FIN, 'year'))
  .filter(ee.Filter.calendarRange(MES_SECO_INI, MES_SECO_FIN, 'month'))
  .map(prepararLandsat);

print('N° de imágenes Landsat usadas (época seca)', landsat.size());


// 3. ÁREA DEL LAGO POR AÑO
// Para cada año mediana del MNDWI de jun–set - agua = MNDWI  0
// - área = suma de (píxel agua × área del píxel)

var anios = ee.List.sequence(ANIO_INICIO, ANIO_FIN);

function areaDelAnio(anio) {
  anio = ee.Number(anio);
  var col = landsat.filter(ee.Filter.calendarRange(anio, anio, 'year'));
  var mndwi = col.select('MNDWI').median();
  var agua = mndwi.gt(UMBRAL_MNDWI).and(zonaLago);
  var areaKm2 = agua.multiply(ee.Image.pixelArea()).divide(1e6)
    .reduceRegion({
      reducer ee.Reducer.sum(),
      geometry caja,
      scale 30,
      maxPixels 1e10
    }).get('MNDWI');
  return ee.Feature(null, {
    'anio' anio,
    'area_km2' areaKm2,
    'n_imagenes' col.size()
  });
}

var tablaArea = ee.FeatureCollection(anios.map(areaDelAnio));
print('Área del lago por año (km²)', tablaArea);

print(ui.Chart.feature.byFeature(tablaArea, 'anio', 'area_km2')
  .setChartType('LineChart')
  .setOptions({
    title 'Chinchaycocha · Área de agua en época seca (jun–set)',
    hAxis {title 'Año', format '####'},
    vAxis {title 'Área (km²)'},
    pointSize 5, lineWidth 2, colors ['#05B4EB'], legend {position 'none'}
  }));



 // 4. COMO SE VE EL LAGO  (ANTES VS DESPUES)

function mapaAnio(anio) {
  var col = landsat.filter(ee.Filter.calendarRange(anio, anio, 'year'));
  var rgb = col.select(['SR_B4', 'SR_B3', 'SR_B2']).median().clip(caja);
  var agua = col.select('MNDWI').median().gt(UMBRAL_MNDWI).and(zonaLago).selfMask();
  Map.addLayer(rgb, {min 0, max 0.25}, 'Color real ' + anio, false);
  Map.addLayer(agua, {palette ['#05B4EB']}, 'Agua ' + anio, false);
}
mapaAnio(ANIO_INICIO);
mapaAnio(2023);
mapaAnio(ANIO_FIN);


 // 5.  CLIMA LLUVIA (CHIRPS) Y TEMPERATURA  EVAPORACIÓN (ERA5-Land)
 // Año hidrológico = setiembre (año−1) a agosto (año)

var cuenca = ee.FeatureCollection('WWFHydroSHEDSv1Basinshybas_7')
  .filterBounds(centroLago).geometry();
Map.addLayer(cuenca, {color 'orange'}, 'Cuenca aportante', false);
print('Área de la cuenca (km²)', cuenca.area(1000).divide(1e6));

var chirps = ee.ImageCollection('UCSB-CHGCHIRPSDAILY');
var era5   = ee.ImageCollection('ECMWFERA5_LANDMONTHLY_AGGR');

function climaDelAnio(anio) {
  anio = ee.Number(anio);
  var ini = ee.Date.fromYMD(anio.subtract(1), 9, 1);
  var fin = ee.Date.fromYMD(anio, 9, 1);

  // Lluvia total del año hidrológico (mm) – promedio en la cuenca
  var pp = chirps.filterDate(ini, fin).select('precipitation').sum()
    .reduceRegion({reducer ee.Reducer.mean(), geometry cuenca, scale 5000})
    .get('precipitation');

   // Temperatura media (°C) y evaporación (mm) de ERA5-Land sobre la cuenca
  var e5 = era5.filterDate(ini, fin);
  var temp = e5.select('temperature_2m').mean().subtract(273.15)
    .reduceRegion({reducer ee.Reducer.mean(), geometry cuenca, scale 11000})
    .get('temperature_2m');
  // En ERA5 la evaporación sale negativa (es agua que SE VA) la volteamos
  var evap = e5.select('total_evaporation_sum').sum().multiply(-1000)
    .reduceRegion({reducer ee.Reducer.mean(), geometry cuenca, scale 11000})
    .get('total_evaporation_sum');

  return ee.Feature(null, {
    'anio_hidrologico' anio,
    'lluvia_mm' pp,
    'temp_C' temp,
    'evap_mm' evap,
    'meses_ERA5' e5.size()    si es  12, el año aún está incompleto
  });
}

var tablaClima = ee.FeatureCollection(anios.map(climaDelAnio));
print('Clima por año hidrológico', tablaClima);

print(ui.Chart.feature.byFeature(tablaClima, 'anio_hidrologico', ['lluvia_mm', 'evap_mm'])
  .setChartType('ColumnChart')
  .setOptions({
    title 'Entradas vs salidas del ciclo hidrológico (cuenca)',
    hAxis {title 'Año hidrológico', format '####'},
    vAxis {title 'mm  año'},
    colors ['#05B4EB', '#FBBC05']
  }));

print(ui.Chart.feature.byFeature(tablaClima, 'anio_hidrologico', 'temp_C')
  .setChartType('LineChart')
  .setOptions({
    title 'Temperatura media del aire (ERA5-Land)',
    hAxis {title 'Año hidrológico', format '####'},
    vAxis {title '°C'},
    pointSize 5, colors ['#E5533D'], legend {position 'none'}
  }));



 // 6. DESCARGAR LAS TABLAS 

Export.table.toDrive({
  collection tablaArea,
  description 'Chinchaycocha_area_2016_2026',
  fileFormat 'CSV'
});
Export.table.toDrive({
  collection tablaClima,
  description 'Chinchaycocha_clima_2016_2026',
  fileFormat 'CSV'
});
