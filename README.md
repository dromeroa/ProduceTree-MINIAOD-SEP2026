# ProduceTree-MINIAOD-SEP2026

En este espacio vamos a indicar los pasos para producir los trees con los archivos MINI-AOD previamente descargados

En la nube dirigirse a

/hepdata/David/cms_open_data_python

Luego hacer:(Para ver los contenedores activos)

docker ps -a

Aqui indetificar que se encuentre:  analisis_ul16v2

Si esta, hacemos:

docker stop analisis_ul16v2

Luego:

docker rm analisis_ul16v2

Despues:

docker run -it --name analisis_ul16v2   -v /hepdata/David/cms_open_data_work/CMSSW_10_6_30:/code   -u $(id -u):$(id -g)   cmsopendata/cmssw_10_6_30-slc7_amd64_gcc700 /bin/bash
Setting up CMSSW_10_6_30

Ahi ingresamos a:

/code/CMSSW_10_6_30/src/Analisis/MyAnalyzer/python






