# CAMxtools

## Descriptions
This CAMxtools package is a collection of following tools,

1. combine : Linearly combine variables in CAMx (UAM IV) or CMAQ (IOAPI) files and create either format file. For example, if you want to add 0.5x"O3" variable in file 1 and 0.5x"O3" in file 2. Following spec_def file is needed,

```
/new species    ,units     ,expression
O3              ,ppmv      ,0.5\*O3[1]+0.5\*O3[2]
```

2. combine\_timezones : When you have two CAMx/CMAQ files, metrics in one file is calculated based on Pacific Standard Time (PST) and another file is calculated based on Mountain Standard Time (MST). The map you want to show includes area in both time zones. This tool selects a value from either files for one grid cell based on the time zone that the grid cell belongs to.

3. extcells : Extract values from specific grid cells which are defined in a reference file. The reference cell can be grid indexes or lon/lat.

4. MATS : MATS is an EPA Windows program which projects observed pollutant concentrations using model results. Two model results are needed; One is basecase and another is for future year. Using the two results relative fashion, observation is projected. CAMxtools/MATS provides MATS input from CAMx or CMAQ results. O3 and PM can be projected using this tool.

5. metrics : Following metrics are frequently asked to evaluate model results. This tool processes metrics below.
   - Annual or monthly average
   - 1st or 4th highest of running 8 hour O3 over a year
   - 1st or 4th highest of 1 hour O3 over a year
   - 8th highest daily PM2.5

6. naaqs : When the metrics are calculated, this tool provides contributions from other species or sources at the time when a sum of pollutant hits the metric. For example, at a certain grid cell, total (sum of all sources) O3 reaches an annual maximum at 11:00 on 8/1. This tool lists contributions of sources at the time.

7. psd : It is similar to naaqs. But the source contribution is calculated independently. From the above example, instead of getting contribution of a source when the total O3 reaches the annual maximum. This calculates annual maximum of the source.

8. regrid : Convert map projection for CAMx or CMAQ files

9. vis : Calculates visibility from CAMx or CMAQ output files

10. W126: W126 is set by EPA and a metric which multiplies O3 mixing ratio and time exposed. Additional weightening factor is applied based on O3 mixing ratio. W126 is used to estimate O3 pollution impact to vegitation. This W126 tool calculates this metric from CAMx or CMAQ results.

## Install python and modules
1. Install Anaconda3:
  1.1 Go to https://www.continuum.io/downloads
  1.2 Select Linux and download Python 3.x
2. Install Pseudonetcdf: https://github.com/barronh/pseudonetcdf
  2.1 pip install http://github.com/barronh/pseudonetcdf/archive/master.zip
3. Install pandasql
  3.1 pip install pandasql
4. Install xarray
  4.1 pip install xarray
5. Install pytzwhere for the timezone processing
  5.1  pip install https://github.com/pegler/pytzwhere/archive/master.zip

