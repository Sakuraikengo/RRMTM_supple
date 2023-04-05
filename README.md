# RRMTM_supple

Data which was used in the "Evaluating drought tolerance stability in soybean by the response of irrigation change captured from time-series multispectral data" was uploaded here.

If you want to see the raw data which is used in this manuscript, please click the "supplementary.zip". Then, Download the zip file. 

___
* scripts
    * 0.0.soilMoisture.R
Visualizing the soil moisture in each treatment
    * 0.1.makeAmat.R
Making additive genetic matrix from SNP data
    * 0.2.spectralValue.R
Calculating various VIs (including NDVI and NDRE)
    * 0.3.dataBind.R
Binding the data of VIs and fresh weight
    * 0.4.timeSeriesMS.R
Changing each day VI as time-series data
    * 1.functionCode.R
    Registering the function which was used in the scripts
    * 1.MTM_2.R
    Registering the function of "MTM_2"
    * 2.0.freshWeightCV.R
Calculating the coefficient of variation (CV)
    * 2.1.genomicPrediction.R
Predicting the CV using only genome data
    * 3.0.RRM.R
    Building random regression models (RRMs) with the various parameters (nr = 0, 1, 2)
    * 3.1.RRMlinear.R
    Building a RRM using nr = 1
    * 4.0.MTMfitting.R
Building multi-trait model (MTM) using time-series multi-spectral (MS) data and calculating the genetic correlation
    * 4.1.MTM.R
Predicting CV in Case1 using time-series MS data
    * 4.2.MTMyear.R
Predicting CV in Case3 using time-series MS data
    * 4.3.MTMsmallData.R
Predicting CV in Case3 using time-series MS data
    * 5.0.RRMMTMfitting.R
Building MTM using the parameters calculated in RRM and calculating the genetic correlation
    * 5.1.RRMMTM.R
Predicting CV in Case1 using the parameters calculated in RRM
    * 5.2.RRMMTMyear.R
Predicting CV in Case3 using the parameters calculated in RRM
    * 5.3.RRMMTMsmallData.R
Predicting CV in Case2 using the parameters calculated in RRM

 
