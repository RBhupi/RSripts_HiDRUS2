#!---------------------------------------------------------------------------
# @author Bhupendra Raut (www.baraut.info)
# @brief This R script reads a netCDF file containing recurrent period for the region.
# The TS data from nc2nc_ents programm are read and bias corrected using the ARI for 20 years.
# the data is written back to netCDF file with suffix, "_bc20ARI".
# The data is obtained from BoM and the contact person is Alan Seed or Mark Curtis.

#@todo :
#==========================================================================================
library(ncdf4)


#=============================================================================== Function Definions Starts
#--------------------------------------------------------------------get.fPath()
#get the full path of a file using pattern
get.fPath<-function(path, pat){
  path<-list.files(path=path, pattern = pat, recursive = T, full.names = T)
  print(path)
  if(length(path)!=1){
    stop("more/less than one file found.")
  }
  print(paste(length(path), "file(s) found."))
  return(path)
}

#--------------------------------------------------------------get.start_count()
#get start and count vectors for given lat-lon values.
#It also adds the time dim for default reading.
get.start_count<-function(ncfile, lons, lats, lat_name, lon_name){
  lon_ind <- start_count(ncfile, lons, lon_name)
  lat_ind <- start_count(ncfile, lats, lat_name)
  out<-list(c(lon_ind[1], lat_ind[1], 5), c(lon_ind[2], lat_ind[2], 1))
  return(out)
}

#this function finds start and count indices for lat and lon
start_count<-function(ncfile, reqVals, dim_name){
  inVals<-ncvar_get(ncfile, varid = dim_name) #dim values in the file
  ind1<-which(abs(inVals-reqVals[1])==min(abs(inVals-reqVals[1])))
  ind2<-which(abs(inVals-reqVals[2])==min(abs(inVals-reqVals[2])))
  start<-min(ind1, ind2)
  end<-max(ind1, ind2)
  count<-end-start+1
  return(c(start, count))
}
#========================================================================

#read TS data
tspath <- "/Users/bhupendra1/data/MLB/Hist/ts/"
location <- c("MRO", "Mlt", "SCr", "NWN", "Tlg")
loc_name <- c("Melbourne Regional Office", "Melton", "Strath Creek", 
              "Narre Warren North", "Toolangi")
model_name <- c("MIROC5", "CSIRO-Mk3")
city <- c("MLB", "SYD")
fname_pat<- paste("raints_", model_name[1], ".*", city[1], ".*", location[2], ".*", sep="")
ts_fName<-get.fPath(tspath, fname_pat)
ts_ncFile<-nc_open(ts_fName)
ts <- ncvar_get(ts_ncFile, varid="rain")
time <- ncvar_get(ts_ncFile, varid = "time")

# read the ARI data for given loaction
ari_fPath <-"/Users/bhupendra1/Dropbox/Bhupendra_shared/data/65_2_19700101_000000.ari_base.nc"
ari_ncfile <- nc_open(fPath)
indVec<- get.start_count(ari_ncfile, lons =c(144.55, 144.55), lats = c(44.45, 44.45), 
                lat_name = "y", lon_name = "x")

ari_loc <- ncvar_get(ari_ncfile, varid="ari_base", start = indVec[[1]], count = indVec[[2]])

#now read all ARI the data for computing bias
ari_map <- ncvar_get(ari_ncfile, varid = "ari_base", start=c(1, 1, 5), count = c(-1, -1, 1))
ari_mean <-mean(ari_map, na.rm = T)

local_bias <- as.vector(ari_loc/ari_mean)

#multiply with the bias factor
ts_bcARI20 <- local_bias * ts


#Create the output file with .nc extension
ofileName <- str_replace(ts_fName, pattern = ".nc", "_bc20ARI.nc")
ofile<-nc_create(ofileName, vars =list(e_date, e_rainAccum, e_rainMax,
                                                             e_rainDur, ts_prms))

id <- ncdim_def("id", units = "", vals = c(1:1), unlim = T, 
                longname = "identity number for ensemble member") 

time <- ncdim_def("time", units = "seconds since", vals = c(1:1), unlim = T, 
                longname = "identity number for ensemble member") 

e_rainAccum <- ncvar_def("e_rainMean", units = "mm", dim = id, 
                         longname = "total rainfall for the event" )










