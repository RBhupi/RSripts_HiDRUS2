#!---------------------------------------------------------------------------
# @author Bhupendra Raut
# @brief We are making a parameter file for HiDRUS-2
#1. This programm reads a reanalysis/GCM/Radar daily fldmean rainfall data.
#2. Reads events and Syn library files and search for the closest event for the mean accumulation.
#3. Writes down all the required paramters for the HiDRUS-2 simulation.

#@todo :
# 1. Needs optimisation for the large array operations.
# 2. Add clause to select the day that is closest to the date in GCM to keep seasonality.
#               (Not needed when selected by clusters.)
#==========================================================================================
library(ncdf4)
library(Hmisc) #for find.match
source("/home/bhupendra/Dropbox/Bhupendra_shared/Scripts/RScripts/myFUNC.R")
options("scipen"=100) #force to print fixed numbers not scientific notations.

#=================================================================================
#----------------------set file paths.------------------------#
#keep all the files in the same directory for convinience.
setwd("/home/bhupendra/sim_real/erai-sim/pardir/era_parm_2000-2007")
gcm_rainFile <- "../erain1989-2010_50x50km_daysum_fmean_mlb.nc" #daily-rain file name
gcm_synFile <- "../ERA-I_MSLP_UV_1979-2014_0.75deg_MLB.nc" #daily-rain file name

rainLib_File <- "../RainLib_2008-14_MLB.nc"      #Rain Library file name
synLib_File<- "../SynLib_2008-14_MLB.nc"
outParm_File_prefix<-"h2parms_erai_2000-2007_" #no extension

ens_id<-seq(1:100)

min_rain=0.2
#-----------------set START and END dates of the simulation.------#
s_date=as.Date("2000-01-01", format = "%Y-%m-%d")  # start date of the simulation
e_date=as.Date("2007-12-31", format = "%Y-%m-%d")
#=================================================================================

#--------------------read the daily rain file----------------------------------#
gcm_rainNC<-nc_open(gcm_rainFile)

#read time first
gcm_time<-ncvar_get(nc=gcm_rainNC, varid = "time")
#convert hours to seconds and then convert it to the calendar time
gcm_date<- as.Date(as.POSIXct(gcm_time*60*60, origin="1900-01-01 00:00:00", tz = "UTC"))


#crop GCM rain for given period
s_index<-which(gcm_date==s_date)
e_index<-which(gcm_date==e_date)
gcm_date<-gcm_date[s_index:e_index]

#Now read the data and subset only for the selected period.
gcm_rain<-ncvar_get(nc = gcm_rainNC, varid = "rain")
gcm_rain<-gcm_rain[s_index:e_index]

#plot it
plot(gcm_date, gcm_rain, type = "h", ylab="Daily Rainfall (mm)", main="input time series of rainfall")

#--------------------read the daily GCM UV MSLP file----------------------------------#
gcm_synNC<-nc_open(gcm_synFile)

#read time first
gcmSyn_time<-ncvar_get(nc=gcm_synNC, varid = "time")

#convert hours to seconds and then convert it to the calendar time
gcmSyn_date<- as.Date(as.POSIXct(gcmSyn_time*60*60*24, origin="1800-01-01 00:00:00", tz = "UTC"))

#crop GCM sunoptic data for given period
s_index<-which(gcmSyn_date==s_date)
e_index<-which(gcmSyn_date==e_date)
gcmSyn_date<-gcmSyn_date[s_index:e_index]

#Now read the data and subset only for the selected period.
gcm_u<-ncvar_get(gcm_synNC, varid="u10")
gcm_u<-gcm_u[s_index:e_index]

gcm_v<-ncvar_get(gcm_synNC, varid="v10")
gcm_v<-gcm_v[s_index:e_index]

gcm_mslp<-ncvar_get(gcm_synNC, varid="msl")
gcm_mslp<-gcm_mslp[s_index:e_index]

#--------------------- Read the Rain and Syn library files --------------------#
rainLib_NC<-nc_open(rainLib_File)
event_mean<-ncvar_get(rainLib_NC, varid = "e_rainMean")

#check the data for integrity
nOverflow <- length(gcm_rain[gcm_rain>max(event_mean)])
if(nOverflow>0) warning(paste(nOverflow, "event(s), having more rainfall than the Radar MaxAccum. --B"))


synLib_NC<-nc_open(filename = synLib_File)
event_u<-ncvar_get(synLib_NC, varid = "u10")
event_v<-ncvar_get(synLib_NC, varid = "v10")
event_p<-ncvar_get(synLib_NC, varid = "mslp")

#================ Select the events================#
tolFactor<- c(0.02, 0.05, 0.1, 0.15, 0.20, 0.5)  #tolerenace factor %
ids<-rep(0L, length(gcm_rain))    # to store match ids


#Here we will make a loop for ensemble of input parameters
for(ens in ens_id){
    print(paste("making parametres for ensemble no.", ens))
    outParm_File<-paste(outParm_File_prefix, ens, ".txt", sep = "")

# We have two criteria 1. total accum and 2. UV and MSLP
for(day in 1:length(gcm_rain)){
    if(gcm_rain[day] < min_rain) next

    #1. find more than 3 days with comparable rainfall or increase the tolerence
    for (i in 1:length(tolFactor)){
        matches<-find.matches(event_mean, gcm_rain[day], tol =gcm_rain[day]*tolFactor[i])
        matchVec<-which(matches$matches==1)
        if(length(matchVec>3)) break
    }

    #2. Find closest day in U-V-P space. This is not normalised hence mslp will have more impact.
    syn_dist<-sqrt((event_u[matchVec]-gcm_u[day])^2 +
                       (event_v[matchVec]-gcm_v[day])^2 + ((event_p[matchVec]-gcm_mslp[day])/100)^2)
    syn_dist<-round(syn_dist)

    matchVec=matchVec[which(syn_dist==min(syn_dist))]

    if(length(matchVec)==1){
        ids[day]=matchVec
        } else {
        ids[day]=sample(x = matchVec, size = 1)
        }
}

# Total timesteps = number of days * number of tsteps each day
ntstep<-length(ncvar_get(rainLib_NC, varid = "tstep"))
tot.tstep<-length(ids)* ntstep

nparm<-length(ncvar_get(rainLib_NC, varid = "prms")) #total parameters in library
prms<-ncvar_get(rainLib_NC, varid = "ts_prms") #read all the parameters


#make a matrix for output parameters
oparm<-matrix(0, nrow =tot.tstep, ncol = nparm+1) # +1 for pMatch
tstep=1  #initialise

#read all the paramters and assign them a day
for(i in 1:length(ids)){
    if(ids[i]==0) next #no id means no rain
    index.dayEnd <- i * ntstep
    index.dayStart <- index.dayEnd+1-ntstep
    #put all the paramters for the day into the ouput array at the end
    oparm[index.dayStart:index.dayEnd, 2:(nparm+1)]=t(prms[, , ids[i]])
}

#set precision for parameters only, not date columns.
oparm[, 3:nparm]<-signif(oparm[, 3:nparm], digits = 3)


#Now we have all the parameters in our output matrix except time
s_time<-as.integer(as.POSIXct(s_date))


intvl<-(86400/ntstep)  #time interval
oparm[, 1]<-s_time+(intvl * 0:(tot.tstep-1))

#get the paramters name to write into the file
pname<-ncatt_get(rainLib_NC, varid = "ts_prms", attname = "param_names")
pname<-unlist(strsplit(pname$value, " "))
pname[1]<-"pMatch"
pname<-c("Time", pname)

#Save this to an ASCII file
write.table(x = signif(oparm, digit=10), file = outParm_File, append = F, quote = F, sep = "\t",
            row.names = F, col.names = pname)

} #Ensemble loop finished







