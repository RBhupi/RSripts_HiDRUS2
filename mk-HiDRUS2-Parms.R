#!---------------------------------------------------------------------------
# @author Bhupendra Raut
# @brief We are making a parameter file for HiDRUS-2
#1. This programm reads a reanalysis/GCM/Radar daily fldmean rainfall data.
#2. Reads events library file and search for the closest event for the mean accumulation.
#3. Writes down all the required paramters for the HiDRUS-2 simulation.

#@todo :
# 1. Needs optimisation for the large array operations.
# 2. Add clause to select the day that is closest to the date in GCM to keep seasonality.
#               (Not needed when selected by clusters.)
# use cdo fldmean -sellonlatbox,143.5,146,-36.5,-39 ifile.nc ofile.nc
#==========================================================================================
library(ncdf4)
library(Hmisc) #for find.match

options("scipen"=100) #force to print fixed numbers not scientific notations.

#keep all the files in the same directory for convinience.
setwd("/home/bhupendra/SimTests/hidrus2")
dfName <- "erain1989-2010_50x50km_daysum_fmean_mlb.nc" #daily-rain file name
lfName <- "RainLib_MLB.nc"      #Rain Library file name

ofName<-"inParms.txt"


#-----------------set START and END dates of the simulation.------
s_date=as.Date("2010-06-01", format = "%Y-%m-%d")  # start date of the simulation
e_date=as.Date("2010-07-31", format = "%Y-%m-%d")

#read the daily rain file
dfile<-nc_open(dfName, readunlim = T)

#read time first
dtime<-ncvar_get(nc=dfile, varid = "time")
#convert hours to seconds and then convert it to the calendar time
date<- as.Date(as.POSIXct(dtime*60*60, origin="1900-01-01 00:00:00", tz = "UTC"))
s_index<-which(date==s_date)
e_index<-which(date==e_date)
date<-date[s_index:e_index]

#Now read the data and subset only for the selected period.
drain<-ncvar_get(nc = dfile, varid = "rain")
drain<-drain[s_index:e_index]

#plot it for the viewing purpose
plot(date, drain, type = "h", ylab="Daily Rainfall (mm)", main="input time series of rainfall")

#------------- Read the Rainfall library files ------------------#
lfile<-nc_open(lfName, readunlim = F)
rmean<-ncvar_get(lfile, varid = "e_rainMean")

#check the data for integrity
nOverflow <- length(drain[drain>max(rmean)])
if(nOverflow>0) warning(paste(nOverflow, "event(s), having more rainfall than the Radar MaxAccum. --B"))

#================ Select the events================#
tolFactor<- c(0.02, 0.05, 0.1, 0.15, 0.20, 0.5)  #tolerenace factor
ids<-rep(0L, length(drain))# to store match ids

# We have only one criterion
for(day in 1:length(drain)){
    if(drain[day] < min(rmean)) next #skip if rain is less than the library min
    #find at least 3 matching events for the day. increase the tolerence
    for (i in 1:length(tolFactor)){
        matches<-find.matches(rmean, drain[day], tol =drain[day]*tolFactor[i])
        matchVec<-which(matches$matches==1)
        if(length(matchVec>=3)) break
    }
    ids[day]=sample(x = matchVec, size = 1)
}

# Total timesteps = number of days * number of tsteps each day
ntstep<-length(ncvar_get(lfile, varid = "tstep"))
tot.tstep<-length(ids)* ntstep

nparm<-length(ncvar_get(lfile, varid = "prms")) #total parameters in library
prms<-ncvar_get(lfile, varid = "ts_prms") #read all the parameters


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
pname<-ncatt_get(lfile, varid = "ts_prms", attname = "param_names")
pname<-unlist(strsplit(pname$value, " "))
pname[1]<-"pMatch"
pname<-c("Time", pname)

#Save this to an ASCII file
write.table(x = signif(oparm, digit=10), file = ofName, append = F, quote = F, sep = "\t",
            row.names = F, col.names = pname)









