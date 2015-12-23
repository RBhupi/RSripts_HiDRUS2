#!------------------------------------------------------------------------------
# @author Bhupendra Raut
# @brief We are making a parameter file for HiDRUS-2 using GCM data.
#1. This programm reads a GCM daily rainfall and SLP, U-V data from the respective folders.
#2. Reads events and Syn library + Extended Library files and search for the closest event for the mean accumulation.
#3. Writes down all the required paramters for the HiDRUS-2 simulation.

#@todo:
# 1. Needs optimisation for the large array operations.
# 2. add 360days calendar conversion
#===============================================================================
library(ncdf4)
library(Hmisc) #for find.match
library(plyr)
options("scipen"=100) #force to print fixed numbers not scientific notations.

#=============================================================================== Function Definions Starts
#--------------------------------------------------------------------get.fPath()
#get the full path of a file using pattern
get.fPath<-function(path, pat){
    path<-list.files(path=path, pattern = pat, recursive = TRUE, full.names = TRUE)
    print(path)
    if(length(path)!=1){
        stop("more/less than one file found.")
    }
    return(path)
}

#--------------------------------------------------------------get.start_count()
#get start and count vectors for given lat-lon values.
#It also adds the time dim for default reading.
get.start_count<-function(ncfile, lons, lats, lat_name, lon_name){
    lon_ind <- start_count(ncfile, lons, lon_name)
    lat_ind <- start_count(ncfile, lats, lat_name)
    out<-list(c(lon_ind[1], lat_ind[1], 1), c(lon_ind[2], lat_ind[2], -1))
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

#----------------------------------------------------------Lib reading functions
# 1step reading of ncvar. suitable for apply functions
read.NCvar<-function(ncName, varName){
    ncfile<-nc_open(ncName)
    var<-ncvar_get(ncfile, varid = varName)
    nc_close(ncfile)
    invisible(var)
}

#reads 1d variable from list of synLib and rainLib files
read.1dVar<-function(ncName.list, varName){
    vlist<-llply(ncName.list, read.NCvar, varName)
    var<-as.vector(unlist(vlist))
    invisible(var)
}

#reads parameters from LibRain files and returns it as an array
read.prms<-function(ncName.list, varName){
    prms.list<-llply(rainLib_File, read.NCvar, varName)
    dims<-as.vector(unlist(lapply(prms.list, dim))) #get list of all the dims
    num_ids<-sum(dims[seq(3, length(dims), by=3)]) #sum of only 3rd dim
    npar <- dims[1]
    nsteps<-dims[2]
    prms <- array(unlist(prms.list), c(npar,  nsteps, num_ids))
    invisible(prms)
}

ncvar.length<-function(ncName, varName){
    ncfile<-nc_open(ncName)
    ncvar_length<-length(ncvar_get(ncfile, varid = varName))
    nc_close(ncfile)
    return(ncvar_length)
}

make_outParNames<-function(ncName, varName){
    ncfile<-nc_open(ncName)
    #get the paramters name to write into the file
    pname<-ncatt_get(ncfile, varid = varName, attname = "param_names")
    pname<-unlist(strsplit(pname$value, " "))
    pname[1]<-"pMatch"
    pname<-c("Time", pname)
    nc_close(ncfile)
    return(pname)
}



#=============================================================================== Program Starts
#-----------------------------------------------------------------set file paths.
setwd("/home/bhupendra/data/cmip5/temp_pardir/")
gcm_dir <- "/home/bhupendra/data/cmip5" #directory to search GCM files
model <- "CSIRO-Mk3" #"MIROC5"  # as in file names
scenario <- "rcp85"   # c("Hist", "rcp45", "rcp85")
city <- "MLB"        # c("MLB", "SYD", "ADL", "BRI")

# This location will be use to read U-V and SLP data
domain_lats <- c(-38, -38)
domain_lons <- c(146, 146)

gcm_rainFile<- get.fPath(gcm_dir, paste("Precip", city, model, tolower(scenario), "nc", sep = ".*"))
gcm_slpFile <- get.fPath(gcm_dir, paste("psl", model, toupper(scenario), ".nc", sep = ".*"))
gcm_uFile <- get.fPath(gcm_dir, paste("uas", model, toupper(scenario), ".nc", sep = ".*"))
gcm_vFile <- get.fPath(gcm_dir, paste("vas", model, toupper(scenario), ".nc", sep = ".*"))

#get list of all the rainlib and synLib files
rainLib_File <- Sys.glob("~/sim_real/SimLibs/RainLib*2008-14_MLB.nc")
synLib_File<- Sys.glob("~/sim_real/SimLibs/SynLib*2008-14_MLB.nc")

ens_id<-seq(1:100)

min_rain=0.2
#--------------------------------------set START and END dates of the simulation
s_date <- as.Date("2040-01-01", format = "%Y-%m-%d", tz = "UTC")
#s_date <- as.PCICt(as.POSIXct(s_date), cal = "360_day") #change it to model calendar

e_date <- as.Date("2049-12-31", format = "%Y-%m-%d", tz = "UTC")
#e_date <- as.PCICt(as.POSIXct(e_date), cal = "360_day")

outFilePrefix <- paste("h2xprm", city, as.numeric(format(s_date, "%Y")),
                       as.numeric(format(e_date, "%Y")), model,
                       scenario, sep = "_")

#-------------------------------------------------------read the daily rain file
gcm_rainNC <- nc_open(gcm_rainFile)

#read time first
gcm_time <- ncvar_get(nc=gcm_rainNC, varid = "time")
#convert hours to seconds and then convert it to the calendar time
seconds.per.day <- 86400
gcm_date <- as.Date(as.POSIXct(gcm_time*seconds.per.day, origin="1800-01-01 00:00:00", tz = "UTC"))


#crop GCM rain for given period
s_index <- which(gcm_date==s_date)
e_index <- which(gcm_date==e_date)
if(length(s_index)==0 || length(e_index)==0)
    stop("start and/or end date(s) out of range.")
gcm_date<-gcm_date[s_index:e_index]

#Now read the data and subset only for the selected period.
gcm_rain<-ncvar_get(nc = gcm_rainNC, varid = "rain")

#------------------------------------------------------- take area mean and crop
gcm_rain<-apply(gcm_rain, MARGIN = 3, FUN = mean)
gcm_rain<-gcm_rain[s_index:e_index]    #crop for the period
gcm_rain[which(gcm_rain<min_rain)]=0.0 #make light rain zero
#----------------------------------------------------------------------- plot it
pdf("gcm_input.pdf", width=6, height=4)
plot(gcm_date, gcm_rain, type = "h", ylab="Daily Rainfall (mm)", main="input time series of rainfall")


#---------------------------------------------------read daily GCM UV MSLP files
gcm_slpNC<-nc_open(gcm_slpFile)
gcm_uNC<-nc_open(gcm_uFile)
gcm_vNC<-nc_open(gcm_vFile)


gcmslp_time<-ncvar_get(nc=gcm_slpNC, varid = "time")
#convert hours to seconds and then convert it to the calendar time
gcmslp_date<- as.Date(as.POSIXct(gcmslp_time*60*60*24, origin="1800-01-01 00:00:00", tz = "UTC"))

#crop GCM sunoptic data for given period
s_index<-which(gcmslp_date==s_date)
e_index<-which(gcmslp_date==e_date)
gcmslp_date<-gcmslp_date[s_index:e_index]

sc<-get.start_count(ncfile = gcm_slpNC, lons = domain_lons, lats = domain_lats, lon_name="lon", lat_name="lat")

#Now read the data and subset only for the selected period.
gcm_mslp<-ncvar_get(gcm_slpNC, varid="psl", start = sc[[1]], count = sc[[2]])
gcm_mslp<-gcm_mslp[s_index:e_index]

gcm_u<-ncvar_get(gcm_uNC, varid="uas", start = sc[[1]], count = sc[[2]])
gcm_u<-gcm_u[s_index:e_index]

gcm_v<-ncvar_get(gcm_vNC, varid="vas", start = sc[[1]], count = sc[[2]])
gcm_v<-gcm_v[s_index:e_index]

#------------------------------------------------------------------ plot it
plot(gcm_date, gcm_mslp, type = "l", ylab="SLP [hPa]", main="SLP")
plot(gcm_mslp, gcm_u)
plot(gcm_mslp, gcm_v)
dev.off()
#------------------------------------------ Read the Rain and Syn library files
event_mean<-read.1dVar(rainLib_File, "e_rainMean")

#check the data for integrity
nOverflow <- length(gcm_rain[gcm_rain>max(event_mean)])
if(nOverflow>0) warning(paste(nOverflow, "event(s), having more rainfall than the Radar MaxAccum. --B"))

ntstep<-ncvar.length(rainLib_File[1], "tstep")
nparm<-ncvar.length(rainLib_File[1], "prms")

pname<-make_outParNames(rainLib_File[1], "ts_prms")


prms<-read.prms(rainLib_File, "ts_prms") #read all the parameters

event_u<-read.1dVar(synLib_File, "u10")
event_v<-read.1dVar(synLib_File, "v10")
event_p<-read.1dVar(synLib_File, "mslp")

#================ Select the events================#
tolFactor<- c(0.02, 0.05, 0.1, 0.15, 0.20, 0.5)  #tolerenace factor %
ids<-rep(0L, length(gcm_rain))    # to store match ids


#Here we will make a loop for ensemble of input parameters
for(ens in ens_id){
    print(paste("making parameters for ensemble no.", ens, "of", length(ens_id)))
    outParm_File<-paste(outFilePrefix, "_", ens, ".txt", sep="")

    # We have two criteria 1. total accum and 2. UV and MSLP
    for(day in 1:length(gcm_rain)){
        if(gcm_rain[day] < min_rain) next

        event_found=FALSE
        #1. find more than 3 days with comparable rainfall or increase the tolerence
        for (i in 1:length(tolFactor)){
            matches<-find.matches(event_mean, gcm_rain[day], tol =gcm_rain[day]*tolFactor[i])
            matchVec<-which(matches$matches==1)
            if(length(matchVec>3)) {
                event_found=TRUE
                break
            }
        }

        if(!event_found){
            print(paste("no matching event found for gcm day", gcmslp_date[day]))
            print(paste("total GCM rain is ", gcm_rain[day], "| Available event mean is ", max(event_mean)))
            print("assigning the largest event to this day.")
            ids[day] <- which(event_mean==max(event_mean))
            next
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
    tot.tstep<-length(ids)* ntstep


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
    s_time<-as.double(as.POSIXct(s_date))

    intvl<-(86400/ntstep)  #time interval
    oparm[, 1]<-s_time+(intvl * 0:(tot.tstep-1))

    #Save this to an ASCII file
    write.table(x = signif(oparm, digit=10), file = outParm_File, append = F, quote = F, sep = "\t",
                row.names = F, col.names = pname)

} #Ensemble loop finished

#comput the required time for simulations on NCI
tm<-oparm[,9]
rainSteps <- length(tm[tm>0])/1000 #rainy slices in 1000s
nci_time <- (0.14 * rainSteps)-5
dsize <- (0.04 * rainSteps)+0.5

warning(paste("Required CPU time on NCI machine ~", ceiling(nci_time),
              " hours. \nEstimated disc space ~", ceiling(dsize), "GB per file", sep=""))
