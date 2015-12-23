#!---------------------------------------------------------------------------
# @author Bhupendra Raut
# @brief  This script reads all the parameters (.prm) files from Cascade
# program and makes one netcdf file from it.
# The netcdf file contains several 24-hour time series of rainfall
# and parameters, sampled from the .prm files. It also saves the exact time
# of the parameters as variable pm.
# This will be used with HiDRUS-2 to generate cascade simulations.
# The 'X2'and 'X3' means ts was repeated two/three times
#===========================================================================

# R libraries used
library(plyr)
library(ncdf4)

#set working directory to parameters directory and give output file name.
setwd("/home/bhupendra/radar_2008-2013/mlb/c_parms_dbz")
ofileName<-"RainLib-X2_2008-14_MLB"  #no extension

#--------- Some initialisation settings----------------#
window <- 240           #number of point in sampled ts
cut_off=0.1             #mean area rain mm/hr
duration=10              #duration (number of tstep) of rain to be a valid ts
nprms=65                # number of paramters in the .prm file
fn_pat <- "???[!5]*_2.prm"    #file name pattern for .prm files
id_number=10000  #start identity number for the event
case_number=0 # actual number of the event.
extd_mean<-5  # only event for more rain than value will be  selected

#-------------------read the prm data-----------------#
flist<-Sys.glob(paths = fn_pat)       #get the file names
prms<-ldply(flist, read.table, header=T, sep=",")  #read .prm
prms<-prms[,1:nprms]                            #crop the prm data to remove last column

#change column name for the first colum and convert time to POSIXct.
colnames(prms)[1]<- "Time"
prms$Time<- as.POSIXct(as.character(prms$Time), format="%Y%m%d%H%M", tz = "UTC")

prms.dim<-dim(prms)                                   # get the length of prm data

#------------create netcdf file---------------------#
#define nc dimensions
id <- ncdim_def("id", units = "", vals = c(1:1), unlim = T, longname = "identity number for the event")  #identity number as dimension
e_prms<-ncdim_def("prms", units="", val=seq(1:nprms)) #parameters count as dimensions
tstep<-ncdim_def("tstep", units="", val=seq(1:window), longname = "time steps") #timestep as dimension

#define nc variables
e_rainAccum <- ncvar_def("e_rainMean", units = "mm", dim = id, longname = "total rainfall for the event" )
e_rainMax <- ncvar_def("e_rainMax", units = "mm/hr", dim = id , longname = "Max area average rain for the event")
e_rainDur <- ncvar_def("e_rainDur", units = "tstep", dim = id, longname="Duration of the event" )
e_date <- ncvar_def("e_date", units = "days since 1970-01-01 00:00:0.0 UTC", dim = id, longname = "date of the central point of the event" )
ts_prms<-ncvar_def("ts_prms", units="", dim=list(e_prms, tstep, id),
                   longname = "cascade parameters", prec = "double")

#Create the output file with .nc extension
ofile<-nc_create(paste(ofileName, ".nc", sep=""), vars =list(e_date, e_rainAccum, e_rainMax,
                                                             e_rainDur, ts_prms))
#write correct parameter names to the file
c_names<-paste(colnames(prms), sep="", collapse = " ")

ncatt_put(nc = ofile, varid = ts_prms, attname = "param_names", attval = c_names, prec = "text")
print(paste("Output netcdf file", ofileName, "created.", sep=" "))
#-----------------------------------------------------------

# Open a pdf device for ploting sampled ts
pdf(paste(ofileName, ".pdf", sep=""), width=8, height=6)
print(paste("Output pdf file", ofileName, "created.", sep=" "))
par(mfrow=c(3, 3))

isFirst=TRUE    #always set to true
window2 <- window/2  #a window for searching events

#starting the search loop with double size window
for (i in 1:(prms.dim[1]-window2))
{
    #subset the rainmean for 'ts' of the size of given window.
    ts1<-prms$rainMean[i:(i+window2-1)]
    ts1[ts1<cut_off] <- 0.0           #set smaller values to '0'.

    # A valid 'ts' should satisfy following criteria.
    # 1. the first and the last point of the ts shouldn't be heavily raining.
    # 2. there should be enough raining points
    #and 3. mean daily rain should be more than threshold.
    if(ts1[1] >= cut_off || ts1[length(ts1)] >= cut_off || length(ts1[ts1>=cut_off]) < duration ||  mean(ts1*24) < extd_mean )
    {
        next #skip this and go to the next iteration
    }


    if(isFirst)   #if this is the first valid 'ts' we have found then save it.
    {
        isFirst=FALSE      #set to false
        ts2 <- ts1     #save it as ts2 for future comparison.

        #save all its parms in a dataframe, twice so that we will have double rainfall
        df_out <- rbind(prms[i:(i+window2-1),], prms[i:(i+window2-1),], deparse.level = 1)

        #comput basic stats to write in to the file
        id_number <- id_number + 1 #give it an identity number id=1
        case_number <- case_number +1
        date1 <-  as.Date(df_out$Time[window/2]) #date of the event
        rain_acc <- mean(df_out$rainMean*24)  #event mean rainfall
        rain_max <- max(df_out$rainMean)      #event Max rainfall
        rain_dur <- length(df_out$rainMean[df_out$rainMean>cut_off]) #duration in timestep

        #Save calendar date as nc date prepare parameters for output
        df_out$Time <- as.integer(df_out$Time)
        outParms<-as.vector(unlist(t(df_out)))  #transpose, unlist and vectorise the prm data

        #Write all the event's data in the output file
        ncvar_put(nc = ofile, varid = id, vals = id_number, start = case_number, count = 1)
        ncvar_put(nc = ofile, varid = e_date, vals = as.double(date1), start = case_number, count = 1)
        ncvar_put(nc = ofile, varid = e_rainAccum, vals = rain_acc, start = case_number, count = 1)
        ncvar_put(nc = ofile, varid = e_rainMax, vals = rain_max, start = case_number, count = 1)
        ncvar_put(nc = ofile, varid = e_rainDur, vals = rain_dur, start = case_number, count = 1)
        ncvar_put(nc=ofile, varid=ts_prms, vals = outParms,
                  start = c(1, 1, case_number), count = c(nprms, window, 1))

        #and plot ts for viewing purpose
        plot(1:240, df_out$rainMean, type="l", xlab="time steps", ylab="mean area rain")
        mtext(paste("id=", id_number, ", date:", date1, ", duration =", rain_dur, sep=""), side = 3, cex = 0.5)
        grid()
        next
    }


    # else, when there is earlier ts, check if the ts1 and ts2
    # are actually different rain events. They should have different
    # number of rainy points. NOTE: other criterion can be added here.
    if(length(ts1[ts1>cut_off])==length(ts2[ts2>cut_off]))
    {
        next
    }

    # We have found another event. Lets write it into the nc file.
    df_out <- rbind(prms[i:(i+window2-1),], prms[i:(i+window2-1),], deparse.level = 1)


    id_number=id_number+1              #give it an identity number
    case_number <- case_number +1
    date1 <-  as.Date(df_out$Time[window/2]) #date of the event
    rain_acc=mean(df_out$rainMean*24)  #event mean rainfall
    rain_max=max(df_out$rainMean)      # event max rain
    rain_dur=length(df_out$rainMean[df_out$rainMean>cut_off]) #duration in timestep

    #Save calendar date as nc date prepare parameters for output
    df_out$Time <- as.integer(df_out$Time)
    outParms<-as.vector(unlist(t(df_out)))  #transpose, unlist and vectorise the prm data

    #Write all the event's data in the output file
    ncvar_put(nc = ofile, varid = id, vals = id_number, start = case_number, count = 1)
    ncvar_put(nc = ofile, varid = e_date, vals = as.integer(date1), start = case_number, count = 1)
    ncvar_put(nc = ofile, varid = e_rainAccum, vals = rain_acc, start = case_number, count = 1)
    ncvar_put(nc = ofile, varid = e_rainMax, vals = rain_max, start = case_number, count = 1)
    ncvar_put(nc = ofile, varid = e_rainDur, vals = rain_dur, start = case_number, count = 1)
    ncvar_put(nc=ofile, varid=ts_prms, vals = outParms,
              start = c(1, 1, case_number), count = c(nprms, window, 1))

    #plot it
    plot(1:240, df_out$rainMean, type="l", xlab="time steps", ylab="mean area rain")
    mtext(paste("id=", id_number, ", date:", date1, ", duration =", rain_dur, sep=""), side = 3, cex = 0.5)
    grid()



    #save this event as an old event
    ts2<-ts1
}#for loop ends

print(paste(case_number, "valid ts generated!!!"))

nc_close(ofile) #close the nc file
dev.off()       #close the pdf device
