#!---------------------------------------------------------------------------
# @author Bhupendra Raut
# @brief  This script reads all the parameters (.prm) files from Cascade
# program and makes one .prm file for hidrus-2.
# This will be used with HiDRUS-2 to generate cascade simulations of given radar data.
#===========================================================================

# R libraries used
library(plyr)

#set working directory to parameters directory and give output file name.
setwd("/home/bhupendra/radar_2008-2013/mlb/c_parms_dbz/")
ofName<-"inParm-radar-200801-201508.txt"

options("scipen"=100) #force to print fixed numbers not scientific notations.

#--------- Some initialisation settings----------------#
fn_pat <- "_2.prm"    #file name pattern for .prm files
nprms=65                # number of paramters in the .prm file
#-------------------read the prm data-----------------#
flist<-list.files("./", pattern = fn_pat)       #get the file names
prms<-ldply(flist, read.table, header=T, sep=",")  #read .prm
prms<-prms[,1:nprms]                            #crop the prm data to remove last column

#change column name for the first colum and convert time to POSIXct.
colnames(prms)[1]<- "Time"
prms$Time<- as.POSIXct(as.character(prms$Time), format="%Y%m%d%H%M", tz = "UTC")

prm.dim<-dim(prms)
oparm<-matrix(0, nrow = prm.dim[1], ncol = prm.dim[2]+1)
oparm[, 1]<-prms$Time
oparm[, 2] <-prms$Time
oparm[, 3:(prm.dim[2]+1)] <-as.matrix(prms[, 2:prm.dim[2]])

#set precision for parameters only, not date columns.
oparm[, 3:nprms]<-signif(oparm[, 3:nprms], digits = 3)

#get parameters names
pname<-colnames(prms)
pname<-c("Time", pname)
pname[2]<-"pMatch"

#Save this to an ASCII file
write.table(x = signif(oparm, digit=10), file = ofName, append = F, quote = F, sep = "\t",
            row.names = F, col.names = pname)


