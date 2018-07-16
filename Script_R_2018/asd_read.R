asd_read <- function(f)
{
# Function to read .asd files.
# Converted from asd_read.m (Matlab) by MEC

# Test if file is a .asd file
if (substring(f,nchar(f)-3,nchar(f)) != ".asd")
{
	print(" Use .asd files only !!!!!!!!!!! (espece de con)")
	return()
}


################################################################
#  support functions.

get_char <- function(fi,num_chars)
{
if (num_chars > 0)
	{
    output_string <- readChar(fi,num_chars,useBytes=TRUE)  #rawToChar(readBin(fi,"raw",num_chars))
	}else
	{
    # The first bytes of referenceheaderSpectrumDescription is the number of character to read
    output_string <- '';
    nchar <- as.numeric(readBin(fi,"raw",1,1))
    rien<- as.numeric(readBin(fi,"raw",1,1))
    if (nchar>0) { output_string <- readChar(fi,nchar,useBytes=TRUE) }  #rawToChar(readBin(fi,"raw",nchar))
	}
return(output_string)
}
################################################################

result <- data.frame(id_spectres=1)

fi <- file(f, "rb")

# Tous les noms sont sans majuscules et sans point car postgresql a du mal avec.

#result <- data.frame(spectrumheaderco=readChar(fi,3,useBytes=TRUE))
result$spectrumheaderco <- readChar(fi,3,useBytes=TRUE)
result$spectrumheadercomments <- readChar(fi,157,useBytes=TRUE) # rawToChar(readBin(fi,"raw",157))

#load up the struct tm object.
result$spectrumheaderwhensec <- readBin(fi, integer(),1,2)
result$spectrumheaderwhenmin <- readBin(fi, integer(),1,2)
result$spectrumheaderwhenhour <- readBin(fi, integer(),1,2)
result$spectrumheaderwhenmday <- readBin(fi, integer(),1,2)
result$spectrumheaderwhenmon <- readBin(fi, integer(),1,2)
result$spectrumheaderwhenyear <- readBin(fi, integer(),1,2)+1900
result$spectrumheaderwhenwday <- readBin(fi, integer(),1,2)
result$spectrumheaderwhenyday <- readBin(fi, integer(),1,2)
result$spectrumheaderwhenisdst <- readBin(fi, integer(),1,2)

#load up the file parameters.
result$spectrumheaderprogram_version <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderfile_version <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderitime <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderdc_corr <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderdc_time <- readBin(fi, integer(),1,4)
result$spectrumheaderdata_type <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderref_time <- readBin(fi, integer(),1,4)
result$spectrumheaderch1_wavel <- readBin(fi, numeric(),1,4)
result$spectrumheaderwavel_step <- readBin(fi, numeric(),1,4)
result$spectrumheaderdata_format <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderold_dc_count <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderold_ref_count <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderold_sample_count <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderapplication <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderwchannels <- readBin(fi, integer(),1,2)

#load up the app_data
result$spectrumheaderapp_data <- readChar(fi,128,useBytes=TRUE)  #rawToChar(readBin(fi,"raw",128))

#load up the GPS data.  NOTE: not used!
result$spectrumheadergps_data <- readChar(fi,56,useBytes=TRUE)  #rawToChar(readBin(fi,"raw",56))

#load up more header fields.
result$spectrumheaderit <- readBin(fi, integer(),1,4)
result$spectrumheaderfo <- readBin(fi, integer(),1,2)
result$spectrumheaderdcc <- readBin(fi, integer(),1,2)
result$spectrumheadercalibration <- readBin(fi, integer(),1,2)
result$spectrumheaderinstrument_num <- readBin(fi, integer(),1,2)
result$spectrumheaderymin <- readBin(fi, numeric(),1,4)
result$spectrumheaderymax <- readBin(fi, numeric(),1,4)
result$spectrumheaderxmin <- readBin(fi, numeric(),1,4)
result$spectrumheaderxmax <- readBin(fi, numeric(),1,4)
result$spectrumheaderip_numbits <- readBin(fi, integer(),1,2)
result$spectrumheaderxmode <- as.numeric(readBin(fi,"raw",1,1))
result$spectrumheaderflags <- paste(as.character(readBin(fi,"raw",4)),collapse=" ")  #readChar(fi,4,useBytes=TRUE)   rawToChar(readBin(fi,"raw",4))
result$spectrumheaderdc_count <- readBin(fi, integer(),1,2)
result$spectrumheaderref_count <- readBin(fi, integer(),1,2)
result$spectrumheadersample_count <- readBin(fi, integer(),1,2)
result$spectrumheaderinstrument <- as.numeric(readBin(fi,"raw",1,1)) #readChar(fi,1,useBytes=TRUE)  #rawToChar(readBin(fi,"raw",1))
result$spectrumheaderbulb <- readBin(fi, integer(),1,4)
result$spectrumheaderswir1_gain <- readBin(fi, integer(),1,2)
result$spectrumheaderswir2_gain <- readBin(fi, integer(),1,2)
result$spectrumheaderswir1_offset <- readBin(fi, integer(),1,2)
result$spectrumheaderswir2_offset <- readBin(fi, integer(),1,2)
result$spectrumheadersplice1_wavelength <- readBin(fi, numeric(),1,4)
result$spectrumheadersplice2_wavelength <- readBin(fi, numeric(),1,4)
result$spectrumheaderwhen_in_ms <- readChar(fi,12,useBytes=TRUE)  #rawToChar(readBin(fi,"raw",12))
result$spectrumheaderspare <- readChar(fi,20,useBytes=TRUE)  #rawToChar(readBin(fi,"raw",20))


#########################################################
# load up the spectrum.
#spectrum <- readBin(fi, double(),result$spectrumheaderwchannels)
spectrum <- t(readBin(fi, double(),result$spectrumheaderwchannels))

#########################################################
# load up the reference header.
result$referenceheaderreferenceflag <- paste(as.character(readBin(fi,"raw",2)),collapse=" ") #readBin(fi, integer(),1,2)
result$referenceheaderreferencetime <- readBin(fi, integer(),1,8)
result$referenceheaderspectrumtime <- readBin(fi, integer(),1,8)
result$referenceheaderspectrumdescription <- get_char(fi,-1);

#########################################################
# load up the reference data.
#reference <- readBin(fi, double(),result$spectrumheaderwchannels)
reference <- t(readBin(fi, double(),result$spectrumheaderwchannels))


close(fi)

return(list(info=result,spectrum=spectrum,reference=reference))
}
