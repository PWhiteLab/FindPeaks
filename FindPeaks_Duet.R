########
# THIS IS FOR ABR FILE EXPORTED FROM IHS DUET.
# 
# This script finds peaks for individual animal or a group of animals.
# It recognizes the waves in the ABR traces, and returns the latency.
# and amplitude of first FIVE waves.
# It takes ASCII files from IHS SmartEP as input, and exports two files for each
# animal: csv file with wave amplitude and latency, pdf file with the waveform 
# in which peak(red) and trough(blue) for each wave are labelled.
# Important helper files: Time.csv with recording time information of IHS system,
# Info.csv with animal information, col 1 = ID, col 2 = Sex, col 3 = Genotype.
# 
# 
#'[Important: revising the results is highly recommended, pay attention to three]
#'[things - the first peak recognized, the trough of Wave IV, and if there are ] 
#'[other mislabeled peaks in the first five ones recognized.]
#'
#'[Note: When exporting ASCII files, only traces should be present on the page,]
#'[labels should be removed. Export one frequency at a time]
#
# source("FindPeaks.R")
# Example usage 1:
# FindPeaks_group("folder_name") # find peaks for a group of animals
# 
# Example usage 2:
# FindPeaks_single("PATH/file_name") # find peaks for one animal
#
# Example usage 3:
# See_trace("PATH/file_name") # export traces of a specific animal, find latency and amplitude interactively
# 
# Example usage 4:
# Compile("folder_name") # put all (revised) csv outcomes together, requires "Info.csv"
#
# algorithm is adapted from William A. Huber (http://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset).
# author: Daxiang Na (daxiang_na@urmc.rochester.edu)
#
# ver 2022.04.14
########

library(tidyverse)
library(plotly)
library(zoo)

# Define ABRs that begin (min) at 1 (ms) and end (max) at 8 (ms), amplitude threshold > 0.1 (μV). 
# Change those parameters as you need.

min <- 1 # beginning of waves
max <- 8 # end of waves
threshold <- 0.1 # amplitude threshold of waves

# Parameters for peak/trough labelling: w is the half-width of the window used to compute the local maximum (to give you a sense, in the IHS Smartbox system, 18 = 0.56 ms); span determines how much you want to smooth the ABR trace, should be between 0 and 1, higher -> smoother.
#'[w = 2, span = 0.05 yield best result in Click-ABR peak labelling so far.]
w <- 2
span <- 0.05

# helper function to find the peaks within defined range
argmax <- function(x, y, w=1, ...) {
        require(zoo)
        n <- length(y)
        y.smooth <- loess(y ~ x, ...)$fitted
        y.max <- rollapply(zoo(y.smooth), 2*w+1, max, 
                           align="center")
        delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
        i.max <- which(delta <= 0) + w
        df <- data.frame(x=x[i.max], i=i.max)
        df <- df[df$x >= min & df$x <= max,]
        list(x=df$x, i=df$i, y.hat=y.smooth)
}

# helper function to find the troughs within defined range
argmin <- function(x, y, w=1, ...) {
        require(zoo)
        n <- length(y)
        y.smooth <- loess(y ~ x, ...)$fitted
        y.min <- rollapply(zoo(y.smooth), 2*w+1, min, align="center")
        delta <- y.min - y.smooth[-c(1:w, n+1-1:w)]
        i.min <- which(delta >= 0) + w
        df <- data.frame(x=x[i.min], i=i.min)
        df <- df[df$x >= min & df$x <= max,]
        list(x=df$x, i=df$i, y.hat=y.smooth)
}

# ASCII_extract <- function(file) {
#         data <- read.csv(file, row.names = 1, header = TRUE)
#         # Clean data frame
#         df <- data[, str_which(data[1,],"LA")+1]
#         colnames(df) <- data[1, str_which(data[1,],"LA")]
#         # Extract numbers from trace file names
#         colnames(df) <- parse_number(colnames(df))
#         df <- df[c(21:nrow(df)),]
#         df <- df[,order(ncol(df):1)]
#         Time <- read.csv("Time.csv")
#         Waveform <- cbind(Time, df) %>% 
#                 filter(Data_Pnt_ms >= 0 & Data_Pnt_ms <= 10) %>%
#                 mutate_all(as.numeric)
#         return(Waveform)
# }

ASCII_extract <- function(file) {
        data <- read.csv(file, row.names = 1, header = TRUE)
        df <- data.frame(data[, str_which(data[1,],"LA")+1])
        colnames(df) <- data[5, str_which(data[1,],"LA")]
        df <- df[c(22:nrow(df)),]
        # if (is.null(ncol(df)) == FALSE) {
        #         df <- df[,order(ncol(df):1)]  
        # }
        Time <- read.csv("Time.csv")
        Waveform <- cbind(Time, df) %>% 
                filter(Data_Pnt_ms >= 0 & Data_Pnt_ms <= 10) %>%
                mutate_all(as.numeric)
        return(Waveform)
}

Sort_Peak <- function(Wave_data,...) {
        ls <- list()
        for (i in 1:length(which(Wave_data$Type == "Peak"))) {
                Wave <- c(Amp = Wave_data[which(Wave_data$Type == "Peak")[i],2] - 
                                  Wave_data[which(Wave_data$Type == "Peak")[i]+1,2],
                          Lat = Wave_data[which(Wave_data$Type == "Peak")[i],1], 
                          peak_id = Wave_data[which(Wave_data$Type == "Peak")[i],3],
                          trough_id = Wave_data[which(Wave_data$Type == "Peak")[i]+1,3])
                ls[[i]] <- Wave
        }
        
        list <- c()
        for (k in 1:length(ls)) if (is.na(ls[[k]]["Amp"])) {
                list <- c(list, k)
        } else if (ls[[k]]["Amp"] < threshold){
                list <- c(list, k)
        }
        
        if (is.null(list) == FALSE) ls <- ls[-list]
        
        if (length(ls) < 5) {
                for (i in 1:(5-length(ls))) {
                        ls[[5-(i-1)]] <- c(Amp = NA, Lat = NA, peak_id = NA, trough_ID = NA)
                }
        }
        
        return(ls)
}

See_trace <- function(file) {
        animalID <- tools::file_path_sans_ext(basename(file))
        Waveform <- ASCII_extract(file)
        ls <- list()
        annotations <- list()
        for (j in c(3:length(names(Waveform)))) {
                ls[[j-2]] <- plot_ly(x = Waveform$Data_Pnt_ms, y = Waveform[,j], type = 'scatter', mode = 'lines')%>%
                        add_trace()%>%
                        layout(showlegend = F,
                               annotations = list(text = paste("Sound Level = ", colnames(Waveform)[j]," dB", sep = "")))
        }
        fig <- subplot(ls, nrows = length(ls))
        htmlwidgets::saveWidget(as_widget(fig), paste(animalID, ".html", sep = ""))
}

# Function to find waves for individual animal
FindPeaks_single <- function(file) {
        animalID <- tools::file_path_sans_ext(basename(file))
        Waveform <- ASCII_extract(file)
        individual_final <- data.frame()
        par(mar=c(1,1,1,1))
        pdf(paste(animalID,".pdf",sep = ""),10,20)
        par(mfrow=c(8,2))
        # read data for each sound level
        for (j in c(3:length(names(Waveform)))) {
                min <- min
                max <- max
                threshold <- threshold
                index <- j
                x <- Waveform$Data_Pnt_ms
                y <- Waveform[,j]
                Peaks <- argmax(x, y, w = w, span = span)
                Troughs <- argmin(x, y, w = w, span = span)
                Peaks_data <- Waveform[Peaks[["i"]], c(2,index)] %>% mutate(i = Peaks[["i"]], Type = "Peak")
                Troughs_data <- Waveform[Troughs[["i"]], c(2,index)] %>% mutate(i = Troughs[["i"]], Type = "Trough")
                Wave_data <- rbind(Peaks_data, Troughs_data) %>% arrange(Data_Pnt_ms)
                ls <- Sort_Peak(Wave_data)
                # make data frame
                result <- data.frame(Sound_Level = colnames(Waveform)[j], 
                                     WaveILat = ls[[1]]["Lat"], WaveIAmp = ls[[1]]["Amp"], 
                                     WaveIILat = ls[[2]]["Lat"], WaveIIAmp = ls[[2]]["Amp"], 
                                     WaveIIILat = ls[[3]]["Lat"], WaveIIIAmp = ls[[3]]["Amp"], 
                                     WaveIVLat = ls[[4]]["Lat"], WaveIVAmp = ls[[4]]["Amp"],
                                     WaveVLat = ls[[5]]["Lat"], WaveVAmp = ls[[5]]["Amp"],
                                     ID = animalID)
                ls <- do.call(rbind,ls) %>% as.data.frame()
                plot(x, y, lwd = 2, col="black", main=paste("Sound Level = ", colnames(Waveform)[j]," dB", sep = ""), type = 'l')
                if (all(is.na(ls)) == FALSE) {
                        points(x[ls$peak_id], y[ls$peak_id], col="Red", pch=19, cex=1)
                        points(x[ls$trough_id], y[ls$trough_id], col="Blue", pch=19, cex=1)
                        text(x[ls$peak_id], y[ls$peak_id], labels = c(1:length(ls$peak_id)), cex = 1.5, adj = -0.5)
                        text(x[ls$trough_id], y[ls$trough_id], labels = c(1:length(ls$trough_id)), cex = 1.5, adj = -0.5)
                }
                individual_final <- rbind(individual_final, result)
        }
        dev.off()
        write.csv(individual_final, paste(animalID, ".csv", sep = ""), row.names = FALSE)
}

# Function to find waves for a group of animals
FindPeaks_group <- function(directory) {
        file_list <- list.files(directory, full.names = TRUE)
        for (i in 1:length(file_list)) {
                file <- file_list[i]
                FindPeaks_single(file)
        }
}

# helper function to put revised csv files into one master sheet
Compile <- function(directory) {
        file_list <- list.files(directory,full.names = TRUE)
        info <- read.csv("Info.csv")
        df <- data.frame()
        for (i in 1:length(file_list)) {
                data <- read.csv(file_list[i], header = TRUE)
                animalID <- data[1,"ID"]
                animalInfo <- info[info[,1]==animalID,][,2:3] # 1 = ID, 2 = Sex, 3 = Genotype
                process <- data.frame(WaveIItoILatDif = data$WaveIILat - data$WaveILat, 
                                      WaveIItoIAmpRat = data$WaveIIAmp / data$WaveIAmp, 
                                      WaveIIItoILatDif = data$WaveIIILat - data$WaveILat, 
                                      WaveIIItoIAmpRat = data$WaveIIIAmp / data$WaveIAmp, 
                                      WaveIVtoILatDif = data$WaveIVLat - data$WaveILat, 
                                      WaveIVtoIAmpRat = data$WaveIVAmp / data$WaveIAmp,
                                      WaveVtoILatDif = data$WaveVLat - data$WaveILat, 
                                      WaveVtoIAmpRat = data$WaveVAmp / data$WaveIAmp)
                data <- cbind(data, animalInfo, process, row.names = NULL)
                df <- rbind(df, data)
        }
        write.csv(df, paste(directory, ".csv", sep = ""), row.names = FALSE)
}
