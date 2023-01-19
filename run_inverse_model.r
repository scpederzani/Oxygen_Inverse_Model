# Running the passey models
# load libraries
library("matrixcalc")
library("dplyr")
library("ggplot2")
library("beepr")

#### load functions from external script ####
# working directory should be set to source file location

source("functions/mSolv1.R")
source("functions/Emeas.R")

# specify project name for model reports
projectname <- "Oxygen Isotope Inverse Model Example" 

#### read in isotope data ####

# specify tooth to analyse
tooth_name = "tooth1"

# read data input for tooth
tooth_data <- read.csv(paste0("input/", tooth_name, "_inverse_model_input.csv"))
tooth_data

#### arrange from distal to proximal ####
# note: The distal end of the tooth is the occlusal surface
# the proximal end is the cervical margin

tooth_data <- tooth_data %>%
  arrange(desc(dist_ERJ))

# transform sample length data into 10*mm units for model input

tooth_data <- tooth_data %>%
  mutate(sample_length_mmx10 = sample_length*10)

#### define some empty objects ####

numtrials <- 100 # number of model solutions
# some empty objects for code to work
Edist <- rep(0,numtrials)
allTrials <- rep(0,numtrials)

#### running the Emeas function ####

## input variables:

# numtrials - number of trials to run iteration
# numsam - number of isotope samples/values
# length - distance between each sample and the next sample in mm*10, occlusal sample to go first
# dMeas - d18O values, occlusal sample to go first
# r1 - (real) isotope 1 sigma reproducibility
# r2 - (integer) sample length 1 sigma reproducibility (in mm * 10)
# r3 - (integer) sample depth 1-sigma reproducibility (in mm * 10)
# la - (integer) length of apposition (in mm * 10) --> get species specific data from literature
# here la data for Bison is used

r1 <- 0.3
r2 <- 1
r3 <- 2
la <- 15

## run the Emeas function to generate Emeas estimates

Emeasout <- Emeasfun(numtrials = numtrials, 
                     numsam = length(tooth_data$mean_d18O), 
                     length = tooth_data$sample_length_mmx10, 
                     dMeas = round(tooth_data$mean_d18O, 1), 
                     r1 = r1, 
                     r2 = r2, 
                     r3 = r3, 
                     la = la)

# write input parameters into a list to print in model reports
Emeas_params <- list(numtrials = numtrials, r1 = r1, r2 = r2, r3 = r3, la = la)

#### plot output ####  
ggplot(Emeasout, aes(x = totallength, y = allTrials))+
  theme_classic()+
  xlim(500,0)+
  geom_line()+
  geom_point(shape = 3)+
  geom_point(aes(x = totallength, y = dMeas), colour = "black", shape = 0)+
  geom_line(aes(x = totallength, y = dMeas), colour = "black")+
  geom_line(aes(x = totallength, y = dMeasError), colour = "black")+
  geom_point(aes(x = totallength, y = dMeasError), colour = "black", shape = 1)

#### evaluate Emeas to adjust df in mSolv code ####
# the distribution of Emeas values is stored in an object called 'Edist'
# values of Emeas/Edist should be close to DPE values from mSolv code
# adjust df stepwise to match Emeas/Edist to DPE

hist(Edist)
mean_edist <- mean(Edist)
mean_edist


#### running the mSolv1_1 function ######


#### input parameters ####

nsolxns <- 200 # number of solutions to be computed
dMeas <- round(tooth_data$mean_d18O, 1) # isotope data input
numsam <- length(dMeas) # number of samples
openindx <- 1 # degree of open-endedness, less than lm, openended (profile mature) --> index = 1; 
# close ended (enamel immature) index = lm
length <- tooth_data$sample_length_mmx10 # length in mm*10
avelength <- round(mean(tooth_data$sample_length_mmx10), digits = 0) # average sample length
# lm = length of maturation in mm*10 (from literature, species specific)
lm <- 250
#lamm = length of apposition in mm
lamm <- la/10
# enter sample depth as fraction of la in mm*10
# in the example data depth is given as measurement in mm
enamel_thickness_mm <- 1 # enamel thickness in mm of your tooth specimen (measured here at occlusal surface)
sample_depth_fraction <- tooth_data$sample_depth/enamel_thickness_mm # calculate fraction of enamel thickness
depth <- round(sample_depth_fraction*lamm*10,0) # multiply fraction with la in mmm, then * 10 for mm*10 units

# if depth is given as fraction of enamel thickness, use this option instead:
# depth <- round(tooth_data$sample_depth*lamm, 1)*10


finit <- 0.25 # initial mineral %weight (from literature, species specific)

# some derived numbers used during modelling

numbefore <- ceiling(la/avelength)
numafter <- ceiling((lm-openindx)/avelength)+1
# empty object for DPE to enable writing to global environment from inside function
MEST = matrix(0, nrow = (numsam+numbefore+numafter-1), ncol = nsolxns)
DPE = rep(0,nsolxns)
S = matrix(0, nrow = numsam, ncol = nsolxns)

# input parameters of the reference vector

maxlength <- 40 # approximately maximum sample length in data (mm * 10)
minlength <- 3 # approximately minimum sample length in data (mm * 10)
mindepth <- 5 # approximately minimum sample depth in data (mm * 10)
maxratio <- 17.5  # maxratio - max value of the generated reference vector (choose close to average of actual isotope value)
minratio <- 15.5 # minratio - min value of the generated reference vector (choose close to average of actual isotope value)
stdev <- 0.5 # stdev - sd for random draws that produce reference vector
df <- 0.001 # damping factor


## after adjusting df, the model needs to be tested for sensitivity to the reference vector
## --> check whether more oscillating (max and min further apart) generates drastically different solution
# df - damping factor. Needs to be chosen to minimize difference between the estimated measurement error
# (E~meas~) and the prediction error (E~pred~)

# write parameter input into a list for printing the model report
mSolvparams <- list(nsolxns = nsolxns, openindx = openindx, lm = lm, la = la, finit = finit, maxlength = maxlength, minlength = minlength, 
                    mindepth = mindepth, r1 = r1, r2 = r2, r3 = r3, maxratio = maxratio, minratio = minratio, stdev = stdev, df = df)

### run mSolv function ###

lfsolv <- PasseyInverse(Length = length,dMeas = dMeas,depth = depth,la = la,lm = lm, maxlength = maxlength, minlength = minlength, 
                      mindepth = mindepth, df = df, nsolxns = nsolxns, finit = finit, openindx = openindx, avelength = avelength, 
                      r1 = r1, r2 = 2, r3 = r3, maxratio = maxratio, minratio = minratio, stdev = stdev, numsam = numsam)
beep(sound = 2)
Sys.sleep(1)

# plot some example solutions of mSolv
# reference vector in black
# original data with numafter extension in pink


ggplot(solvout, aes(x = totallength, y = dMeasd))+
  theme_bw()+
  geom_point(colour = "black")+
  geom_line(colour = "black")+
  geom_line(aes(x = totallength, y = dMeas), colour = "pink", lwd = 2)+
  geom_line(aes(x = totallength, y = trial1), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial2), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial3), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial4), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial5), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial6), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial7), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial8), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial9), colour = "lightblue")+
  geom_line(aes(x = totallength, y = trial10), colour = "lightblue")

#### evaluate DPE (= Epred, prediction error) vs Edist (= Emeas, estimated measurement error) ####

# DPE should be close to Edist
# adjust df to achieve this
# increase df to increase DPE

hist(DPE, breaks = 10)
mean(DPE)
# compare to measurement error calculated in the first part of the model
mean_edist

#### generate report of model parameters ####

rmarkdown::render("make_model_report.Rmd", params = list(
  project = projectname, 
  tooth = tooth_name, 
  Emeasparams = Emeas_params, 
  Edist = Edist, 
  
  mSolvparams = mSolvparams, 
  DPE = DPE
), output_file = paste0("output/",tooth_name, "_inverse_model_report.html"))

#### extract output with mean and CI and write to csv ####

# extract only trial output

tdata <- solvout %>%
  select(-totallength)

# transpose
tdatatrans <- t(tdata)
tdatatrans

nl <- 1:(numsam + numbefore + numafter - 1)
lengths <- paste("l",nl,sep="")

colnames(tdatatrans) <- lengths

# calculate mean and confidence interval (95%)

tdata_lci <- apply(tdatatrans, MARGIN = 2, quantile, prob = 0.025, na.rm = T) # calculate lower boundary of Conf Interval 
tdata_uci <- apply(tdatatrans, MARGIN = 2, quantile, prob = 0.975, na.rm = T) # calculate upper boundary of Conf Interval 
tdata_mean <- apply(tdatatrans, MARGIN = 2,mean) # calculate mean (essentially the most likely model solution)

# bind into data frame and clean up names
tdataci <- cbind(solvout$totallength, tdata_lci, tdata_uci, tdata_mean)
tdataci
colnames(tdataci) <- c("ci_length","lower_CI","upper_CI","mean")

# convert to data frame
tdataci_d <- as.data.frame(tdataci)

# convert lengths back into mm
#tdataci_d$ci_length <- tdataci_d$ci_length/10

# combine all output
all_out <- cbind(solvout, tdataci_d)

all_out <- all_out %>%
  mutate(tooth = tooth_name) %>%
  mutate(layer = tooth_data$layer[1])

# write output to csv
outputfile <- paste("output/",tooth_name, "_inverse_model_output.csv", sep = "")
write.csv(all_out, file = outputfile)

#### plot model output ####

# crop dMeas to original length, by removing numafter rows added in the mSolv script
# also add record columns to label measured d18O data and reference vector in the plot

dMeas_cropped <- solvout %>%
  select(totallength, dMeas) %>%
  slice_head(n = nrow(solvout)-numafter) %>%
  mutate(record = "measured d18O")

all_out <- all_out %>%
  mutate(record = "reference vector")

# axis label expressions and plot colors
ylab_name_o <- expression("Enamel "*delta^18*'O (\u2030 VPDB)' )
plotcolors <- c("#51B2C3", "#C04B2C")

# plot

model_plot <- ggplot(all_out, aes(x = ci_length, y = dMeas))+
  theme_classic()+
  theme(legend.position = "bottom", 
        strip.background = element_rect(colour = "white", fill = "grey95", size = 1), 
        panel.background=element_blank(), 
        panel.border=element_rect(colour="grey95",size = 0.1, fill = NA), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        strip.text = element_text(size = 12))+
  ylab(ylab_name_o)+
  xlab("model total length (mm)")+
  geom_ribbon(aes(x = ci_length, ymin = lower_CI, ymax = upper_CI, linetype = "solid"), 
              fill = "#094860", alpha = 0.1)+
  geom_line(aes(x = ci_length, y = mean, linetype = "solid"), lwd = 1, 
            colour = "#094860", alpha = 0.3)+
  scale_linetype_manual(name = "", values = c("solid"),labels = c("estimated input"))+
  #geom_line(aes(x = ci_length, y = dMeas, colour = tooth), lwd = 2)+
  geom_line(aes(x = ci_length, y = dMeasd, color = record), lwd = 1, alpha = 0.5)+
  geom_point(aes(x = ci_length, y = dMeasd, fill = record), 
             colour = "black", shape = 21, size = 2, alpha = 0.5)+
  geom_line(data = dMeas_cropped, aes(x = totallength, y = dMeas, color = record), lwd = 2)+
  geom_point(data = dMeas_cropped, aes(x = totallength, y = dMeas, fill = record), 
             colour = "black", shape = 21, size = 4)+
  scale_color_manual(values = plotcolors, 
                     labels = c(expression('measured '*italic(delta)^18*'O  '), 
                                "reference vector"), 
                     name = "")+
  scale_fill_manual(values = plotcolors, 
                    labels = c(expression('measured '*italic(delta)^18*'O  '), 
                               "reference vector"), 
                    name = "")+
  guides(fill = guide_legend(order = 1, override.aes = list(shape = 21)))+
  guides(colour = guide_legend(order = 1))+
  ggtitle(label = tooth_name)+
  NULL

ggsave(model_plot, filename = paste("output/",tooth_name, "_inverse_model_plot.png", sep = ""), 
       device = "png", width = 20, height = 15, units = "cm", dpi = 150)



