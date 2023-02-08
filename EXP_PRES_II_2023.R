# //////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# https://www.tec.ac.cr
# Session: MEDICION DE PRESION

# M.Sc. Eng. Maikel Mendez M
# Water Resources + GIS + DataScience
# Instituto Tecnologico de Costa Rica
# https://www.tec.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://github.com/maikelonu
# //////////////////////////////////////////////////////////////////////////////////

# References: Lawrence Lin, A. S Hedayat, Bikas Sinha, Min Yang. 
# Statistical Methods in Assessing Agreement:
# Models, Issues, and Tools
# Journal of the American Statistical Association. 
# March 1, 2002, 97(457): 257-270.

# //////////////////////////////////////////////////////////////////////////////////
# INFO
# Pruebas_Concordancia (Accuracy, Precision, CCC, TDI, CP)
# t-test
# Shapiro_test
# Estadisticas_descriptivas
# Correlacion_Pearson

# //////////////////////////////////////////////////////////////////////////////////
# Q= flow (GPM)
# YB= bourdon gauge (kgf/cm2)
# XD= digital gauge (kgf/cm2)

# 1 kgf/cm2 = 10 mca
# 0.1 kgf/cm2 = 1 mca
# 0.01 kgf/cm2 = 0.1 mca = 10 cm water, 
# As uncertainty for YB=  0.025 kgf/cm2 (25 cm water), what's your target???
# Can we replace X with Y?
# //////////////////////////////////////////////////////////////////////////////////

# Workspace is cleared
rm(list = ls())

# Working directory is selected
setwd("/mnt/BACKUP/R_ITC/R_LABHYD/EXP_PRES")

# CRAN libraries are loaded
require(Agreement)
require(DescTools)
require(ggplot2)
require(MASS)
require(pastecs)
require(reshape)

# /////////////////////////////////////////////////////////////
# BLOCK: Custom function, round data.frame to specif digits
# /////////////////////////////////////////////////////////////
round_df <- function(df, digits) {
  options(scipen = 0)
  options(scipen = -3)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return(df)
}

# Agreement function source is loaded
options(scipen=999)
source("agreement_2020.R")

# Input data is loaded and a data.frame is created
df.base <- read.table("test_manometros.txt", header = TRUE)

# Desc {DescTools} function is requested
Desc(df.base, plotit = TRUE)

# names function is requested
names(df.base)

# a two sample t-test allows us to test whether the MEANS of 
# two independent groups are different. 
# H0 (null hypothesis): The true probability of success is equal to what we expected
# H1 (alternative hypothesis): The true probability of success is not equal to what we expected
# If (p-value > 0.05) the null hypothesis CANNOT be reject
t.test(df.base$XD, df.base$YB)

# A simple deviation (SD) column is created
df.base$SD <- (df.base$XD - df.base$YB)

# shapiro.test {stats} Normality Test is applied to df.base$SD
# if p-value > 0.05 normality stands true, meaning that
# the comparison is parametric
shapiro.test(df.base$SD)

# Mean absolute deviation (MAD) column is created
df.base$MAD <- abs(df.base$XD - df.base$YB)

# Mean squared deviation (MSD) column is created
df.base$MSD <- (df.base$XD - df.base$YB)^2

# A standard R summary is requested
summary(df.base)

# Descriptive statistics are requested and round to 5 decimals
df.base.desc <- round(stat.desc(df.base), 5)

# mean-MAD is requested from df.base data.frame
df.base.desc[9,5]

# mean-MSD is requested from df.base data.frame
df.base.desc[9,6]

# A sq(mean-MSD) is calculated
(df.base.desc[9,6])^0.5

# A correlation test is executed for X vs Y
cor.test(df.base$XD,df.base$YB)

# An agreement {Agreement} function is executed
test1 <- agreement(x = df.base$XD,
                   y = df.base$YB,
                   error = "const",
                   target = "fixed",
                   CCC_a = 0.975, # no more of 2.5% should be accepted for instruments comparison
                   TDI_a = 0.025, # an absolute difference when error="const" (units-interval)
                   alpha = 0.05, # 100(1-alpha)
                   CP_a = 0.90, # mirrored of TDI usually taken as 0.90 
                   H_label = "Digital Gauge",
                   V_label = "Bourdon Gauge")

# An agreement-list summary is requested
summary.agreement.Allowance(test1)
summary.agreement.Estimate(test1)
summary.agreement.Conf_Limit(test1)

# Accuracy:
# Accuracy refers to the closeness of a measured value to a known value (or standard),
# where 0 represents no agreement and 1 represents perfect agreement.
#
# Precision:
# Precision refers to the closeness of two or more measurements to each other,
# where 0 represents no agreement and 1 represents perfect agreement.
# For example, if on average, your measurements for a given substance are close 
# to the known value, but the measurements are far from each other, 
# then you have accuracy without precision.
#
# Concordance Correlation Coefficient (CCC): 
# The CCC can be written as the product of the accuracy and the 
# precision coefficients. It measures the agreement along the identity line,
# in which a value of 1 represents a perfect agreement, a value of -1 represents
# a perfect disagreement and a value of 0 represents no agreement.
#
# Total Deviation Index (TDI) + Coverage probability (CP):
# An intuitively clear measurement of agreement is a measure that captures
# a large proportion of data within a predetermined boundary from target values.
# For example, we may want to capture at least 90% of individual observations that are
# within 10 mmHg of their target values (blood pressure). In other words, 
# one can assume that the manual instrument can be replaced by the automatic
# device if a large proportion of paired-measurement differences are within 
# a boundary of 10 mmHg. The CP describes the proportion captured within a 
# pre-specified boundary of the absolute paired-measurement differences from two devices.

# A ggplot2 object is created
fg01 <- ggplot(aes(x = XD,y = YB),data=df.base) +
  geom_abline(intercept = 0, slope = 1, colour = "gray", size = 1) +
  geom_point(shape = 19, colour = "red", size = 4.0, alpha = 0.95) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Bourdon vs digital gauge comparison") +
  xlab("Digital gauge pressure (kgf/cm2)") +
  ylab("Bourdon gauge pressure (kgf/cm2)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg01

# A ggplot2 object is created with xy log10 scales
fg02 <- ggplot(aes(x = XD,y = YB),data=df.base) +
  geom_abline(intercept = 0, slope = 1, colour = "gray", size = 1) +
  geom_point(shape = 19, colour = "red", size = 4.0, alpha = 0.95) +
  scale_x_log10(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_y_log10(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Bourdon vs digital gauge comparison") +
  xlab("Digital gauge pressure (kgf/cm2)") +
  ylab("Bourdon gauge pressure (kgf/cm2)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg02

# A new "CLASS" character column is created in data.frame data=df.base
df.base$CLASS <- c("Deviation")

# A ggplot2 object is created
fg03 <- ggplot(aes(x = CLASS,y = SD),data=df.base) +
  geom_boxplot(size = 0.55,alpha = 0.70,outlier.colour = '#ff0033',
               outlier.size = 4.5) +
  geom_point(colour = '#0000ff',size = 3.5,
             position = position_jitter(width = 0.05)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Boxplot of simple deviation (SD)") +
  xlab("Clase") +
  ylab("simple deviation (SD)") +
  theme_bw()

# A ggplot2 object is requested
fg03

# A subset data.frame is created
df.melt <- df.base[c("YB","XD")]

# melt {reshape} function is requested to convert from wide to long format
df.melt <- melt(df.melt)

# A ggplot2 object is created
fg04 <- ggplot(aes(y = value,x = variable,colour = variable),data=df.melt) +
  geom_boxplot(size = 0.75,alpha = 0.95) +
  geom_point(size = 4.0,alpha = 0.95) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Instruments Boxplot Comparison") +
  xlab("Instrument") +
  ylab("Pressure (kgf/cm2)") +
  theme_bw(base_size = 20.0)

# A ggplot2 object is requested
fg04

# round_df function is applied to relevant data.frames
df.base.desc <- round_df(df=df.base.desc, digits=3)

# Objects to export
# Desc(None)
# t.test
# fg01, fg02, fg03, df.base.desc, shapiro.test(df.base$SD), mean-MAD, sq(mean-MSD)
# cor.test(df.base$XF,df.base$YR), summary(test1)
write.csv(df.base.desc, file = "df.base.desc.csv")

# /////////////////////////////////////////////////////////////
# END OF SCRIPT
# /////////////////////////////////////////////////////////////
