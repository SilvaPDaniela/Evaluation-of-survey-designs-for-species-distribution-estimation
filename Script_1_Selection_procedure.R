#**********************
# Paper Adapting the sampling design of research surveys to improve the biomass estimation of non-target species -  the case study of Raja clavata #
# Evaluation of alternative survey designs: Selection procedure                                                                                    #
#*****************************
# Daniela Silva        #
#*****************************
rm(list=ls()) # Clean workspace
library(DT)
library(fields)


# Read dataset -----------------------------------------------------------------

# First dataset corresponds to the dataset resulting from the estimation of abundance 
# indicators of the species.
# It should contain the coordinates (eg. longitude and latitude) of a thinned grid,
# and all variables that can be important for the selection procedure of stations 
# as the predicted abundance for all species and the respective error measure 
# resulting from the prediction process.

# In our case, the variables are: longitude, latitude, biomass of Raja clavata and 
# the corresponding standard error of spatial effects, abundance of Merluccius 
# merluccius and the corresponding standard error of spatial effects, stratum and 
# sector (See section 2.3 of the paper).

estimated_grid<-read.csv2("./Github/Data/estimated_grid.csv")

# The survey design to be proposed is defined by considering the standard survey 
# protocol, which includes a total of 65 stations, from which 54 were previously 
# identified as fixed stations and 11 stations selected.

# Read the predicted data for fixed locations. This dataset should contain the 
# coordinates, predicted abundance indicators and variables that can help to define
# the constraint. This dataset is unnecessary if there is no constraints regarding
# to distribution of the stations in the study region.

# In our case, the variables are: longitude, latitude, biomass of Raja clavata and 
# abundance of Merluccius merluccius and stratum.
fixed_grid<-read.csv2("./Github/Data/fixed_grid.csv")


# Selection procedure ---------------------------------------------------------

# 11 stations are selected according to the sampling procedure considering the 
# survey condition (at least two stations by stratum).

# Define the number of sampling units (locations) to consider in the survey design
n_s=65

# Define the minimum number of sampling units by stratum
n_r<-2

# Define the constraint: all strata should contain at least n_r by stratum
constraint<-data.frame(stratum=unique(estimated_grid$stratum),min=n_r)

# Compute the number of fixed stations by strata in order to identify which strata
# have missing stations
units.by.stratum<-data.frame(t(table(fixed_grid$stratum)))[,-1]
names(units.by.stratum)[1]<-"stratum"
units.by.stratum<-merge(units.by.stratum,constraint,all=T)
units.by.stratum$Freq[which(is.na(units.by.stratum$Freq))]<-rep(0,sum(is.na(units.by.stratum$Freq)))
units.by.stratum$missing<-ifelse(units.by.stratum$min-units.by.stratum$Freq<0,0,units.by.stratum$min-units.by.stratum$Freq)

# Identify the strata with missing stations
missing.units<-units.by.stratum[which(units.by.stratum$missing>0),]
estimated_stratum<-unique(estimated_grid$stratum)
missing_stratum<-unique(missing.units$stratum)

# Calculate the number of missing locations after considering the constraint
n_m<-n_s-nrow(fixed_grid)-sum(missing.units$missing)

# Compute the minimum distance between two locations in order to avoid to close locations.
# In our case, we defined it as the minimum distance between two fixed locations.
dist.matrix.fixed.grid<-rdist(fixed_grid[,c("long","lat")],fixed_grid[,c("long","lat")])
diag(dist.matrix.fixed.grid)<-rep(100000000,nrow(dist.matrix.fixed.grid))
dist.vector<-as.vector(dist.matrix.fixed.grid)
min.dist<-min(dist.matrix.fixed.grid)


# Compute the maximum of predicted biomass/abundance and predicted standard deviations of
# spatial effects
max.grid.specie1<-max(estimated_grid$abund_specie1)
max.grid.specie2<-max(estimated_grid$abund_specie2)
max.grid.sd.specie1<-max(estimated_grid$sd_specie1)
max.grid.sd.specie2<-max(estimated_grid$sd_specie2)

# Compute the sector weights (identified as q_s in Section 2.3.1. of the paper)
q_by_stratum<-tapply(estimated_grid$abund_specie1,estimated_grid$sector,mean)
q_by_stratum<-data.frame(sector=names(q_by_stratum),q=q_by_stratum)
q_by_stratum$q<-10-rank(q_by_stratum$q)
estimated_grid<-merge(estimated_grid,q_by_stratum,by="sector",all.x=T)

# Define the measures to select locations for new station (presented in Table 2 of the
# paper). The functions are in the order presented in the table.
points.selection.measures<-list(function(x){(x[,"sd_specie2"]/max.grid.sd.specie2)*(x[,"sd_specie1"]/max.grid.sd.specie1)},
                                function(x){((x[,"sd_specie1"]/max.grid.sd.specie1)*(1-(x[,"abund_specie1"]/max.grid.specie1)))*((x[,"sd_specie2"]/max.grid.sd.specie2)*(1-(x[,"abund_specie2"]/max.grid.specie2)))},
                                function(x){((x[,"sd_specie1"]/max.grid.sd.specie1)*(1-(x["abund_specie1"]/max.grid.specie1)))*(x[,"sd_specie2"]/max.grid.sd.specie2)},
                                function(x){(x[,"sd_specie1"]/max.grid.sd.specie1)*(x[,"sd_specie2"]/max.grid.sd.specie2)*x[,"q"]},
                                function(x){(x[,"sd_specie1"])},
                                function(x){(x[,"sd_specie1"])*(1-(x["abund_specie1"]/max.grid.specie1))},
                                function(x){(x[,"sd_specie1"]*x[,"q"])},
                                function(x){(sqrt((sqrt(x["abund_specie1"]/max.grid.specie1))-sqrt(x[,"abund_specie2"]/max.grid.specie2))^2)})

# Create a clean list that it will be used to save the designs (set of locations)
# corresponding to each measure above detailed.
samples=vector(mode = "list", length = length(points.selection.measures))


# The points.selection.function allows to selected the locations according to the 
# measures above mentioned.
# As inputs, it is necessary a thinned estimated grid, the fixed (if the case) and
# the clean list to save the results of this function.
points.selection.function<-function(estimated_grid,fixed_grid,samples){
  for(s in 1:length(samples)){
    samp.data1<-data.frame()
    # The first step is given by the selection of locations in order to guarantee
    # the standard survey protocol, which, in our case, it is to consider a 
    # minimum number of sampling units by stratum
    points.setor<-fixed_grid[,c("long","lat","stratum")]
    for (i in 1:length(missing.units$stratum)){
      pos<-estimated_grid[which(estimated_grid$stratum==missing_stratum[i]),]
      for(j in 1:(missing.units[i,"missing"])){
        dist.df<-rdist(pos[,c("long","lat")],points.setor[,c("long","lat")])
        mat.dist<-dist.df>=min.dist
        pos.dist<-which(apply(mat.dist,1,sum)==nrow(points.setor))
        if(length(pos.dist)==0){
          available.points.setor<-pos
        }
        if(length(pos.dist)>0){
          available.points.setor<-pos[pos.dist,]
        }
        min<-order(points.selection.measures[[s]](available.points.setor))[1]
        choose.points<-data.frame(long=available.points.setor[min,"long"],
                                  lat=available.points.setor[min,"lat"],
                                  stratum=available.points.setor[min,"stratum"],
                                  specie1_pred=available.points.setor[min,"abund_specie1"],
                                  specie2_pred=available.points.setor[min,"abund_specie2"])
        samp.data1<-rbind(samp.data1,choose.points)
        points.setor<-rbind(points.setor,choose.points[,1:3])
      }
    }
    
    samp.data1.1<-data.frame(long=c(fixed_grid$long,samp.data1$long),
                             lat=c(fixed_grid$lat,samp.data1$lat),
                             stratum=c(fixed_grid$stratum,as.character(samp.data1$stratum)),
                             abund_specie1=c(fixed_grid$abund_specie1,samp.data1$specie1_pred),
                             abund_specie2=c(fixed_grid$abund_specie2,samp.data1$specie2_pred))
    
    # Select the remaining points considering the entire study region
    estimated_grid2<-estimated_grid
    for(k in 1:n_m){
      for (l in 1:nrow(samp.data1.1)){
        dist.df<-rdist(estimated_grid2[,c("long","lat")],samp.data1.1[l,c("long","lat")])
        estimated_grid2<-estimated_grid2[which(dist.df>=min.dist),]
      }
      min<-order(points.selection.measures[[s]](estimated_grid2))[1]
      choose.points<-data.frame(long=estimated_grid2[min,"long"],
                                lat=estimated_grid2[min,"lat"],
                                stratum=estimated_grid2[min,"stratum"],
                                abund_specie1=estimated_grid2[min,"abund_specie1"],
                                abund_specie2=estimated_grid2[min,"abund_specie2"])
      samp.data1.1<-rbind(samp.data1.1,choose.points)
    }
    samples[[s]]<-cbind(samp.data1.1,sample=s)
  }
  do.call("rbind", samples)
}

sample.surveys<-points.selection.function(estimated_grid = estimated_grid,fixed_grid = fixed_grid,samples=samples)

# Save the results
write.csv(sample.surveys,"./Github/samples4survey.csv")
