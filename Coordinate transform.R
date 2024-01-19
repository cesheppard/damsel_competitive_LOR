library(tidyverse)

###SET WORKING DIRECTORY###
setwd("")

# Transect 1
# SVS 3D point data

T1_SVS <- read.csv("T1_SVS.csv", header = T, stringsAsFactors = F)
T1_fish <- T1_SVS %>% filter(Genus != "Reference")
T1_refs <- T1_SVS %>% filter(Genus == "Reference")

# Reference points
T1_all_refs <- read.csv("T1_all_references.csv", header = T, stringsAsFactors = F)

# Fetch reference numbers for 3D point data
for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$ref1[i] <- T1_refs$Number[T1_refs$OpCode == T1_fish$OpCode[i]
                               & T1_refs$Time == T1_fish$Time[i]][1]
}

for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$ref2[i] <- T1_refs$Number[T1_refs$OpCode == T1_fish$OpCode[i]
                              & T1_refs$Time == T1_fish$Time[i]][2]
}

# Fetch z and x coordinates for 3D reference points (note z will be x, and x will be y on LOR)
for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$x1[i] <- T1_refs$Z[T1_refs$OpCode == T1_fish$OpCode[i]
                            & T1_refs$Number == T1_fish$ref1[i]
                             & T1_refs$Time == T1_fish$Time[i]]
}

for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$y1[i] <- T1_refs$X[T1_refs$OpCode == T1_fish$OpCode[i]
                       & T1_refs$Number == T1_fish$ref1[i]
                       & T1_refs$Time == T1_fish$Time[i]]
}

for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$x2[i] <- T1_refs$Z[T1_refs$OpCode == T1_fish$OpCode[i]
                       & T1_refs$Number == T1_fish$ref2[i]
                       & T1_refs$Time == T1_fish$Time[i]]
}

for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$y2[i] <- T1_refs$X[T1_refs$OpCode == T1_fish$OpCode[i]
                       & T1_refs$Number == T1_fish$ref2[i]
                       & T1_refs$Time == T1_fish$Time[i]]
}

# Fetch x and y LOR reference point coordinates
for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$x1b[i] <- T1_all_refs$x[T1_all_refs$OpCode == T1_fish$OpCode[i]
                            & T1_all_refs$ID == T1_fish$ref1[i]]
}

for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$y1b[i] <- T1_all_refs$y[T1_all_refs$OpCode == T1_fish$OpCode[i]
                            & T1_all_refs$ID == T1_fish$ref1[i]]
}

for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$x2b[i] <- T1_all_refs$x[T1_all_refs$OpCode == T1_fish$OpCode[i]
                            & T1_all_refs$ID == T1_fish$ref2[i]]
}

for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$y2b[i] <- T1_all_refs$y[T1_all_refs$OpCode == T1_fish$OpCode[i]
                            & T1_all_refs$ID == T1_fish$ref2[i]]
}

# Transformation equation
T1_fish <- T1_fish %>% mutate(p = ((Z-x1)*(x2-x1))+((X-y1)*(y2-y1)), 
                        q = -((Z-x1)*(y2-y1))+((X-y1)*(x2-x1)))

T1_fish <- T1_fish %>% mutate(
  x3 = x1b + ((1/((x2-x1)^2+(y2-y1)^2))*((p*(x2b-x1b))-(q*(y2b-y1b))))
)

T1_fish <- T1_fish %>% mutate(
  y3 = y1b + ((1/((x2-x1)^2+(y2-y1)^2))*((p*(y2b-y1b))+(q*(x2b-x1b))))
)

# TRANSECT 2
# SVS 3D point data

T2_SVS <- read.csv("T2_SVS.csv", header = T, stringsAsFactors = F)
T2_fish <- T2_SVS %>% filter(Genus != "Reference")
T2_refs <- T2_SVS %>% filter(Genus == "Reference")

# Reference points
T2_all_refs <- read.csv("T2_all_references.csv", header = T, stringsAsFactors = F)

# Fetch reference numbers for 3D point data
for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$ref1[i] <- T2_refs$Number[T2_refs$Time == T2_fish$Time[i]][1]
}

for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$ref2[i] <- T2_refs$Number[T2_refs$Time == T2_fish$Time[i]][2]
}

# Fetch z and x coordinates for 3D reference points (note z will be x, and x will be y on LOF)
for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$x1[i] <- T2_refs$Z[T2_refs$Number == T2_fish$ref1[i]
                             & T2_refs$Time == T2_fish$Time[i]]
}

for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$y1[i] <- T2_refs$X[T2_refs$Number == T2_fish$ref1[i]
                             & T2_refs$Time == T2_fish$Time[i]]
}

for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$x2[i] <- T2_refs$Z[T2_refs$Number == T2_fish$ref2[i]
                             & T2_refs$Time == T2_fish$Time[i]]
}

for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$y2[i] <- T2_refs$X[T2_refs$Number == T2_fish$ref2[i]
                             & T2_refs$Time == T2_fish$Time[i]]
}

# Fetch x and y reference point coordinates
for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$x1b[i] <- T2_all_refs$x[T2_all_refs$id == T2_fish$ref1[i]]
}

for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$y1b[i] <- T2_all_refs$y[T2_all_refs$id == T2_fish$ref1[i]]
}

for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$x2b[i] <- T2_all_refs$x[T2_all_refs$id == T2_fish$ref2[i]]
}

for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$y2b[i] <- T2_all_refs$y[T2_all_refs$id == T2_fish$ref2[i]]
}

# Transformation equation
T2_fish <- T2_fish %>% mutate(p = ((Z-x1)*(x2-x1))+((X-y1)*(y2-y1)), 
                              q = -((Z-x1)*(y2-y1))+((X-y1)*(x2-x1)))

T2_fish <- T2_fish %>% mutate(
  x3 = x1b + (((p*(x2b-x1b))-(q*(y2b-y1b)))/((x2-x1)^2+(y2-y1)^2))
)

T2_fish <- T2_fish %>% mutate(
  y3 = y1b + (((p*(y2b-y1b))+(q*(x2b-x1b)))/((x2-x1)^2+(y2-y1)^2))
)

# dietary data
dietary <- read.csv("dietary.csv", header = T, stringsAsFactors = F)

for (i in 1:length(T1_fish$OpCode)) {
  T1_fish$dietary[i] <- dietary$Diet[dietary$Genus == T1_fish$Genus[i]
                                        & dietary$Species == T1_fish$Species[i]]
}

for (i in 1:length(T2_fish$OpCode)) {
  T2_fish$dietary[i] <- dietary$Diet[dietary$Genus == T2_fish$Genus[i]
                                        & dietary$Species == T2_fish$Species[i]]
}

write.csv(T1_fish, "T1_fish_transformed.csv")
write.csv(T2_fish, "T2_fish_transformed.csv")
