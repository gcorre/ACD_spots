## analysis of ACD staining
## july 21th, 2021

# CORRE Guillaume @ GENETHON, Evry, France

##
##
##
##
# Load libraries ---------------------------------------------------------------

library(tidyverse)
library(data.table)
library(RImageJROI);
library(sp);
library(spatstat);


# configuration -----------------------------------------------------------

channels <- data.frame(channel = as.character(1:3), staining = c("DAPI","Dystrophin","ACD"))




# load datasets -----------------------------------------------------------

files <- list.files(pattern = "nuclei.csv|spots.csv|fibre.csv") 


raw = lapply(files,fread)
names(raw) <- str_remove(files,pattern = ".csv")


## extract metadata, annotate channels, 
raw_df <- raw %>%
  bind_rows(.id = "sample") %>% 
  separate(col = "sample", into = c("sample","object"), sep = "_(?=[[:alnum:]]+$)") %>%
  filter(str_starts(string = Label, pattern = "RGB")) %>% 
  separate(Label, into = c("image","ROI","channel"), sep = ":") %>%           # extract metadata from LABEL columns
  mutate(animal = str_pad(str_match(sample,pattern = "_([0-9]+)-.+$")[,2],width = 3,side = "left",pad = "0"),
         slice = str_match(sample,pattern = "-([0-9]+)$")[,2]) %>% 
  left_join(channels, by = "channel") %>%
  select(-V1,-image,-XM,-YM,-channel) %>% pivot_wider(names_from  ="staining", values_from = c(Mean:Max,Median))


raw_list <- split(raw_df, f = raw_df$object)

#change column names and split list to dataframe
raw_list <- lapply(raw_list, function(x){
  obj <- unique(x$object)
  x %>% 
    select(sample,animal, slice, ROI,everything(),-object) %>% 
    rename_at(vars(ROI:last_col()), ~paste(.x,obj,sep = "."))
  
})

for(i in 1:length(raw_list)){
  
  assign(names(raw_list)[i], raw_list[[i]])
  
}


rm(raw_list); rm(raw); rm(raw_df); rm(files);rm(i);rm(obj);rm(channels)





# Assign spot and nuclei to their parent fiber ----------------------------

images <- unique(fibre$sample)
spots_in_fiber_ls <-list()
nuclei_in_fiber_ls <- list()
fibre_bordertouch_ls <- list()
polygon_rois <- list()

for(i in images){
  
  cat("....reading",i,"\n")
  rois <- read.ijzip(paste(i,"_fibre.roi.zip",sep = ""), names = T)
  
  ## convert to polygon
  p = lapply(rois, function(y){
    Polygon(y$coords)
  })
  
  ## convert polygone to spatial objects
  p <- lapply(seq_along(p), function(j) Polygons(list(p[[j]]), ID = names(p)[j]))
  p <- SpatialPolygons(p)
  
  polygon_rois[[i]] <- p
  
  
  ## determine if the fibers are touching the image border
  status = rep(NA,length(p))
  for(j in 1:length(p)){
    
    coords <- p@polygons[[j]]@Polygons[[1]]@coords
    if(min(coords[,1])==0 | max(coords[,1]) == 2048 | min(coords[,2])==0 | max(coords[,2])==2048){
      status[j] <- "border"
    } else {
      status[j] <- "internal"
    }
  }
  
  
  fibre_bordertouch_ls[[i]] <-  data.frame(sample = i,ROI_fibre = names(p), "status" = status)
  
  
  
  
  
  ## Assign object to parent fiber
  spots_i <- spots %>% filter(sample==i)
  spots.pp <- SpatialPoints(coords = spots_i[,c("X.spots","Y.spots")])
  
  spots_in_fiber <-  sp::over(x = spots.pp, y = p, returnList = F)
  spots_in_fiber_df <- data.frame(sample = i,
                                  ROI.spots = spots_i$ROI.spots,
                                  ROI.fibre = names(p)[spots_in_fiber]) %>% 
    filter(!is.na(ROI.fibre))
  
  spots_in_fiber_ls[[i]] <- spots_in_fiber_df
  
  
  
  nuclei_i <- nuclei %>% filter(sample == i)
  nuclei.pp <- SpatialPoints(coords = nuclei_i[,c("X.nuclei","Y.nuclei")])
  nuclei_in_fiber <-  sp::over(x = nuclei.pp, y = p, returnList = F)
  nuclei_in_fiber_df <- data.frame(sample = i,
                                  ROI.nuclei = nuclei_i$ROI.nuclei,
                                  ROI.fibre = names(p)[nuclei_in_fiber]) %>% 
    filter(!is.na(ROI.fibre))
  
  nuclei_in_fiber_ls[[i]] <- nuclei_in_fiber_df
  
}

fibre_bordertouch <- fibre_bordertouch_ls %>% bind_rows()

nuclei_in_fiber <- nuclei_in_fiber_ls %>% bind_rows()

spots_in_fiber <- spots_in_fiber_ls %>% bind_rows()


fibre <- fibre %>% left_join(fibre_bordertouch, by = c("sample","ROI.fibre"="ROI_fibre"))
nuclei <- nuclei %>% left_join(nuclei_in_fiber, by = c("sample","ROI.nuclei"))
spots <- spots %>% left_join(spots_in_fiber, by = c("sample","ROI.spots"))


rm(list =setdiff(ls(),c("fibre","nuclei","spots","polygon_rois")))


save(list = ls(),file = "ACD_20210721.rdata")



