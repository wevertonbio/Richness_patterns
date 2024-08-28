#### Select species and get unique lifeform ####
#Last run: 202 August
library(dplyr)
library(data.table)
library(tidyverse)
library(stringr)
library(pbapply)
library(future.apply) #Parallell processing
library(progressr) #Progress bar
library(tidyverse)
library(florabr)
devtools::load_all("../florabr")

# Species in AF according to florabr #
#Import data
#Import brazilian flora
my_dir <- "C:/Users/wever/Downloads/dataflorabr"
# get_florabr(output_dir = my_dir, data_version = "393.370")
bf <- load_florabr(data_dir = "C:/Users/wever/Downloads/dataflorabr", 
                   data_version = "393.370", type = "complete")

#Solve incongruencies
bf <- solve_discrepancies(bf)
#Select species in Atlantic Forest
d_af <- select_species(bf, group = c("Angiosperms","Gymnosperms"),
                       include_subspecies = FALSE,
                       include_variety = FALSE,
                       biome = "Atlantic_Forest", origin = "Native",
                       taxonomicStatus  = "Accepted")
#Create columns with source of occurrence
d_af$source_occurrence <- "florabr"

# Species in AF according to Misiones checklist #
spp_arg_par <- fread("Data/Misiones_checklist_native_endemic_species.csv", encoding = "Latin-1")
spp_arg_par$species <- get_binomial(paste(spp_arg_par$GENERO,spp_arg_par$`EPITETO ESPECIFICO`, " "))

spp <- unique(spp_arg_par$species)

#Check names with florabr
spp_check <- check_names(data = bf, species = spp)
table(spp_check$taxonomicStatus)
#Get only accepted names and species not found in Brazilian Flora
spp_accepted <- na.omit(spp_check$acceptedName) %>% unique()
spp_not_found <- spp_check$input_name[which(is.na(spp_check$acceptedName))] %>% 
  unique()
spp2 <- c(spp_accepted, spp_not_found) %>% unique()

#Check unique species in Argentinian-Paraguayan flora
spp_only <- setdiff(spp2, d_af$species)

#Subset misiones species from florabr
m_bf <- spp_check %>% filter(input_name %in% spp_only,
                             taxonomicStatus == "Accepted", Spelling == "Correct")
m_bf <- subset_species(data = bf, species = m_bf$acceptedName) %>% 
  filter(taxonomicStatus == "Accepted")
table(m_bf$origin)
#Create columns with source of occurrence
m_bf$source_occurrence <- "misiones_checklist"
#Merge data
data <- bind_rows(d_af, m_bf)
#Check species duplicates
which(duplicated(data$species))

#Create dataframe with other species
d_other <- setdiff(spp_only, m_bf$species)
#Subset from misiones
d_other <- spp_arg_par %>% filter(species %in% d_other) %>% 
  select(species, lifeForm = HABITO, family = FAMILIA, origin = STATUS)
#Create habitat
d_other$habitat <- NA
d_other$habitat[grepl("terrestre", d_other$lifeForm)] <- "Terrestrial"
d_other$habitat[grepl("epífita", d_other$lifeForm)] <- "Epiphytic"
#Fix lifeForm
d_other$lifeForm <- gsub(" terrestre o rupícola| suculento| suculenta| terrestre", ";", d_other$lifeForm)
d_other$lifeForm <- gsub(" o ", ";", d_other$lifeForm)
#Translate
d_other$lifeForm %>% unique()
d_other$lifeForm <- gsub("Árbol", "Tree", d_other$lifeForm)
d_other$lifeForm <- gsub("árbol", "Tree", d_other$lifeForm)
d_other$lifeForm <- gsub("Hierba", "Herb", d_other$lifeForm)
d_other$lifeForm <- gsub("Subarbusto", "Subshrub", d_other$lifeForm)
d_other$lifeForm <- gsub("subarbusto", "Subshrub", d_other$lifeForm)
d_other$lifeForm <- gsub("enredadera", "Liana/scandent/vine", d_other$lifeForm)
d_other$lifeForm <- gsub("Enredadera", "Liana/scandent/vine", d_other$lifeForm)
d_other$lifeForm <- gsub("Arbusto", "Shrub", d_other$lifeForm)
d_other$lifeForm <- gsub("arbusto", "Shrub", d_other$lifeForm)
d_other$lifeForm <- gsub("arbolito", "Tree", d_other$lifeForm)
d_other$lifeForm <- gsub("Arbolito", "Tree", d_other$lifeForm)
d_other$lifeForm[d_other$lifeForm == "Herb;"] <- "Herb"
d_other$lifeForm[d_other$lifeForm == "Herb epífita"] <- "Herb"
d_other$lifeForm[d_other$lifeForm == "Subshrub;"] <- "Subshrub"
d_other$lifeForm[d_other$lifeForm == "Tree;Tree"] <- "Tree"
d_other$lifeForm[d_other$lifeForm == "Liana" ] <- "Liana/scandent/vine"

#Get endemism
d_other$endemism <- "Non-endemic"
d_other$endemism[d_other$origin == "Endémica"] <- "Endemic"
d_other$origin <- "Native"
#Create columns with source of occurrence
d_other$source_occurrence <- "misiones_checklist"

#Merge all data
data2 <- data %>% select(colnames(d_other)) %>% 
  rbind(d_other) %>% distinct()
which(duplicated(data2$species))

#Remove naturalized and cultived species
data2$origin %>% table()
data2 <- data2 %>% filter(!(origin %in% c("Cultivated", "Naturalized")))

#### Get unique lifeforms ####
d_af <- data2
#Convert unknown and dracaenoid lifeforms to NA
d_af$lifeForm[which(d_af$lifeForm %in% c("Unknown", "Dracaenoid", "", 
                                         "Not_found_in_brazil"))] <- NA
#Creta columns with genus
d_af$genus <- stringr::word(d_af$species, 1)


####Get species without lifeform####
lf.na <- d_af %>% filter(is.na(lifeForm) | lifeForm == "")
#Get species with multiple lifeforms
lf.more <- lengths(strsplit(d_af$lifeForm, ';'))
lf.more <- cbind(d_af, "lf.more" = lf.more) %>% filter(lf.more > 1) %>% select(-lf.more)
#Merge lifeforms to fix (get only one lifeform by species)
lf.prob <- rbind(lf.na, lf.more)
#Get lifeforms ok to check
lf.ok <- d_af %>% filter(!(species %in% lf.prob$species))
table(lf.ok$lifeForm)

#Get consensus on lifeform based on the paper:
#"A plant growth form dataset for the New World"
#https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.1569#support-information-section
lf.data <- fread("Data/GrowthForm_Final.txt")
colnames(lf.data)
lf.data <- lf.data %>% dplyr::select(species = SPECIES_STD,
                                     Lf_consensus = GROWTHFORM_STD)

#Merge data with problematic lifeforms
lf.j <- left_join(lf.prob, lf.data, by = "species")
lf.j$lifeForm[which(lf.j$lifeForm == "")] <- NA
#Fix names to match brazilian flora
table(lf.j$Lf_consensus)
lf.j$Lf_consensus[which(lf.j$Lf_consensus == "Liana")] <- "Liana/scandent/vine"
lf.j$Lf_consensus[which(lf.j$Lf_consensus == "Non-woody epiphyte")] <- "Herb"
lf.j$Lf_consensus[which(lf.j$Lf_consensus == "Vine")] <- "Liana/scandent/vine"

#Create new column to fill with new lifeform
lf.j$new.lifeForm <- NA
#Fill NAs
lf.j$new.lifeForm[which(is.na(lf.j$lifeForm))] <- lf.j$Lf_consensus[which(is.na(lf.j$lifeForm))]
#Get consensus
#Tree
lf.j$new.lifeForm[which(grepl("Tree", lf.j$lifeForm) &
                          lf.j$Lf_consensus == "Tree")] <- "Tree"
#Shrub
lf.j$new.lifeForm[which(grepl("Shrub", lf.j$lifeForm) &
                          lf.j$Lf_consensus == "Shrub")] <- "Shrub"
#Herb
lf.j$new.lifeForm[which(grepl("Herb", lf.j$lifeForm) &
                          lf.j$Lf_consensus == "Herb")] <- "Herb"
#Liana/scandent/vine
lf.j$new.lifeForm[which(grepl("Liana/scandent/vine", lf.j$lifeForm) &
                          lf.j$Lf_consensus == "Liana/scandent/vine")] <- "Liana/scandent/vine"
#For species with lifeform OK, merge with lf.ok
lf.ok2 <- lf.j %>% filter(!is.na(new.lifeForm)) %>% 
  dplyr::select(species, lifeForm = new.lifeForm, genus, family) %>% bind_rows(lf.ok)
table(lf.ok2$lifeForm)
lf.prob2 <- lf.prob %>% filter(!species %in% lf.ok2$species)

#For species with more than one lifeform, get most common lifeform in the genus or family
#Get all native species in Brazil
d_br <- select_species(bf, group = c("Angiosperms","Gymnosperms"), include_subspecies = FALSE,
                       include_variety = FALSE,
                       origin = "Native",
                       taxonomicStatus  = "Accepted")
d_br <- d_br %>% dplyr::select(species, genus, family, lifeForm)
#Expand lifeforms
d_br_exp <- d_br %>% 
  tidyr::separate_rows(lifeForm, sep = ";")
#Get only genus with 3 or more species
genus3 <- table(d_br$genus) %>% as.data.frame() %>% filter(Freq >= 3) %>%
  pull(Var1) %>% as.character()
#Count lifeform by genus
genus_lf <- d_br_exp %>%
  filter(genus %in% genus3) %>% 
  group_by(genus, lifeForm) %>%
  summarise(count = n()) %>%
  ungroup() %>% 
  filter(lifeForm != "Unknown",
         lifeForm != "",
         lifeForm != "Succulent")
unique(genus_lf$lifeForm)
head(genus_lf)
#Most common lifeform by genus
genus_final <- genus_lf %>% group_by(genus) %>%
  slice(which.max(count)) %>%
  ungroup()
table(genus_final$lifeForm)
#Second most common lifeform by genus
genus_final2 <- genus_lf %>%
  arrange(genus, desc(count)) %>% 
  group_by(genus) %>%
  slice(2) %>%
  ungroup() %>% 
  dplyr::rename(lifeForm2 = lifeForm) %>% dplyr::select(-count)
genus_final <- left_join(genus_final, genus_final2) %>%
  dplyr::select(-count) %>% 
  dplyr::rename(lifeForm1 = lifeForm)
#Join data
fix_genus <- left_join(lf.prob2, genus_final, by = "genus")
#Get genus fixed 1
genus_fixed <- fix_genus %>% filter(is.na(lifeForm) | lifeForm == "") %>% 
  filter(!is.na(lifeForm1)) %>% 
  dplyr::mutate(lifeForm = lifeForm1) %>% 
  dplyr::select(species, lifeForm, genus, family)
#Get genus to fix
fix_genus2 <- fix_genus %>% filter(!(species %in% genus_fixed$species),
                                   !is.na(lifeForm1))
#For genus with more than one lifeform, get the first or second option
first_lf_genus <- pbsapply(1:nrow(fix_genus2), function(i){
  fg_i <- fix_genus2[i,]
  lf_i <- strsplit(fg_i$lifeForm, ";") %>% unlist()
  new_i1 <- fg_i$lifeForm1
  new_i2 <- fg_i$lifeForm2
  if(new_i1 %in% lf_i) {
    lf <- new_i1
  }
  if(!(new_i1 %in% lf_i) & new_i2 %in% lf_i) {
    lf <- new_i2
  }
  if(!(new_i1 %in% lf_i) & !(new_i2 %in% lf_i)) {
    lf <- NA
  }
  return(lf)
})
fix_genus2$new_lf <- first_lf_genus
genus_fixed2 <- fix_genus2 %>% filter(!is.na(new_lf)) %>% 
  dplyr::mutate(lifeForm = new_lf) %>% 
  dplyr::select(species, lifeForm, genus, family)
table(genus_fixed2$lifeForm)
#Join data
lf.ok3 <- bind_rows(genus_fixed, genus_fixed2, lf.ok2)
#Get genus with problems
lf.prob3 <- lf.prob %>% filter(!species %in% lf.ok3$species)

#Do the same thing with family now
#Get only family with 3 or more species
family3 <- table(d_br$family) %>% as.data.frame() %>% filter(Freq >= 3) %>%
  pull(Var1) %>% as.character()
#Count lifeform by family
family_lf <- d_br_exp %>%
  filter(family %in% family3) %>% 
  group_by(family, lifeForm) %>%
  summarise(count = n()) %>%
  ungroup() %>% 
  filter(lifeForm != "Unknown",
         lifeForm != "",
         lifeForm != "Succulent")
unique(family_lf$lifeForm)
head(family_lf)
#Most common lifeform by family
family_final <- family_lf %>% group_by(family) %>%
  slice(which.max(count)) %>%
  ungroup()
table(family_final$lifeForm)
#Second most common lifeform by family
family_final2 <- family_lf %>%
  arrange(family, desc(count)) %>% 
  group_by(family) %>%
  slice(2) %>%
  ungroup() %>% 
  dplyr::rename(lifeForm2 = lifeForm) %>% dplyr::select(-count)
family_final <- left_join(family_final, family_final2) %>%
  dplyr::select(-count) %>% 
  dplyr::rename(lifeForm1 = lifeForm)
#Join data
fix_family <- left_join(lf.prob3, family_final, by = "family")
#Get family fixed 1
family_fixed <- fix_family %>% filter(is.na(lifeForm) | lifeForm == "") %>% 
  filter(!is.na(lifeForm1)) %>% 
  dplyr::mutate(lifeForm = lifeForm1) %>% 
  dplyr::select(species, lifeForm, family, family)
#Get family to fix
fix_family2 <- fix_family %>% filter(!(species %in% family_fixed$species),
                                     !is.na(lifeForm1))
#For family with more than one lifeform, get the first or second option
first_lf_family <- pbsapply(1:nrow(fix_family2), function(i){
  fg_i <- fix_family2[i,]
  lf_i <- strsplit(fg_i$lifeForm, ";") %>% unlist()
  new_i1 <- fg_i$lifeForm1
  new_i2 <- fg_i$lifeForm2
  if(new_i1 %in% lf_i) {
    lf <- new_i1
  }
  if(!(new_i1 %in% lf_i) & new_i2 %in% lf_i) {
    lf <- new_i2
  }
  if(!(new_i1 %in% lf_i) & !(new_i2 %in% lf_i)) {
    lf <- NA
  }
  return(lf)
})
fix_family2$new_lf <- first_lf_family
family_fixed2 <- fix_family2 %>% filter(!is.na(new_lf)) %>% 
  dplyr::mutate(lifeForm = new_lf) %>% 
  dplyr::select(species, lifeForm, genus, family)
table(family_fixed2$lifeForm)
#Join data
lf.ok4 <- bind_rows(family_fixed, family_fixed2, lf.ok3)
#Get genus with problems
lf.prob4 <- lf.prob %>% filter(!species %in% lf.ok4$species) %>% 
  left_join(genus_final)
#For species with consensus, use consensus by genus 1
lf.ok5 <- lf.prob4 %>% filter(!is.na(lifeForm1)) %>% 
  mutate(lifeForm = lifeForm1) %>%
  dplyr::select(species, lifeForm, genus, family) %>% 
  bind_rows(lf.ok4)
table(lf.ok5$lifeForm)
#Get family with problems
lf.prob5 <- lf.prob %>% filter(!species %in% lf.ok5$species) %>% 
  left_join(family_final)
#For species with consensus, use consensus by genus 1
lf.ok6 <- lf.prob5 %>% filter(!is.na(lifeForm1)) %>% 
  mutate(lifeForm = lifeForm1) %>%
  dplyr::select(species, lifeForm, genus, family) %>% 
  bind_rows(lf.ok5)
table(lf.ok5$lifeForm)

#Last fix
lf.prob6 <- lf.prob %>% filter(!species %in% lf.ok6$species) %>% 
  left_join(family_final)
#For these species, use first option
lf7 <- strsplit(lf.prob6$lifeForm, ";")
lf7 <- sapply(lf7, function(i) i[1])
lf.ok7 <- lf.prob6 %>% mutate(lifeForm = lf7) %>% 
  dplyr::select(species, lifeForm, genus, family) %>% 
  bind_rows(lf.ok6)
table(lf.ok7$lifeForm)
table(lf.ok7$lifeForm) %>% sum()

#Fix columns data
spp_ok_lf <- lf.ok7 %>% select(species, lifeForm)
#Merge
data_final <- d_af %>% select(-lifeForm) %>% left_join(spp_ok_lf)
table(data_final$lifeForm)

#Fix lifeforms
data_final$lifeForm[data_final$lifeForm == "Liana/scandent/vine"] <- "Liana"

####Split herbs in terrestrial and Epiphytic (and remove aquatic)####
data_final <- data_final %>% filter(habitat != "Aquatic")
#Herbs
herbs <- data_final %>% filter(lifeForm == "Herb")
table(herbs$habitat)
herbs$habitat[grepl("Epiphytic|Hemiepiphyte|Hemiparasite", herbs$habitat)] <- "Epiphytic"
herbs$habitat[grepl("Terrestrial", herbs$habitat)] <- "Terrestrial"
herbs$habitat[!grepl("Terrestrial|Epiphytic", herbs$habitat)] <- "Terrestrial"

herbs$lifeForm[herbs$lifeForm == "Herb" & herbs$habitat == "Terrestrial"] <- "Terrestrial_herb"
herbs$lifeForm[herbs$lifeForm == "Herb" & herbs$habitat == "Epiphytic"] <- "Epiphytic_herb"

table(herbs$lifeForm)

data_final2 <- data_final %>% filter(lifeForm != "Herb") %>% rbind(herbs)
table(data_final2$lifeForm)

#Save final results
fwrite(data_final2, "Data/SpeciesData.gz", compress = "gzip", row.names = FALSE)
