#!/opt/common/CentOS_7/R/R-4.2.0/bin/Rscript 
suppressPackageStartupMessages({
    library(varhandle)
    library(argparse)
    library(dplyr)
    library(ggplot2)
    library(tibble)
    library(tidyr)
    library(stringr)
    library(stringi)
    library(ggpubr)
    library(fmsb)
    
})

args = commandArgs(TRUE)
print(args)


print(args[2])
parser = ArgumentParser(description = 'Intake grouped tax_out files and generates alphadiv plots and oddsratio plots for comparisons derived from specified csv file.')

parser$add_argument('-v', '--verbose', action="store_true", default = TRUE,
                    help = 'Print run info')
parser$add_argument('-f', '--awk-file', required = TRUE,
                    help = 'text file output of awk output from microbiome pipeline')
parser$add_argument('-m', '--manifest', required = FALSE,
                    help = '.csv file with first column with case id matching microbiome files with classes labeled in columns')                    
parser$add_argument('-o', '--output', required = TRUE,
                    help = 'Enter the General Tumor Type')
parser$add_argument('-g', '--generaltumortype', required = TRUE,
                    help = 'Enter the General Tumor Type')
parser$add_argument('-d', '--detailedtumortype', required = TRUE,
                    help = 'Enter the detailed Tumor Type')
parser$add_argument('-D', '--directory', required = TRUE,
                    help = 'Output directory to which all output files are written to')
parser$add_argument('-A', '--div_fig_dir', required = FALSE,
                    help = 'Output directory to which alpha diveristy table and plot output files are written to')
parser$add_argument('-O', '--odds_fig_dir', required = FALSE,
                    help = 'Output directory to which odds ratio plots and table output files are written to')


args = parser$parse_args()

print(args)

print(args$copynumber_file)

finaljoin <- read.delim("/lila/data/vanderbilt/microbiome_files/finaljoin.txt", sep="\t", header = TRUE)
better_taxonomy <- read.delim("/lila/data/vanderbilt/microbiome_files/better_taxonomy.txt", sep="\t", header = TRUE)


input <- read.table(args$awk_file, quote="#", stringsAsFactors=FALSE)


casetax <-  input %>%
  mutate(escore = 1.1) %>%
  select(V3, V1, V2 , escore) %>% 
  rename(readcount = V1, tax_id = V2, case = V3)




casetax <- transform(casetax, tax_id = as.numeric(tax_id))
casetax <- transform(casetax, readcount = as.numeric(readcount))



casetax$tax_id[casetax$tax_id == 330] <- 301
casetax$tax_id[casetax$tax_id == 9478] <- 1868482


trialjoin <- casetax %>% 
  left_join(finaljoin, c("tax_id" = "tax_id")) %>% 
  filter(tax_name != "NA")

#Filtering column

d <- data.frame()


for(i in 1:nrow(trialjoin)){ifelse(trialjoin$rank0[i] %in% c("species", "species_group"),  d <- rbind(d, trialjoin$tax_id[i]),  
                                   ifelse(trialjoin$rank1[i] %in% c("species", "species_group"), d <- rbind(d, trialjoin$parent_id1[i]), 
                                          ifelse(trialjoin$rank2[i] %in% c("species", "species_group"), d <- rbind(d, trialjoin$parent_id2[i]), 
                                                 ifelse(trialjoin$rank3[i] %in% c("species", "species_group"), d <- rbind(d, trialjoin$parent_id3[i]),
                                                        ifelse(trialjoin$rank4[i] %in% c("species", "species_group"), d <- rbind(d, trialjoin$parent_id4[i]),
                                                               ifelse(trialjoin$rank5[i] %in% c("species", "species_group"), d <- rbind(d, trialjoin$parent_id5[i]),
                                                                      ifelse(trialjoin$rank6[i] %in% c("species", "species_group"), d <- rbind(d, trialjoin$parent_id6[i]),
                                                                             ifelse(trialjoin$rank7[i] %in% c("species", "species_group"), d <- rbind(d, trialjoin$parent_id7[i]),
                                                                                    ifelse(trialjoin$rank8[i] %in% c("species", "species_group"), d <- rbind(d, trialjoin$parent_id8[i]), d <- rbind(d, 0))))))))))}

g <- data.frame()

for(i in 1:nrow(trialjoin)){ifelse(trialjoin$rank0[i] %in% "genus",  g <- rbind(g, trialjoin$tax_id[i]),  
                                   ifelse(trialjoin$rank1[i] %in% "genus", g <- rbind(g, trialjoin$parent_id1[i]), 
                                          ifelse(trialjoin$rank2[i] %in% "genus", g <- rbind(g, trialjoin$parent_id2[i]), 
                                                 ifelse(trialjoin$rank3[i] %in% "genus", g <- rbind(g, trialjoin$parent_id3[i]),
                                                        ifelse(trialjoin$rank4[i] %in% "genus", g <- rbind(g, trialjoin$parent_id4[i]),
                                                               ifelse(trialjoin$rank5[i] %in% "genus", g <- rbind(g, trialjoin$parent_id5[i]),
                                                                      ifelse(trialjoin$rank6[i] %in% "genus", g <- rbind(g, trialjoin$parent_id6[i]),
                                                                             ifelse(trialjoin$rank7[i] %in% "genus", g <- rbind(g, trialjoin$parent_id7[i]),
                                                                                    ifelse(trialjoin$rank8[i] %in% "genus", g <- rbind(g, trialjoin$parent_id8[i]), g <- rbind(g, 0))))))))))}


newmicroall <- data.frame(d, g, trialjoin) %>% 
  select(case, readcount,X0, X0.1,  tax_id, escore) %>% 
  rename(species = X0, genus =X0.1) %>% 
  left_join(better_taxonomy, c("species" = "tax_id")) %>%
  select(case, readcount, tax_name, genus,  tax_id, escore) %>% 
  rename(species_name= tax_name) %>% 
  left_join(better_taxonomy, c("genus" = "tax_id")) %>%
  select(case, readcount, species_name, tax_name,  tax_id, escore) %>% 
  rename(genus_name= tax_name) %>% 
  mutate(GeneralTumorType = args$generaltumortype) %>% 
  mutate(DetailedTumorType = args$detailedtumortype) %>% 
  mutate(MSI_Class = "Do not Report", SampleType = "null") %>% 
  left_join(better_taxonomy, c("tax_id" = "tax_id")) %>% 
  rename(DMP_ASSAY_ID = case) %>% 
  select(DMP_ASSAY_ID, readcount, species_name, genus_name, GeneralTumorType, DetailedTumorType, tax_name, MSI_Class, SampleType)

microall_bacteria <- newmicroall %>% unfactor(.) %>% 
  left_join(finaljoin, c("tax_name" ="tax_name")) %>%
  rename(tax_id = tax_id) %>% 
  mutate(isbacteria = ifelse(tax_id==2,1,ifelse(parent_id1==2,1,ifelse(parent_id2==2,1,ifelse(parent_id3==2,1,ifelse(parent_id4==2,1,ifelse(parent_id5==2,1,ifelse(parent_id6==2,1,ifelse(parent_id7==2,1,ifelse(parent_id8==2,1,ifelse(parent_id9==2,1,ifelse(parent_id10==2,1,ifelse(parent_id11==2,1,ifelse(parent_id12==2,1,0)))))))))))))) %>% 
  filter(isbacteria==1) %>% select(1:9) %>% mutate(type="bacteria")


microall_fungi <- newmicroall %>% unfactor(.) %>% 
  left_join(finaljoin, c("tax_name" ="tax_name")) %>%
  rename(tax_id = tax_id) %>% 
  mutate(isfungi = ifelse(tax_id==4751,1,ifelse(parent_id1==4751,1,ifelse(parent_id2==4751,1,ifelse(parent_id3==4751,1,ifelse(parent_id4==4751,1,ifelse(parent_id5==4751,1,ifelse(parent_id6==4751,1,ifelse(parent_id7==4751,1,ifelse(parent_id8==4751,1,ifelse(parent_id9==4751,1,ifelse(parent_id10==4751,1,ifelse(parent_id11==4751,1,ifelse(parent_id12==4751,1,0)))))))))))))) %>% 
  filter(isfungi==1) %>% select(1:9) %>% mutate(type="fungi")

microall_virus <- newmicroall %>% unfactor(.) %>% 
  left_join(finaljoin, c("tax_name" ="tax_name")) %>%
  rename(tax_id = tax_id) %>% 
  mutate(isvirus = ifelse(tax_id==10239,1,ifelse(parent_id1==10239,1,ifelse(parent_id2==10239,1,ifelse(parent_id3==10239,1,ifelse(parent_id4==10239,1,ifelse(parent_id5==10239,1,ifelse(parent_id6==10239,1,ifelse(parent_id7==10239,1,ifelse(parent_id8==10239,1,ifelse(parent_id9==10239,1,ifelse(parent_id10==10239,1,ifelse(parent_id11==10239,1,ifelse(parent_id12==10239,1,0)))))))))))))) %>% 
  filter(isvirus==1) %>% select(1:9) %>% mutate(type="virus") 

species_file <- paste0(args$directory, args$output, "_DB_species.txt")
genus_file <- paste0(args$directory, args$output, "_DB_genus.txt")

filtered_virus_bacteria_fungus_output_species <- 
  rbind(microall_bacteria, microall_fungi, microall_virus) %>% 
  select(DMP_ASSAY_ID, GeneralTumorType, DetailedTumorType, species_name, readcount) %>%
  rename(Taxonomy_ID_Label = species_name) %>% 
  group_by(DMP_ASSAY_ID, Taxonomy_ID_Label) %>% 
  summarise(readcount = sum(readcount)) 
filtered_virus_bacteria_fungus_output_species %>% 
  write.table(species_file, sep = ",", row.names = FALSE)

filtered_virus_bacteria_fungus_output_genus <- 
  rbind(microall_bacteria, microall_fungi, microall_virus) %>% 
  select(DMP_ASSAY_ID, GeneralTumorType, DetailedTumorType, genus_name, readcount) %>%
  rename(Taxonomy_ID_Label = genus_name) %>% 
  group_by(DMP_ASSAY_ID, Taxonomy_ID_Label) %>% 
  summarise(readcount = sum(readcount)) 
filtered_virus_bacteria_fungus_output_genus %>% 
  write.table(genus_file, sep = ",", row.names = FALSE)

  
#functions

# x: Species count vector
shannon <- function(x) {
  
  # Ignore zeroes
  x <- x[x > 0]
  
  # Species richness (number of species)
  S <- length(x)
  
  # Relative abundances
  p <- x/sum(x)
  
  # Shannon index
  (-sum(p * log(p)))
  
}

inverse_simpson <- function(x) {
  
  # Simpson index
  lambda <- simpson_index(x)
  
  # Inverse Simpson diversity
  (1/lambda)
  
}

# x: Species count vector
gini_simpson <- function(x) {
  
  # Simpson index
  lambda <- simpson_index(x)
  
  # Gini-Simpson diversity
  1 - lambda
  
}

simpson_index <- function(x) {
  
  # Relative abundances
  p <- x/sum(x)
  
  # Simpson index
  lambda <- sum(p^2)
  
  lambda
  
}



# x: Species count vector
shannon <- function(x) {
  
  # Ignore zeroes
  x <- x[x > 0]
  
  # Species richness (number of species)
  S <- length(x)
  
  # Relative abundances
  p <- x/sum(x)
  
  # Shannon index
  (-sum(p * log(p)))
  
}

# x: Species count vector
observed <- function(x) {
  
  # Ignore zeroes
  x <- x[x > 0]
  
  # Species richness (number of species)
  (length(x))
  
}






Generate_alphadiv <- 
  function(ClassifcationTask, DB, figfile,  title_fig) {
    #ClassifcationTask <- tmp2
    #DB <- DB_All_Gastric_species
    
    ## Merge dataframe with classification file.  datafram must include column with header 'DMP_ASSAY_ID' and 'Classification'.  Classification must be labeled '1' and '2' as negative and postive classes.
    MergedClassification <-  ClassifcationTask %>% 
      left_join(DB  , c("DMP_ASSAY_ID"="DMP_ASSAY_ID")) %>% 
      filter(is.na(Taxonomy_ID_Label)==FALSE) %>% 
      filter(readcount >1)
    
   
    ## Generate separate dateframe for each class
    #finaljoin <- unfactor(finaljoin)
    MergedClassification1 <- MergedClassification %>%
      #filter(Classification == 1) %>%
      left_join(finaljoin, c("Taxonomy_ID_Label" ="tax_name")) %>%
      mutate(isbacteria = ifelse(tax_id==2,1,ifelse(parent_id1==2,1,ifelse(parent_id2==2,1,ifelse(parent_id3==2,1,ifelse(parent_id4==2,1,ifelse(parent_id5==2,1,ifelse(parent_id5==2,1,ifelse(parent_id7==2,1,ifelse(parent_id8==2,1,ifelse(parent_id9==2,1,0))))))))))) %>% 
      filter(isbacteria==1) %>% 
      filter(parent_id1 != 136841) %>% 
      filter(parent_id1 != 1232139) %>%
      filter(parent_id1 != 136845) %>%
      filter(parent_id1 != 995085) %>%
      group_by(DMP_ASSAY_ID, Taxonomy_ID_Label) %>% summarise(readcount = sum(readcount)) %>% 
      select(Taxonomy_ID_Label, DMP_ASSAY_ID, readcount) %>% 
      spread(key = DMP_ASSAY_ID, value = readcount, fill = 0)
    
    MergedClassification1 <- as.data.frame(MergedClassification1)
    rownames(MergedClassification1) <- MergedClassification1[,1]
    
    matshannon <- as.matrix(MergedClassification1[,2:ncol(MergedClassification1)])
    
    #col1 <- shannon(matshannon[,1])
    s<-NULL
    S<-NULL
    for(i in 1:ncol(matshannon)){s <- rbind(s,shannon(matshannon[,i]))}
    for(i in 1:ncol(matshannon)){S <- rbind(S,observed(matshannon[,i]))}
    si <- NULL
    for(i in 1:ncol(matshannon)){si <- rbind(si,simpson_index(matshannon[,i]))}
    gsi <- NULL
    for(i in 1:ncol(matshannon)){gsi <- rbind(gsi, gini_simpson(matshannon[,i]))}
    invs <- NULL
    for(i in 1:ncol(matshannon)){invs <- rbind(invs,inverse_simpson(matshannon[,i]))}
    #fisha <- rep(0,ncol(matshannon))
    #cover <- rep(0,ncol(matshannon))
    alphascoresmc1 <- data.frame(colnames(matshannon), invs, gsi, s, S)
    colnames(alphascoresmc1) <- c("DMP_ASSAY_ID", "invers_simpson", "gini_simpson", "shannon", "observed")
    
    
    
    MergedClassification2 <- MergedClassification %>% 
      filter(Classification == 2) %>%
      left_join(finaljoin, c("Taxonomy_ID_Label" ="tax_name")) %>%
      mutate(isbacteria = ifelse(tax_id==2,1,ifelse(parent_id1==2,1,ifelse(parent_id2==2,1,ifelse(parent_id3==2,1,ifelse(parent_id4==2,1,ifelse(parent_id5==2,1,ifelse(parent_id5==2,1,ifelse(parent_id7==2,1,ifelse(parent_id8==2,1,ifelse(parent_id9==2,1,0))))))))))) %>% 
      filter(isbacteria==1)%>%
      filter(parent_id1 != 136841) %>% 
      filter(parent_id1 != 1232139) %>%
      filter(parent_id1 != 136845) %>%
      filter(parent_id1 != 995085) %>%
      group_by(DMP_ASSAY_ID, Taxonomy_ID_Label) %>%  summarise(readcount = sum(readcount)) %>% 
      select(Taxonomy_ID_Label, DMP_ASSAY_ID, readcount) %>% 
      spread(key = DMP_ASSAY_ID, value = readcount, fill = 0)
    MergedClassification2 <- as.data.frame(MergedClassification2)
    rownames(MergedClassification2) <- MergedClassification2[,1]
    
    matshannon <- as.matrix(MergedClassification2[,2:ncol(MergedClassification2)])
    
    col1 <- shannon(matshannon[,1])
    s<-NULL
    S<-NULL
    for(i in 1:ncol(matshannon)){s <- rbind(s,shannon(matshannon[,i]))}
    for(i in 1:ncol(matshannon)){S <- rbind(S,observed(matshannon[,i]))}
    si <- NULL
    for(i in 1:ncol(matshannon)){si <- rbind(si,simpson_index(matshannon[,i]))}
    gsi <- NULL
    for(i in 1:ncol(matshannon)){gsi <- rbind(gsi, gini_simpson(matshannon[,i]))}
    invs <- NULL
    for(i in 1:ncol(matshannon)){invs <- rbind(invs,inverse_simpson(matshannon[,i]))}
    fisha <- rep(0,ncol(matshannon))
    cover <- rep(0,ncol(matshannon))
    alphascoresmc2 <- data.frame(colnames(matshannon), invs, gsi, s, S)
    colnames(alphascoresmc2) <- c("DMP_ASSAY_ID", "invers_simpson", "gini_simpson", "shannon", "observed")
    
    ClassifcationTask <- unfactor(ClassifcationTask)
    alphamargedall<-NULL
    alphamargedall <- rbind(alphascoresmc1, alphascoresmc2) 
    
    alphamargedall$DMP_ASSAY_ID <- unfactor(alphamargedall$DMP_ASSAY_ID) 
    
    alphamargedall <- alphamargedall %>%
      right_join(ClassifcationTask, c("DMP_ASSAY_ID" = "DMP_ASSAY_ID")) %>% 
      distinct(DMP_ASSAY_ID, .keep_all = TRUE)
    merge_args <- paste0(args$div_fig_dir,"%s", "_alpha_diversity.csv")
    write.table(alphamargedall, file = sprintf(merge_args, figfile),sep = ",", row.names = FALSE )
    
    alphadivplot <- alphamargedall %>% 
      distinct(DMP_ASSAY_ID, .keep_all = TRUE) %>% filter(shannon>0) %>%  
      mutate(MSI = ifelse(Classification==2,
                          paste0(figfile, "is_TRUE", sep = "_"), ifelse(Classification==1,paste0(figfile, "is_FALSE", sep = "_"),3))) %>% 
      ggplot(aes(x = MSI, y = shannon)) + 
      coord_flip() + 
      #geom_violin(width=0.6, fill= "blue") + 
      geom_boxplot(width=0.8) + geom_jitter() +
      labs(y = "Alpha Diversity (Shannon)", x = "Group") + 
      ggtitle(label = paste0(title_fig, "Shannon Alpha Diversity", sep = " ")) +
      theme_classic(base_size = 25) +
      stat_compare_means(method = "t.test", 
                         label.x = 0.5, label.y = 0,
                         inherit.aes = TRUE) +
      theme(axis.text = element_text(size = 20, color = "black"),  
            axis.title = element_text(size = 25), 
            legend.text = element_text(size = 20), 
            legend.title = element_text(size = 20))
    
    merge_args <- paste0(args$div_fig_dir,"%s", "_shannon_alpha_diversity.pdf")
    file <- sprintf(merge_args, figfile)
    ggsave(file, plot = alphadivplot, dpi = 100, units = "cm", width = 60, height = 30)
    
    alphadivplot <- alphamargedall %>% 
      distinct(DMP_ASSAY_ID, .keep_all = TRUE) %>%
      filter(invers_simpson>0) %>%  
      mutate(MSI = ifelse(Classification==2,
                          paste0(figfile, "is_TRUE", sep = "_"), ifelse(Classification==1,paste0(figfile, "is_FALSE", sep = "_"),3))) %>% 
      ggplot(aes(x = MSI, y = invers_simpson)) + 
      coord_flip() + 
      #geom_violin(width=0.6, fill= "blue") + 
      geom_boxplot(width=0.8) + geom_jitter() +
      labs(y = "Alpha Diversity (Inv Simpson)", x = "Group") + 
      ggtitle(label = paste0(title_fig, "Inverse Simpson Alpha Diversity", sep = " ")) +
      theme_classic(base_size = 25) +
      stat_compare_means(method = "t.test", 
                         label.x = 0.5, label.y = 0,
                         inherit.aes = TRUE) +
      theme(axis.text = element_text(size = 20, color = "black"),  
            axis.title = element_text(size = 25), 
            legend.text = element_text(size = 20), 
            legend.title = element_text(size = 20))
    
    merge_args <- paste0(args$div_fig_dir,"%s", "_inv_simpson_alpha_diversity.pdf")
    file <- sprintf(merge_args, figfile)
    ggsave(file, plot = alphadivplot, dpi = 100, units = "cm", width = 60, height = 30)
    
    alphadivplot <- alphamargedall %>% 
      distinct(DMP_ASSAY_ID, .keep_all = TRUE) %>%
      filter(gini_simpson>0) %>%  
      mutate(MSI = ifelse(Classification==2,
                          paste0(figfile, "is_TRUE", sep = "_"), ifelse(Classification==1,paste0(figfile, "is_FALSE", sep = "_"),3)))%>% 
      ggplot(aes(x = MSI, y = gini_simpson)) + 
      coord_flip() + 
      #geom_violin(width=0.6, fill= "blue") + 
      geom_boxplot(width=0.8) + geom_jitter() +
      labs(y = "Alpha Diversity (Gini Simpson)", x = "Group") + 
      ggtitle(label = paste0(title_fig, "Gini Simpson Alpha Diversity", sep = " ")) +
      theme_classic(base_size = 25) +
      stat_compare_means(method = "t.test", 
                         label.x = 0.5, label.y = 0,
                         inherit.aes = TRUE) +
      theme(axis.text = element_text(size = 20, color = "black"),  
            axis.title = element_text(size = 25), 
            legend.text = element_text(size = 20), 
            legend.title = element_text(size = 20))
    
    merge_args <- paste0(args$div_fig_dir,"%s", "_gini_simpson_alpha_diversity.pdf")
    file <- sprintf(merge_args, figfile)
    ggsave(file, plot = alphadivplot, dpi = 100, units = "cm", width = 60, height = 30)
    
    alphadivplot <- alphamargedall %>% 
      distinct(DMP_ASSAY_ID, .keep_all = TRUE) %>%
      filter(observed>0) %>%  
      mutate(MSI = ifelse(Classification==2,
                          paste0(figfile, "is_TRUE", sep = "_"), ifelse(Classification==1,paste0(figfile, "is_FALSE", sep = "_"),3))) %>% 
      ggplot(aes(x = MSI, y = observed)) + 
      coord_flip() + 
      #geom_violin(width=0.6, fill= "blue") + 
      geom_boxplot(width=0.8) + geom_jitter() +
      labs(y = "Alpha Diversity (Inv Simpson)", x = "Group") + 
      ggtitle(label = paste0(title_fig, "Observed Alpha Diversity", sep = " ")) +
      theme_classic(base_size = 25) +
      stat_compare_means(method = "t.test", 
                         label.x = 0.5, label.y = 0,
                         inherit.aes = TRUE) +
      theme(axis.text = element_text(size = 20, color = "black"),  
            axis.title = element_text(size = 25), 
            legend.text = element_text(size = 20), 
            legend.title = element_text(size = 20))
    
    merge_args <- paste0(args$div_fig_dir,"%s", "_observed_alpha_diversity.pdf")
    file <- sprintf(merge_args, figfile)
    ggsave(file, plot = alphadivplot, dpi = 100, units = "cm", width = 60, height = 30)
    
  }


Generate_odds <- 
  function(ClassifcationTask, DB, figfile, txtfile, title_fig) {
 
  #DB <- Gall_bladder_DB
  
  ## Merge dataframe with classification file.  datafram must include column with header 'DMP_ASSAY_ID' and 'Classification'.  Classification must be labeled '1' and '2' as negative and postive classes.
  MergedClassification <-  ClassifcationTask %>% 
    left_join(DB  , c("DMP_ASSAY_ID"="DMP_ASSAY_ID")) %>% 
    filter(is.na(Taxonomy_ID_Label)==FALSE) %>% 
    filter(readcount > 1)
  
  Query1 <- MergedClassification %>% 
    distinct(DMP_ASSAY_ID,.keep_all = TRUE) %>% 
    filter(Classification == 1)
  Query1 <- length(Query1$Classification)
  
  Query2 <- MergedClassification %>% 
    distinct(DMP_ASSAY_ID,.keep_all = TRUE) %>% 
    filter(Classification == 2)
  Query2 <- length(Query2$Classification)
  
  
  ## Generate separate dateframe for each class
  #finaljoin <- unfactor(finaljoin)
  MergedClassification1 <- MergedClassification %>% 
    filter(Classification == 1) 
  
  MergedClassification2 <- MergedClassification %>% 
    filter(Classification == 2) 
  
  ## Rename columns
  MergedClassification1 <- MergedClassification1 %>% 
    group_by(Taxonomy_ID_Label, Classification) %>% 
    tally(sort = FALSE)
  names(MergedClassification1) <- c("Taxonomy_ID_Label", "NegativeClass", "NegativeClassPositive")
  
  MergedClassification2 <- MergedClassification2 %>% group_by(Taxonomy_ID_Label, Classification) %>% tally(sort = FALSE)
  names(MergedClassification2) <- c("Taxonomy_ID_Label", "PositiveClass", "PositiveClassPositive")
  
  totaltested <- MergedClassification2 %>% filter(PositiveClassPositive >2) 
  totaltested <- length(totaltested$Taxonomy_ID_Label)
  
  ## Merge by species for single row for each taxonomyID
  
  OddsMerged <- NULL
  
  OddsMerged <- full_join(MergedClassification1, MergedClassification2, c("Taxonomy_ID_Label" = "Taxonomy_ID_Label"))
  ## Haldane-Anscombe correction
  
  OddsMerged$NegativeClassPositive[is.na(OddsMerged[,"NegativeClassPositive"])] <- 0.5
  OddsMerged$PositiveClassPositive[is.na(OddsMerged[,"PositiveClassPositive"])] <- 0.5
  OddsMerged$PositiveClass[is.na(OddsMerged[,"PositiveClass"])] <- 2
  OddsMerged$NegativeClass[is.na(OddsMerged[,"NegativeClass"])] <- 1
  
  
  #OddsMerged <- OddsMerged %>% 
   # mutate(NegativeClassPositive = as.numeric(NegativeClassPositive), PositiveClassPositive = as.numeric(PositiveClassPositive), PositiveClass = is.numeric(PositiveClass), NegativeClass = is.numeric(NegativeClass)) %>% 
    #replace_na(list(NegativeClassPositive = 0.5, PositiveClassPositive = 0.5, PositiveClass = 2, NegativeClass = 1))
  ## Calculate Negative value for each class
  
  
  OddsMerged <- OddsMerged %>% mutate(NegativeClassNegative = Query1-NegativeClassPositive)
  OddsMerged <- OddsMerged %>% mutate(PositiveClassNegative = Query2-PositiveClassPositive)
  
  
  ## Rename columns
  
  OddsMerged <- OddsMerged %>% select(Taxonomy_ID_Label, NegativeClass, PositiveClass, PositiveClassPositive, PositiveClassNegative, NegativeClassPositive, NegativeClassNegative)
  
  
  
  ## calculate oddsratio and generate dataframe containing calculations
  ##with bonferronni
  k <- oddsratio(OddsMerged$PositiveClassPositive, OddsMerged$PositiveClassNegative, OddsMerged$NegativeClassPositive, OddsMerged$NegativeClassNegative, conf.level = 0.95, p.calc.by.independence = FALSE)
  ##without bonferroni
  #k <- oddsratio(OddsMerged$PositiveClassPositive, OddsMerged$PositiveClassNegative, OddsMerged$NegativeClassPositive, OddsMerged$NegativeClassNegative, conf.level = 0.95, p.calc.by.independence = FALSE)
  #k <- oddsratio(OddsMerged$PositiveClassPositive, 
          #       OddsMerged$PositiveClassNegative, 
           #      OddsMerged$NegativeClassPositive, 
            #     OddsMerged$NegativeClassNegative, conf.level = (0.996), p.calc.by.independence = FALSE)
  
  OddsToPlot <- data.frame(OddsMerged[,1:7], k$estimate, k$conf.int[1:(length(k$conf.int)/2)], k$conf.int[(((length(k$conf.int))/2)+1):length(k$conf.int)], k$p.value)
  
  
  ## Rename collumns
  
  names(OddsToPlot) <- c("Taxonomy_ID_Label", "NegativeClass", "PositiveClass", "PositiveClassPositive", "PositiveClassNegative", "NegativeClassPositive", "NegativeClassNegative", "OddsRatio", "OddsRatio95CIlow", "OddsRatio95CIhigh", "pvalue")
  
  ## Sort for most significant TaxonomyIDS and filter for top performers
  
  OddsToPlot <- OddsToPlot %>% 
    filter(Taxonomy_ID_Label!="NA") %>% 
    dplyr::rename(newcol = Taxonomy_ID_Label) %>% 
    mutate(Taxonomy_ID_Label = paste(newcol,
                                      PositiveClassPositive,
                                      PositiveClassNegative,
                                      NegativeClassPositive, 
                                      NegativeClassNegative, round(pvalue,6), sep = "_")) %>% 
    #select(-newcol) %>% 
    arrange(desc(OddsRatio95CIlow))  
  ## Plot only top X species
  OddstoPlot10 <- OddsToPlot[1:25,]
  #OddstoPlot10 <-tmp
  ## Plot by generating dataframe for GGplot
  
  df10 <- data.frame(yAxis = length(OddstoPlot10$Taxonomy_ID_Label):1, 
                     boxOdds = OddstoPlot10$OddsRatio, 
                     boxCILow = OddstoPlot10$OddsRatio95CIlow, 
                     boxCIHigh = OddstoPlot10$OddsRatio95CIhigh)
  
  K <- ggplot(df10, aes(x = boxOdds, y = yAxis))
  ## Bells & whistles
  Autoplot <- K + geom_vline(aes(xintercept = 1.01), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 1, height = .8, color = "black") + 
    geom_point(shape=20, size = 15, color = "#33bfff") + 
    ggtitle(label =  title_fig) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank())+
    theme(title = element_text(size = 40,hjust = 0),
      axis.text.x = element_text(size = 25, face = "bold",color = "black"), 
      axis.text.y = element_text(size = 30, face = "italic",color = "black"), 
      axis.title.x = element_text(size = 20,face = "bold"))  + 
    scale_x_continuous(breaks = c(0.01,0.05, 0.1, 0.50, 5, 10, 50)) + 
    coord_trans(x = "log10") + ylab("") + xlab("Odds ratio (log scale)") + 
    scale_y_continuous(labels = OddstoPlot10$Taxonomy_ID_Label, breaks = length(OddstoPlot10$Taxonomy_ID_Label):1)
  
  ## plot visualize
  Autoplot       
  ## plot save
  merge_args <- paste0(args$odds_fig_dir,"%s",".pdf")
  file <- sprintf(merge_args, figfile)
  ggsave(file, plot = Autoplot, dpi = 100, units = "cm", width = 60, height = 30)
  merge_args <- paste0(args$odds_fig_dir,"%s",".txt")
  file2 <- sprintf(merge_args, txtfile)
  write.table(file =  file2, OddsToPlot, sep = "\t", row.names = FALSE)
}





tmp2 <- read.delim(args$manifest, sep = ",", header = T)

for (i in 2:ncol(tmp2)) {
  ClassifcationTask <- tmp2 %>%  
    select(DMP_ASSAY_ID, i) %>% 
    rename(new = 2) %>% 
    mutate(Classification = ifelse(new==2, 2, ifelse(new==1, 1, 3)) ) %>%
    select(-new) %>% 
    filter(Classification<3)
  Generate_odds(ClassifcationTask = ClassifcationTask,
                DB = filtered_virus_bacteria_fungus_output_genus %>% ungroup(),
                figfile = paste("Enriched_if_TRUE", colnames(tmp2[i]), sep = "_"),
                txtfile = paste("odds_ratios", colnames(tmp2[i]), sep = "_"),
                title_fig = paste("Enriched_if", colnames(tmp2[i]), "is_TRUE", sep = "_"))
  Generate_alphadiv(ClassifcationTask = ClassifcationTask,
                    DB = filtered_virus_bacteria_fungus_output_genus %>% ungroup(),
                    figfile = colnames(tmp2[i]),
                    title_fig = colnames(tmp2[i]))
  ClassifcationTask <- tmp2 %>%  select(DMP_ASSAY_ID, i) %>% rename(new = 2) %>% 
    mutate(Classification = ifelse(new==2, 1, ifelse(new==1, 2, 3)) ) %>% 
    select(-new) %>% 
    filter(Classification<3)
  Generate_odds(ClassifcationTask = ClassifcationTask, 
                DB = filtered_virus_bacteria_fungus_output_genus %>% ungroup(),
                figfile = paste("Enriched_if_FALSE", colnames(tmp2[i]), sep = "_"),
                txtfile = "tmp_ovewrite", 
                title_fig = paste("Enriched_if", colnames(tmp2[i]), "is_FALSE", sep = "_")  )
}

