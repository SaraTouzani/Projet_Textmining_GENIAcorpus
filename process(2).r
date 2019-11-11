library(xml2)
library(doParallel)

############ Classes

molecule = c(

  "G#organic"	 ,
  "G#inorganic",
  "G#atom"	,
  "G#carbohydrate"	,
  "G#lipid",
  "G#amino_acid_monomer",
  "G#peptide",
  "G#protein_N/A"	,
  "G#protein_complex"	,
  "G#protein_domain_or_region"	,
  "G#protein_family_or_group"	,
  "G#protein_molecule"	,
  "G#protein_substructure"	,
  "G#protein_subunit",
  "G#other_organic_compound"

)
gene = c (
  "G#nucleotide"	,
  "G#polynucleotide"	,
  "G#DNA_N/A"	,
  "G#DNA_domain_or_region"	,
  "G#DNA_family_or_group"	,
  "G#DNA_molecule"	,
  "G#DNA_substructure"	,
  "G#RNA_N/A"	,
  "G#RNA_domain_or_region"	,
  "G#RNA_family_or_group"	,
  "G#RNA_molecule"	,
  "G#RNA_substructure"
  
)

Cell = c("G#virus"	,
          "G#mono_cell"	,
          "G#tissue",
          "G#multi_cell"	,
          "G#body_part")

Cell_term = c("G#cell_type"	,
              "G#cell_component"	,
              "G#cell_line"	,
              "G#other_artificial_source")

other_name = "G#other_name"
cordinated= c("G#cordinated")

#######################

cl <-makeCluster(3)
registerDoParallel(cl)
clusterExport(cl, c( 'proteine','gene','Cell','Cell_term','other_name','cordinated'))


##################################### extract futures


molecules <-
  foreach::foreach(i = 1:length(molecule), .combine = rbind) %dopar% {
    Geniacorpus.data <- 	data.frame(therme_type = NULL,
                                    therme_key = NULL,
                                    therme_value = NULL)
    doc = XML::xmlParse("~/GENIAcorpus3.02.xml")
    tag = sprintf("//cons[@sem=\'%s\']\n", molecule[i])
    node_withtags <- XML::getNodeSet(doc, tag)
    values = sapply(node_withtags, XML::xmlValue)
    if(length( values)>0){
    for (j in 1:length(values)) {
      Geniacorpus.newdata <- 	data.frame(therme_type = "moleculle",
                                         therme_key = molecule[i],
                                         therme_value = values[j])
      Geniacorpus.data  <-
        rbind(Geniacorpus.data , Geniacorpus.newdata)
    }}
    else{
      Geniacorpus.data <- 	data.frame(therme_type = "moleculle",
                                     therme_key = molecule[i] ,
                                     therme_value = "NONE")
    }
    
    return(Geniacorpus.data)
  }


genes <-
  foreach::foreach(i = 1:length(gene), .combine = rbind) %dopar% {
    Geniacorpus.data <- 	data.frame(therme_type = NULL,
                                    therme_key = NULL,
                                    therme_value = NULL)
    doc = XML::xmlParse("~/GENIAcorpus3.02.xml")
    tag = sprintf("//cons[@sem=\'%s\']\n", gene[i])
    node_withtags <- XML::getNodeSet(doc, tag)
    values = sapply(node_withtags, XML::xmlValue)
    for (j in 1:length(values)) {
      Geniacorpus.newdata <- 	data.frame(therme_type = "gene",
                                         therme_key = gene[i],
                                         therme_value = values[j])
      Geniacorpus.data  <-
        rbind(Geniacorpus.data , Geniacorpus.newdata)
    }
    return(Geniacorpus.data)
  }

Cells <-
  foreach::foreach(i = 1:length(Cell), .combine = rbind) %dopar% {
    Geniacorpus.data <- 	data.frame(therme_type = NULL,
                                    therme_key = NULL,
                                    therme_value = NULL)
    doc = XML::xmlParse("~/GENIAcorpus3.02.xml")
    tag = sprintf("//cons[@sem=\'%s\']\n", Cell[i])
    node_withtags <- XML::getNodeSet(doc, tag)
    values = sapply(node_withtags, XML::xmlValue)
    for (j in 1:length(values)) {
      Geniacorpus.newdata <- 	data.frame(therme_type = "Cell",
                                         therme_key = Cell[i],
                                         therme_value = values[j])
      Geniacorpus.data  <-
        rbind(Geniacorpus.data , Geniacorpus.newdata)
    }
    return(Geniacorpus.data)
  }


Cells_term <-
  foreach::foreach(i = 1:length(Cell_term), .combine = rbind) %dopar% {
    Geniacorpus.data <- 	data.frame(therme_type = NULL,
                                    therme_key = NULL,
                                    therme_value = NULL)
    doc = XML::xmlParse("~/GENIAcorpus3.02.xml")
    tag = sprintf("//cons[@sem=\'%s\']\n", Cell_Term[i])
    node_withtags <- XML::getNodeSet(doc, tag)
    values = sapply(node_withtags, XML::xmlValue)
    for (j in 1:length(values)) {
      Geniacorpus.newdata <- 	data.frame(therme_type = "Cell_term",
                                         therme_key = Cell_term[i],
                                         therme_value = values[j])
      Geniacorpus.data  <-
        rbind(Geniacorpus.data , Geniacorpus.newdata)
    }
    return(Geniacorpus.data)
  }

other_names <-
  foreach::foreach(i = 1:length(other_name), .combine = rbind) %dopar% {
    Geniacorpus.data <- 	data.frame(therme_type = NULL,
                                    therme_key = NULL,
                                    therme_value = NULL)
    doc = XML::xmlParse("~/GENIAcorpus3.02.xml")
    tag = sprintf("//cons[@sem=\'%s\']\n", other_name[i])
    node_withtags <- XML::getNodeSet(doc, tag)
    values = sapply(node_withtags, XML::xmlValue)
    for (j in 1:length(values)) {
      Geniacorpus.newdata <- 	data.frame(therme_type = "othernames",
                                         therme_key = other_name[i],
                                         therme_value = values[j])
      Geniacorpus.data  <-
        rbind(Geniacorpus.data , Geniacorpus.newdata)
    }
    return(Geniacorpus.data)
  }



cordinateds <-
  foreach::foreach(i = 1:length(cordinated) , .combine = rbind ) %dopar% {
    Geniacorpus.data <- 	data.frame(therme_type = NULL,
                                    therme_key = NULL,
                                    therme_value = NULL)
    doc = XML::xmlParse("~/GENIAcorpus3.02.xml")
    tag = sprintf("//cons[@sem=\'%s\']\n", cordinated[i])
    node_withtags <- XML::getNodeSet(doc, tag)
    values = sapply(node_withtags, XML::xmlValue)
    if(length( values)>0){
    for (j in 1:length(values)) {
 
        Geniacorpus.newdata <- 	data.frame(therme_type = "cordinated",
                                           therme_key =cordinated[i],
                                           therme_value = values[j])
        
        Geniacorpus.data  <-
          rbind(Geniacorpus.data , Geniacorpus.newdata)


        
      
    }}
    
    else {
      Geniacorpus.newdata <- 	data.frame(therme_type = "cordinated",
                                      therme_key =cordinated[i],
                                      therme_value = "NONE")
      Geniacorpus.data  <-
        rbind(Geniacorpus.data , Geniacorpus.newdata)
    }

    return(Geniacorpus.data)
  } ## ne pas executer
################"
data<-data.frame()
data<-rbind(data,cordinateds)
data<-rbind(data,other_names)
data<-rbind(data,Cells_term)
data<-rbind(data,Cells)
data<-rbind(data,genes)
data<-rbind(data,molecules)


Uniq <- unique(data[, "therme_value"])

others<-data.frame()





##############
library(xml2)
library(tm)
library(stringr)
library(stringi)

pg <- read_xml("~/GENIAcorpus3.02.xml")

# get all the <record>s
recs <- xml_find_all(pg, "//article/abstract/sentence")

# extract and clean all the columns
vals <- trimws(xml_text(recs))
c<-paste(vals[1], collapse = " ")
for (i in 2:length(vals)) {
   c<-paste(c,vals[i], collapse = " ")
}


## remove ponctuation
p<-gsub(pattern = "\\W", replacement = " ", c)
## remove numbers
n<-gsub(pattern = "\\d", replacement = " ", p)
## llower case
l<-tolower(n)
## remove stopwords
w<-removeWords(l, stopwords())
## remove lettrs of lenght one  remaining from the cleaning done before
one<-gsub(pattern = "\\b[A-z]\\b{1}", replacement = " ", w)
## stripwhitespace
s<-stripWhitespace(one)
#extract text bag
sp<-str_split(s, pattern = "\\s")
textbag<-unlist(sp)
dic<-unique(unlist(sp))

uniq_char<-tolower(as.character(Uniq))

textbag2<-textbag

for (i in 1:length(uniq_char)) {
 for (j in 1:length(dic)) {
   if(stri_cmp(uniq_char[i],dic[j])==0){
     textbag2<-removeWords(textbag2,dic[j])
   }
 }
  
}


#######
Geniacorpus.data <- 	data.frame(therme_type = NULL,
                                therme_key = NULL,
                                therme_value = NULL)
for (i in 1:length(textbag2) ) {
  if(textbag2[i]!=""){
    Geniacorpus.newdata <- 	data.frame(therme_type = "others",
                                       therme_key ="others",
                                       therme_value = textbag2[i])
    Geniacorpus.data  <-rbind(Geniacorpus.data , Geniacorpus.newdata)    
  }

  
}
data<-rbind(data,Geniacorpus.data)
write.csv(data, file = "MyData.csv")


