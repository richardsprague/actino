# test-phyloseq.r
# test that it's a proper phyloseq object

library(phyloseq)
library(actino)
context("Prove it's a proper phyloseq object")

DATA_DIR <- system.file("extdata", package = "actino") # "../../inst/extdata"

data("kombucha.genus")
p <- kombucha.genus
# save some time during testing by using the following pre-loaded data
pj <- phyloseq_from_JSON_at_rank(just_json_files_in(DATA_DIR),paste0(DATA_DIR,"/kombucha-mapfile.xlsx"))
#save(pj,file=paste0(DATA_DIR,"/pj.RData"))
#load(paste0(DATA_DIR,"/pj.RData"))

d <- paste0(DATA_DIR,"/dummyData/2017-01-10-dummy-755190.json")
d.j <- jsonlite::fromJSON(d)[["ubiome_bacteriacounts"]]
d.p <- phyloseq_from_JSON(d, paste0(DATA_DIR,"/dummyData/Dummy_User_Mapfile.xlsx"))

test_that("dummy data works",{
  expect_equal(ntaxa(d.p),12)
  expect_equal(as.character(tax_table(d.p)[4,"Genus"]), "Bifidobacterium")
  expect_equal(as.numeric(otu_table(d.p)["Root",1]),75924) # total counts in dummmy file

})

test_that(paste("Current directory=",getwd()),{
  #print("running phyloseq tests now")
  #print(getwd())
  #expect_equal(TRUE,FALSE)
})

test_that("phyloseq_from_JSON_at_rank generates a correct date value",{
  expect_equal(as.character(sample_data(pj)$Date[2]),"2016-07-27")
  expect_equal(class(sample_data(pj)$Date),c("POSIXct","POSIXt"))
})



test_that("Phyloseq object has correct OTU value",{
  expect_equal(as.numeric(otu_table(p)[10,3]),2786) # just an arbitrary element
})

test_that("Phyloseq object has correct number of samples",{
  expect_equal(nsamples(p),19)
})

test_that("Phyloseq tax_table has correct rank names",{

  expect_equal(colnames(tax_table(pj)), "Genus")
  #expect_equal(colnames(tax_table(pj)), c("Root","Genus"))
  expect_equal(colnames(tax_table(p)), "Genus")

})


test_that("Phyloseq object has correct rank names",{
  expect_equal(rownames(tax_table(pj))[1:3], c("Flavobacterium","Kluyvera","Bacteroides"))
})


test_that("create Phyloseq object from a directory of json files",{
  expect_equal(length(just_json_files_in(DATA_DIR)),19)
  pj <- phyloseq_from_JSON_at_rank(just_json_files_in(DATA_DIR),paste0(DATA_DIR,"/kombucha-mapfile.xlsx"))
  expect_equal(nsamples(pj),nsamples(p))
  expect_equivalent(p,pj)

})

test_that("create full Phyloseq object from a directory of json files",{
  jsonFileList <- just_json_files_in(DATA_DIR)[1:3]
  expect_equal(length(jsonFileList),3)
  pj <- phyloseq_from_JSON(jsonFileList,paste0(DATA_DIR,"/kombucha-mapfile.xlsx"))
  expect_equal(get_taxa_unique(pj,taxonomic.rank = "Rank1"),"Root")
  firmicutes <- subset_taxa(pj, Phylum=="Firmicutes")
  blautia <- subset_taxa(firmicutes, Genus == "Blautia")
  expect_equal(as.numeric(otu_table(blautia)[2,2]),5522)
  expect_equal(length(otu_table(firmicutes)[,2]),109)
  expect_equal(as.numeric(otu_table(firmicutes)[1,2]),4300) # abundance of Roseburia in 2nd sample (ssr=73837)

})

test_that("correct metadata is loaded into phyloseq object",{
  #pj <- phyloseq_from_JSON_at_rank(just_json_files_in(DATA_DIR),paste0(DATA_DIR,"/kombucha-mapfile.xlsx"))
  expect_equal(as.character(sample_data(kombucha.genus)[9,5]),"pizza")
})

test_that("tax_abbrev returns correct shortened abbrev for a tax_rank",{
  expect_equal(tax_abbrev("genus"),"g")
  expect_equal(tax_abbrev("family"),"f")
  expect_equal(tax_abbrev("s",inverse=TRUE),"species")
  expect_equal(tax_abbrev("r",inverse=TRUE),"root")
  expect_equal(tax_abbrev("i",inverse=TRUE),"subgenus")
  expect_equal(tax_abbrev("k",inverse=TRUE),"superkingdom")
  expect_equal(tax_abbrev("subphylum"),"3")
})

test_that("full_taxa returns correct qiime-formated name for a taxa",{
 expect_equal(full_taxa(kombucha.csv,kombucha.csv[3,]), "Root;n__cellular organisms;k__Bacteria;p__Proteobacteria;3__delta/epsilon subdivisions;c__Epsilonproteobacteria;o__Campylobacterales;f__Campylobacteraceae;g__Campylobacter")
 #expect_equal(full_taxa(kombucha.csv,as.list(kombucha.csv[3,])),"Root;n__cellular organisms;k__Bacteria;p__Proteobacteria;3__delta/epsilon subdivisions;c__Epsilonproteobacteria;o__Campylobacterales;f__Campylobacteraceae;g__Campylobacter")

})

test_that("Unpacking a qiime taxa returns correct tax_name and rank",{
  sample_taxon = "Root;n__cellular organisms;k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Alcaligenaceae"
  sample_taxon2 = "Root;n__cellular organisms;k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Kluyvera"
  sample_taxon3 = "Root;n__cellular organisms;k__Bacteria;p__Actinobacteria;c__Actinobacteria;k__Actinobacteridae;o__Actinomycetales;r__Corynebacterineae;f__Corynebacteriaceae;g__Corynebacterium;s__Corynebacterium sp. NML 97-0186"
  expect_equal(tax_rank_of_full_taxa(sample_taxon2),list("genus","Kluyvera"))
  expect_equal(tax_rank_of_full_taxa(sample_taxon),
                                     list("family","Alcaligenaceae"))
  expect_equal(tax_rank_of_full_taxa("Root"),list("Root","Root"))
  expect_equal(tax_rank_of_full_taxa(sample_taxon3),list("species","Corynebacterium sp. NML 97-0186"))
}
)

