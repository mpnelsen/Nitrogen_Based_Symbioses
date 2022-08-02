#Remove duplicates, unnamed, ambiguously named (cf, aff) and subspecies, so only left with species
#Then find out whether genus is ECM or not

require(corHMM)

#Trees and Species Name Match from from Varga et al. were downloaded from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.gc2k9r9


#Read one tree to test
tr<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_8.tree2")

#read in all trees, and make multiphylo
tr8<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_8.tree2")
tr12<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_12.tree2")
tr16<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_16.tree2")
tr53<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_53.tree2")
tr88<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_88.tree2")
tr165<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_165.tree2")
tr170<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_170.tree2")
tr188<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_188.tree2")
tr216<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_216.tree2")
tr222<-read.tree(file=".../PhyloBayes_FastDate/FastDate_analysis/fastdate_kronogram_222.tree2")
trs<-c(tr8,tr12,tr16,tr53,tr88,tr165,tr170,tr188,tr216,tr222)
class(trs)<-"multiPhylo"

#read species match file
conv<-read.csv(file="...Species_name_match.csv",stringsAsFactors=FALSE,sep="\t")

###############RECONCILING TRAIT DATA FROM VARGA ET AL W PHYLOGENIES
##read  trait data from Varga et al. and add column for Speceis_names_final to it
#stuff<-read.csv(file=".../41559_2019_834_MOESM3_ESM.csv",stringsAsFactors=FALSE)
#stuff<-stuff[!stuff$Species %in% "",]
#nrow(stuff)
#
##conv$Species_names_final[!conv$Species_names_final %in% stuff$Species]
##stuff$Species[!stuff$Species %in% conv$Species_names_final]
##"Tricholoma_sp_G1208"
#
##tr$tip.label[!tr$tip.label %in% stuff$Species]
##tr$tip.label[!tr$tip.label %in% conv$BayesTrait_names]
##tr$tip.label[!tr$tip.label %in% conv$Alignment_names]
##tr$tip.label[!tr$tip.label %in% conv$Species_names_final]
#
##conv$BayesTrait_names[conv$BayesTrait_names %in% conv$Species_names_final]
#
##stuff$Species_names_final<-NA
##stuff$BayesTrait_names<-NA
##for(x in 1:nrow(stuff)){
##	if(length(conv$BayesTrait_names[conv$BayesTrait_names %in% stuff$Species[x]])==1){
##		stuff$BayesTrait_names[x]<-conv$BayesTrait_names[conv$BayesTrait_names %in% stuff$Species[x]]
##		stuff$Species_names_final[x]<-conv$Species_names_final[conv$BayesTrait_names %in% stuff$Species[x]]
##	}
##}
##sum(is.na(stuff$Species_names_final))
##many of the names that would have been used in BayesTraits analysis for angio/gymno do not match
#
#table(stuff[,4])
##strange...no "?"'s...figure 1 in Varga et al. has angio, gymno, generalist and missing and is discussed that way in M&M, but file has angio, gymno or uncertain (presumably is generalist+missing)...since missing seem like a small fraction from Fig 1, it might make sense to code as generalists.
#
##far more match from Species_names_final
#stuff$Species_names_final<-NA
#for(x in 1:nrow(stuff)){
#	if(length(conv$Species_names_final[conv$Species_names_final %in% stuff$Species[x]])==1){
#		stuff$Species_names_final[x]<-conv$Species_names_final[conv$Species_names_final %in% stuff$Species[x]]
#	}
#}
#sum(is.na(stuff$Species_names_final))
#stuff.missing<-sort(stuff$Species[is.na(stuff$Species_names_final)])
##these are duplicate species...not sure these can be linked to specific taxa, at least the duplicates...
#stuff.missing
#
##so...do same as above, but this time check if length if duplicates if length is greater than one and put that in.
#stuff$Species_names_final<-NA
#for(x in 1:nrow(stuff)){
#	stuff$Species_names_final[x]<-unique(conv$Species_names_final[conv$Species_names_final %in% stuff$Species[x]])
#}
#sum(is.na(stuff$Species_names_final))
###############
###############
###############


#some species in tree are not in Alignment_names or BayesTrait_names
tr$tip.label[!tr$tip.label %in% conv$Alignment_names]
conv$Alignment_names[!conv$Alignment_names %in% tr$tip.label]

tr$tip.label[!tr$tip.label %in% conv$BayesTrait_names]
conv$BayesTrait_names[!conv$BayesTrait_names %in% tr$tip.label]

#tr$tip.label[!tr$tip.label %in% conv$Species_names_final]

#first adjust the tree names to match either the alignment or bayestraits.
#> tr$tip.label[!tr$tip.label %in% conv$Alignment_names]
#[1] "Melanoleuca_sp._blue_stipe_G0291_NL52"
#[2] "Galerella_nigeriensis_strain"         
#[3] "Eonema_pyriforme_lett_Agaricales_G106"
#> conv$Alignment_names[!conv$Alignment_names %in% tr$tip.label]
#[1] "Eonema_pyriforme_lett_Agaricales!_G106"  
#[2] "GAlerella_nigeriensis_strain_CNF1/5859_J"
#[3] "Melanoleuca_sp._blue_stipe_!_G0291_NL52" 
#[4] NA                                        
#> tr$tip.label[!tr$tip.label %in% conv$BayesTrait_names]
#[1] "Entoloma_sp._thujae_G0016_NL5095"     
#[2] "Melanoleuca_sp._blue_stipe_G0291_NL52"
#[3] "Eonema_pyriforme_lett_Agaricales_G106"
#> conv$BayesTrait_names[!conv$BayesTrait_names %in% tr$tip.label]
#[1] "Entoloma_sp._thujae_!_G0016_NL5095"     
#[2] "Eonema_pyriforme_lett_Agaricales!_G106" 
#[3] "Melanoleuca_sp._blue_stipe_!_G0291_NL52"
#[4] NA               

###############rename for individual test tree and then for all trees
tr$tip.label[tr$tip.label %in% "Entoloma_sp._thujae_G0016_NL5095"]<-"Entoloma_sp._thujae_!_G0016_NL5095"
tr$tip.label[tr$tip.label %in% "Melanoleuca_sp._blue_stipe_G0291_NL52"]<-"Melanoleuca_sp._blue_stipe_!_G0291_NL52"
tr$tip.label[tr$tip.label %in% "Eonema_pyriforme_lett_Agaricales_G106"]<-"Eonema_pyriforme_lett_Agaricales!_G106"
tr$tip.label[!tr$tip.label %in% conv$BayesTrait_names]

for(x in 1:length(trs)){
	trs[[x]]$tip.label[trs[[x]]$tip.label %in% "Entoloma_sp._thujae_G0016_NL5095"]<-"Entoloma_sp._thujae_!_G0016_NL5095"
	trs[[x]]$tip.label[trs[[x]]$tip.label %in% "Melanoleuca_sp._blue_stipe_G0291_NL52"]<-"Melanoleuca_sp._blue_stipe_!_G0291_NL52"
	trs[[x]]$tip.label[trs[[x]]$tip.label %in% "Eonema_pyriforme_lett_Agaricales_G106"]<-"Eonema_pyriforme_lett_Agaricales!_G106"
	trs[[x]]$tip.label[!trs[[x]]$tip.label %in% conv$BayesTrait_names]
}

#sort(trs[[x]]$tip.label)


###############

#some of the final names are duplicated, meaning some of the final ones in the tree are also duplicated
dd<-conv$Species_names_final[duplicated(conv$Species_names_final)]
sort(conv$BayesTrait_names[conv$Species_names_final %in% dd])

#drop the duplicated one
drop.dup<-conv$BayesTrait_names[duplicated(conv$Species_names_final)]
drop.dup

#save this for when when dropping taxa in future...this is important, because other trees have tips in different order, so which one is considered the "duplicate" may vary across trees
write.csv(drop.dup,file="...drop.dup.csv")

#read back in and make to vector
#drop.dup<-read.csv(file="...drop.dup.csv",stringsAsFactors=FALSE)
#drop.dup<-drop.dup[,2]

###############drop for single tree and then for all trees
tr
tr<-drop.tip(tr,drop.dup)
tr

#drop for all trees
for(x in 1:length(trs)){
	trs[[x]]<-drop.tip(trs[[x]],drop.dup)
}
trs[[1]]

#################
###############
#########
################then adjust names in tree to be the final names for test tree and all trees
for(x in 1:length(tr$tip.label)){
	tr$tip.label[x]<-conv$Species_names_final[conv$BayesTrait_names %in% tr$tip.label[x]]
}
tr

for(x in 1:length(trs)){
	for(z in 1:length(trs[[x]]$tip.label)){
		trs[[x]]$tip.label[z]<-conv$Species_names_final[conv$BayesTrait_names %in% trs[[x]]$tip.label[z]]
	}
}
trs[[1]]
###############
###############
###############
######

#remove empty rows at end
tax<-tr$tip.label

drops<-c("_var_","_subsp_","_f_","_aff$","_aff_","_sp_","_sp$","_cf_","_cf$","Cortinarius_sp. 'pyrrhomarmarus'_G1174","_sp[0-9]")

droptaxa<-NULL
for(x in 1:length(drops)){
	newdrops<-tax[grep(drops[x],tax)]
	droptaxa<-c(droptaxa,newdrops)
}

#first remove those that are duplicated...NONE ARE DUPLICATED
taxred<-tax[!tax %in% droptaxa]
length(taxred)
taxred<-taxred[!duplicated(taxred)]
length(taxred)

#remove endings
taxred.mod<-sub("^([^_]*_[^_]*).*", "\\1", taxred)
length(taxred.mod)
#get original name prior to stripping endings
drops.taxred<-taxred[duplicated(taxred.mod)]
length(drops.taxred)

#save this for when when dropping taxa in future...this is important, because other trees have tips in different order, so which one is considered the "duplicate" may vary across trees
write.csv(drops.taxred,file="...drops.taxred.csv")

#read back in and make to vector
#drops.taxred<-read.csv(file=".../drops.taxred.csv",stringsAsFactors=FALSE)
#drops.taxred<-drops.taxred[,2]

length(taxred)
length(taxred.mod)
taxred<-taxred[!taxred %in% drops.taxred]
length(taxred)
#now clean up again
taxred<-sub("^([^_]*_[^_]*).*", "\\1", taxred)
taxred[duplicated(taxred)]

tax<-as.data.frame(as.matrix(taxred),stringsAsFactors=FALSE)
colnames(tax)<-"Taxon"
tax$Genus<-NA
tax$Species<-NA

for(x in 1:nrow(tax)){
	tax$Genus[x]<-strsplit(tax$Taxon[x],split="_")[[1]][1]
	tax$Species[x]<-strsplit(tax$Taxon[x],split="_")[[1]][2]
}

tax.gen<-sort(unique(tax$Genus))

##
###############
#READ IN TRAIT DATA
#https://plutof.ut.ee/#/doi/10.15156/BIO/807446
#from Tedersoo et al. 2020 Frontiers in Microbiology
#FungalTraits beta version 2020-03-01 (Põlme et al. Unpublished).
dat<-read.csv(file="...ft_downloaded_5oct2020.csv",stringsAsFactors=FALSE)

#see if anything in secondary lifestyle column, and add it in if is present
colnames(dat)[colnames(dat) %in% "Secondary_lifestyle"]<-"secondary_lifestyle"
colnames(dat)[colnames(dat) %in% "Genus"]<-"GENUS"

table(dat$secondary_lifestyle)

tax$primary_lifestyle<-NA
tax$secondary_lifestyle<-NA

for(x in 1:nrow(tax)){
	if(tax$Genus[x] %in% dat$GENUS){
		if(length(dat[dat$GENUS %in% tax$Genus[x],"primary_lifestyle"])==1){
			tax$primary_lifestyle[x]<-dat[dat$GENUS %in% tax$Genus[x],"primary_lifestyle"]
			if(length(dat[dat$GENUS %in% tax$Genus[x],"secondary_lifestyle"])==1){
				tax$secondary_lifestyle[x]<-dat[dat$GENUS %in% tax$Genus[x],"secondary_lifestyle"]
			}
		}
		#608 Brunneocorticium has a primary_lifestyle of "wood_saprotroph" "wood_saprotroph" and causes problems 
		#(says there are two answers)...take first
		if(length(dat[dat$GENUS %in% tax$Genus[x],"primary_lifestyle"])>1){
			if(all(dat[dat$GENUS %in% tax$Genus[x],"primary_lifestyle"]==c("wood_saprotroph","wood_saprotroph"))){
				tax$primary_lifestyle[x]<-"wood_saprotroph"
				if(length(dat[dat$GENUS %in% tax$Genus[x],"secondary_lifestyle"])==1){
					tax$secondary_lifestyle[x]<-dat[dat$GENUS %in% tax$Genus[x],"secondary_lifestyle"]
				}
			}	
		}
	}
}


#check if anything in secondary_lifestyle
table(tax$secondary_lifestyle)

#summarize what is in primary_lifestyle
table(tax$primary_lifestyle)

#############
#ladderize
tr<-ladderize(tr,right=FALSE)

#identify drops and drop them
tree.drops<-tr$tip.label[!tr$tip.label %in% tax$Taxon]
tr.clean<-drop.tip(tr,tree.drops)

trs.clean<-NULL
for(x in 1:length(trs)){
	trs[[x]]<-ladderize(trs[[x]],right=FALSE)
	trs.clean[[x]]<-drop.tip(trs[[x]],tree.drops)
}
class(trs.clean)<-"multiPhylo"
###############

#reduce tax to only those in tree
tax.clean<-tax[tax$Taxon %in% tr.clean$tip.label,]

#double check no tips are missing from tax
tr.clean$tip.label[!tr.clean$tip.label %in% tax.clean$Taxon]
trs.clean[[1]]$tip.label[!trs.clean[[1]]$tip.label %in% tax.clean$Taxon]


#double-check no tips are duplicated
tr.clean$tip.label[duplicated(tr.clean$tip.label)]
trs.clean[[1]]$tip.label[duplicated(trs.clean[[1]]$tip.label)]

#check which taxa are missing data
tax.clean[tax.clean$primary_lifestyle %in% "",]

#Tsuchiyaea_wingfieldii - Tremellaceae
######Laeticorticium_roseocarneum IS MISSING IN A LATER VERSION BUT NOT FUNGALTRAITS v USED HERE Laeticorticium_roseocarneum - Dendrocorticium roseocarneum looks like currently accepted name.

#get rid of taxa lacking any ecological info ("", "unspecified") and those that are "unspecificed_symbiotroph" since this could theoretically be an ECM?
nas<-tax.clean$Taxon[is.na(tax.clean$primary_lifestyle)]
nas
# [1] "Battarraea_laciniata"          "Pachylepirium_fulvidula"      
# [3] "Torrendia_pulchella"           "Laeticorticium_roseocarneum"  
# [5] "Pleurogala_panuoides"          "Pseudotrametes_gibbosa"       
# [7] "Pachytospora_tuberculosa"      "Pacychytospora_tuberculosa"   
# [9] "Physosporinus_sanguinolentus"  "Phaeoporus_ulmicola"          
#[11] "Phaeoporus_pruinosus"          "Pyrenogaster_pityophilus"     
#[13] "Afrocantharellus_platyphyllus" "Afrocantharellus_symoensii"   
#[15] "Guehomyces_pullulans"          "Asterotremella_humicola" 

abs<-tax.clean$Taxon[tax.clean$primary_lifestyle %in% ""]
abs
#[1] "Tsuchiyaea_wingfieldii"

uns<-tax.clean$Taxon[tax.clean$primary_lifestyle %in% "unspecified"]
uns
#[1] "Albomagister_subaustralis"     "Gloeotulasnella_cystidiophora"

uns.sym<-tax.clean$Taxon[tax.clean$primary_lifestyle %in% "unspecified_symbiotroph"]
uns.sym
#[1] "Phlebopus_marginatus"     "Phlebopus_beniensis"     
#[3] "Phlebopus_sudanicus"      "Phlebopus_portentosus"   
#[5] "Phlebopus_spongiosus"     "Boletinellus_rompelii"   
#[7] "Boletinellus_merulioides"

#merge
more.drops<-c(nas,abs,uns,uns.sym)

###################
#drop from test tree and all trees
tr.cleaner<-drop.tip(tr.clean,more.drops)
tr.cleaner

trs.cleaner<-NULL
for(x in 1:length(trs.clean)){
	trs.cleaner[[x]]<-drop.tip(trs.clean[[x]],more.drops)
}
class(trs.cleaner)<-"multiPhylo"
trs.cleaner[[1]]
###################

#drop from tax
tax.cleaner<-tax.clean[!tax.clean$Taxon %in% more.drops,]
nrow(tax.cleaner)

#double check no tips are missing from tax
tr.cleaner$tip.label[!tr.cleaner$tip.label %in% tax.cleaner$Taxon]
trs.cleaner[[1]]$tip.label[!trs.cleaner[[1]]$tip.label %in% tax.cleaner$Taxon]

#double check no tax are missing from tips
tax.cleaner$Taxon[!tax.cleaner$Taxon %in% tr.cleaner$tip.label]
tax.cleaner$Taxon[!tax.cleaner$Taxon %in% trs.cleaner[[1]]$tip.label]

#double check no taxa duplicated
tax.cleaner$Taxon[duplicated(tax.cleaner$Taxon)]

sort(unique(tax.cleaner$primary_lifestyle))

tax.cleaner$Ecto<-1
tax.cleaner$Ecto[tax.cleaner$primary_lifestyle %in% "ectomycorrhizal"]<-2
tax.cleaner.ecto<-tax.cleaner[,c("Taxon","Ecto"),]

#write file
write.csv(tax.cleaner.ecto,file="...Ecto_File.csv",row.names=FALSE)

#write trees
write.tree(trs.cleaner,file="...Ecto_Trees_Cleaned.tre")


#Modify data file a little
tax.cleaner.ecto<-read.csv(file="...Ecto_File.csv",stringsAsFactors=FALSE)
tax.cleaner.ecto$Ecto[tax.cleaner.ecto$Ecto==1]<-"Absent"
tax.cleaner.ecto$Ecto[tax.cleaner.ecto$Ecto==2]<-"Present"
colnames(tax.cleaner.ecto)<-c("Taxon","State")
write.csv(tax.cleaner.ecto,"...EctoStates_File.csv",row.names=FALSE)



#Stochastic Mapping
#Pi constrained to be non-ECM, Q="empirical"
require(phytools)
#read in trait data
tax<-read.csv(file="...EctoStates_File.csv")
dat<-tax$State
names(dat)<-tax$Taxon
head(dat)

#read in trees
trs.cleaner<-read.tree(file="...Ecto_Trees_Cleaned.tre")

#set prior on root node
pi<-setNames(c(1.0,0.0),c("Absent","Present"))
pi

for(z in 1:length(trs.cleaner)){
	cat(paste("Beginning tree ", z, "\n",sep=""))
	tr<-NULL
	tr<-trs.cleaner[[z]]
	tr.fit<-NULL
	tr.fit<-make.simmap(tree=tr,x=dat,model="ARD",pi=pi,Q="empirical",nsim=1000)
	save(tr.fit,dat,tr,file=paste("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_Tree_",z,".RData",sep=""))
	cat(paste(z, "_ARD finished\n",sep=""))
	tr.fit<-NULL
	tr.fit<-make.simmap(tree=tr,x=dat,model="ER",pi=pi,Q="empirical",nsim=1000)
	save(tr.fit,dat,tr,file=paste("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ER_Tree_",z,".RData",sep=""))
	cat(paste(z, "_ER finished\n",sep=""))
}


#make density maps for all and save
require(phytools)
for(z in 1:10){
	tr.fit<-NULL
	dat<-NULL
	tr<-NULL
	obj<-NULL
	load(paste("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_Tree_",z,".RData",sep=""))
	obj<-densityMap(tr.fit,states=levels(dat)[1:2],plot=FALSE)
	save(tr.fit,dat,tr,obj,file=paste("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_Tree_",z,".RData",sep=""))
}


#summarize all 
require(phytools)
for(z in 1:10){
	tr.fit<-NULL
	dat<-NULL
	tr<-NULL
	obj<-NULL
	load(paste("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_Tree_",z,".RData",sep=""))
	sums<-describe.simmap(tr.fit)
	save(tr.fit,dat,tr,obj,sums,file=paste("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_",z,".RData",sep=""))
}


#basic plots of individual trees with tip labels
require(phytools)

#tree numbers from Varga et al.
alt.nos<-c(8,12,16,53,88,165,170,188,216,222)

for(m in 1:10){
	tr.fit<-NULL
	dat<-NULL
	tr<-NULL
	obj<-NULL
	load(paste("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_",m,".RData",sep=""))
	obj<-setMap(obj,c("gray48","green4"))
	cols<-setNames(c("gray48","green4"),c("non-ECM","ECM"))
	obj$states<-c("non-ECM","ECM")
	pdf(paste("...tree",m,"_",alt.nos[m],"_tip.labels.pdf",sep=""))
	plot(obj,ftype="i",type="fan",lwd=0.05,fsize=0.05,colors=cols,legend=0)
	axisPhylo(cex=0.5)
	legend(x="bottomright",legend=names(cols),pt.cex=2.5,pch=16,col=cols,bty="n")
	title(main=paste("Tree ",m," (",alt.nos[m],")",sep=""),font.main=2,line=-0.5)
	dev.off()
}


################################
################################
################################
################################
################################
####LECANOROMYCETES ANY CYANOS
################################
#Data matrix from Nelsen et al. (2020 PNAS) downloaded from: https://github.com/mpnelsen/Lecanoromycetes_megaphylogeny/blob/master/HiSSE/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_anycyano_forcorhmm.csv
tax<-read.csv(file="...ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_anycyano_forcorhmm.csv",stringsAsFactors=FALSE)
colnames(tax)<-c("Taxon","State")
rownames(tax)<-tax$Taxon
tax$State[tax$State==0]<-"Absent"
tax$State[tax$State==1]<-"Present"
write.csv(tax,file="...Lecanoromycetes_anycyano.csv",row.names=TRUE)

#Topologies from Nelsen et al. (2020 PNAS) downloaded from: https://github.com/mpnelsen/Lecanoromycetes_megaphylogeny/tree/master/Trees
require(phytools)
#read in ML tree
trs.cleaner<-read.tree(file="...Lecanoromycetes_ML_Tree_Timescaled")

#read in boots
boots<-read.tree(file="...Lecanoromycetes_Timescaled_Bootstrap_Topologies_100.tre")

#randomly select 10 trees from boots
keeps<-sort(sample(1:100,10,FALSE))

#> keeps
# [1]  1 19 30 35 62 68 84 88 90 92
boots.mod<-NULL
for(x in 1:length(keeps)){
	print(keeps[x])
	if(x==1){
		boots.mod<-boots[[keeps[x]]]
	} else {
		boots.mod<-c(boots.mod,boots[[keeps[x]]])
		class(boots.mod)<-"multiPhylo"
	}
}

#now add ML tree to end as tree 11
boots.mod<-c(boots.mod,trs.cleaner)
class(boots.mod)<-"multiPhylo"

#save
write.tree(boots.mod,file="...Lecanoromycetes_10boots_1ML.tre")


#Stochastic Mapping
#Pi constrained to be non-ECM, Q="empirical"
require(phytools)
#read in trait data
tax<-read.csv(file="...Lecanoromycetes_anycyano.csv")
dat<-tax$State
names(dat)<-tax$Taxon
head(dat)

#read in ML tree
trs.cleaner<-read.tree(file="...Lecanoromycetes_10boots_1ML.tre")

#set prior on root node
pi<-setNames(c(1.0,0.0),c("Absent","Present"))
pi

for(z in 1:length(trs.cleaner)){
	cat(paste("Beginning tree ", z, "\n",sep=""))
	tr<-NULL
	tr<-trs.cleaner[[z]]
	tr.fit<-NULL
	tr.fit<-make.simmap(tree=tr,x=dat,model="ARD",pi=pi,Q="empirical",nsim=1000)
	save(tr.fit,dat,tr,file=paste("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_Tree_",z,".RData",sep=""))
	cat(paste(z, "_ARD finished\n",sep=""))
}


#make density maps for all and save
require(phytools)
for(z in 1:11){
	tr.fit<-NULL
	dat<-NULL
	tr<-NULL
	obj<-NULL
	load(paste("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_Tree_",z,".RData",sep=""))
	obj<-densityMap(tr.fit,states=levels(dat)[1:2],plot=FALSE)
	save(tr.fit,dat,tr,obj,file=paste("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_Tree_",z,".RData",sep=""))
}


#summarize all 
require(phytools)
for(z in 1:11){
	tr.fit<-NULL
	dat<-NULL
	tr<-NULL
	obj<-NULL
	load(paste("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_Tree_",z,".RData",sep=""))
	sums<-describe.simmap(tr.fit)
	save(tr.fit,dat,tr,obj,sums,file=paste("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_",z,".RData",sep=""))
}


#basic plots of individual trees with tip labels
require(phytools)

#tree numbers 
alt.nos<-c(1,19,30,35,62,68,84,88,90,92,"ML")

for(m in 1:11){
	tr.fit<-NULL
	dat<-NULL
	tr<-NULL
	obj<-NULL
	load(paste("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_",m,".RData",sep=""))
	obj<-setMap(obj,c("gray48","lightsteelblue2"))
	cols<-setNames(c("gray48","lightsteelblue2"),c("non-Cyanobacterial","Cyanobacterial"))
	obj$states<-c("non-ECM","ECM")
	pdf(paste("...tree",m,"_",alt.nos[m],"_tip.labels.pdf",sep=""))
	plot(obj,ftype="i",type="fan",lwd=0.05,fsize=0.05,colors=cols,legend=0)
	axisPhylo(cex=0.5)
	legend(x="bottomright",legend=names(cols),pt.cex=2.5,pch=16,col=cols,bty="n")
	title(main=paste("Tree ",m," (",alt.nos[m],")",sep=""),font.main=2,line=-0.5)
	dev.off()
}



#########################
#########################
######Jointly plot two trees
#########################
#plot first tree nicely
require(phytools)
require(circlize)

#ML Tree
load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_11.RData")
obj<-setMap(obj,c("gray48","darkgoldenrod1"))
cols<-setNames(c("gray48","darkgoldenrod1"),c("non-Cyanobacterial","Cyanobacterial"))
obj$states<-c("non-Cyanobacterial","Cyanobacterial")

#get age of root
tr.max<-max(nodeHeights(obj$tree))

#read timescale in
timescale_ics2020<-read.csv(file="...timescale_ics2020.csv")

#reduce to just periods of interest
timescale_ics2020<-timescale_ics2020[timescale_ics2020$Type %in% "Period",]
timescale_ics2020<-timescale_ics2020[timescale_ics2020$End <= tr.max,]

#make rgb color
for(x in 1:nrow(timescale_ics2020)){
	timescale_ics2020$RGB[x]<-rgb(timescale_ics2020$Col_R[x]/255,timescale_ics2020$Col_G[x]/255,timescale_ics2020$Col_B[x]/255,alpha=0.2)
}

pdf("...both.pdf",width=7,height=3.5)
par(mfrow=c(1,2))
plot(obj,ftype="off",type="fan",lwd=0.05,fsize=0.05,colors=cols,legend=0)
par(new = T)
circos.par(gap.degree=0,start.degree = 90,canvas.xlim=c(-1,1),canvas.ylim=c(-1,1))
circos.initialize("a", xlim = c(0, 1),sector.width=1) # 'a` just means there is one sector 
for(q in 1:nrow(timescale_ics2020)){
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(tr.max-timescale_ics2020$End[q])/tr.max,rou2=(tr.max-timescale_ics2020$Start[q])/tr.max,col=timescale_ics2020$RGB[q],border=NA)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(tr.max-100)/tr.max,col=NA,border="gray77",lwd=0.5,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(tr.max-200)/tr.max,col=NA,border="gray77",lwd=0.5,lty=3)
}
circos.clear()
par(new = T)
plot(obj,ftype="off",type="fan",lwd=0.25,fsize=0.05,colors=cols,legend=0,add=TRUE)
legend(x=-52,y=-15,legend=names(cols),pt.cex=1.5,pch=16,col=cols,box.lty=0,cex=0.75,bg="white")
for(q in 3:nrow(timescale_ics2020)){
	text(x=tr.max-rev(timescale_ics2020$Midpoint[q]),y=seq(from=-8,to=-9,by=-1/8)[q],labels=timescale_ics2020$Abbrev[q],cex=0.5,col="black")
}
title(main="Lecanoromycetes",cex.main=1,font.main=1,line=-0.4)




load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_4.RData")
obj<-setMap(obj,c("gray48","darkgoldenrod"))
cols<-setNames(c("gray48","darkgoldenrod"),c("non-ECM","ECM"))
obj$states<-c("non-ECM","ECM")

#get age of root
tr.max<-max(nodeHeights(obj$tree))

#read timescale in
timescale_ics2020<-read.csv(file="...timescale_ics2020.csv")

#reduce to just periods of interest
timescale_ics2020<-timescale_ics2020[timescale_ics2020$Type %in% "Period",]
timescale_ics2020<-timescale_ics2020[timescale_ics2020$End <= tr.max,]

#make rgb color
for(x in 1:nrow(timescale_ics2020)){
	timescale_ics2020$RGB[x]<-rgb(timescale_ics2020$Col_R[x]/255,timescale_ics2020$Col_G[x]/255,timescale_ics2020$Col_B[x]/255,alpha=0.2)
}

plot(obj,ftype="off",type="fan",lwd=0.05,fsize=0.05,colors=cols,legend=0)
par(new = T)
circos.par(gap.degree=0,start.degree = 90,canvas.xlim=c(-1,1),canvas.ylim=c(-1,1))
circos.initialize("a", xlim = c(0, 1),sector.width=1) # 'a` just means there is one sector 
for(q in 1:nrow(timescale_ics2020)){
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(tr.max-timescale_ics2020$End[q])/tr.max,rou2=(tr.max-timescale_ics2020$Start[q])/tr.max,col=timescale_ics2020$RGB[q],border=NA)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(tr.max-100)/tr.max,col=NA,border="gray77",lwd=0.5,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(tr.max-200)/tr.max,col=NA,border="gray77",lwd=0.5,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(tr.max-300)/tr.max,col=NA,border="gray77",lwd=0.5,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(tr.max-400)/tr.max,col=NA,border="gray77",lwd=0.5,lty=3)
}
circos.clear()
par(new = T)
plot(obj,ftype="off",type="fan",lwd=0.25,fsize=0.05,colors=cols,legend=0,add=TRUE)
legend(x=-120,y=-75,legend=names(cols),pt.cex=1.5,pch=16,col=cols,box.lty=0,cex=0.75,bg="white")
for(q in 4:nrow(timescale_ics2020)){
	text(x=tr.max-rev(timescale_ics2020$Midpoint[q]),y=seq(from=-8,to=-9-(1/8),by=-1/8)[q],labels=timescale_ics2020$Abbrev[q],cex=0.5,col="black")
}
title(main="Agaricomycotina",cex.main=1,font.main=1,line=-0.4)
dev.off()







#########################
#########################
####Jointly plot rates/changes in Lecanoromycetes ML tree and exemplar Agaricomycotina tree
#########################

require(ggplot2)
require(gridExtra)
require(phytools)

#read timescale in
timescale_ics2020<-read.csv(file="...timescale_ics2020.csv")

#reduce to just periods of interest
timescale_ics2020<-timescale_ics2020[timescale_ics2020$Type %in% "Period",]
timey<-timescale_ics2020[timescale_ics2020$Start<550,]
timey$Name<-c(NA,"Ng","Pg","K","J","T","P","C","D","S","O","Ca")

#make rgb color
for(x in 1:nrow(timey)){
	timey$RGB[x]<-rgb(timey$Col_R[x]/255,timey$Col_G[x]/255,timey$Col_B[x]/255)
}

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_11.RData")
ctt.to.plot.cyano<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")

load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_4.RData")
ctt.to.plot.ecm<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")

ch.cyano<-as.data.frame(cbind(max(ctt.to.plot.cyano$segments)-as.vector(t(ctt.to.plot.cyano$segments)),as.vector(rbind(ctt.to.plot.cyano$nchanges,ctt.to.plot.cyano$nchanges))))
colnames(ch.cyano)<-c("Time","Change")
ch.cyano$Kind<-"Cyanobacteria"

ch.r.cyano<-as.data.frame(cbind(max(ctt.to.plot.cyano$segments)-as.vector(t(ctt.to.plot.cyano$segments)),as.vector(rbind(ctt.to.plot.cyano$nchanges/ctt.to.plot.cyano$edge.length,ctt.to.plot.cyano$nchanges/ctt.to.plot.cyano$edge.length))))
colnames(ch.r.cyano)<-c("Time","Change")
ch.r.cyano$Kind<-"Cyanobacteria"

ch.ecm<-as.data.frame(cbind(max(ctt.to.plot.ecm$segments)-as.vector(t(ctt.to.plot.ecm$segments)),as.vector(rbind(ctt.to.plot.ecm$nchanges,ctt.to.plot.ecm$nchanges))))
colnames(ch.ecm)<-c("Time","Change")
ch.ecm$Kind<-"ECM"

ch.r.ecm<-as.data.frame(cbind(max(ctt.to.plot.ecm$segments)-as.vector(t(ctt.to.plot.ecm$segments)),as.vector(rbind(ctt.to.plot.ecm$nchanges/ctt.to.plot.ecm$edge.length,ctt.to.plot.ecm$nchanges/ctt.to.plot.ecm$edge.length))))
colnames(ch.r.ecm)<-c("Time","Change")
ch.r.ecm$Kind<-"ECM"

ch<-rbind(ch.cyano,ch.ecm)
ch.r<-rbind(ch.r.cyano,ch.r.ecm)

ch.comb<-ch
colnames(ch.comb)<-c("Time","ChangeNo","Kind")
ch.comb<-cbind(ch.comb,ch.r$Change)
colnames(ch.comb)<-c("Time","ChangeNo","Kind","ChangeRate")

mycols<-c("darkgoldenrod","darkgoldenrod1")

timeplot<-ggplot(ch.comb,aes(x=Time,ymin=-0.5))+scale_x_reverse(expand=c(0,0),limits=c(500,0), breaks=c(500,0))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x.top=element_blank(), axis.title.x.bottom=element_blank(), axis.line.x=element_blank(), panel.background=element_blank(), axis.text.x=element_text(color="black"), axis.text.y=element_text(color="black"),axis.ticks.x = element_blank(), axis.ticks.y = element_line(color="black")) + annotate(geom='segment',y=-Inf,yend=Inf,x=Inf,xend=Inf) + annotate(geom='segment',y=-Inf,yend=Inf,x=-Inf,xend=-Inf) +annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=-0.25,fill=timey$RGB[1],color="black",alpha=1)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=-0.25,fill=timey$RGB[2],color="black",alpha=1)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=-0.25,fill=timey$RGB[3],color="black",alpha=1)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=-0.25,fill=timey$RGB[4],color="black",alpha=1)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=-0.25,fill=timey$RGB[5],color="black",alpha=1)+annotate("rect",xmin=timey$End[6],xmax=timey$Start[6],ymin=-Inf,ymax=-0.25,fill=timey$RGB[6],color="black",alpha=1)+annotate("rect",xmin=timey$End[7],xmax=timey$Start[7],ymin=-Inf,ymax=-0.25,fill=timey$RGB[7],color="black",alpha=1)+annotate("rect",xmin=timey$End[8],xmax=timey$Start[8],ymin=-Inf,ymax=-0.25,fill=timey$RGB[8],color="black",alpha=1)+annotate("rect",xmin=timey$End[9],xmax=timey$Start[9],ymin=-Inf,ymax=-0.25,fill=timey$RGB[9],color="black",alpha=1)+annotate("rect",xmin=timey$End[10],xmax=timey$Start[10],ymin=-Inf,ymax=-0.25,fill=timey$RGB[10],color="black",alpha=1)+annotate("rect",xmin=timey$End[11],xmax=timey$Start[11],ymin=-Inf,ymax=-0.25,fill=timey$RGB[11],color="black",alpha=1)+annotate("rect",xmin=timey$End[12],xmax=500,ymin=-Inf,ymax=-0.25,fill=timey$RGB[12],color="black",alpha=1)+annotate("text",x=timey$Midpoint[11],y=-0.475,label=timey$Name[11],size=3,colour="black")+annotate("text",x=timey$Midpoint[10],y=-0.475,label=timey$Name[10],size=3,colour="black")+annotate("text",x=timey$Midpoint[9],y=-0.475,label=timey$Name[9],size=3,colour="black")+annotate("text",x=timey$Midpoint[8],y=-0.475,label=timey$Name[8],size=3,colour="black")+annotate("text",x=timey$Midpoint[7],y=-0.475,label=timey$Name[7],size=3,colour="black")+annotate("text",x=timey$Midpoint[6],y=-0.475,label=timey$Name[6],size=3,colour="black")+annotate("text",x=timey$Midpoint[5],y=-0.475,label=timey$Name[5],size=3,colour="black")+annotate("text",x=timey$Midpoint[4],y=-0.475,label=timey$Name[4],size=3,colour="black")+annotate("text",x=timey$Midpoint[3],y=-0.475,label=timey$Name[3],size=3,colour="black")+annotate("text",x=timey$Midpoint[2],y=-0.475,label=timey$Name[2],size=3,colour="black")+annotate("text",x=1.3,y=-0.475,label=timey$Name[1],size=3,colour="black") + annotate(geom='segment',y=3.2,yend=3.2,x=450,xend=425,color="darkgoldenrod",linetype=1) + annotate(geom='segment',y=2.8,yend=2.8,x=450,xend=425,color="darkgoldenrod1",linetype=1)+ annotate(geom='segment',y=1.8,yend=1.8,x=450,xend=425,color="darkgoldenrod",linetype=3) + annotate(geom='segment',y=1.4,yend=1.4,x=450,xend=425,color="darkgoldenrod1",linetype=3) + annotate(geom="text", x=450, y=3.6, label="Number of Gains",color="black",hjust=0,size=2.5,fontface=2) + annotate(geom="text", x=420, y=3.2, label="ECM Associations (Agaricomycotina)",color="black",hjust=0,size=2.5) + annotate(geom="text", x=420, y=2.8, label="Cyanobacterial Associations (Lecanoromycetes)",color="black",hjust=0,size=2.5) + annotate(geom="text", x=450, y=2.2, label="Gain Rate",color="black",hjust=0,size=2.5,fontface=2)+annotate(geom="text", x=420, y=1.8, label="ECM Associations (Agaricomycotina)",color="black",hjust=0,size=2.5) + annotate(geom="text", x=420, y=1.4, label="Cyanobacterial Associations (Lecanoromycetes)",color="black",hjust=0,size=2.5)
timeplot.comb<-timeplot+geom_step(data=subset(ch.comb[ch.comb$Kind=="ECM",]),linetype=1,alpha=1,color=mycols[1],aes(y=ChangeNo))
timeplot.combb<-timeplot.comb+geom_step(data=subset(ch.comb[ch.comb$Kind=="Cyanobacteria",]),linetype=1,alpha=1,color=mycols[2],aes(y=ChangeNo))
timeplot.combbb<-timeplot.combb+geom_step(data=subset(ch.comb[ch.comb$Kind=="ECM",]),linetype=3,alpha=1,color=mycols[1],aes(y=ChangeRate*666))
timeplotc<-timeplot.combbb+geom_step(data=subset(ch.comb[ch.comb$Kind=="Cyanobacteria",]),linetype=3,alpha=1,color=mycols[2],aes(y=ChangeRate*666))
timeplotd<-timeplotc+scale_y_continuous(name="Number", sec.axis = sec_axis(~./666, name = "Rate"))
timeplote<-timeplotd
pdf(file="... testmixed_CTT_10my.pdf",width=7,height=2)
timeplote
dev.off()


####Jointly PLOT CTT's TOGETHER ON ONE PAGE NICELY

require(ggplot2)
require(gridExtra)
require(phytools)
require(scico)

#Get files
load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_1.RData")
ctt.to.plot.mod.lec.1<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_2.RData")
ctt.to.plot.mod.lec.2<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_3.RData")
ctt.to.plot.mod.lec.3<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_4.RData")
ctt.to.plot.mod.lec.4<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_5.RData")
ctt.to.plot.mod.lec.5<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_6.RData")
ctt.to.plot.mod.lec.6<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_7.RData")
ctt.to.plot.mod.lec.7<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_8.RData")
ctt.to.plot.mod.lec.8<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_9.RData")
ctt.to.plot.mod.lec.9<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_10.RData")
ctt.to.plot.mod.lec.10<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonCyano_empirical_Lecanoromycetes_Cyano_ARD_wDensityMap_summary_Tree_11.RData")
ctt.to.plot.mod.lec.11<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL


#put into a data frame for plotting
ch.lec.1<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.1$segments)-as.vector(t(ctt.to.plot.mod.lec.1$segments)),as.vector(rbind(ctt.to.plot.mod.lec.1$nchanges,ctt.to.plot.mod.lec.1$nchanges))))
colnames(ch.lec.1)<-c("Time","Change")
ch.lec.1$Kind<-"1"

ch.r.lec.1<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.1$segments)-as.vector(t(ctt.to.plot.mod.lec.1$segments)),as.vector(rbind(ctt.to.plot.mod.lec.1$nchanges/ctt.to.plot.mod.lec.1$edge.length,ctt.to.plot.mod.lec.1$nchanges/ctt.to.plot.mod.lec.1$edge.length))))
colnames(ch.r.lec.1)<-c("Time","Change")
ch.r.lec.1$Kind<-"1"

ch.lec.2<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.2$segments)-as.vector(t(ctt.to.plot.mod.lec.2$segments)),as.vector(rbind(ctt.to.plot.mod.lec.2$nchanges,ctt.to.plot.mod.lec.2$nchanges))))
colnames(ch.lec.2)<-c("Time","Change")
ch.lec.2$Kind<-"19"

ch.r.lec.2<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.2$segments)-as.vector(t(ctt.to.plot.mod.lec.2$segments)),as.vector(rbind(ctt.to.plot.mod.lec.2$nchanges/ctt.to.plot.mod.lec.2$edge.length,ctt.to.plot.mod.lec.2$nchanges/ctt.to.plot.mod.lec.2$edge.length))))
colnames(ch.r.lec.2)<-c("Time","Change")
ch.r.lec.2$Kind<-"19"

ch.lec.3<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.3$segments)-as.vector(t(ctt.to.plot.mod.lec.3$segments)),as.vector(rbind(ctt.to.plot.mod.lec.3$nchanges,ctt.to.plot.mod.lec.3$nchanges))))
colnames(ch.lec.3)<-c("Time","Change")
ch.lec.3$Kind<-"30"

ch.r.lec.3<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.3$segments)-as.vector(t(ctt.to.plot.mod.lec.3$segments)),as.vector(rbind(ctt.to.plot.mod.lec.3$nchanges/ctt.to.plot.mod.lec.3$edge.length,ctt.to.plot.mod.lec.3$nchanges/ctt.to.plot.mod.lec.3$edge.length))))
colnames(ch.r.lec.3)<-c("Time","Change")
ch.r.lec.3$Kind<-"30"

ch.lec.4<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.4$segments)-as.vector(t(ctt.to.plot.mod.lec.4$segments)),as.vector(rbind(ctt.to.plot.mod.lec.4$nchanges,ctt.to.plot.mod.lec.4$nchanges))))
colnames(ch.lec.4)<-c("Time","Change")
ch.lec.4$Kind<-"35"

ch.r.lec.4<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.4$segments)-as.vector(t(ctt.to.plot.mod.lec.4$segments)),as.vector(rbind(ctt.to.plot.mod.lec.4$nchanges/ctt.to.plot.mod.lec.4$edge.length,ctt.to.plot.mod.lec.4$nchanges/ctt.to.plot.mod.lec.4$edge.length))))
colnames(ch.r.lec.4)<-c("Time","Change")
ch.r.lec.4$Kind<-"35"

ch.lec.5<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.5$segments)-as.vector(t(ctt.to.plot.mod.lec.5$segments)),as.vector(rbind(ctt.to.plot.mod.lec.5$nchanges,ctt.to.plot.mod.lec.5$nchanges))))
colnames(ch.lec.5)<-c("Time","Change")
ch.lec.5$Kind<-"62"

ch.r.lec.5<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.5$segments)-as.vector(t(ctt.to.plot.mod.lec.5$segments)),as.vector(rbind(ctt.to.plot.mod.lec.5$nchanges/ctt.to.plot.mod.lec.5$edge.length,ctt.to.plot.mod.lec.5$nchanges/ctt.to.plot.mod.lec.5$edge.length))))
colnames(ch.r.lec.5)<-c("Time","Change")
ch.r.lec.5$Kind<-"62"

ch.lec.6<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.6$segments)-as.vector(t(ctt.to.plot.mod.lec.6$segments)),as.vector(rbind(ctt.to.plot.mod.lec.6$nchanges,ctt.to.plot.mod.lec.6$nchanges))))
colnames(ch.lec.6)<-c("Time","Change")
ch.lec.6$Kind<-"68"

ch.r.lec.6<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.6$segments)-as.vector(t(ctt.to.plot.mod.lec.6$segments)),as.vector(rbind(ctt.to.plot.mod.lec.6$nchanges/ctt.to.plot.mod.lec.6$edge.length,ctt.to.plot.mod.lec.6$nchanges/ctt.to.plot.mod.lec.6$edge.length))))
colnames(ch.r.lec.6)<-c("Time","Change")
ch.r.lec.6$Kind<-"68"

ch.lec.7<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.7$segments)-as.vector(t(ctt.to.plot.mod.lec.7$segments)),as.vector(rbind(ctt.to.plot.mod.lec.7$nchanges,ctt.to.plot.mod.lec.7$nchanges))))
colnames(ch.lec.7)<-c("Time","Change")
ch.lec.7$Kind<-"84"

ch.r.lec.7<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.7$segments)-as.vector(t(ctt.to.plot.mod.lec.7$segments)),as.vector(rbind(ctt.to.plot.mod.lec.7$nchanges/ctt.to.plot.mod.lec.7$edge.length,ctt.to.plot.mod.lec.7$nchanges/ctt.to.plot.mod.lec.7$edge.length))))
colnames(ch.r.lec.7)<-c("Time","Change")
ch.r.lec.7$Kind<-"84"

ch.lec.8<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.8$segments)-as.vector(t(ctt.to.plot.mod.lec.8$segments)),as.vector(rbind(ctt.to.plot.mod.lec.8$nchanges,ctt.to.plot.mod.lec.8$nchanges))))
colnames(ch.lec.8)<-c("Time","Change")
ch.lec.8$Kind<-"88"

ch.r.lec.8<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.8$segments)-as.vector(t(ctt.to.plot.mod.lec.8$segments)),as.vector(rbind(ctt.to.plot.mod.lec.8$nchanges/ctt.to.plot.mod.lec.8$edge.length,ctt.to.plot.mod.lec.8$nchanges/ctt.to.plot.mod.lec.8$edge.length))))
colnames(ch.r.lec.8)<-c("Time","Change")
ch.r.lec.8$Kind<-"88"

ch.lec.9<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.9$segments)-as.vector(t(ctt.to.plot.mod.lec.9$segments)),as.vector(rbind(ctt.to.plot.mod.lec.9$nchanges,ctt.to.plot.mod.lec.9$nchanges))))
colnames(ch.lec.9)<-c("Time","Change")
ch.lec.9$Kind<-"90"

ch.r.lec.9<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.9$segments)-as.vector(t(ctt.to.plot.mod.lec.9$segments)),as.vector(rbind(ctt.to.plot.mod.lec.9$nchanges/ctt.to.plot.mod.lec.9$edge.length,ctt.to.plot.mod.lec.9$nchanges/ctt.to.plot.mod.lec.9$edge.length))))
colnames(ch.r.lec.9)<-c("Time","Change")
ch.r.lec.9$Kind<-"90"

ch.lec.10<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.10$segments)-as.vector(t(ctt.to.plot.mod.lec.10$segments)),as.vector(rbind(ctt.to.plot.mod.lec.10$nchanges,ctt.to.plot.mod.lec.10$nchanges))))
colnames(ch.lec.10)<-c("Time","Change")
ch.lec.10$Kind<-"92"

ch.r.lec.10<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.10$segments)-as.vector(t(ctt.to.plot.mod.lec.10$segments)),as.vector(rbind(ctt.to.plot.mod.lec.10$nchanges/ctt.to.plot.mod.lec.10$edge.length,ctt.to.plot.mod.lec.10$nchanges/ctt.to.plot.mod.lec.10$edge.length))))
colnames(ch.r.lec.10)<-c("Time","Change")
ch.r.lec.10$Kind<-"92"

ch.lec.11<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.11$segments)-as.vector(t(ctt.to.plot.mod.lec.11$segments)),as.vector(rbind(ctt.to.plot.mod.lec.11$nchanges,ctt.to.plot.mod.lec.11$nchanges))))
colnames(ch.lec.11)<-c("Time","Change")
ch.lec.11$Kind<-"ML"

ch.r.lec.11<-as.data.frame(cbind(max(ctt.to.plot.mod.lec.11$segments)-as.vector(t(ctt.to.plot.mod.lec.11$segments)),as.vector(rbind(ctt.to.plot.mod.lec.11$nchanges/ctt.to.plot.mod.lec.11$edge.length,ctt.to.plot.mod.lec.11$nchanges/ctt.to.plot.mod.lec.11$edge.length))))
colnames(ch.r.lec.11)<-c("Time","Change")
ch.r.lec.11$Kind<-"ML"


ch.lec<-rbind(ch.lec.1,ch.lec.2,ch.lec.3,ch.lec.4,ch.lec.5,ch.lec.6,ch.lec.7,ch.lec.8,ch.lec.9,ch.lec.10,ch.lec.11)
ch.lec.r<-rbind(ch.r.lec.1,ch.r.lec.2,ch.r.lec.3,ch.r.lec.4,ch.r.lec.5,ch.r.lec.6,ch.r.lec.7,ch.r.lec.8,ch.r.lec.9,ch.r.lec.10,ch.r.lec.11)

ch.comb.lec<-ch.lec
colnames(ch.comb.lec)<-c("Time","ChangeNo","Kind")
ch.comb.lec<-cbind(ch.comb.lec,ch.lec.r$Change)
colnames(ch.comb.lec)<-c("Time","ChangeNo","Kind","ChangeRate")

mycols.lec<-c(scico(10, palette = "nuuk"),"black")

#read timescale in
timescale_ics2020<-read.csv(file="...timescale_ics2020.csv")

#reduce to just periods of interest
timescale_ics2020<-timescale_ics2020[timescale_ics2020$Type %in% "Period",]
timey<-timescale_ics2020[timescale_ics2020$Start<550,]
timey$Name<-c(NA,"Ng","Pg","K","J","T","P","C","D","S","O","Ca")

#make rgb color
for(x in 1:nrow(timey)){
	timey$RGB[x]<-rgb(timey$Col_R[x]/255,timey$Col_G[x]/255,timey$Col_B[x]/255)
}


timeplot.lec<-ggplot(ch.comb.lec,aes(x=Time,ymin=-0.5))+scale_x_reverse(expand=c(0,0),limits=c(500,0), breaks=c(500,0))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x.top=element_blank(), axis.title.x.bottom=element_blank(), axis.line.x=element_blank(), panel.background=element_blank(), axis.text.x=element_text(color="black"), axis.text.y=element_text(color="black"),axis.ticks.x = element_blank(), axis.ticks.y = element_line(color="black")) + annotate(geom='segment',y=-Inf,yend=Inf,x=Inf,xend=Inf) + annotate(geom='segment',y=-Inf,yend=Inf,x=-Inf,xend=-Inf) +annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=-0.25,fill=timey$RGB[1],color="black",alpha=1)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=-0.25,fill=timey$RGB[2],color="black",alpha=1)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=-0.25,fill=timey$RGB[3],color="black",alpha=1)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=-0.25,fill=timey$RGB[4],color="black",alpha=1)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=-0.25,fill=timey$RGB[5],color="black",alpha=1)+annotate("rect",xmin=timey$End[6],xmax=timey$Start[6],ymin=-Inf,ymax=-0.25,fill=timey$RGB[6],color="black",alpha=1)+annotate("rect",xmin=timey$End[7],xmax=timey$Start[7],ymin=-Inf,ymax=-0.25,fill=timey$RGB[7],color="black",alpha=1)+annotate("rect",xmin=timey$End[8],xmax=timey$Start[8],ymin=-Inf,ymax=-0.25,fill=timey$RGB[8],color="black",alpha=1)+annotate("rect",xmin=timey$End[9],xmax=timey$Start[9],ymin=-Inf,ymax=-0.25,fill=timey$RGB[9],color="black",alpha=1)+annotate("rect",xmin=timey$End[10],xmax=timey$Start[10],ymin=-Inf,ymax=-0.25,fill=timey$RGB[10],color="black",alpha=1)+annotate("rect",xmin=timey$End[11],xmax=timey$Start[11],ymin=-Inf,ymax=-0.25,fill=timey$RGB[11],color="black",alpha=1)+annotate("rect",xmin=timey$End[12],xmax=500,ymin=-Inf,ymax=-0.25,fill=timey$RGB[12],color="black",alpha=1)+annotate("text",x=timey$Midpoint[11],y=-0.475,label=timey$Name[11],size=3,colour="black")+annotate("text",x=timey$Midpoint[10],y=-0.475,label=timey$Name[10],size=3,colour="black")+annotate("text",x=timey$Midpoint[9],y=-0.475,label=timey$Name[9],size=3,colour="black")+annotate("text",x=timey$Midpoint[8],y=-0.475,label=timey$Name[8],size=3,colour="black")+annotate("text",x=timey$Midpoint[7],y=-0.475,label=timey$Name[7],size=3,colour="black")+annotate("text",x=timey$Midpoint[6],y=-0.475,label=timey$Name[6],size=3,colour="black")+annotate("text",x=timey$Midpoint[5],y=-0.475,label=timey$Name[5],size=3,colour="black")+annotate("text",x=timey$Midpoint[4],y=-0.475,label=timey$Name[4],size=3,colour="black")+annotate("text",x=timey$Midpoint[3],y=-0.475,label=timey$Name[3],size=3,colour="black")+annotate("text",x=timey$Midpoint[2],y=-0.475,label=timey$Name[2],size=3,colour="black")+annotate("text",x=1.3,y=-0.475,label=timey$Name[1],size=3,colour="black") 
timeplot.lec<-timeplot.lec+ annotate(geom='segment',y=3.2,yend=3.2,x=475,xend=450,color=mycols.lec[1],linetype=1) + annotate(geom='segment',y=2.9,yend=2.9,x=475,xend=450,color=mycols.lec[2],linetype=1)+ annotate(geom='segment',y=2.6,yend=2.6,x=475,xend=450,color=mycols.lec[3],linetype=1)+ annotate(geom='segment',y=2.3,yend=2.3,x=475,xend=450,color=mycols.lec[4],linetype=1)+ annotate(geom='segment',y=2.0,yend=2.0,x=475,xend=450,color=mycols.lec[5],linetype=1)+ annotate(geom='segment',y=1.7,yend=1.7,x=475,xend=450,color=mycols.lec[6],linetype=1)+ annotate(geom='segment',y=1.4,yend=1.4,x=475,xend=450,color=mycols.lec[7],linetype=1)+ annotate(geom='segment',y=1.1,yend=1.1,x=475,xend=450,color=mycols.lec[8],linetype=1)+ annotate(geom='segment',y=0.8,yend=0.8,x=475,xend=450,color=mycols.lec[9],linetype=1)+ annotate(geom='segment',y=0.5,yend=0.5,x=475,xend=450,color=mycols.lec[10],linetype=1)+ annotate(geom='segment',y=0.2,yend=0.2,x=475,xend=450,color=mycols.lec[11],linetype=1) + annotate(geom='segment',y=3.2,yend=3.2,x=440,xend=415,color=mycols.lec[1],linetype=3) + annotate(geom='segment',y=2.9,yend=2.9,x=440,xend=415,color=mycols.lec[2],linetype=3)+ annotate(geom='segment',y=2.6,yend=2.6,x=440,xend=415,color=mycols.lec[3],linetype=3)+ annotate(geom='segment',y=2.3,yend=2.3,x=440,xend=415,color=mycols.lec[4],linetype=3)+ annotate(geom='segment',y=2.0,yend=2.0,x=440,xend=415,color=mycols.lec[5],linetype=3)+ annotate(geom='segment',y=1.7,yend=1.7,x=440,xend=415,color=mycols.lec[6],linetype=3)+ annotate(geom='segment',y=1.4,yend=1.4,x=440,xend=415,color=mycols.lec[7],linetype=3)+ annotate(geom='segment',y=1.1,yend=1.1,x=440,xend=415,color=mycols.lec[8],linetype=3)+ annotate(geom='segment',y=0.8,yend=0.8,x=440,xend=415,color=mycols.lec[9],linetype=3)+ annotate(geom='segment',y=0.5,yend=0.5,x=440,xend=415,color=mycols.lec[10],linetype=3)+ annotate(geom='segment',y=0.2,yend=0.2,x=440,xend=415,color=mycols.lec[11],linetype=3) + annotate(geom="text", x=405, y=3.2, label="1",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=2.9, label="19",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=2.6, label="30",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=2.3, label="35",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=2.0, label="62",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=1.7, label="68",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=1.4, label="84",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=1.1, label="88",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=0.8, label="90",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=0.5, label="92",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=0.2, label="ML",color="black",hjust=0,size=2.5) + annotate(geom="text", x=475, y=3.6, label="Number",color="black",hjust=0,size=2.5,fontface=2)+ annotate(geom="text", x=440, y=3.6, label="Rate",color="black",hjust=0,size=2.5,fontface=2)+ annotate(geom="text", x=405, y=3.6, label="Tree",color="black",hjust=0,size=2.5,fontface=2)+ annotate(geom="text", x=410, y=4.5, label="Evolution of cyanobacterial associations in Lecanoromycetes",color="black",hjust=0,size=3.5,fontface=2)
timeplot.lec.comb<-timeplot.lec+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="1",]),linetype=1,alpha=1,color=mycols.lec[1],aes(y=ChangeNo))
timeplot.lec.combb<-timeplot.lec.comb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="19",]),linetype=1,alpha=1,color=mycols.lec[2],aes(y=ChangeNo))
timeplot.lec.combb<-timeplot.lec.combb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="30",]),linetype=1,alpha=1,color=mycols.lec[3],aes(y=ChangeNo))
timeplot.lec.combb<-timeplot.lec.combb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="35",]),linetype=1,alpha=1,color=mycols.lec[4],aes(y=ChangeNo))
timeplot.lec.combb<-timeplot.lec.combb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="62",]),linetype=1,alpha=1,color=mycols.lec[5],aes(y=ChangeNo))
timeplot.lec.combb<-timeplot.lec.combb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="68",]),linetype=1,alpha=1,color=mycols.lec[6],aes(y=ChangeNo))
timeplot.lec.combb<-timeplot.lec.combb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="84",]),linetype=1,alpha=1,color=mycols.lec[7],aes(y=ChangeNo))
timeplot.lec.combb<-timeplot.lec.combb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="88",]),linetype=1,alpha=1,color=mycols.lec[8],aes(y=ChangeNo))
timeplot.lec.combb<-timeplot.lec.combb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="90",]),linetype=1,alpha=1,color=mycols.lec[9],aes(y=ChangeNo))
timeplot.lec.combb<-timeplot.lec.combb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="92",]),linetype=1,alpha=1,color=mycols.lec[10],aes(y=ChangeNo))
timeplot.lec.combb<-timeplot.lec.combb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="ML",]),linetype=1,alpha=1,color=mycols.lec[11],aes(y=ChangeNo))

timeplot.lec.combbb<-timeplot.lec.combb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="1",]),linetype=3,alpha=1,color=mycols.lec[1],aes(y=ChangeRate*666))
timeplot.lecc<-timeplot.lec.combbb+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="19",]),linetype=3,alpha=1,color=mycols.lec[2],aes(y=ChangeRate*666))
timeplot.lecc<-timeplot.lecc+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="30",]),linetype=3,alpha=1,color=mycols.lec[3],aes(y=ChangeRate*666))
timeplot.lecc<-timeplot.lecc+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="35",]),linetype=3,alpha=1,color=mycols.lec[4],aes(y=ChangeRate*666))
timeplot.lecc<-timeplot.lecc+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="62",]),linetype=3,alpha=1,color=mycols.lec[5],aes(y=ChangeRate*666))
timeplot.lecc<-timeplot.lecc+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="68",]),linetype=3,alpha=1,color=mycols.lec[6],aes(y=ChangeRate*666))
timeplot.lecc<-timeplot.lecc+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="84",]),linetype=3,alpha=1,color=mycols.lec[7],aes(y=ChangeRate*666))
timeplot.lecc<-timeplot.lecc+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="88",]),linetype=3,alpha=1,color=mycols.lec[8],aes(y=ChangeRate*666))
timeplot.lecc<-timeplot.lecc+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="90",]),linetype=3,alpha=1,color=mycols.lec[9],aes(y=ChangeRate*666))
timeplot.lecc<-timeplot.lecc+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="92",]),linetype=3,alpha=1,color=mycols.lec[10],aes(y=ChangeRate*666))
timeplot.lecc<-timeplot.lecc+geom_step(data=subset(ch.comb.lec[ch.comb.lec$Kind=="ML",]),linetype=3,alpha=1,color=mycols.lec[11],aes(y=ChangeRate*666))
timeplot.lecd<-timeplot.lecc+scale_y_continuous(name="Mean number of gains /\n5 my time bin", sec.axis = sec_axis(~./666, name = "Mean number of gains /\n total branch length / 5 my time bin"))+theme(axis.title.y = element_text(size=8))
timeplot.lece<-timeplot.lecd
#pdf(file="~/lecanoromycetes_cyano_stochmap/bsrt_ml.pdf",width=7,height=2)
#timeplot.lece
#dev.off()


#Get files
load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_1.RData")
ctt.to.plot.mod.agarico.1<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_2.RData")
ctt.to.plot.mod.agarico.2<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_3.RData")
ctt.to.plot.mod.agarico.3<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_4.RData")
ctt.to.plot.mod.agarico.4<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_5.RData")
ctt.to.plot.mod.agarico.5<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_6.RData")
ctt.to.plot.mod.agarico.6<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_7.RData")
ctt.to.plot.mod.agarico.7<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_8.RData")
ctt.to.plot.mod.agarico.8<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_9.RData")
ctt.to.plot.mod.agarico.9<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL

load("...SIMMAP_1K_RootNonECM_empirical_Agaricomycotina_ECM_ARD_wDensityMap_summary_Tree_10.RData")
ctt.to.plot.mod.agarico.10<-ctt.mod2(tr.fit,timebin=5,direction="Absent->Present")
tr.fit<-NULL



#put into a data frame for plotting
ch.agarico.1<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.1$segments)-as.vector(t(ctt.to.plot.mod.agarico.1$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.1$nchanges,ctt.to.plot.mod.agarico.1$nchanges))))
colnames(ch.agarico.1)<-c("Time","Change")
ch.agarico.1$Kind<-"Tree 8"

ch.r.agarico.1<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.1$segments)-as.vector(t(ctt.to.plot.mod.agarico.1$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.1$nchanges/ctt.to.plot.mod.agarico.1$edge.length,ctt.to.plot.mod.agarico.1$nchanges/ctt.to.plot.mod.agarico.1$edge.length))))
colnames(ch.r.agarico.1)<-c("Time","Change")
ch.r.agarico.1$Kind<-"Tree 8"

ch.agarico.2<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.2$segments)-as.vector(t(ctt.to.plot.mod.agarico.2$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.2$nchanges,ctt.to.plot.mod.agarico.2$nchanges))))
colnames(ch.agarico.2)<-c("Time","Change")
ch.agarico.2$Kind<-"Tree 12"

ch.r.agarico.2<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.2$segments)-as.vector(t(ctt.to.plot.mod.agarico.2$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.2$nchanges/ctt.to.plot.mod.agarico.2$edge.length,ctt.to.plot.mod.agarico.2$nchanges/ctt.to.plot.mod.agarico.2$edge.length))))
colnames(ch.r.agarico.2)<-c("Time","Change")
ch.r.agarico.2$Kind<-"Tree 12"

ch.agarico.3<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.3$segments)-as.vector(t(ctt.to.plot.mod.agarico.3$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.3$nchanges,ctt.to.plot.mod.agarico.3$nchanges))))
colnames(ch.agarico.3)<-c("Time","Change")
ch.agarico.3$Kind<-"Tree 16"

ch.r.agarico.3<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.3$segments)-as.vector(t(ctt.to.plot.mod.agarico.3$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.3$nchanges/ctt.to.plot.mod.agarico.3$edge.length,ctt.to.plot.mod.agarico.3$nchanges/ctt.to.plot.mod.agarico.3$edge.length))))
colnames(ch.r.agarico.3)<-c("Time","Change")
ch.r.agarico.3$Kind<-"Tree 16"

ch.agarico.4<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.4$segments)-as.vector(t(ctt.to.plot.mod.agarico.4$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.4$nchanges,ctt.to.plot.mod.agarico.4$nchanges))))
colnames(ch.agarico.4)<-c("Time","Change")
ch.agarico.4$Kind<-"Tree 53"

ch.r.agarico.4<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.4$segments)-as.vector(t(ctt.to.plot.mod.agarico.4$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.4$nchanges/ctt.to.plot.mod.agarico.4$edge.length,ctt.to.plot.mod.agarico.4$nchanges/ctt.to.plot.mod.agarico.4$edge.length))))
colnames(ch.r.agarico.4)<-c("Time","Change")
ch.r.agarico.4$Kind<-"Tree 53"

ch.agarico.5<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.5$segments)-as.vector(t(ctt.to.plot.mod.agarico.5$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.5$nchanges,ctt.to.plot.mod.agarico.5$nchanges))))
colnames(ch.agarico.5)<-c("Time","Change")
ch.agarico.5$Kind<-"Tree 88"

ch.r.agarico.5<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.5$segments)-as.vector(t(ctt.to.plot.mod.agarico.5$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.5$nchanges/ctt.to.plot.mod.agarico.5$edge.length,ctt.to.plot.mod.agarico.5$nchanges/ctt.to.plot.mod.agarico.5$edge.length))))
colnames(ch.r.agarico.5)<-c("Time","Change")
ch.r.agarico.5$Kind<-"Tree 88"

ch.agarico.6<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.6$segments)-as.vector(t(ctt.to.plot.mod.agarico.6$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.6$nchanges,ctt.to.plot.mod.agarico.6$nchanges))))
colnames(ch.agarico.6)<-c("Time","Change")
ch.agarico.6$Kind<-"Tree 165"

ch.r.agarico.6<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.6$segments)-as.vector(t(ctt.to.plot.mod.agarico.6$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.6$nchanges/ctt.to.plot.mod.agarico.6$edge.length,ctt.to.plot.mod.agarico.6$nchanges/ctt.to.plot.mod.agarico.6$edge.length))))
colnames(ch.r.agarico.6)<-c("Time","Change")
ch.r.agarico.6$Kind<-"Tree 165"

ch.agarico.7<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.7$segments)-as.vector(t(ctt.to.plot.mod.agarico.7$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.7$nchanges,ctt.to.plot.mod.agarico.7$nchanges))))
colnames(ch.agarico.7)<-c("Time","Change")
ch.agarico.7$Kind<-"Tree 170"

ch.r.agarico.7<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.7$segments)-as.vector(t(ctt.to.plot.mod.agarico.7$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.7$nchanges/ctt.to.plot.mod.agarico.7$edge.length,ctt.to.plot.mod.agarico.7$nchanges/ctt.to.plot.mod.agarico.7$edge.length))))
colnames(ch.r.agarico.7)<-c("Time","Change")
ch.r.agarico.7$Kind<-"Tree 170"

ch.agarico.8<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.8$segments)-as.vector(t(ctt.to.plot.mod.agarico.8$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.8$nchanges,ctt.to.plot.mod.agarico.8$nchanges))))
colnames(ch.agarico.8)<-c("Time","Change")
ch.agarico.8$Kind<-"Tree 188"

ch.r.agarico.8<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.8$segments)-as.vector(t(ctt.to.plot.mod.agarico.8$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.8$nchanges/ctt.to.plot.mod.agarico.8$edge.length,ctt.to.plot.mod.agarico.8$nchanges/ctt.to.plot.mod.agarico.8$edge.length))))
colnames(ch.r.agarico.8)<-c("Time","Change")
ch.r.agarico.8$Kind<-"Tree 188"

ch.agarico.9<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.9$segments)-as.vector(t(ctt.to.plot.mod.agarico.9$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.9$nchanges,ctt.to.plot.mod.agarico.9$nchanges))))
colnames(ch.agarico.9)<-c("Time","Change")
ch.agarico.9$Kind<-"Tree 216"

ch.r.agarico.9<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.9$segments)-as.vector(t(ctt.to.plot.mod.agarico.9$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.9$nchanges/ctt.to.plot.mod.agarico.9$edge.length,ctt.to.plot.mod.agarico.9$nchanges/ctt.to.plot.mod.agarico.9$edge.length))))
colnames(ch.r.agarico.9)<-c("Time","Change")
ch.r.agarico.9$Kind<-"Tree 216"

ch.agarico.10<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.10$segments)-as.vector(t(ctt.to.plot.mod.agarico.10$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.10$nchanges,ctt.to.plot.mod.agarico.10$nchanges))))
colnames(ch.agarico.10)<-c("Time","Change")
ch.agarico.10$Kind<-"Tree 222"

ch.r.agarico.10<-as.data.frame(cbind(max(ctt.to.plot.mod.agarico.10$segments)-as.vector(t(ctt.to.plot.mod.agarico.10$segments)),as.vector(rbind(ctt.to.plot.mod.agarico.10$nchanges/ctt.to.plot.mod.agarico.10$edge.length,ctt.to.plot.mod.agarico.10$nchanges/ctt.to.plot.mod.agarico.10$edge.length))))
colnames(ch.r.agarico.10)<-c("Time","Change")
ch.r.agarico.10$Kind<-"Tree 222"

ch.agarico<-rbind(ch.agarico.1,ch.agarico.2,ch.agarico.3,ch.agarico.4,ch.agarico.5,ch.agarico.6,ch.agarico.7,ch.agarico.8,ch.agarico.9,ch.agarico.10)
ch.agarico.r<-rbind(ch.r.agarico.1,ch.r.agarico.2,ch.r.agarico.3,ch.r.agarico.4,ch.r.agarico.5,ch.r.agarico.6,ch.r.agarico.7,ch.r.agarico.8,ch.r.agarico.9,ch.r.agarico.10)

ch.comb.agarico<-ch.agarico
colnames(ch.comb.agarico)<-c("Time","ChangeNo","Kind")
ch.comb.agarico<-cbind(ch.comb.agarico,ch.agarico.r$Change)
colnames(ch.comb.agarico)<-c("Time","ChangeNo","Kind","ChangeRate")

require(viridis)
mycols.agarico<-inferno(10)


timeplot.agarico<-ggplot(ch.comb.agarico,aes(x=Time,ymin=-0.5))+scale_x_reverse(expand=c(0,0),limits=c(500,0), breaks=c(500,0))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x.top=element_blank(), axis.title.x.bottom=element_blank(), axis.line.x=element_blank(), panel.background=element_blank(), axis.text.x=element_text(color="black"), axis.text.y=element_text(color="black"),axis.ticks.x = element_blank(), axis.ticks.y = element_line(color="black")) + annotate(geom='segment',y=-Inf,yend=Inf,x=Inf,xend=Inf) + annotate(geom='segment',y=-Inf,yend=Inf,x=-Inf,xend=-Inf) +annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=-0.25,fill=timey$RGB[1],color="black",alpha=1)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=-0.25,fill=timey$RGB[2],color="black",alpha=1)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=-0.25,fill=timey$RGB[3],color="black",alpha=1)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=-0.25,fill=timey$RGB[4],color="black",alpha=1)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=-0.25,fill=timey$RGB[5],color="black",alpha=1)+annotate("rect",xmin=timey$End[6],xmax=timey$Start[6],ymin=-Inf,ymax=-0.25,fill=timey$RGB[6],color="black",alpha=1)+annotate("rect",xmin=timey$End[7],xmax=timey$Start[7],ymin=-Inf,ymax=-0.25,fill=timey$RGB[7],color="black",alpha=1)+annotate("rect",xmin=timey$End[8],xmax=timey$Start[8],ymin=-Inf,ymax=-0.25,fill=timey$RGB[8],color="black",alpha=1)+annotate("rect",xmin=timey$End[9],xmax=timey$Start[9],ymin=-Inf,ymax=-0.25,fill=timey$RGB[9],color="black",alpha=1)+annotate("rect",xmin=timey$End[10],xmax=timey$Start[10],ymin=-Inf,ymax=-0.25,fill=timey$RGB[10],color="black",alpha=1)+annotate("rect",xmin=timey$End[11],xmax=timey$Start[11],ymin=-Inf,ymax=-0.25,fill=timey$RGB[11],color="black",alpha=1)+annotate("rect",xmin=timey$End[12],xmax=500,ymin=-Inf,ymax=-0.25,fill=timey$RGB[12],color="black",alpha=1)+annotate("text",x=timey$Midpoint[11],y=-0.475,label=timey$Name[11],size=3,colour="black")+annotate("text",x=timey$Midpoint[10],y=-0.475,label=timey$Name[10],size=3,colour="black")+annotate("text",x=timey$Midpoint[9],y=-0.475,label=timey$Name[9],size=3,colour="black")+annotate("text",x=timey$Midpoint[8],y=-0.475,label=timey$Name[8],size=3,colour="black")+annotate("text",x=timey$Midpoint[7],y=-0.475,label=timey$Name[7],size=3,colour="black")+annotate("text",x=timey$Midpoint[6],y=-0.475,label=timey$Name[6],size=3,colour="black")+annotate("text",x=timey$Midpoint[5],y=-0.475,label=timey$Name[5],size=3,colour="black")+annotate("text",x=timey$Midpoint[4],y=-0.475,label=timey$Name[4],size=3,colour="black")+annotate("text",x=timey$Midpoint[3],y=-0.475,label=timey$Name[3],size=3,colour="black")+annotate("text",x=timey$Midpoint[2],y=-0.475,label=timey$Name[2],size=3,colour="black")+annotate("text",x=1.3,y=-0.475,label=timey$Name[1],size=3,colour="black") 
timeplot.agarico<-timeplot.agarico+ annotate(geom='segment',y=3.2,yend=3.2,x=475,xend=450,color=mycols.agarico[1],linetype=1) + annotate(geom='segment',y=2.9,yend=2.9,x=475,xend=450,color=mycols.agarico[2],linetype=1)+ annotate(geom='segment',y=2.6,yend=2.6,x=475,xend=450,color=mycols.agarico[3],linetype=1)+ annotate(geom='segment',y=2.3,yend=2.3,x=475,xend=450,color=mycols.agarico[4],linetype=1)+ annotate(geom='segment',y=2.0,yend=2.0,x=475,xend=450,color=mycols.agarico[5],linetype=1)+ annotate(geom='segment',y=1.7,yend=1.7,x=475,xend=450,color=mycols.agarico[6],linetype=1)+ annotate(geom='segment',y=1.4,yend=1.4,x=475,xend=450,color=mycols.agarico[7],linetype=1)+ annotate(geom='segment',y=1.1,yend=1.1,x=475,xend=450,color=mycols.agarico[8],linetype=1)+ annotate(geom='segment',y=0.8,yend=0.8,x=475,xend=450,color=mycols.agarico[9],linetype=1)+ annotate(geom='segment',y=0.5,yend=0.5,x=475,xend=450,color=mycols.agarico[10],linetype=1) + annotate(geom='segment',y=3.2,yend=3.2,x=440,xend=415,color=mycols.agarico[1],linetype=3) + annotate(geom='segment',y=2.9,yend=2.9,x=440,xend=415,color=mycols.agarico[2],linetype=3)+ annotate(geom='segment',y=2.6,yend=2.6,x=440,xend=415,color=mycols.agarico[3],linetype=3)+ annotate(geom='segment',y=2.3,yend=2.3,x=440,xend=415,color=mycols.agarico[4],linetype=3)+ annotate(geom='segment',y=2.0,yend=2.0,x=440,xend=415,color=mycols.agarico[5],linetype=3)+ annotate(geom='segment',y=1.7,yend=1.7,x=440,xend=415,color=mycols.agarico[6],linetype=3)+ annotate(geom='segment',y=1.4,yend=1.4,x=440,xend=415,color=mycols.agarico[7],linetype=3)+ annotate(geom='segment',y=1.1,yend=1.1,x=440,xend=415,color=mycols.agarico[8],linetype=3)+ annotate(geom='segment',y=0.8,yend=0.8,x=440,xend=415,color=mycols.agarico[9],linetype=3)+ annotate(geom='segment',y=0.5,yend=0.5,x=440,xend=415,color=mycols.agarico[10],linetype=3) + annotate(geom="text", x=405, y=3.2, label="8",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=2.9, label="12",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=2.6, label="16",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=2.3, label="53",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=2.0, label="88",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=1.7, label="165",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=1.4, label="170",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=1.1, label="188",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=0.8, label="216",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=0.5, label="222",color="black",hjust=0,size=2.5) + annotate(geom="text", x=475, y=3.6, label="Number",color="black",hjust=0,size=2.5,fontface=2)+ annotate(geom="text", x=440, y=3.6, label="Rate",color="black",hjust=0,size=2.5,fontface=2)+ annotate(geom="text", x=405, y=3.6, label="Tree",color="black",hjust=0,size=2.5,fontface=2)+ annotate(geom="text", x=410, y=4.5, label="Evolution of ectomycorrhizal associations in Agaricomycotina",color="black",hjust=0,size=3.5,fontface=2)
timeplot.agarico.comb<-timeplot.agarico+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 8",]),linetype=1,alpha=1,color=mycols.agarico[1],aes(y=ChangeNo))
timeplot.agarico.combb<-timeplot.agarico.comb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 12",]),linetype=1,alpha=1,color=mycols.agarico[2],aes(y=ChangeNo))
timeplot.agarico.combb<-timeplot.agarico.combb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 16",]),linetype=1,alpha=1,color=mycols.agarico[3],aes(y=ChangeNo))
timeplot.agarico.combb<-timeplot.agarico.combb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 53",]),linetype=1,alpha=1,color=mycols.agarico[4],aes(y=ChangeNo))
timeplot.agarico.combb<-timeplot.agarico.combb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 88",]),linetype=1,alpha=1,color=mycols.agarico[5],aes(y=ChangeNo))
timeplot.agarico.combb<-timeplot.agarico.combb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 165",]),linetype=1,alpha=1,color=mycols.agarico[6],aes(y=ChangeNo))
timeplot.agarico.combb<-timeplot.agarico.combb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 170",]),linetype=1,alpha=1,color=mycols.agarico[7],aes(y=ChangeNo))
timeplot.agarico.combb<-timeplot.agarico.combb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 188",]),linetype=1,alpha=1,color=mycols.agarico[8],aes(y=ChangeNo))
timeplot.agarico.combb<-timeplot.agarico.combb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 216",]),linetype=1,alpha=1,color=mycols.agarico[9],aes(y=ChangeNo))
timeplot.agarico.combb<-timeplot.agarico.combb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 222",]),linetype=1,alpha=1,color=mycols.agarico[10],aes(y=ChangeNo))

timeplot.agarico.combbb<-timeplot.agarico.combb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 8",]),linetype=3,alpha=1,color=mycols.agarico[1],aes(y=ChangeRate*666))
timeplot.agaricoc<-timeplot.agarico.combbb+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 12",]),linetype=3,alpha=1,color=mycols.agarico[2],aes(y=ChangeRate*666))
timeplot.agaricoc<-timeplot.agaricoc+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 16",]),linetype=3,alpha=1,color=mycols.agarico[3],aes(y=ChangeRate*666))
timeplot.agaricoc<-timeplot.agaricoc+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 53",]),linetype=3,alpha=1,color=mycols.agarico[4],aes(y=ChangeRate*666))
timeplot.agaricoc<-timeplot.agaricoc+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 88",]),linetype=3,alpha=1,color=mycols.agarico[5],aes(y=ChangeRate*666))
timeplot.agaricoc<-timeplot.agaricoc+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 165",]),linetype=3,alpha=1,color=mycols.agarico[6],aes(y=ChangeRate*666))
timeplot.agaricoc<-timeplot.agaricoc+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 170",]),linetype=3,alpha=1,color=mycols.agarico[7],aes(y=ChangeRate*666))
timeplot.agaricoc<-timeplot.agaricoc+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 188",]),linetype=3,alpha=1,color=mycols.agarico[8],aes(y=ChangeRate*666))
timeplot.agaricoc<-timeplot.agaricoc+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 216",]),linetype=3,alpha=1,color=mycols.agarico[9],aes(y=ChangeRate*666))
timeplot.agaricoc<-timeplot.agaricoc+geom_step(data=subset(ch.comb.agarico[ch.comb.agarico$Kind=="Tree 222",]),linetype=3,alpha=1,color=mycols.agarico[10],aes(y=ChangeRate*666))
timeplot.agaricod<-timeplot.agaricoc+scale_y_continuous(name="Mean number of gains /\n5 my time bin", sec.axis = sec_axis(~./666, name = "Mean number of gains / \ntotal branch length / 5 my time bin"))+theme(axis.title.y = element_text(size=8))
timeplot.agaricoe<-timeplot.agaricod
#pdf(file="~/agaricomycotina_ecm_vortex/ctt_trees.pdf",width=7,height=2)
#timeplot.agaricoe
#dev.off()

require(gridExtra)
pdf(file="...ctt_combined_no_names.pdf",width=7,height=4)
grid.arrange(timeplot.lece, timeplot.agaricoe, nrow = 2)
dev.off()


#########################
####EXTRA FUNCTIONS
#########################

#modified so user can specify timebin duration - truncates oldest one near root.
#standardizing timebins allows users to plot trees of different heights together with standardized time windows.
#also modified so use can specify which changes specifically to look at.
ctt.mod2<-function (trees, timebin = 5,direction=NULL){
	#modified from phytools 0.7.70
    if (!(inherits(trees, "multiSimmap"))) 
        stop("trees should be an object of class \"multiSimmap\".")
    tree <- as.phylo(trees[[1]])
    changes <- sapply(trees, phytools:::getChanges)
    #if wanting to focus on specific changes, keep changes from the changes list that are the desired changes
    if(!is.null(direction)){
		#make an empty changes list
		changes.mod <- vector(mode = "list", length = length(changes))
		for(chng in 1:length(changes)){
			changes.mod[[chng]]<-changes[[chng]][names(changes[[chng]]) %in% direction]
		}
		changes<-changes.mod	    
    }
    h <- max(nodeHeights(tree))
    b <- ceiling(h/timebin)
    segs <- cbind(rev(c(seq(h,0,-timebin),0))[1:b], rev(c(seq(h,0,-timebin),0))[1:b+1])
    #first cell needs to be partial because 0 is actually root in this situation
    #segs <- cbind(seq(0, h, timebin), seq(timebin, timebin*b, timebin))
    #change end of max segment to max height of tree
    #segs[nrow(segs),2]<-h    
    nchanges <- rep(0, b)
    for (i in 1:length(changes)) {
        for (j in 1:length(changes[[i]])) {
            ind <- which((changes[[i]][j] > segs[, 1]) + (changes[[i]][j] <= 
                segs[, 2]) == 2)
            nchanges[ind] <- nchanges[ind] + 1/length(changes)
        }
    }
    LTT <- ltt(tree, plot = FALSE)
    LTT <- cbind(LTT$ltt[2:(length(LTT$ltt) - 1)], LTT$times[2:(length(LTT$ltt) - 
        1)], LTT$times[3:length(LTT$ltt)])
    ii <- 1
    edge.length <- rep(0, b)
    for (i in 1:nrow(segs)) {
        done.seg <- FALSE
        while (LTT[ii, 2] <= segs[i, 2] && done.seg == FALSE) {
            edge.length[i] <- edge.length[i] + LTT[ii, 1] * (min(segs[i, 
                2], LTT[ii, 3]) - max(segs[i, 1], LTT[ii, 2]))
            if (LTT[ii, 3] >= segs[i, 2]) 
                done.seg <- TRUE
            if (LTT[ii, 3] <= segs[i, 2]) 
                ii <- if (ii < nrow(LTT)) 
                  ii + 1
                else ii
        }
    }
    object <- list(segments = segs, nchanges = nchanges, edge.length = edge.length, 
        tree = tree)
    class(object) <- "ctt"
    object
}
