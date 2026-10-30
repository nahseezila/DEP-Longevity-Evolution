library(ape)
library(phytools)
library(caper)
library(picante)

#read in tree
Flytree<-read.nexus("../FlyTree_correctednames.3.0.nex")
#will fail the ultrametric test due to numerical imprescision
ult.Flytree<-force.ultrametric(Flytree)

#read in data
raw.data<-read.table(file="DEP_Lifespan_data_10Dec2025.txt",header=T,row.names=1,sep = "\t")

#Combine tree and traits
Combo<-match.phylo.data(ult.Flytree,raw.data)

#Individual traits
AgeReproMatF<-as.numeric(Combo$data$Age_reproductive_maturity_F)
AgeReproMatM<-as.numeric(Combo$data$Age_reproductive_maturity_M)
e0_f<-as.numeric(Combo$data$e0_f)
e0_m<-as.numeric(Combo$data$e0_m)
ls90_f<-as.numeric(Combo$data$ls90_f)
ls90_m<-as.numeric(Combo$data$ls90_m)
Entropy_Ketfitz_f<-as.numeric(Combo$data$Keyfitz_Entropy_f)
Entropy_Ketfitz_m<-as.numeric(Combo$data$Keyfitz_Entropy_m)
Lifespan_Gini_f<-as.numeric(Combo$data$Lifespan_Gini_f)
Lifespan_Gini_m<-as.numeric(Combo$data$Lifespan_Gini_m)
rate_f<-as.numeric(Combo$data$rate_f)
rate_m<-as.numeric(Combo$data$rate_m)
netreproout<-as.numeric(Combo$data$net_reproductive_output)
Senescence<-as.numeric(Combo$data$Senescence_Slope)
Shape_rep<-as.numeric(Combo$data$Shape_Rep)
Peakfecundity<-as.numeric(Combo$data$Peak_Fecundity)

names(AgeReproMatF)<-Combo$phy$tip.label
names(AgeReproMatM)<-Combo$phy$tip.label
names(e0_f)<-Combo$phy$tip.label
names(e0_m)<-Combo$phy$tip.label
names(ls90_f)<-Combo$phy$tip.label
names(ls90_m)<-Combo$phy$tip.label
names(Entropy_Ketfitz_f)<-Combo$phy$tip.label
names(Entropy_Ketfitz_m)<-Combo$phy$tip.label
names(Lifespan_Gini_f)<-Combo$phy$tip.label
names(Lifespan_Gini_m)<-Combo$phy$tip.label
names(rate_f)<-Combo$phy$tip.label
names(rate_m)<-Combo$phy$tip.label
names(netreproout)<-Combo$phy$tip.label
names(Senescence)<-Combo$phy$tip.label
names(Shape_rep)<-Combo$phy$tip.label
names(Peakfecundity)<-Combo$phy$tip.label

#Model test
#e0_f
fitContinuous(Combo$phy,dat = e0_f,model = "BM")
fitContinuous(Combo$phy,dat = e0_f,model = "OU")
fitContinuous(Combo$phy,dat = e0_f,model = "EB")
fitContinuous(Combo$phy,dat = e0_f,model = "rate_trend")
fitContinuous(Combo$phy,dat = e0_f,model = "lambda")
AICs.e0_f<-c(1253.68,1188.06,1255.68,1244.80,1174.39)
names(AICs.e0_f)<-c("BM","OU","EB","Rate-trend","lambda")
aic.w(AICs.e0_f)
aicw(AICs.e0_f)

#e0_m
fitContinuous(Combo$phy,dat = e0_m,model = "BM")
fitContinuous(Combo$phy,dat = e0_m,model = "OU")
fitContinuous(Combo$phy,dat = e0_m,model = "EB")
fitContinuous(Combo$phy,dat = e0_m,model = "rate_trend")
fitContinuous(Combo$phy,dat = e0_m,model = "lambda")
AICs.e0_m<-c(1279.58,1228.43,1281.58,1271.51,1223.94)
names(AICs.e0_m)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.e0_m)
aic.w(AICs.e0_m)

#ls_90_f
fitContinuous(Combo$phy,dat = ls90_f,model = "BM")
fitContinuous(Combo$phy,dat = ls90_f,model = "OU")
fitContinuous(Combo$phy,dat = ls90_f,model = "EB")
fitContinuous(Combo$phy,dat = ls90_f,model = "rate_trend")
fitContinuous(Combo$phy,dat = ls90_f,model = "lambda")
AICs.ls90_f<-c(1346.54,1302.61,1348.54,1338.79,1292.31)
names(AICs.ls90_f)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.ls90_f)
aic.w(AICs.ls90_f)

#ls_90_m
fitContinuous(Combo$phy,dat = ls90_m,model = "BM")
fitContinuous(Combo$phy,dat = ls90_m,model = "OU")
fitContinuous(Combo$phy,dat = ls90_m,model = "EB")
fitContinuous(Combo$phy,dat = ls90_m,model = "rate_trend")
fitContinuous(Combo$phy,dat = ls90_m,model = "lambda")
AICs.ls90_m<-c(1376.57,1341.96,1378.57,1369.53,1339.64)
names(AICs.ls90_m)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.ls90_m)
aic.w(AICs.ls90_m)

#rate_f
fitContinuous(Combo$phy,dat = rate_f,model = "BM")
fitContinuous(Combo$phy,dat = rate_f,model = "OU")
fitContinuous(Combo$phy,dat = rate_f,model = "EB")
fitContinuous(Combo$phy,dat = rate_f,model = "rate_trend")
fitContinuous(Combo$phy,dat = rate_f,model = "lambda")
AICs.rate_f<-c(614.58,513.65,616.58,604.98,510.38)
names(AICs.rate_f)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.rate_f)
aic.w(AICs.rate_f)

#rate_m
fitContinuous(Combo$phy,dat = rate_m,model = "BM")
fitContinuous(Combo$phy,dat = rate_m,model = "OU")
fitContinuous(Combo$phy,dat = rate_m,model = "EB")
fitContinuous(Combo$phy,dat = rate_m,model = "rate_trend")
fitContinuous(Combo$phy,dat = rate_m,model = "lambda")
AICs.rate_m<-c(622.54,525.01,624.54,612.95,521.01)
names(AICs.rate_m)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.rate_m)
aic.w(AICs.rate_m)

#Gini_f
fitContinuous(Combo$phy,dat = Lifespan_Gini_f,model = "BM")
fitContinuous(Combo$phy,dat = Lifespan_Gini_f,model = "OU")
fitContinuous(Combo$phy,dat = Lifespan_Gini_f,model = "EB")
fitContinuous(Combo$phy,dat = Lifespan_Gini_f,model = "rate_trend")
fitContinuous(Combo$phy,dat = Lifespan_Gini_f,model = "lambda")
AICs.Lifespan_Gini_f<-c()
names(AICs.Lifespan_Gini_f)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.Lifespan_Gini_f)
aic.w(AICs.Lifespan_Gini_f)

#Gini_m
fitContinuous(Combo$phy,dat = Lifespan_Gini_m,model = "BM")
fitContinuous(Combo$phy,dat = Lifespan_Gini_m,model = "OU")
fitContinuous(Combo$phy,dat = Lifespan_Gini_m,model = "EB")
fitContinuous(Combo$phy,dat = Lifespan_Gini_m,model = "rate_trend")
fitContinuous(Combo$phy,dat = Lifespan_Gini_m,model = "lambda")
AICs.Lifespan_Gini_m<-c(-182.05,-281.15,-180.05,-191.43,-285.95)
names(AICs.Lifespan_Gini_m)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.Lifespan_Gini_m)
aic.w(AICs.Lifespan_Gini_m)

#Keyfitz-entropy_f
fitContinuous(Combo$phy,dat = Entropy_Ketfitz_f,model = "BM")
fitContinuous(Combo$phy,dat = Entropy_Ketfitz_f,model = "OU")
fitContinuous(Combo$phy,dat = Entropy_Ketfitz_f,model = "EB")
fitContinuous(Combo$phy,dat = Entropy_Ketfitz_f,model = "rate_trend")
fitContinuous(Combo$phy,dat = Entropy_Ketfitz_f,model = "lambda")
AICs.Entropy_Ketfitz_f<-c(-54.58,-167.63,-52.58,-64.47,-175.65)
names(AICs.Entropy_Ketfitz_f)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.Entropy_Ketfitz_f)
aic.w(AICs.Entropy_Ketfitz_f)

#Keyfitz-entropy_m
fitContinuous(Combo$phy,dat = Entropy_Ketfitz_m,model = "BM")
fitContinuous(Combo$phy,dat = Entropy_Ketfitz_m,model = "OU")
fitContinuous(Combo$phy,dat = Entropy_Ketfitz_m,model = "EB")
fitContinuous(Combo$phy,dat = Entropy_Ketfitz_m,model = "rate_trend")
fitContinuous(Combo$phy,dat = Entropy_Ketfitz_m,model = "lambda")
AICs.Entropy_Ketfitz_m<-c(23.78,-73.10,25.78,14.48,-76.97)
names(AICs.Entropy_Ketfitz_m)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.Entropy_Ketfitz_m)
aic.w(AICs.Entropy_Ketfitz_m)

#Net_reproductive_output
fitContinuous(Combo$phy,dat = netreproout,model = "BM")
fitContinuous(Combo$phy,dat = netreproout,model = "OU")
fitContinuous(Combo$phy,dat = netreproout,model = "EB")
fitContinuous(Combo$phy,dat = netreproout,model = "rate_trend")
fitContinuous(Combo$phy,dat = netreproout,model = "lambda")
fitContinuous(Combo$phy,dat = netreproout,model = "delta")
AICs.netreproout<-c(2029.84,1938.05,2031.85,2020.28,1936.05)
names(AICs.netreproout)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.netreproout)
aic.w(AICs.netreproout)

#Peak fecundity
fitContinuous(Combo$phy,dat = Peakfecundity,model = "BM")
fitContinuous(Combo$phy,dat = Peakfecundity,model = "OU")
fitContinuous(Combo$phy,dat = Peakfecundity,model = "EB")
fitContinuous(Combo$phy,dat = Peakfecundity,model = "rate_trend")
fitContinuous(Combo$phy,dat = Peakfecundity,model = "lambda")
fitContinuous(Combo$phy,dat = Peakfecundity,model = "delta")
AICs.Peakfecundity<-c()
names(AICs.Peakfecundity)<-c("BM","OU","EB","Rate-trend","lambda","delta")
aicw(AICs.Peakfecundity)
aic.w(AICs.Peakfecundity)

#Shape_rep
fitContinuous(Combo$phy,dat = Shape_rep,model = "BM")
fitContinuous(Combo$phy,dat = Shape_rep,model = "OU")
fitContinuous(Combo$phy,dat = Shape_rep,model = "EB")
fitContinuous(Combo$phy,dat = Shape_rep,model = "rate_trend")
fitContinuous(Combo$phy,dat = Shape_rep,model = "lambda")
fitContinuous(Combo$phy,dat = Shape_rep,model = "delta")
AICs.Shape_rep<-c()
names(AICs.Shape_rep)<-c("BM","OU","EB","Rate-trend","lambda","delta")
aicw(AICs.Shape_rep)
aic.w(AICs.Shape_rep)

#Senescence_slope
fitContinuous(Combo$phy,dat = Senescence,model = "BM")
fitContinuous(Combo$phy,dat = Senescence,model = "OU")
fitContinuous(Combo$phy,dat = Senescence,model = "EB")
fitContinuous(Combo$phy,dat = Senescence,model = "rate_trend")
fitContinuous(Combo$phy,dat = Senescence,model = "lambda")
fitContinuous(Combo$phy,dat = Senescence,model = "delta")
AICs.Senescence<-c(-408.12,-517.27,-406.12,-417.48,-518.28)
names(AICs.Senescence)<-c("BM","OU","EB","Rate-trend","lambda")
aicw(AICs.Senescence)
aic.w(AICs.Senescence)

phylosig(Combo$phy,x = e0_f,method = "lambda")
phylosig(Combo$phy,x = e0_f,method = "K")

phylosig(Combo$phy,x = e0_m,method = "lambda")
phylosig(Combo$phy,x = e0_m,method = "K")

phylosig(Combo$phy,x = ls90_f,method = "lambda")
phylosig(Combo$phy,x = ls90_f,method = "K")

phylosig(Combo$phy,x = ls90_m,method = "lambda")
phylosig(Combo$phy,x = ls90_m,method = "K")

phylosig(Combo$phy,x = rate_f,method = "lambda")
phylosig(Combo$phy,x = rate_f,method = "K")

phylosig(Combo$phy,x = rate_m,method = "lambda")
phylosig(Combo$phy,x = rate_m,method = "K")

phylosig(Combo$phy,x = Lifespan_Gini_f,method = "lambda")
phylosig(Combo$phy,x = Lifespan_Gini_f,method = "K")

phylosig(Combo$phy,x = Lifespan_Gini_m,method = "lambda")
phylosig(Combo$phy,x = Lifespan_Gini_m,method = "K")

phylosig(Combo$phy,x = Entropy_Ketfitz_f,method = "lambda")
phylosig(Combo$phy,x = Entropy_Ketfitz_f,method = "K")

phylosig(Combo$phy,x = Entropy_Ketfitz_m,method = "lambda")
phylosig(Combo$phy,x = Entropy_Ketfitz_m,method = "K")

phylosig(Combo$phy,x = netreproout,method = "lambda")
phylosig(Combo$phy,x = netreproout,method = "K")

phylosig(Combo$phy,x = Peakfecundity,method = "lambda")
phylosig(Combo$phy,x = Peakfecundity,method = "K")

phylosig(Combo$phy,x = Shape_rep,method = "lambda")
phylosig(Combo$phy,x = Shape_rep,method = "K")

phylosig(Combo$phy,x = Senescence,method = "lambda")
phylosig(Combo$phy,x = Senescence,method = "K")

#Focus in on species groups
#Do species group explain some variationa?
groups<-Combo$data$Group
names(groups)<-Combo$phy$tip.label

PrunedtoGroups<-drop.tip(Combo$phy,c("S_lebanonensis","C_procnemis"))

#Prune Combo to match
PrunedCombo<-match.phylo.data(PrunedtoGroups,Combo$data)

groups.pruned<-PrunedCombo$data$Group
names(groups.pruned)<-PrunedCombo$phy$tip.label

AgeReproMatF.p<-as.numeric(PrunedCombo$data$Age_reproductive_maturity_F)
AgeReproMatM.p<-as.numeric(PrunedCombo$data$Age_reproductive_maturity_M)
e0_f.p<-as.numeric(PrunedCombo$data$e0_f)
e0_m.p<-as.numeric(PrunedCombo$data$e0_m)
ls90_f.p<-as.numeric(PrunedCombo$data$ls90_f)
ls90_m.p<-as.numeric(PrunedCombo$data$ls90_m)
Entropy_Ketfitz_f.p<-as.numeric(PrunedCombo$data$Keyfitz_Entropy_f)
Entropy_Ketfitz_m.p<-as.numeric(PrunedCombo$data$Keyfitz_Entropy_m)
Lifespan_Gini_f.p<-as.numeric(PrunedCombo$data$Lifespan_Gini_f)
Lifespan_Gini_m.p<-as.numeric(PrunedCombo$data$Lifespan_Gini_m)
rate_f.p<-as.numeric(PrunedCombo$data$rate_f)
rate_m.p<-as.numeric(PrunedCombo$data$rate_m)
netreproout.p<-as.numeric(PrunedCombo$data$net_reproductive_output)
Senescence.p<-as.numeric(PrunedCombo$data$Senescence_Slope)
Shape_rep.p<-as.numeric(PrunedCombo$data$Shape_Rep)
Peakfecundity.p<-as.numeric(PrunedCombo$data$Peak_Fecundity)

names(AgeReproMatF.p)<-PrunedCombo$phy$tip.label
names(AgeReproMatM.p)<-PrunedCombo$phy$tip.label
names(e0_f.p)<-PrunedCombo$phy$tip.label
names(e0_m.p)<-PrunedCombo$phy$tip.label
names(ls90_f.p)<-PrunedCombo$phy$tip.label
names(ls90_m.p)<-PrunedCombo$phy$tip.label
names(Entropy_Ketfitz_f.p)<-PrunedCombo$phy$tip.label
names(Entropy_Ketfitz_m.p)<-PrunedCombo$phy$tip.label
names(Lifespan_Gini_f.p)<-PrunedCombo$phy$tip.label
names(Lifespan_Gini_m.p)<-PrunedCombo$phy$tip.label
names(rate_f.p)<-PrunedCombo$phy$tip.label
names(rate_m.p)<-PrunedCombo$phy$tip.label
names(netreproout.p)<-PrunedCombo$phy$tip.label
names(Senescence.p)<-PrunedCombo$phy$tip.label
names(Shape_rep.p)<-PrunedCombo$phy$tip.label
names(Peakfecundity.p)<-PrunedCombo$phy$tip.label

#Phylogenetic ANOVAs
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=AgeReproMatF.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=AgeReproMatM.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=e0_f.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=e0_m.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=ls90_f.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=ls90_m.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=Entropy_Ketfitz_f.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=Entropy_Ketfitz_m.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=Lifespan_Gini_f.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=Lifespan_Gini_m.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=rate_f.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=rate_m.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=netreproout.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=Senescence.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=Shape_rep.p)
phylANOVA(PrunedCombo$phy,x=groups.pruned,y=Peakfecundity.p)

phylosig(PrunedCombo$phy,AgeReproMatF.p,method="lambda")
phylosig(PrunedCombo$phy,AgeReproMatM.p,method="lambda")
phylosig(PrunedCombo$phy,e0_f.p,method="lambda")
phylosig(PrunedCombo$phy,e0_m.p,method="lambda")
phylosig(PrunedCombo$phy,ls90_f.p,method="lambda")
phylosig(PrunedCombo$phy,ls90_m.p,method="lambda")
phylosig(PrunedCombo$phy,Entropy_Ketfitz_f.p,method="lambda")
phylosig(PrunedCombo$phy,Entropy_Ketfitz_m.p,method="lambda")
phylosig(PrunedCombo$phy,Lifespan_Gini_f.p,method="lambda")
phylosig(PrunedCombo$phy,Lifespan_Gini_m.p,method="lambda")
phylosig(PrunedCombo$phy,rate_f.p,method="lambda")
phylosig(PrunedCombo$phy,rate_m.p,method="lambda")
phylosig(PrunedCombo$phy,netreproout.p,method="lambda")
phylosig(PrunedCombo$phy,Senescence.p,method="lambda")
phylosig(PrunedCombo$phy,Shape_rep.p,method="lambda")
phylosig(PrunedCombo$phy,Peakfecundity.p,method="lambda")

#USING broader clades/groups
broadgroup<-Combo$data$Broad_clade
names(broadgroup)<-Combo$phy$tip.label
strings_to_remove<-c("Miscellaneous")
updated_vector <- broadgroup[!broadgroup %in% strings_to_remove]
broadgroup<-updated_vector
BroadGroupTreeMatch<-match.phylo.data(Combo$phy,broadgroup)
BroadGroupDataMatch<-match.phylo.data(BroadGroupTreeMatch$phy,Combo$data)

#extract "broad" group vector
broadgroup <- BroadGroupDataMatch$data$Broad_clade
names(broadgroup)<-BroadGroupDataMatch$phy$tip.label

AgeReproMatF.bg<-as.numeric(BroadGroupDataMatch$data$Age_reproductive_maturity_F)
AgeReproMatM.bg<-as.numeric(BroadGroupDataMatch$data$Age_reproductive_maturity_M)
e0_f.bg<-as.numeric(BroadGroupDataMatch$data$e0_f)
e0_m.bg<-as.numeric(BroadGroupDataMatch$data$e0_m)
ls90_f.bg<-as.numeric(BroadGroupDataMatch$data$ls90_f)
ls90_m.bg<-as.numeric(BroadGroupDataMatch$data$ls90_m)
Entropy_Ketfitz_f.bg<-as.numeric(BroadGroupDataMatch$data$Keyfitz_Entropy_f)
Entropy_Ketfitz_m.bg<-as.numeric(BroadGroupDataMatch$data$Keyfitz_Entropy_m)
Lifespan_Gini_f.bg<-as.numeric(BroadGroupDataMatch$data$Lifespan_Gini_f)
Lifespan_Gini_m.bg<-as.numeric(BroadGroupDataMatch$data$Lifespan_Gini_m)
rate_f.bg<-as.numeric(BroadGroupDataMatch$data$rate_f)
rate_m.bg<-as.numeric(BroadGroupDataMatch$data$rate_m)
netreproout.bg<-as.numeric(BroadGroupDataMatch$data$net_reproductive_output)
Senescence.bg<-as.numeric(BroadGroupDataMatch$data$Senescence_Slope)
Shape_rep.bg<-as.numeric(BroadGroupDataMatch$data$Shape_Rep)
Peakfecundity.bg<-as.numeric(BroadGroupDataMatch$data$Peak_Fecundity)

names(AgeReproMatF.bg)<-BroadGroupDataMatch$phy$tip.label
names(AgeReproMatM.bg)<-BroadGroupDataMatch$phy$tip.label
names(e0_f.bg)<-BroadGroupDataMatch$phy$tip.label
names(e0_m.bg)<-BroadGroupDataMatch$phy$tip.label
names(ls90_f.bg)<-BroadGroupDataMatch$phy$tip.label
names(ls90_m.bg)<-BroadGroupDataMatch$phy$tip.label
names(Entropy_Ketfitz_f.bg)<-BroadGroupDataMatch$phy$tip.label
names(Entropy_Ketfitz_m.bg)<-BroadGroupDataMatch$phy$tip.label
names(Lifespan_Gini_f.bg)<-BroadGroupDataMatch$phy$tip.label
names(Lifespan_Gini_m.bg)<-BroadGroupDataMatch$phy$tip.label
names(rate_f.bg)<-BroadGroupDataMatch$phy$tip.label
names(rate_m.bg)<-BroadGroupDataMatch$phy$tip.label
names(netreproout.bg)<-BroadGroupDataMatch$phy$tip.label
names(Senescence.bg)<-BroadGroupDataMatch$phy$tip.label
names(Shape_rep.bg)<-BroadGroupDataMatch$phy$tip.label
names(Peakfecundity.bg)<-BroadGroupDataMatch$phy$tip.label

#Phylogenetic ANOVAs using "broader" groups
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=AgeReproMatF.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=AgeReproMatM.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=e0_f.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=ae0_m.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=ls90_f.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=ls90_m.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=Entropy_Ketfitz_f.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=Entropy_Ketfitz_m.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=Lifespan_Gini_f.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=Lifespan_Gini_m.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=rate_f.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=rate_m.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=netreproout.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=Senescence.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=Shape_rep.bg)
phylANOVA(BroadGroupDataMatch$phy,x=broadgroup,y=Peakfecundity.bg)

#Fan Tree with traits
colors<-list(
  c("steelblue","white","red"),  c("steelblue","white","red"),  c("steelblue","white","red"),  c("steelblue","white","red"),
  terrain.colors(n=10),terrain.colors(n=10), terrain.colors(n=10),terrain.colors(n=10), terrain.colors(n=10),terrain.colors(n=10),c("black","gray","white"),c("black","gray","white"),c("black","gray","white"),c("black","gray","white"))

cols<-plotFanTree.wTraits(ladderize(Combo$phy,right=T),m_numeric, right=FALSE,
                          part=0.5,arc_height=1,spacer=0.1,fsize=0.7,colors=colors)
h<-max(nodeHeights(Combo$phy))
arc_height<-1
text(x=h+seq(0.05,by=0.07,length.out=ncol(m_numeric))*h+
       arc_height*h,y=rep(-0.5,ncol(m_numeric)),
     c(paste("",colnames(m_numeric)[1:14],")",sep="")),
     srt=-90,adj=c(0, 0.25),cex=0.5)
text(x=-h-seq(0.05,by=0.07,length.out=ncol(m_numeric))*h-
       arc_height*h,y=rep(-0.5,ncol(m_numeric)),
     c(paste("",colnames(m_numeric)[1:14],")",sep="")),
     srt=90,adj=c(1, 0.25),cex=0.5)
add.color.bar(1.5*h,cols[[1]],title="",lims=NULL,digits=0,
              subtitle="",x=-0.75*h,y=0,prompt=FALSE,outline=FALSE,
              lwd=8)
text(x=-0.75*h,y=-0.25,"smaller",adj=c(0,1))
text(x=0.75*h,y=-0.25,"larger",adj=c(1,1))
add.color.bar(1.5*h,cols[[5]],title="",lims=NULL,digits=0,
              subtitle="",x=-0.75*h,y=10,prompt=FALSE,outline=FALSE,
              lwd=8)
add.color.bar(1.5*h,cols[[14]],title="",lims=NULL,digits=0,
              subtitle="",x=-0.75*h,y=20,prompt=FALSE,outline=FALSE,
              lwd=8)
