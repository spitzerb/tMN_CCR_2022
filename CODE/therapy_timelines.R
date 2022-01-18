##########
#  
#  Generate timelines for all patients (Fig S5 and parts in main Figure 3)
#  3 rows per patient: treatment, mutations, structural variants/copy number changes
#     Output should go in directory "FIGURES" in the repositry root directory
#
#  Copyright (C) 2021 Barbara Spitzer - spitzerb@mskcc.org
#  This program is free software: you can redistribute it and/or modify 
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or 
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
##########


library("shape")
library("seqinr")
library("dplyr")
library("ggplot2")

samplepink <- "#E1B7B9"
lightpurple <- "lavenderblush3"

tx.colors <- c("lightcyan3","darkorange","thistle4","darkorange") #RT and RI should have the same color for publication
cg.colors <- c("lightcyan3","darkorange","thistle4","salmon","#6092d9","#bc4862", "#aaa23f")
colors<- c("paleturquoise4","darkorange","thistle4","olivedrab3","lightcyan3", "lavenderblush4","lightsalmon","cornflowerblue")

basey=20
c=1


######
# Plot chemo, RT, immunotherapy, RI (radioimmunotherapy)
#   For the purposes of the final plot, left RI color same as RT, but could distinguish in tx.color vec if desired
# wd is for lwd in segment width
######
plot.therapy=function(tt, wd=18){
  for(therapy in c("chemo","RT","immuno","RI")){
    for(i in 1:nrow(tt)) {
      if(!is.na(tt[i, paste0(therapy,".TTS")]))
        
        segments(x0=tt[i,paste0(therapy,".TTS")], x1=tt[i,paste0(therapy,".TTE")]+1, y0=basey, lty=1, lwd=wd, lend="square", col=tx.colors[c])
                                                                                #fudge factor +1 to account for single day treatments/prevent diamonds
    }
    if(basey==80){
      basey=50
    }
    else{
    basey=basey+30
    }
    c=c+1
  }
}

plot.SCharvest <- function(tt) {
  for(i in 1:nrow(tt)) {
    pt <- tt$patient[i] 
    ystart <- (pt-1)*200
    if(!is.na(tt[i, "SC_harvest_TTS"])) {
      segments(x0=tt[i,"SC_harvest_TTS"], y0=1, y1=200, lty=2, lwd=1.5, lend="square", col="gray70")
    }
  }
}

plot.auto <- function(tt) {
  for(i in 1:nrow(tt)) {
    pt <- tt$patient[i] 
    ystart <- (pt-1)*200
    if(!is.na(tt[i, "auto_trans_TTS"])) {
      segments(x0=tt[i,"auto_trans_TTS"], y0=1, y1=200, lty=2, lwd=1.5, lend="square")
    }
  }
}

plot.allo <- function(tt) {
  for(i in 1:nrow(tt)) {
    pt <- tt$patient[i] 
    ystart <- (pt-1)*200
    if(!is.na(tt[i, "allo_trans_TTS"])) {
      segments(x0=tt[i,"allo_trans_TTS"], y0=1, y1=200, lty=1, lwd=1.5, lend="square")
    }
  }
}

######
# Plot available samples
# This function can distinguish between those "already sent" and those "tbd" (to be sent) and those not intended to be sent (diff colors)
#   but ultimately only used "already sent" (allsamps=FALSE)
# default allsamps=FALSE=plot only those sent - TRUE prints all samples in the table
# cgannot: TRUE changes the shape based on whether cytogenetics were normal, abnormal, or not done.  default=FALSE, all circles
######
plot.NKCTB <- function(tt, allsamps=FALSE, cgannot=FALSE, cx=5) {
  if(!allsamps)
    tt <- tt[!tt[,"NKC_TB.sent."]=="n",]
  for(i in 1:nrow(tt)) {
    if(!is.na(tt[i, "NKC_TB.TTS"])) {
      cg.nl <- tt[i, "CG_nl"]
      is.nb <- tt[i, "is_NB"]
      
      color <- if(tt[i, "NKC_TB.sent."]=="y" | tt[i, "NKC_TB.sent."]=="tbd" ) lightpurple else if(tt[i, "NKC_TB.sent."]=="n") purple else brown
      if(cgannot) {
        pt <- if(cg.nl=="yes")  21  
          else if(cg.nl=="no") 22  
          else if(cg.nl=="unk") 23 
          else 13
      } else pt <- 21
      
      ## prior versions had * on samples that were "tbd", but ultimately samples were either sent or not sent.. depricated..
      ## prior versions had different border colors to distinguish based on the presence or absence of NB in the sample.
      ## left this logic for the ability to change as desired, but as written below, all borders are black
      if(is.na(is.nb)) {
        points(x=tt[i,"NKC_TB.TTS"], y=135, pch=pt, cex=cx, col=purple, bg=color)
        if(tt[i, "NKC_TB.sent."]=="tbd") points(x=tt[i,"NKC_TB.TTS"], y=155, pch="*", col="black")
      }
      else {
        border <- if(is.nb) "black" else if(!is.nb) "black" else "black"
        points(x=tt[i,"NKC_TB.TTS"], y=135, pch=pt, cex=cx, bg=color, col=border)
        if(tt[i, "NKC_TB.sent."]=="tbd") points(x=tt[i,"NKC_TB.TTS"], y=155, pch="*", col="black")
      }
    }
  }
}

plot.LFU <- function(tt,w=2) {
  for(i in 1:nrow(tt)) {
    if(!is.na(tt[i, "TT_LFU"])) {
      if(tt[i, "alive_LFU"]=='y') arr="curved" else if(tt[i, "alive_LFU"]=='n') arr="T" else arr="ellipse"
      Arrows(x0=tt[i,"TT_LFU"]-90, x1=tt[i,"TT_LFU"], y0=105, y1=105, code=2, 
             lty=1, lwd=3, arr.length=0.2, arr.width=0.3, arr.type=arr, col="gray50")
    }
  }
}

plot.AML <- function(tt, cx=2.7) {
  for(i in 1:nrow(tt)) {
    if(!is.na(tt[i, "TT_AML"])) {
      points(x=tt[i,"TT_AML"], y=105, pch='+', cex=cx) 
    }
  }
}

#################
# PLOT MOLECULAR VARIANTS FOR A GIVEN PATIENT
#
# `varlist` is a table containing the variant information for all patients, and includes (at least) the following columns:
##  PID (patient ID)
##  VAR_TYPE (sub, indel, fusion)
##  ANNOTATION (somatic | oncogenic, somatic | likely, somatic | unknown.  remove also acceptable but not in our annotations)
##  VARIANT (here: GENE_p.change, e.g. TP53_R175H, but any text will work)
##  sample.time..mo (sample time in months)
##  TARGET_VAF (VAF at given timepoint)
#   
#  maxtime is the maximum time (in days) needed for the individual patient (to set up plot area up front)
#  leg.cx cx for the legend of variants
#  mainfigs (defaults TRUE) prevents "somatic | unknown" variants from plotting (as in for main figure)
##   when mainfigs=FALSE, "somatic | unknown" variants plot with a dashed line
#
#################
plot.vars <- function(varlist, maxtime, leg.cx=1.5, mainfigs=T) { #to plot a list of one patient's variants
  maxt <- maxtime/30.42
  mos <- vector() 
  for(i in 0:ceiling(maxt/12)) {
    mos[i+1] <- i*12
  }
  if(mainfigs){
    varlist.gm <- varlist[varlist$VAR_TYPE %in% c("sub","indel") & !(varlist$ANNOTATION %in% c("remove","somatic | unknown")),]
  }
  else{
    varlist.gm <- varlist[varlist$VAR_TYPE %in% c("sub","indel") & !(varlist$ANNOTATION %in% c("remove")),]
  }
  vars <- unique(varlist.gm$VARIANT)
  pt <- varlist.gm$PID[1]
  par(mar=c(0,5,0,3))
  plot(x=c(0,maxt), y=c(0,0.55), type='n', xaxt='n', xaxs='i', yaxs='i', yaxt='n', ylab="VAF", xlab="", cex.lab=2, font=1)
  axis(2, at=seq(from=0,to=1,by=0.2), cex.axis=2, font=1)
  sapply(mos, FUN=function(x) {abline(v=x, col="grey90")})
  for(i in 1:length(vars)) {
    if(mainfigs){
      varlist.sub <- varlist.gm[varlist.gm$VARIANT==vars[i] & !(varlist.gm$ANNOTATION %in% c("remove","somatic | unknown")),]
    }
    else{
      varlist.sub <- varlist.gm[varlist.gm$VARIANT==vars[i] & !(varlist.gm$ANNOTATION %in% c("remove")),]
    }
    varvaf <- varlist.sub[,c("sample.time..mo.", "TARGET_VAF")]
    color=colors[i] 
    varvaf=varvaf[order(varvaf),]
        #unknowns plotted with dashed lines
    if(varlist.sub$ANNOTATION[1]=="somatic | unknown" | length(varlist.sub$ANNOTATION)==0){
        points(x=varvaf$sample.time..mo., y=varvaf$TARGET_VAF, col=color, lwd=2, cex=2,pch=16)
        lines(x=varvaf$sample.time..mo., y=varvaf$TARGET_VAF, col=color, lty=2, lwd=3,cex=2)
    }
    else{
      points(x=varvaf$sample.time..mo., y=varvaf$TARGET_VAF, col=color, lwd=2, cex=2,pch=16)
      lines(x=varvaf$sample.time..mo., y=varvaf$TARGET_VAF, col=color, lwd=5, cex=2)
    }
  }
  par(font=2)
  
  # pass through for lty to preserve dashed for unknowns
  if(mainfigs) {
    ly <- varlist.gm %>% distinct(VARIANT, .keep_all=TRUE) %>% filter(ANNOTATION != "somatic | unknown") %>% 
      select(ANNOTATION) %>% unlist %>% { ifelse(.=="somatic | unknown", 2, 1) }
  }
  else
    ly <- varlist.gm %>% distinct(VARIANT, .keep_all=TRUE) %>% select(ANNOTATION) %>% unlist() %>% { ifelse(.=="somatic | unknown", 2, 1) }
  if(length(vars)>0)	legend("topright", legend=vars, col=colors[1:length(vars)],lty=ly, lwd=4, cex=leg.cx, seg.len=2.5)
}


###################
# PLOT FUSIONS AND COPY NUMBER CHANGES - either classic/clinical CG or computationally derived ("facets") or Archer (fusions only)
#
# varlist table of "variants" (fusions, CNAs) - processed from raw file outside this function, see below
##  variant any/each fusion/CNA
##  sample_time in months
##  source = (clinical, facets, archer, both (>1))
##  CLONAL (T/F) - false if did not meet criteria for clonality by ISCN guidelines
## note: KMT2A rearrangements were noted as "11q23 rearr" and plotted first for consistency. 
# maxtime above (plot.vars)
# 
###################
plot.cg.all <- function(varlist, maxtime) {
  maxt <- maxtime/30.42
  mos <- vector() 
  for(i in 0:ceiling(maxt/12)) {
    mos[i+1] <- i*12
  }
  par(mar=c(3,5,0,3))
  plot(x=c(0,maxt), y=c(0,1), type='n', xaxt='n', xaxs='i', xlab="", yaxs='i', yaxt='n', ylab="CG", cex.lab=2, font=1)
  axis(1, at=mos, labels=mos, cex.axis=2, font=1)
  mtext("months", side=1, line=0, cex=1, outer=TRUE, font=1)
  sapply(mos, FUN=function(x) {abline(v=x, col="grey90")})
  if("11q23 rearr" %in% varlist$variant)  #always put 11q23 first if it is present, otherwise don't display it
    varlist$variant <- factor(varlist$variant, levels=unique(c("11q23 rearr",varlist$variant)))
  else
    varlist$variant <- factor(varlist$variant, levels=unique(varlist$variant))
  varlist <- varlist[order(varlist$variant),]
  vars <- levels(varlist$variant)
  for(i in 1:length(vars)) {
    if(length(vars)==0) next      #if no CNVs, move to next sample
    if(vars[i] %in% varlist$variant) {
      cg.plot <- as.data.frame(varlist[varlist$variant==vars[i],c("sample_time", "source", "CLONAL")])
      color <- NULL
                                                                                                #final color is for "both"
      adj.color <- sapply(cg.plot$source, FUN=function(x) { ifelse(x=="clinical","#F08080", ifelse(x %in% c("facets", "archer"),"#6495ED","#8B5F65"))})

          #if sample time is exactly 0, bump it forward just a bit to be visible from the beginning of the plot
      lft=ifelse(cg.plot$sample_time>0,(cg.plot$sample_time)-0.7,cg.plot$sample_time)
      rgt=ifelse(cg.plot$sample_time>0,(cg.plot$sample_time)+0.7,(cg.plot$sample_time)+2)
      top=rep(0.99-(i-1)/length(vars), nrow(cg.plot))
      bot=rep(1.01-i/length(vars), nrow(cg.plot))
      rect(xleft=lft, ybottom=bot, xright=rgt, ytop=top, col=adj.color,border=NA)
      for(j in 1:nrow(cg.plot)) {
        if(!is.na(cg.plot[j,"CLONAL"])) {
          if(!cg.plot[j,"CLONAL"]){
            points(x=(lft[j]+rgt[j])/2, y=(bot[j]+top[j])/2, pch=8, adj=0.5, cex=1.5)
          } 
        } 
      }
    }
    text(x=maxt, y=1-(i-0.5)/length(vars), pos=2, adj=0.5, labels=vars[i], cex=1.5)
    abline(h=i/length(vars), col="grey90")
  }
}


#######
# legend is missing shapes for vital status - need to add manually
#    alive = arrowhead
#    deceased = bar
######
timeline.legend <- function(horizontal=TRUE){
    par(mar=rep(0.5,4))
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(-1,1), ylim=c(-2,2))

    if(horizontal) {
        legend("center",legend=c("chemotherapy", "radiation",
                                 "immunotherapy", "BM sample","","","","clinical CG" , "molecular CNA or SV", "both", "non-clonal by clinical CG", "","","",
                                 "MDS/AML", "SC harvest", "auto transplant","allo transplant","last follow-up (alive)","last follow-up (deceased)"),
                cex=0.8, pch=c(rep(15, 3), 19,NA,NA,NA,rep(15, 3),8,NA,NA,NA, 3,rep(NA, 2),NA,NA,NA), pt.cex=1.9,
                lty=c(rep(NA, 15),2,2,1,NA,NA), lwd=1.5, ncol=3,
                col=c("lightcyan3","darkorange", "thistle4","lavenderblush3",NA,NA,NA,"#F08080","#6495ED","#8B5F65","black",NA,NA, NA,"black", "gray70", rep("black",4))
            )
    } else { #vertical
        legend("center",legend=c("chemotherapy", "radiation",
                             "immunotherapy", "BM sample","","clinical CG" , "molecular CNA or SV", "both", "non-clonal by clinical CG", "",
                             "MDS/AML", "SC harvest", "auto transplant","allo transplant","last follow-up (alive)","last follow-up (deceased)"),
                cex=0.9,
                pch=c(rep(15, 3), 19,NA,rep(15, 3),8,NA, 3,rep(NA, 2),NA,NA,NA), pt.cex=1.9,
                lty=c(rep(NA, 11),2,2,1,NA,NA), lwd=1.5, ncol=3,
                col=c("lightcyan3","darkorange", "thistle4","lavenderblush3",NA,"#F08080","#6495ED","#8B5F65","black", NA,"black", "gray70", rep("black",4))
        )
    }

}

# finalpts includes all cohorts (including controls)
finalpts <- c("I-H-107649", "I-H-118725", "I-H-118728", "I-H-118734", "I-H-118735", "I-H-118738", 
              "I-H-118741", "I-H-118742", "I-H-118726", "I-H-118727", "I-H-118729", "I-H-118730", 
              "I-H-118731", "I-H-118733", "I-H-118737", "I-H-118732", "I-H-118736", "I-H-107650", 
              "I-H-107651", "I-H-133427", "I-H-133429", "I-H-133430", "I-H-133431", "I-H-133432", 
              "I-H-133433", "I-H-133434", "I-H-133435", "I-H-133436", "I-H-133437", "I-H-133522", 
              "I-H-133523", "I-H-133882", "I-H-133883", "I-H-133884", "I-H-133885", "I-H-133886", 
              "I-H-133887", "I-H-133888", "I-H-133889", "I-H-133890", "I-H-133891", "I-H-133892", 
              "I-H-133893", "I-H-133894", "I-H-133895", "I-H-133896", "I-H-133897", "I-H-133898", 
              "I-H-133899", "I-H-133900", "I-H-133901", "I-H-133902")

ctrlpts <- c("I-H-133882", "I-H-133883", "I-H-133884", "I-H-133885", "I-H-133886", "I-H-133887", "I-H-133888", 
             "I-H-133889", "I-H-133890", "I-H-133891", "I-H-133892", "I-H-133893", "I-H-133894", "I-H-133895", 
             "I-H-133896", "I-H-133897", "I-H-133898", "I-H-133899", "I-H-133900", "I-H-133901", "I-H-133902")

# convert between original pt nums (1,2,3) and leukgen PIDs, extract patient categories (transformation, transient, control)
pid.ptnum <- read.table("../DATA/DATA_FILES/pt_num_PID_cat.txt", check.names=FALSE, stringsAsFactors=FALSE, header=TRUE, sep="\t")
rownames(pid.ptnum)<- pid.ptnum$pt_num
pid.ptnum[pid.ptnum$pt_cat=="other","pt_cat"] <- "transient"  # replace "other" with "transient"


#####
# varmat.sv and varmatpile are raw tables for structural variants and mutations respectively
# processing below to feed into plot.cg.all and plot.var respectively
#####
varmat.sv <- read.csv("../DATA/DATA_FILES/p173_variants_SVs_for_timeplots_CLEAN.txt", check.names=FALSE, sep="\t",stringsAsFactors=FALSE, header=TRUE)
varmatpile <- read.csv("../DATA/DATA_FILES/allvariantstopuppiled_CLEAN.tsv", sep="\t",stringsAsFactors=FALSE, header=TRUE)
varmatpile <- varmatpile[!varmatpile$VARIANT %in% c("TP53_R175G","WT1_R380G"),] #remove 1 variant (TP53) where VAF always <2% (maxvaf=0.0173) - added 08/03/20
                                                                               #and one where there are 2 calls for the same event (WT1 R380G and WT1 R380fs*5) 

# extract patients who had variants (i.e., were included in the varmat table), then remove those not included in the final set
#    (some samples sequenced were from patients not in this dataset)
pts <- intersect(union(varmat.sv$PID, varmatpile$PID), finalpts)
varmat.sv <- varmat.sv %>% filter(PID %in% pts)
varmatpile <- varmatpile %>% filter(PID %in% pts)

## extract the appropriate sample times from the original data list
maindf <- readRDS("../DATA/173_frozen.rds")
sampletimes=maindf %>% select(TARGET_NAME, sample.time..mo.) %>% distinct(TARGET_NAME, .keep_all=TRUE)
rownames(sampletimes) <- sampletimes$TARGET_NAME

varmat.sv$sample_time <- sampletimes[substr(varmat.sv$TARGET_NAME,1,13),"sample.time..mo."]

#processing FACETS findings (computational CNA) for varlist for plot.cg.all
varmatcna=readRDS("../DATA/FACETS_TIMELINES.rds")
varmatcna$sample_time <- sampletimes[varmatcna$sample,"sample.time..mo."]

varmatcna.sub <- varmatcna %>% filter(!(arm %in% c("err", "a"))) %>% rowwise() %>%
  mutate(sign=ifelse(type=="amp","+",ifelse(type=="del","-","LOH"))) %>% 
  mutate(variant=paste(chrom,arm,sign,sep="")) %>% mutate(patient2=paste0("I-H-",patient)) %>% filter(patient2 %in% finalpts)
 varmatcna.sub$source <- "facets"
 varmatcna.sub$CLONAL <- TRUE
 
colnames(varmatcna.sub)[c(5,9)] <- c("chr", "PID")

#processing clinical cytogenetic findings for varlist for plot.cg.all
varlist.cg <- varmat.sv[varmat.sv$Algorithm %in% c("CG_PUB") & (!varmat.sv$CG_REDUCED %in% c("","100% normal")) & 
                       varmat.sv$variant_type %in% c("aneuploidy","fusion"),]
varlist.cg$source <- "clinical"
varlist.cg=varlist.cg %>% merge(sampletimes, by='TARGET_NAME', all.x=TRUE) %>% unique()
varlist.cg <- varlist.cg[!substr(varlist.cg$CG_REDUCED,1,3)=="del",]  #remove to keep clean, not losing anything in either patient (del(ETV6), del(MLL))
varlist.cg.plot <- varlist.cg[,c("PID","sample_time", "chr", "arm", "CG_REDUCED", "source", "CLONAL")]
colnames(varlist.cg.plot) <- c("PID", "sample_time", "chr", "arm", "variant", "source", "CLONAL")

#processing Archer fusion findings for varlist for plot.cg.all
varlist.archer <- varmat.sv[varmat.sv$Algorithm %in% c("Archer") & (!varmat.sv$CG_REDUCED %in% c("","100% normal")),]
varlist.archer$source <- "archer"
varlist.archer$CLONAL <- TRUE

cgvarlist.all <- rbind(varmatcna.sub[,c(PID="PID", "sample_time", "chr", "arm", "variant", "source", "CLONAL")], 
                       varlist.cg.plot[,c("PID", "sample_time", "chr", "arm", "variant", "source", "CLONAL")],
                       varlist.archer[,c("PID", "sample_time", "chr", "arm", "variant", "source", "CLONAL")])
cgvarlist.all$chr <- factor(cgvarlist.all$chr, 
                        levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
cgvarlist.all <- cgvarlist.all[cgvarlist.all$PID %in% finalpts,] %>% arrange(PID, chr, arm, variant, sample_time, source, CLONAL)

#modify source to be 'both' if found by both clin and (cnacs or archer)
srces=cgvarlist.all %>% group_by(PID,chr, arm, sample_time, variant) %>% 
  mutate(n=length(unique(source)),source=unique(source)[1]) %>%
  mutate(source=ifelse(n==1,source,"both"))


tt<- read.csv("../DATA/DATA_FILES/therapy_timetable_CLEAN.csv", header=TRUE)
tt$PID <- pid.ptnum[as.character(tt$patient), "PID"]
tt <- tt[!is.na(tt$PID),]

tt <- tt[tt$PID %in% finalpts,]

#pts are those with any abnormalities to plot (variants, CG, etc - mostly transformation, transient, and the 1 control pt)
#intersect with finalpts to order appropriately - first transformation, transient, control
pts <- intersect(finalpts, pts)

# legend separately to keep the size reasonable
pdf("../FIGURES/timeline_legend.pdf", height=1.7, width=7.4)
timeline.legend()
dev.off()

pdf("../FIGURES/TIMELINE.pdf", width=15, height=6)

par(pin=c(8,5.5))
par(oma=c(0,0,0,0))
layout(matrix(c(1,2,3), nrow=3, ncol=1, byrow=TRUE))

for(i in 1:length(pts)){
  cat("pt=",pts[i],"\t")
  tt.sub <- tt[tt$PID==pts[i],]
  maxti <- max(tt.sub[,"NKC_TB.TTS"], tt.sub[1,"TT_LFU"])+365
  cat(maxti,"\n")
  mons <- vector()
  for(j in 0:ceiling(maxti/12/30.42)) {
    mons[j+1] <- j*12*30.42
  }
  par(mar=c(0,5,2,3))
  plot(x=c(0,maxti), y=c(0,160), type='n', xaxt='n', xaxs='i', yaxs='i', yaxt='n', ylab="", xlab="", cex.axis=2, cex.main=1.9, main=paste0(pts[i], "_", pid.ptnum[pid.ptnum$PID==pts[i], "pt_cat"]))
  sapply(mons, FUN=function(x) {abline(v=x, col="grey90")})
  
  plot.therapy(tt.sub,wd=24)
  plot.SCharvest(tt.sub)
  plot.auto(tt.sub)
  plot.allo(tt.sub)
  plot.NKCTB(tt.sub, allsamps=FALSE)
  plot.LFU(tt.sub, w=3)
  plot.AML(tt.sub, cx=2.5)
  plot.vars(varmatpile[varmatpile$PID==pts[i],], maxti)
  plot.cg.all(srces[srces$PID==pts[i],], maxti)

  }

dev.off()
