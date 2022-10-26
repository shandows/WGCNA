#######此版本适用于组织的表达量的WGCNA
rm(list=ls())
library(WGCNA)
library(reshape2)
library(limma)
library(stringr)
options(stringsAsFactors = FALSE)
exprMat <- "fpkm1.CSV"
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
####################导入数据###############
dataExpr <- read.csv(exprMat,row.names = 1, sep=',',header=T, 
                     quote="", comment="", check.names=F)
index <- rownames(dataExpr)
dataExpr <- as.data.frame(lapply(dataExpr,as.numeric))
sum(is.na(dataExpr))
dataExpr[is.na(dataExpr)] <-0
row.names(dataExpr) <- index
#y2 <- dataExpr
###########消除批次效应##########
batch <- c(rep("DAP8",3),rep("DAP24",3),rep("DAP8",3),rep("DAP24",3),rep("DAP8",3),rep("DAP24",3),rep("DAP8",3),rep("DAP24",3))
y2 <- removeBatchEffect(dataExpr,batch = batch)

####################筛选基因和样本############
#用mad筛选，取前75%
m.mad = apply(y2,1,mad)
dataExpr <- y2[which(m.mad >
                       max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
datExpr <- as.data.frame(t(dataExpr))
###筛选掉异常的样本
A=adjacency(t(datExpr),type="signed") 
k=as.numeric(apply(A,2,sum))-1  
Z.k=scale(k)  
thresholdZ.k=- 2.5  
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black") 
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datExpr=datExpr[!remove.samples,]

nGenes=ncol(datExpr)
nSamples=nrow(datExpr)

sampleTree = hclust(dist(datExpr),method = "average")
sizeGrwindow(12,9)

par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree,main="Sample clustering",sub="",xlab="",cex.lab=1.5,cex.axis=1.5,cex.main=2)

####################获取处理的文件##################
traitData=read.csv("characterD.csv",header = T)

Sample=rownames(datExpr)
traitRows=match(Sample,traitData$Sample)
datTraits=traitData[traitRows,-1]
rownames(datTraits)=traitData[traitRows,1]
collectGarbage()
traitColors= numbers2colors(datTraits,signed=FALSE)
pdf("sample.pdf",width = 25,height = 12)
plotDendroAndColors(sampleTree,traitColors,groupLabels = names(datTraits),main="Sample dendrogram and treat heatmap")
dev.off()

save(datExpr,traitData,datTraits,sampleTree,traitColors,file = "01-data_sample_L.RData")


##筛选软阈值
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector=powers, networkType=type, verbose=5)

pdf('power.pdf',width = 14,height = 9)
par(mfrow = c(1,2))
cex1 = 0.88
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.88,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
dev.off()

power=sft$powerEstimate

if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                     ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                           ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                   ifelse(type == "unsigned", 6, 12))       
                     )
  )
}

##网络构建

cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("TOM.RData"),
                       verbose = 3)

table(net$colors)
cor <- stats::cor
pdf('net.pdf',width = 14,height = 9)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

MEs=net$MEs
geneTree=net$dendrograms[[1]];
save(MEs,moduleLabels,power,moduleColors,geneTree,file="02-network_L-auto.RData")

nGenes=ncol(datExpr)
nSamples=nrow(datExpr)
#datTraits[45,"period"] <- 0
MEs=moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs=orderMEs(MEs)

moduleTraitCor=cor(MEs,datTraits,use = "p")

moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

sizeGrWindow(10,6) #绘图框大小
textMatrix=paste(signif(moduleTraitCor,2)," (",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)


pdf("Module-trait relationships.pdf",width = 20,height = 10)
par(mar=c(9,12,3,3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               xColorWidth = 1,
               setStdMargins = FALSE,#默认参数出图
               #yColorWidth=2,
               
               textMatrix = textMatrix,
               cex.text = 0.75,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()

############ 寻找hubgene ###########

salt=as.data.frame(datTraits$L27_Kernel_DAP8)
names(salt)="L27_Kernel_DAP8" #对CM的列名进行更换，上一步操作后的列名为“datTrait$CM”

modNames=substring(names(MEs),3)#获取模块颜色名称，从ME后开始取
geneModuleMembership=as.data.frame(cor(datExpr,MEs,use = "p"))
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
names(geneModuleMembership)=paste("MM",modNames,sep = "")
names(MMPvalue)=paste("p.MM",modNames,sep="")
geneTraitSignificance=as.data.frame(cor(datExpr,salt,use = "p"))
GSPvalue=as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance)=paste("GS.",names(salt),sep = "")
names(GSPvalue)=paste("p.GS.",names(salt),sep="")

#########输出各模块基因及相关结果#########
probes=colnames(datExpr)
geneInfo = data.frame(geneID=probes,
                      moduleColor=moduleColors,
                      geneTraitSignificance,
                      GSPvalue)
modOrder = order(-abs(cor(MEs,salt,use = "p")))
for(mod in 1:ncol(geneModuleMembership)){
  oldnames=names(geneInfo)
  geneInfo=data.frame(geneInfo,geneModuleMembership[,modOrder[mod]],MMPvalue[,modOrder[mod]])
  names(geneInfo)=c(oldnames,paste("MM.",modNames[modOrder[mod]],sep = ""),paste("p.MM.",modNames[modOrder[mod]],sep = ""))
  print(mod)
}
geneOrder = order(geneInfo$moduleColor,-abs(geneInfo$GS.L27_Kernel_DAP8))
geneInfo=geneInfo[geneOrder,]
write.csv(geneInfo,file = "geneInfo_L27_24.csv")

##################################观察具体模块与处理的相关性##########################
#########循环绘图####
H22_8=list("darkorange","pink","blue","violet","darkolivegreen","midnightblue")
H27_8=list("lightsteelblue1","plum1","black","paleturquoise","brown")
L22_8=list("darkorange2","salmon","green","saddlebrown")
L27_8=list("white","pink","blue","green","saddlebrown","black","yellow")
H22_24=list("orangered4")
H27_24=list("greenyellow","ivory")
L22_24=list("darkgreen","royalblue","cyan","lightgreen","skyblue3","ivory","plum2","red")
L27_24=list("ivory","lightcyan1","lightcyan","darkturquoise")
for (module in L27_24){
  column=match(module,modNames)
  moduleGenes=moduleColors==module
  
  sizeGrWindow(7,7)
  par(mfrow=c(1,1))
  pdf(paste("L27_24_",module,".pdf",sep = ""),width = 10,height = 6)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                     abs(geneTraitSignificance[moduleGenes,1]),
                     xlab = paste("Module Membership in",module,"module"),
                     ylab = "Gene significance for L27_24",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2,
                     cex.lab = 1.2,
                     cex.axis = 1.2,
                     col = module)
  dev.off()
}


####################分步###########
module="white"
column=match(module,modNames)
moduleGenes=moduleColors==module

sizeGrWindow(7,7)
par(mfrow=c(1,1))
pdf("L27_8_white.pdf",width = 10,height = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in",module,"module"),
                   ylab = "Gene significance for L27_8",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2,
                   cex.lab = 1.2,
                   cex.axis = 1.2,
                   col = "darkgrey")
dev.off()

####################################################

probes = names(datExpr)
cyt=exportNetworkToCytoscape(TOM,
                               edgeFile = "CytoscapeInput-edges2.txt",# paste("CytoscapeInput-edges-",paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = "CytoscapeInput-nodes2.txt",#paste("CytoscapeInput-nodes-",paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = probes,

                               nodeAttr = moduleColors)


