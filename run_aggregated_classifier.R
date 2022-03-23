#!/usr/bin/Rscript
library('reshape2')
library('plyr')
library('optparse')
source('utils.R')

options<-list(
  make_option('--scores',action='store',type='character',help="RDS file containg the classfier's input data matrix"),
  make_option('--genelist',action='store',type='character',
              help='.txt file with list of gold standard genes, each gene in a sepearte row'),
  make_option('--outfile',action='store',type='character',
              help='Name of the output file. Note that te output is in csv format'),
  make_option('--fpr_threshold',action='store',type='numeric',default=0.02,
              help="False positive rate for inclusion in the the final gene list.")
)

parser<-OptionParser(option_list=options)
args<-parse_args(parser)

path_to_scores<-args$scores
path_to_genelist<-args$genelist
outfile<-args$outfile
fpr_threshold<-args$fpr_threshold

scores<-readRDS(file=path_to_scores)
gs<-readLines(path_to_genelist)

#Re-arrange the classifier output
scoreslong<-melt(scores,id.vars = 'Gene',variable.name = 'Feature')

is_binary<-grepl('OMIM|Peng|Slabicki|Adamson|Matsukoa',colnames(scores))
screenstr<-c('OMIM|Peng|Slabicki|Adamson|Matsukoa|RAD51.foci|Paulsen|Lord|Bartz|Elia')
isscrreen<-grepl(screenstr,scoreslong$Feature)

scoreslongNonScreens<-scoreslong[!isscrreen,]
scoreslongNonScreens$Feature<-as.character(scoreslongNonScreens$Feature)
featuresplitpoint<-regexpr('\\.',scoreslongNonScreens$Feature)
scoreslongNonScreens$Feature_gene<-substring(scoreslongNonScreens$Feature,1,featuresplitpoint-1)
scoreslongNonScreens$Feature_type<-substring(scoreslongNonScreens$Feature,featuresplitpoint+1)
scoreslongNonScreens$is_self<-apply(scoreslongNonScreens,1,function(x)
  x['Gene']==x['Feature_gene'])
scoreslongNonScreens$value[which(scoreslongNonScreens$is_self==TRUE)]<-0
scoreslongNonScreens$Feature_type<-factor(scoreslongNonScreens$Feature_type)
scoreslongNonScreens$Feature_type<-mapvalues(scoreslongNonScreens$Feature_type,
                                             from=c('CCI.Malovannaya.2011','Genehopper.biological.process',
                                                    'Genehopper.cellular.component','Genehopper.domains',
                                                    'Genehopper.molecular.function',
                                                    'STRING.database','STRING.experimental','STRING.textmining'),
                                             to=c('Complex-Complex Interactions','GO.biological.process',
                                                  'GO.cellular.component','Shared domains','GO.molecular.function',
                                                  'Shared pathways','Protein-Protein Interaction','Literature'))
#Sum evidence score pre feature type for interaction terrms 
nonScreenAggScore<-tapply(scoreslongNonScreens$value,
                          list(scoreslongNonScreens$Gene,scoreslongNonScreens$Feature_type),sum)

scorelongScreens<-scoreslong[isscrreen,]
ScreenScoreWide<-dcast(scorelongScreens,Gene~Feature,value.var = 'value')
rownames(ScreenScoreWide)<-ScreenScoreWide$Gene
ScreenScoreWide<-ScreenScoreWide[,-1]
allScores<-cbind(nonScreenAggScore,ScreenScoreWide)

#Run classifier using all columns in the data matrix
LRlist<-apply(allScores,2,function(x) LR(x,gs))
#Calculate per-gene classifier score
preds<-getPreds(LRlist,gs)
saturation<-make_saturation_table(preds)


#Run classifier withou litetature and annotations
LRlistnolit<-apply(allScores[,!grepl('GO|Lit|path|OMIM',colnames(allScores))],2,function(x) LR(x,gs))
predsolit<-getPreds(LRlistnolit,gs)
saturationnolit<-make_saturation_table(predsolit)


#Merge genelists from preds (classifier with all data) and predsolt(classifier without literture)

sat_merge<-merge(saturation,saturationnolit,by='gene')
sat_merge_minrank<-ddply(sat_merge,.(gene),mutate,minRank=min(Rank.x,Rank.y),Source=ifelse(Rank.x<=Rank.y,'all_data','without_literature'))
#Define for each gene its rank in its merged list as the minimum of the rank assigned to
#it the full data classifier and the classifier without literature
sat_merge_minrank<-sat_merge_minrank[order(sat_merge_minrank$minRank),]
sat_clean<-sat_merge_minrank[,c('gene','minRank','Source')]
#Create the final satutation table
sat_clean$Rank<-1:nrow(sat_clean)
sat_clean$gold.standard<-sat_clean$gene%in%gs
sat_clean$gs_count<-cumsum(sat_clean$gold.standard)
sat_clean$TPR<-sat_clean$gs_count/sum(sat_clean$gold.standard)
sat_clean$FPR<-(sat_clean$Rank-sat_clean$gs_count)/(nrow(sat_clean)-sum(sat_clean$gold.standard))
sat_clean$Recall<-100*sat_clean$TPR
sat_clean$gs_ratio<-sat_clean$Rank/sat_clean$gs_count
#THE NEXT  LINE DETERMINES THE FINAL GENE LIST
sat_clean$Call<-sat_clean$FPR<=fpr_threshold
hitlist<-sat_clean$gene[sat_clean$Call]

#Arrange output

#merge screens represented by two columns (due to two effect directions) into one column by taking the maximum of the two columns
preds$Peng.DE<-apply(preds[,c('Peng2014.Downregulated','Peng2014.Upregulated')],1,max)
preds$RAD51.foci<-apply(preds[,c('RAD51.foci.impaired','RAD51.foci.retained')],1,max)
preds$H2AX.phosphrylation<-apply(preds[,c('Paulsen.2009.gamma.HA2X.increased','Paulsen.2009.gamma.HA2X.decreased')],1,max)

#Remove the  excess columns for these screens
preds<-preds[,!grepl('increased|decreased|retained|impaired|Downregulated|Upregulated',colnames(preds))]
preds$Gene<-rownames(preds)


#Rename and re-order columns
new_labels<-c('Complex-Complex Interactions','Clade.PP',
              'GO biological process','GO cellular component','Shared domains',
              'GO molecular function','Shared pathways','Protein-Protein Interactions',
              'Text mining','Coexpression','OMIM',
              'HRR siRNA screen (Slabicki 2010)',
              'HRR siRNA screen (Adamson 2012)',
              'DSB induced phosphorylation (Matsukoa.2007)',
              'DSB induced acetylation',
              'DSB induced ubiquitination',
              'DSB induced ubiquitination - protesome inihibitor',
              'DSB induced phosphorylation (Elia 2015)',
              'DSB induced ubiquitination - nuclear extracts',
              'PARP inhibitors sensitivity',
              'Cisplatin sensitivity', 'Total classifier score','in HRR gold standard',
              'Differential expression',
              'RAD51 foci formation' ,'H2AX phosphorylation','Gene')
chk_lab<-cbind(colnames(preds),new_labels)
colnames(preds)<-new_labels
intfeatures<-c('Clade.PP','Coexpression','Protein-Protein Interactions','Complex-Complex Interactions',
               'GO biological process','GO molecular function','GO cellular component','Shared pathways','Text mining',
               'Shared domains')
phenofeatures<-c('HRR siRNA screen (Slabicki 2010)',
                'HRR siRNA screen (Adamson 2012)',
                 'RAD51 foci formation',
                 'H2AX phosphorylation','Cisplatin sensitivity','PARP inhibitors sensitivity','OMIM')
activefeatures<-c('DSB induced phosphorylation (Matsukoa.2007)',
                  'DSB induced phosphorylation (Elia 2015)',
                  'DSB induced acetylation',
                  'DSB induced ubiquitination',
                  'DSB induced ubiquitination - nuclear extracts',
                  'DSB induced ubiquitination - protesome inihibitor',
                  'Differential expression')

new_lab_ord<-match(c('Gene','in HRR gold standard','Total classifier score',intfeatures,phenofeatures,activefeatures),
                   colnames(preds))

#Create output table
logpreds<-preds[,new_lab_ord]
#Transform raw classifier score into log scale. 
#The classifier emits log likelihood ratio for each gene, so by summing the logs
#from the different features we obtain the total score
logpreds[,seq(4,ncol(logpreds))]<-signif(log2(logpreds[,seq(4,ncol(logpreds))]),digits = 3)
logpreds$`Total classifier score`<-signif(logpreds$`Total classifier score`,digits = 3)


logpreds_final<-logpreds[hitlist,]
logpreds_final<-logpreds_final[order(-logpreds_final$`Total classifier score`),]

write.csv(logpreds_final,file=outfile,row.names = FALSE,quote = FALSE)

#Write log
outlog<-gsub('\\.csv','.log',outfile)
cat(sprintf('Run classifeir on %s with the following inputs:\n',Sys.time()),file=outlog,sep = '\n')
cat(sprintf('Classifier input matrix: %s\n',path_to_scores),file=outlog,sep='\n',append = TRUE)
cat(sprintf('False positive rate threshold for inclusion in outpput gene list : %s\n',fpr_threshold),file=outlog,sep='\n',append = TRUE)

cat('Gold standard genes: \n',sep='\n',append = TRUE,file = outlog)
cat(gs,append = TRUE,file = outlog,sep='\n')
cat(sprintf('Output location: %s',outfile),append = TRUE,file = outlog)

#Export classifier output for all genens (not just the ones that exceeded the FPR thershold)
#useful for inspection, QC etc.
outfull<-gsub('\\.csv','_all_genes.csv',outfile)
write.csv(logpreds,file=outfull,row.names = FALSE,quote = FALSE)
