# Genome Annotation using Segway

###working using wig file 

#GSM701525_UW.CD3_Primary_Cells.ChromatinAccessibility.CB1_1-3-2011.DS17702.wig.gz
#GSM1058764_UW.CD3_Primary_Cells.H3K27ac.RO_01701.Histone.DS21704.wig.gz
#GSM537613_BI.CD3_Primary_Cells.H3K27me3.CD3_39661.wig.gz
#GSM537666_BI.CD3_Primary_Cells.H3K36me3.CD3_39804.wig.gz
#GSM537636_BI.CD3_Primary_Cells.H3K4me1.CD3_39811.wig.gz
#GSM537633_BI.CD3_Primary_Cells.H3K4me3.CD3_39804.wig.gz
#GSM537645_BI.CD3_Primary_Cells.H3K9me3.CD3_39840.wig.gz

# transform wig to bigwig using wigToBigWig [http://hgdownload.cse.ucsc.edu/admin/exe/]


#hg19.chrom.sizes matrix containing each chromosome and its size in bp

#chr1    249250621
#chr2    243199373
#chr3    198022430
#chr4    191154276
#chr5    180915260

wigToBigWig GSM701525_UW.CD3_Primary_Cells.ChromatinAccessibility.CB1_1-3-2011.DS17702.wig.gz hg19.chrom.sizes DNASE.bw

wigToBigWig GSM1058764_UW.CD3_Primary_Cells.H3K27ac.RO_01701.Histone.DS21704.wig.gz hg19.chrom.sizes H3K27ac.bw

wigToBigWig GSM537613_BI.CD3_Primary_Cells.H3K27me3.CD3_39661.wig.gz hg19.chrom.sizes H3K27me32.bw

wigToBigWig GSM537666_BI.CD3_Primary_Cells.H3K36me3.CD3_39804.wig.gz hg19.chrom.sizes H3K36me3.bw


wigToBigWig GSM537636_BI.CD3_Primary_Cells.H3K4me1.CD3_39811.wig.gz hg19.chrom.sizes H3K4me1.bw

wigToBigWig GSM537633_BI.CD3_Primary_Cells.H3K4me3.CD3_39804.wig.gz hg19.chrom.sizes H3K4me3.bw

wigToBigWig GSM537645_BI.CD3_Primary_Cells.H3K9me3.CD3_39840.wig.gz hg19.chrom.sizes H3K9me3.bw

###bigwig to bedgraph 
bigWigToBedGraph DNASE.bw DNASE.bedGraph

bigWigToBedGraph H3K27ac.bw H3K27ac.bedGraph

bigWigToBedGraph H3K27me3.bw H3K27me3.bedGraph

bigWigToBedGraph H3K36me3.bw H3K36me3.bedGraph

bigWigToBedGraph H3K4me1.bw H3K4me1.bedGraph

bigWigToBedGraph H3K4me3.bw H3K4me3.bedGraph

bigWigToBedGraph H3K9me3.bw H3K9me3.bedGraph


####100bp windows

module load bedtools2/2.24.0/gcc.4.4.7

#hg19_100bpWin.bed file with 100 bp windows
# bedtools makewindows -g hg19.txt -w 100
#chr1    1       101
#chr1    101     201
#chr1    201     301
#chr1    301     401
#chr1    401     501
#chr1    501     601
#chr1    601     701
#chr1    701     801
#chr1    801     901
#chr1    901     1001
#chr1    1001    1101
#chr1    1101    1201
#chr1    1201    1301
#chr1    1301    1401
#chr1    1401    1501
#chr1    1501    1601
#chr1    1601    1701

bedtools map -a hg19_100bpWin.bed -b DNASE.bedGraph -c 4 -o mean  -null "0" > DNASE_100.bedGraph

bedtools map -a hg19_100bpWin.bed -b H3K27.bedGraph -c 4 -o mean  -null "0" > H3K27_100.bedGraph

bedtools map -a hg19_100bpWin.bed -b H3K27me3.bedGraph -c 4 -o mean  -null "0" > H3K27me3_100.bedGraph

bedtools map -a hg19_100bpWin.bed -b H3K36me3.bedGraph -c 4 -o mean  -null "0" > H3K36me3_100.bedGraph

bedtools map -a hg19_100bpWin.bed -b H3K4me1.bedGraph -c 4 -o mean  -null "0" > H3K4me1_100.bedGraph

bedtools map -a hg19_100bpWin.bed -b H3K4me3.bedGraph -c 4 -o mean  -null "0" > H3K4me3_100.bedGraph

bedtools map -a hg19_100bpWin.bed -b H3K9me3.bedGraph -c 4 -o mean  -null "0" > H3K9me3_100.bedGraph

###segway

module load segway/1.3/python.2.7.10 


genomedata-load-assembly genome.data hg19.fa.gz
genomedata-open-data genome.data H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 DNASE1
genomedata-load-data genome.data H3K27ac < H3K27ac_100.bedGraph
genomedata-load-data genome.data H3K27me3 < H3K27me3_100.bedGraph 
genomedata-load-data genome.data H3K36me3 < H3K36me3_100.bedGraph
genomedata-load-data genome.data H3K4me1 < H3K4me1_100.bedGraph
genomedata-load-data genome.data H3K4me3 < H3K4me3_100.bedGraph
genomedata-load-data genome.data H3K9me3 < H3K9me3_100.bedGraph
genomedata-load-data genome.data DNASE1 < DNASE_100.bedGraph
genomedata-close-data genome.data


genomedata=genome.data 
traindir=traindir
identifydir=identifydir 

##train in chr1

#label_lengths_500bpmin_500bpruler.txt to specify min and max length of fragment

#label   len
#0:      500::500  "0:" for all labels ; 500 min size of fragment ; "::" ; no max size : "500" ruler length


segway --num-labels=7 train --num-instances=1 --include-coords=chr1.hg19.bed --ruler-scale=500 --resolution=500 --seg-table=label_lengths_500bpmin_500bpruler.txt $genomedata $traindir 

segway identify $genomedata $traindir $identifydir 

###this create a bed files in identifydir with the annotation of each feature. Each feature represent a mixture of the different histone mark uploaded to define each feature you can use a ratio observed/expected for the different marks by feature. 


### example 


#in R

options(stringsAsFactors=F)

#function to use bedtools in R 

bedTools.2in<-function(functionstring="intersectBed",bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
 
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
 
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}

H3K4me1=H3K4me1_020316.bedGraph
H3K4me3=H3K4me3_020316.bedGraph
H3K9me3=H3K9me3_020316.bedGraph
H3K36me3=H3K36me3_020316.bedGraph
H3K27me3=H3K27me3_020316.bedGraph
H3K27ac=H3K27ac_020316.bedGraph
DNASE=DNASE_020316.bedGraph

FEATURE=segway.bed

y=H3K4me1[,3]-H3K4me1[,2]
H3K4me1=data.frame(H3K4me1,y)
IH3K4me1<-bedTools.2in(bed1=H3K4me1, bed2=FEATURE,opt.string="-wo") #feature overlapping with H3K4me1

split_data<-split(IH3K4me1,(IH3K4me1[,9]))

output0=split_data[[1]][,10]
output1=split_data[[2]][,10]
output2=split_data[[3]][,10]
output3=split_data[[4]][,10]
output4=split_data[[5]][,10]
output5=split_data[[6]][,10]
output6=split_data[[7]][,10]

x=FEATURE[,3]-FEATURE[,2]
FEATURE=data.frame(FEATURE,x)

#number of bp represented by each feature divided by total number of bp in bed file


feat0=sum(as.numeric(FEATURE[which(FEATURE[,4]==0),5]))*100/sum(as.numeric(FEATURE[,5]))
feat1=sum(as.numeric(FEATURE[which(FEATURE[,4]==1),5]))*100/sum(as.numeric(FEATURE[,5]))
feat2=sum(as.numeric(FEATURE[which(FEATURE[,4]==2),5]))*100/sum(as.numeric(FEATURE[,5]))
feat3=sum(as.numeric(FEATURE[which(FEATURE[,4]==3),5]))*100/sum(as.numeric(FEATURE[,5]))
feat4=sum(as.numeric(FEATURE[which(FEATURE[,4]==4),5]))*100/sum(as.numeric(FEATURE[,5]))
feat5=sum(as.numeric(FEATURE[which(FEATURE[,4]==5),5]))*100/sum(as.numeric(FEATURE[,5]))
feat6=sum(as.numeric(FEATURE[which(FEATURE[,4]==6),5]))*100/sum(as.numeric(FEATURE[,5]))


#number of bp overlapping with H3K4me1 by feature

exp0=sum(as.numeric(IH3K4me1[which(IH3K4me1[,4]!=0),10]))*feat0/100 
exp1=sum(as.numeric(IH3K4me1[which(IH3K4me1[,4]!=0),10]))*feat1/100
exp2=sum(as.numeric(IH3K4me1[which(IH3K4me1[,4]!=0),10]))*feat2/100
exp3=sum(as.numeric(IH3K4me1[which(IH3K4me1[,4]!=0),10]))*feat3/100
exp4=sum(as.numeric(IH3K4me1[which(IH3K4me1[,4]!=0),10]))*feat4/100
exp5=sum(as.numeric(IH3K4me1[which(IH3K4me1[,4]!=0),10]))*feat5/100
exp6=sum(as.numeric(IH3K4me1[which(IH3K4me1[,4]!=0),10]))*feat6/100


obs0=sum(as.numeric(output0))
obs1=sum(as.numeric(output1))
obs2=sum(as.numeric(output2))
obs3=sum(as.numeric(output3))
obs4=sum(as.numeric(output4))
obs5=sum(as.numeric(output5))
obs6=sum(as.numeric(output6))



ratio0=obs0/exp0
ratio1=obs1/exp1
ratio2=obs2/exp2
ratio3=obs3/exp3
ratio4=obs4/exp4
ratio5=obs5/exp5
ratio6=obs6/exp6


final=c("H3K4me1",ratio0,ratio1,ratio2,ratio3,ratio4,ratio5,ratio6)
colnames(final)=c("histone","feature0","feature1","feature2","feature3","feature4","feature5","feature6")


##repeat with each mark
