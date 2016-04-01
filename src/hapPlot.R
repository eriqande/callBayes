
FOLDPATH=/home/biopipe/Genetics_Lab_Data/GTseq/satrovirens01092016/map_se/
for i in {1..144}; do perl src/hap_recap.pl -v ${FOLDPATH}/satro144_panel1_flashed_noMNP_noComplex.vcf -s ${FOLDPATH}/satro_flashed_s${i}_aln.sam -i $i > data/sum/satro_panel1_${i}.summary& done

cat data/sum/satro_panel1_*.summary > data/sum/all_satro_panel1.summary

library("ggplot2")
library("plyr")
library("dplyr")

haplo.sum <- read.table("data/sum/all_satro_panel1.summary", stringsAsFactors = FALSE) %>% 
  tbl_df

colnames(haplo.sum) <- c("id", "locus", "haplo", "depth", "logP.call", "logP.miscall")


# clean the file: remove any haplotype with alignment character ("*" and "_")  &
# remove any loci that does less than half of total individual
# post remove any loci with highly variant haplotypes > 40

num.id <- length(unique(haplo.sum$id))

haplo.cleanup <- haplo.sum %>% 
  filter(!grepl("[_|*]", haplo)) %>%
  group_by(locus, id) %>%
  mutate(n.haplo.per.indiv=n()) %>%
  ungroup() %>%
  group_by(locus) %>%
  mutate(n.indiv.per.locus = length(unique(id)), max.uniq.hapl=max(n.haplo.per.indiv)) %>%
  ungroup() %>%
  filter(n.indiv.per.locus > num.id/2, max.uniq.hapl < 40)  %>%
  select(id, locus, haplo, depth, logP.call, logP.miscall)


haplo.add.balance <- haplo.cleanup %>% 
  group_by(locus,id) %>%
  arrange(-depth) %>%
  mutate(allele.balance = depth/depth[1], rank=row_number() ) %>%
  ungroup() 

saveRDS(haplo.add.balance, "hapPlot/satrovirens01092016_panel1_haplo_filter.rds", compress="xz")
saveRDS(haplo.add.balance, "hapPlot/satrovirens02102016_panel2_haplo_filter.rds")


#write.table(haplo.add.balance, "hapPlot/satrovirens02102016_haplo_filter.tbl", sep="\t", quote=FALSE, row.names=F, col.names=F)












# box plot to show the number of haplotype

haplo.cutoff <- haplo.sum %>%
  group_by(locus, id) %>%
  summarise(hapl.three.pl.st = ifelse(length(depth) > 2, 0, 0),
            hapl.three.pl.end = ifelse(length(depth) > 2, sort(depth, decr=T)[3]-1, 0),
            hapl.one.st = ifelse(sum(depth==max(depth))==1 && length(depth) > 1, sort(depth, decr=T)[2], 0),
            hapl.one.end = ifelse(sum(depth==max(depth))==1, max(depth),0))
            

haplo.sample <- haplo.cutoff %>% filter(locus== "tag_id_1511")
ggplot()+ 
  geom_rect(data=haplo.sample, aes(xmin = hapl.one.st, xmax = hapl.one.end, ymin = id, ymax = id+1, fill="1"))+
  geom_rect(data=haplo.sample, aes(xmin = hapl.three.pl.end, xmax = hapl.one.st, ymin = id, ymax = id+1,fill="2") )+
  geom_rect(data=haplo.sample, aes(xmin = hapl.three.pl.st, xmax = hapl.three.pl.end, ymin = id, ymax = id+1, fill="3+"))+
  scale_x_log10()+
  theme_bw()+
  xlab("read coverage cutoff")+
  ylab("Individual ID")+
  scale_fill_manual(name= "Haplotypes:", values=c("1"="light grey","2"= "#4BBA82", "3+"="#A48A82"))+
  theme(legend.position="bottom")
  
#"#F9707D"

ggplot()+ 
  geom_segment(data=haplo.sample, aes(x = hapl.one.st, xend = hapl.one.end, y = id, yend = id, colour= "1"), size=3 )+
  geom_segment(data=haplo.sample, aes(x = hapl.three.pl.end, xend = hapl.one.st, y = id, yend = id, colour="2"), size=3 )+
  geom_segment(data=haplo.sample, aes(x = hapl.three.pl.st, xend = hapl.three.pl.end, y = id, yend = id, colour="3+"), size=2)+
  scale_x_log10()+
  theme_bw()+
  xlab("read coverage cutoff")+
  ylab("Individual ID")+
  scale_color_manual(name= "Haplotypes:", values=c("1"="light grey","2"= "#4BBA82", "3+"="#A48A82"))+
  theme(legend.position="bottom")


ggplot(data=haplo.sample, aes(x=hapl.one.st/hapl.one.end, y=id,size=log(hapl.one.end, 10)))+
  geom_point()+
  scale_size_continuous("Read Depth of \nthe most common \nhaplotype(log 10)")+
  xlab ("Depth Ratio of the second common hapl to first ")



# prototype for plotting the total number of haplotype for each locus
haplo.ct <- haplo.sum %>% 
  filter(depth > 2) %>%
  group_by(locus, id) %>% 
  summarise(tot.hapl = n()) 

haplo.tot.tbl <- haplo.ct %>% 
  group_by(locus, tot.hapl) %>%
  summarise(ct = n()) %>%
  ungroup() %>%
  group_by(locus) %>%
  mutate(frac = ct/sum(ct))

ggplot()+ 
  geom_point(data=haplo.tot.tbl, aes(x = tot.hapl, y = locus, size=frac, color=frac))+
  scale_x_log10()+
  xlab("total number of unique haplotypes in an individual")+
  ylab("Locus ID")+
  scale_color_continuous("fraction")+
  scale_size_continuous(guide=FALSE)+
  theme_bw()+
  theme(legend.position="bottom")
  
  

 # prototype for casting haplotype composition profile
haplo.profile.frac <- haplo.sum %>% filter(depth >1, locus=="tag_id_1498") %>%
  group_by(haplo) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(frac = n/sum(n)) 

haplo.split.profile <- sapply(1:nrow(haplo.profile.frac), function(i) {
  char.split <- strsplit(haplo.profile.frac[i,]$haplo, "")
  n.char <- length(char.split[[1]])
  sapply(1:n.char, function(j) c(i, j, char.split[[1]][j], haplo.profile.frac[i,]$frac)) 
}) %>% 
  matrix(., ncol=4, byrow=T) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  tbl_df()
  
colnames(haplo.split.profile) <- c("group", "pos", "seq", "frac")
haplo.split.profile <- haplo.split.profile %>% mutate(pos=as.numeric(pos),
                               frac=as.numeric(frac),
                               group=as.numeric(group))


ggplot(haplo.split.profile, aes(x=factor(pos), y=seq, group=group, size=frac, color=factor(group))) +
  xlab("relative position")+
  ylab("sequence")+
  geom_point(alpha=0.9)+
  scale_size_continuous(guide=FALSE)+
  theme_bw()+
  theme(legend.position="bottom")

ggplot(haplo.sum %>% filter(depth >1, locus=="tag_id_1652"), aes(x=locus, y=depth)) +
  xlab("locus id")+
  ylab("Read depth")+
  geom_violin()+
  theme_bw()+
  scale_y_log10()


ggplot(haplo.sum %>% filter(depth >1, locus=="tag_id_1652"), aes(x=factor(id), y=depth)) +
  xlab("indiv id")+
  ylab("Read depth")+
  theme_bw()+
  geom_violin()+
  scale_y_log10()
