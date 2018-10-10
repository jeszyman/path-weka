# Setup
#  R setup
#  Install these packages before starting:
#   stringr
#   reshape2
#   e1071
#   data.table

library('stringr')
library('reshape2')
library('e1071')
library('data.table')

#  Import
#   UNCOMMENT AND ENTER PARAMETERS BELOW: 
# source.cbc.rbo=list.files(path="<PATH TO YOUR RBO csvs>", pattern="<pattern matching your RBO CBC files>", full.names=T)
# e.g.:
#  source.cbc.rbo=list.files(path="~/data/phi/cerner/", pattern="cbc.+rbo.+csv", full.names=T)
input.list=lapply(source.cbc.rbo, read.csv)
cat.cbc=do.call(rbind, input.list)
mod.cbc=cat.cbc

# Convert time
mod.cbc$RECD_DTTM=as.character(mod.cbc$RECD_DTTM)
mod.cbc$RECD_DTTM=str_replace_all(mod.cbc$RECD_DTTM, ' 0:',' 00:')
mod.cbc$rec.time=as.POSIXct(mod.cbc$RECD_DTTM, format = "%m/%d/%y %H:%M", tz="America/Chicago")

# DTA long to wide
mod.cbc=mod.cbc[!duplicated(mod.cbc[c(2,10,23)]),]
wide.cbc=dcast(mod.cbc, AGE + ORDER_PROVIDER + SEX + PERSON_ID +ACCESSION_NO + ENCNTR_ID + ORDR_DTTM + rec.time ~ DTA, value.var = c("RESULT"))

# Define anemia
wide.cbc$anemia=ifelse(wide.cbc$SEX=="Female" & wide.cbc$Hgb<12, 'anemic',
                            ifelse(wide.cbc$SEX=="Male" & wide.cbc$Hgb<13, 'anemic', 'normal'))

# Subset to ML populations
#  Sequence
wide.cbc=wide.cbc[with(wide.cbc,order(PERSON_ID, rec.time)),]
wide.cbc=data.table(wide.cbc, key='PERSON_ID')
wide.cbc[,pt.seq:=seq(.N),by='PERSON_ID']

#  Identify not final and not anemic, classify subsAnemia
wide.cbc$subsAnemia=ifelse(wide.cbc$pt.seq < shift(wide.cbc$pt.seq, type=("lead")) & wide.cbc$anemia=='normal', shift(wide.cbc$anemia, type=('lead')), 'non')

# Case control split
cbc=subset(wide.cbc,wide.cbc$subsAnemia!="non")
cbc$match=ifelse(cbc$subsAnemia=="normal", 'cont','case')
m = matchControls(match ~ Hct + Hgb, data=cbc)
cbc$m=m$factor
cbc=subset(cbc, m %in% c('cont','case'))
table(cbc$m)
boxplot(cbc$Hct ~ cbc$m)

# Deidentify
#  Set seed- replace below with a number:
# set.seed(<ANY NUM>)
#  Create randomized columns
#   Person
patientID=data.frame(mskPERSON_ID=sample(1:length(unique(cbc$PERSON_ID)),replace=F),PERSON_ID=unique(cbc$PERSON_ID))
cbc=merge(cbc,patientID,by="PERSON_ID",all=T)
#   Accession
accessionID=data.frame(mskACCESSION=sample(1:length(unique(cbc$ACCESSION_NO)),replace=F),ACCESSION_NO=unique(cbc$ACCESSION_NO))
cbc=merge(cbc,accessionID,by="ACCESSION_NO",all=T)
#   Encounter 
encounterID=data.frame(mskENCOUNTER=sample(1:length(unique(cbc$ENCNTR_ID)),replace=F),ENCNTR_ID=unique(cbc$ENCNTR_ID))
cbc=merge(cbc,encounterID,by="ENCNTR_ID",all=T)
#  Remove column with identifiers
cbc=as.data.frame(cbc)
cbc=cbc[ , -which(names(cbc) %in% c("ENCNTR_ID","PERSON_ID","ACCESSION_NO"))]

write.csv(cbc, "~/cbc.csv")
