#### SAPtyper ####

saptyper<-function(genomes,outfile){
	
	# Define primers
	# According to Tu et al. (2001) doi://10.1128/IAI.69.4.2237–2244.2001
	
	saf01<-'ATGTTAAACAAAACAGATGT'
	sar01<-'ATCAAGATCACTAGCACTA'
	sbf01<-'TTCAGAGCTATTTATAGTTC'
	sbr01<-'TCAACACTACTACTATTACTA'
	
	system('touch saf01.fasta sar01.fasta sbf01.fasta sbr01.fasta')
	
	cat('>SAF01',sep='\n',file='saf01.fasta')
	cat(saf01,sep='\n',file='saf01.fasta',append=T)
	
	cat('>SAR01',sep='\n',file='sar01.fasta')
	cat(sar01,sep='\n',file='sar01.fasta',append=T)
	
	cat('>SBF01',sep='\n',file='sbf01.fasta')
	cat(sbf01,sep='\n',file='sbf01.fasta',append=T)
	
	cat('>SBR01',sep='\n',file='sbr01.fasta')
	cat(sbr01,sep='\n',file='sbr01.fasta',append=T)
	
	result<-NULL
	
	for (g in genomes){
		
		saf01.cmd<-paste("blastn -query saf01.fasta -subject ",g," -outfmt '6 pident qlen length gaps' -out ",saf.blout,sep='')
		sar01.cmd<-paste("blastn -query sar01.fasta -subject ",g," -outfmt '6 pident qlen length gaps' -out ",sar.blout,sep='')
		sbf01.cmd<-paste("blastn -query sbf01.fasta -subject ",g," -outfmt '6 pident qlen length gaps' -out ",sbf.blout,sep='')
		sbr01.cmd<-paste("blastn -query sbr01.fasta -subject ",g," -outfmt '6 pident qlen length gaps' -out ",sbr.blout,sep='')			
		saf<-read.table('saf.blout',sep='\t',header=F)
		sar<-read.table('sar.blout',sep='\t',header=F)
		sbf<-read.table('sbf.blout',sep='\t',header=F)
		sbr<-read.table('sbr.blout',sep='\t',header=F)
		
		colnames(saf)<-c('pid','qlen','len','gaps')
		colnames(sar)<-c('pid','qlen','len','gaps')
		colnames(sbf)<-c('pid','qlen','len','gaps')
		colnames(sbr)<-c('pid','qlen','len','gaps')
		
		saf$ali<-saf$len-saf$gaps
		sar$ali<-sar$len-sar$gaps
		sbf$ali<-sbf$len-sbf$gaps
		sbr$ali<-sbr$len-sbr$gaps
		
		saf.num<-length(which(saf$pid==100 & saf$ali==saf$qlen))
		sar.num<-length(which(sar$pid==100 & sar$ali==sar$qlen))
		sbf.num<-length(which(sbf$pid==100 & sbf$ali==sbf$qlen))
		sbr.num<-length(which(sbr$pid==100 & sbr$ali==sbr$qlen))
		
		if (saf.num>=1 & sar.num>=1 & sbf.num==0 & sbr.num==0){
			
			result<-c(result,'A')
			
		} else if (saf.num==0 & sar.num==0 & sbf.num>=1 & sbr.num>=1){
			
			result<-c(result,'B')
			
		} else if ((saf.num>=1 & sar.num>=1 & sbf.num>=1 & sbr.num>=1){
			
			result<-c(result,'A/B')
			
		} else {
			
			result<-c(result,'U')
		}
		
	}
	
	dfr<-as.data.frame(cbind(genomes,result))
	
	colnames(dfr)<-c('genome','SAP')
	
	write.table(dfr,file=outfile,sep='\t',quotes=F,row.names=F)
	
	system('rm -rf saf*.fasta sar*.fasta sbf*.fasta sbr*.fasta')
	
}