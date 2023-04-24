#wget http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif50.motif


#scanMotifGenomeWide.pl motif50.motif mm10 -p 50 -bed > ctcf.scan.out.bed

grep "-" ctcf.scan.out.bed > ctcf.scan.out.neg.bed
grep "+" ctcf.scan.out.bed > ctcf.scan.out.pos.bed