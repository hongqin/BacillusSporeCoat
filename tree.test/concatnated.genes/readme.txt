
make sure all seqs are in upper case
pick_fasta_records_by_ids.1b.pl -if cc.cds.fas -id _ids.tab -out cc.cdsUPCASE.fas -upc 1
 align_cds_by_protein_bioperl.00.pl -i cc.cdsUPCASE.fas -o cc.cds.nexus -af nexus

$ align_cds_by_protein_bioperl.00.pl -i cc.cds.fas -o cc.cdsaligned.fas -af fasta
$ seqCat.pl -dseqCat.in 

-----
=> try consel
-> protml 
$ muscle -in cc.prot.fas -out cc.prot-aling.cls -clw
$ clus2mol.pl cc.prot-align.clw cc.prot-align.molphy 
$ cc.prot-align.molphy -sd

