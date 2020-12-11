for i in *.sorted; do bname=`echo $i | cut -d. -f 1`; echo $bname; samtools depth $i | awk -v var="$bname" '{sum+=$3}; END {print var "\t" sum/NR}' >> all_genome_coverage.tsv; done
