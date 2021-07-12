for f in *_1.fastq.gz;
do
out=${f/_1.fastq.gz/_S1_L001_I1_001.fastq.gz}
mv $f $out;
done

for f in *_2.fastq.gz;
do
out=${f/_2.fastq.gz/_S1_L001_R1_001.fastq.gz}
mv $f $out;
done

for f in *_3.fastq.gz;
do
out=${f/_3.fastq.gz/_S1_L001_R2_001.fastq.gz}
mv $f $out;
done


