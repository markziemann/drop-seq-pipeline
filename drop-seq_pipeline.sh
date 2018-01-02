#!/bin/bash

#set -x

CWD=$(pwd)
REFDIR=$CWD/refgenome
RAWDATA=$CWD/raw_fastq
DEMUXMAP=$CWD/demux_map

REFTX_URL="ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
REFTX=$REFDIR/$(basename $REFTX_URL .gz)
REFIDX=$REFTX.idx

RRNA=$REFDIR/mouse_rRNA.fa

REF_GENOME_URL="ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
REF_GENOME=$REFDIR/$(basename $REF_GENOME_URL .gz)

REF_GTF_URL="ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz"
REF_GTF=$REFDIR/$(basename $REF_GTF_URL .gz)

RSCRIPT_PATH=$CWD/drop-seq_plots.R

########################################
echo generate transcriptome reference
########################################
mkdir $REFDIR ; cd $REFDIR
if [ ! -r $REFIDX ] ; then
  wget -N $REFTX_URL
  gunzip -f $REFTX.gz
  kallisto index -i $(basename $REFIDX) $(basename $REFTX)
fi

if [ ! -r "SAindex" ] ; then
  wget -N $REF_GENOME_URL
  gunzip -f $REF_GENOME.gz
  wget -N $REF_GTF_URL
  gunzip -f $REF_GTF.gz

  STAR --runMode genomeGenerate \
  --sjdbGTFfile $REF_GTF \
  --genomeDir $REFDIR  \
  --genomeFastaFiles $REF_GENOME \
  --runThreadN $(nproc)
fi

cd ..

# function to count duplicate reads
run_dupqc(){
FQZ1=$1
FQZ2=$(echo $FQZ1 | sed 's/_1.txt.gz/_2.txt.gz/')
TOT=$(zcat $FQZ1 | sed -n '2~4p' | wc -l )
UNIQ=$(paste <(zcat $FQZ1 | sed -n '2~4p' ) <(zcat $FQZ2 | sed -n '2~4p' | cut -c-20) | sort | uniq -u | wc -l)
DUP=$(paste <(zcat $FQZ1 | sed -n '2~4p' ) <(zcat $FQZ2 | sed -n '2~4p' | cut -c-20) | sort | uniq -d | wc -l)
echo $FQZ1 $TOT $UNIQ $DUP
}
export -f run_dupqc


# function for skewer
run_skewer(){
DEMUX_FQZ=$1
skewer -q 20 $DEMUX_FQZ
BASENAME=$(basename $1 .gz)
TRIMMED="${BASENAME}-trimmed.fastq"
fastx_trimmer -Q33 -l 100 -i $TRIMMED > $TRIMMED.tmp
mv $TRIMMED.tmp $TRIMMED
}
export -f run_skewer


# function for kallisto alignment and report human and mouse et counts
run_kal(){
FQ=$1
IDX=$2
TSV=$FQ.abundance.tsv
LOG=$FQ.kal.log
kallisto quant -t3 -i $IDX -o $FQ.kal --single -l 100 -s 20 --fr-stranded $FQ &> $LOG
mv $FQ.kal/abundance.tsv $TSV
BASE=$(echo $TSV | cut -d '_' -f1)
cut -f1,4 $TSV | sed 1d | sed "s/^/${BASE}\t/" > $TSV.tmp
MMU_CNT=$(grep ENSMUS $TSV | cut -f4 | numsum)
TOT_CNT=$(grep processed $LOG | awk '{print $3}' | tr -d ',')
MAP_CNT=$(grep processed $LOG | awk '{print $5}' | tr -d ',')
echo $TSV $MMU_CNT $TOT_CNT $MAP_CNT
}
export -f run_kal

# function to map with STAR aligner not saving bam files then report human and mouse counts
run_star(){
PFX=$1.star.
TAB=$1.star.ReadsPerGene.out.tab
FQ=$1
REFDIR=$2
#fastx_clipper -n -a AAAAAAAAA -i $FQ -l 25 -Q33 | fastx_trimmer -l 80 -Q33 > $FQ.tmp
fastx_clipper -n -a AAAAAAAAA -i $FQ -l 25 -Q33  > $FQ.tmp ; mv $FQ.tmp $FQ
STAR --runThreadN 3 --quantMode GeneCounts --genomeLoad LoadAndKeep \
--outSAMtype None --genomeDir $REFDIR --readFilesIn=$FQ \
--outFileNamePrefix $FQ.star.
MMU_CNT=$(grep ENSMUS $TAB | cut -f3 | numsum)
echo $TAB $MMU_CNT
BASE=$(echo $FQ | cut -d '_' -f1)
cut -f1,3 $TAB | tail -n +5 | sed "s/^/${BASE}\t/" > $TAB.tmp
}
export -f run_star

# function to map reads to rRNA reference to quantify ribo contamination

setup_rrna_quant(){
REFURL="https://drive.google.com/uc?authuser=0&id=0BweZVB8jfIwPQVJaazZSbnRZcVE&export=download"
REFNAME=mouse_rRNA.fa
REFDIR=$1
REF=$REFDIR/$REFNAME

if [ ! -r $REF ] ; then
  wget -O $REF "https://drive.google.com/uc?authuser=0&id=0BweZVB8jfIwPQVJaazZSbnRZcVE&export=download"
  bwa index $REF
fi
}
export -f setup_rrna_quant

rrna_quant(){
FQ=$1
REFNAME=mouse_rRNA.fa
REFDIR=$2
REF=$REFDIR/$REFNAME
BAM=${FQ}.sort.bam
STATS=$BAM.stats

bwa aln -t 30 $REF $FQ \
| bwa samse $REF - $FQ \
| samtools view -uSh - \
| samtools sort -o $BAM -
samtools index $BAM
samtools flagstat $BAM > $STATS
rm $BAM $BAM.bai
TOTAL=$(awk 'NR==1 {print $1}' $STATS)
RRNA=$(awk 'NR==5 {print $1}' $STATS)
echo $FQ $TOTAL $RRNA
}
export -f rrna_quant

########################################
echo barcoding qc
########################################
mkdir $DEMUXMAP ; cd $DEMUXMAP

for FQZ1 in $RAWDATA/*_R1_001.fastq.gz ; do

  FQZ2=$(echo $FQZ1 | sed 's/_R1_/_R2_/')
  ln $FQZ1 $FQZ2 .
  FQZ1=$(basename $FQZ1)
  FQZ2=$(basename $FQZ2)

  pwd
  TABLE=$FQZ1.bc.tsv
  BCTABLE=$FQZ1.bc
  BCCOUNT=$(zcat $FQZ1 | sed -n '2~4p' | cut -c-12 | grep -v NN | sort -u | wc -l)
  HINDEX=$(zcat $FQZ1 | sed -n '2~4p' | cut -c-12 | grep -v NN | sort | uniq -c | sort -k1gr | nl -n ln | tr -s ' ' | tr ' ' '\t' | tr -s '\t' | tee $TABLE | awk '$1>$2 {print $1}' | head -1)
  TOT=$(cut -f2 $TABLE | numsum)
  HALF=$((TOT/2))
  HALFLINE=$(awk '{print $0, (total += $2)}' $TABLE | awk -v H=$HALF '$4>H {print $1} '| head -1)
  NL=$(wc -l < $TABLE)
  GINI=$(echo $HALFLINE $NL | awk '{print ($2-$1)/$2}')
  echo BCCOUNT:$BCCOUNT HINDEX:$HINDEX TOT:$TOT HALF:$HALF HALFLINE:$HALFLINE NumLines:$NL GINI:$GINI > $FQZ1.log

  awk '$2>1000 {print $3"\t"$3}' $TABLE > $BCTABLE
  DEMUX_DIR=$DEMUXMAP/$FQZ1.demux
  mkdir $DEMUX_DIR

  if [ ! -r $DEMUX_DIR/finished ] ; then

    je demultiplex F1=$FQZ1 F2=$FQZ2 BF=$BCTABLE O=$DEMUX_DIR BPOS=READ_1

    touch $DEMUX_DIR/finished

  fi

  cd $DEMUX_DIR

  parallel run_dupqc ::: *_1.txt.gz > ../${FQZ1}.dupl.txt
  parallel -j30 run_skewer ::: *_2.txt.gz
  parallel -j10 run_kal ::: *_2.txt-trimmed.fastq ::: $REFIDX > ../$FQZ1.kalcounts.txt
  STAR --genomeLoad Remove --genomeDir $REFDIR
  STAR --genomeLoad LoadAndExit --genomeDir $REFDIR
  parallel -j3 run_star ::: *_2.txt-trimmed.fastq ::: $REFDIR | grep ReadsPerGene.out.tab > ../$FQZ1.starcounts.txt
  STAR --genomeLoad Remove --genomeDir $REFDIR

  for LOG in *Log.final.out ; do
    TOT=$(grep -i 'Number of input reads' $LOG | awk '{print $NF}')
    UNIQMAPRATE=$(grep -i 'Uniquely mapped reads %' $LOG | awk '{print $NF}' | tr -d '%')
    echo $LOG $TOT $UNIQMAPRATE
  done > ../${FQZ1}.starmapping.txt

  setup_rrna_quant $REFDIR
  parallel rrna_quant ::: *fastq ::: $REFDIR > ../$FQZ1.rrna.counts.txt

  KAL_DB=$FQZ1.kal.3col
  KAL_MX=$FQZ1.kal.mx
  for TMP in *abundance.tsv.tmp ; do
    cat $TMP
  done > ../$KAL_DB

  STAR_DB=$FQZ1.star.3col
  STAR_MX=$FQZ1.star.mx
  for TMP in *star.ReadsPerGene.out.tab.tmp ; do
    cat $TMP
  done > ../$STAR_DB

  Rscript $RSCRIPT_PATH $FQZ1
  cd -
done

