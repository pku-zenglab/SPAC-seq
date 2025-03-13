#!/bin/bash

# Run Stero-seq SAW Mapping pipeline
# Last Modified: 24-04-14
# Changelog:
#   Created the pipeline using SAW.sh
##########################

show_usage() {
    echo "Usage: $0 -i <in1> <in2> -o <out1> <out2> [options]"
    echo "-"
    echo "Required:"
    echo "  -r <referenceGenome>                The reference path for index STAR genome fasta file"
    echo "  -g <referenceSJdbFile>              The reference path for index gtf file"
    echo "  -s <sampleFilePath>                 the txt summerizing the fastq, maskFile path, ABSOLUTE"
    echo "  -o <outPath>                        The path to output, ABSOLUTE"
    echo "Optional:"
    echo "  -@ <threads>                        preferred max threads, default: 8"
    echo "  -d, --doCellBin [Y/N]               do cell binning, default: N"
}

doCellBin=0
threas=4

while [[ -n "$1" ]]; do
    case "$1" in
        -r) reference="$2"
            shift ;;
        -s) sampleFilePath="$2"
            shift ;;
        -o) outputPath="$2"
            shift ;;
        -@) threads="$2"
            shift ;;
        -d) doCellBin"$2"
            shift ;;
    esac
        shift
done

# generate STAR reference base
# input: gtf file with gRNA annotation
#        fasta file with gRNA seq
# output: STAR index directory
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir . --genomeFastaFiles lib-tumor-refine.fa --sjdbGTFfile lib-tumor-refine.gtf --sjdbOverhang 20 --genomeSAindexNbases 6 --genomeChrBinNbits 4

# CID count
# intput: brcode to position mask file
#         species id
#         genome size
# output: a file containing CID count and memory use
/usr/bin/time -v singularity exec ${sif} CIDCount \
        -i ${maskFile} \
        -s ${refName} \
        -g ${GSize} > ${result_00mapping}/CIDCount

# bcSTAR alignment
# input: a list of read1 file
#        
for ((i=0;i<=`expr $(echo $fqNumber) - 1`;i++)); do
    fqname=$(basename ${read1List[i]})
    fqdir=$(dirname ${read1List[i]})
    fqbase=${fqname%%.*}
    fqbases[i]=$fqbase
    bcPara=${result_00mapping}/${fqbase}.bcPara
    barcodeReadsCount=${result_00mapping}/${fqbase}.barcodeReadsCount.txt
    echo  " ~~~ mapping - $fqname ~~~"
    echo "in=${maskFile}" > $bcPara
    echo "in1=${read1List[i]}" >> $bcPara
    echo "in2=${read2List[i]}" >> $bcPara
    echo "encodeRule=ACTG" >> $bcPara
    echo "action=4" >> $bcPara
    echo "barcodeReadsCount=${barcodeReadsCount}" >> $bcPara
    echo "platform=T10" >> $bcPara
    echo "out=${fqbase}" >> $bcPara
    echo "barcodeStart=0" >> $bcPara
    echo "barcodeLen=25" >> $bcPara
    echo "umiStart=25" >> $bcPara
    echo "umiLen=10" >> $bcPara
    echo "umiRead=1" >> $bcPara
    echo "mismatch=1" >> $bcPara
    echo "bcNum=`head -1 ${result_00mapping}/CIDCount`" >> $bcPara
    echo "polyAnum=15" >> $bcPara
    echo "mismatchInPolyA=2" >> $bcPara
    read1DIR=$(dirname ${read1List[i]})
    read2DIR=$(dirname ${read2List[i]})
    export SINGULARITY_BIND=$read1DIR,$read2DIR,$outDir,$maskDIR,$annoDIR,$refDIR
    /usr/bin/time -v singularity exec ${sif} mapping \
        --outSAMattributes spatial \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir ${GDir} \
        --runThreadN ${threads} \
        --outFileNamePrefix ${result_00mapping}/${fqbase}. \
        --sysShell /bin/bash \
        --stParaFile ${bcPara} \
        --readNameSeparator \" \" \
        --limitBAMsortRAM 63168332971 \
        --limitOutSJcollapsed 10000000 \
        --limitIObufferSize=280000000 \
        --outBAMsortingBinsN 50 \
        > ${result_00mapping}/${fqbase}_barcodeMap.stat

    starBam=${result_00mapping}/${fqbase}.Aligned.sortedByCoord.out.bam
    starBams[i]=$starBam
    bcStat[i]=${result_00mapping}/${fqbase}_barcodeMap.stat
    bcFinalOut[i]=${result_00mapping}/${fqbase}.Log.final.out
    bcReadsCounts[i]=$barcodeReadsCount
done

#merge to barcodeReadsCount file
echo "`data` "
/usr/bin/time -v singularity