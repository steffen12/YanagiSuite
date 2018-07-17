#!/bin/bash
#./YanagiPipeline.sh sim_reads txome_ENSG00000100842.fa polyester

outputDir=$1
transcriptomeLib=$2
readsType=$3 #can be polyester, rlsim, or real

if [ "$readsType" = "polyester" ]; then
	echo "Simulating polyester Reads"
	outputDir="polyester_${outputDir}"
	mkdir -p "${outputDir}"
	Rscript PolyesterSimulateData.R "${transcriptomeLib}" "${outputDir}"
	perl fasta_to_fastq.pl "${outputDir}/sample_01.fasta" > "${outputDir}/sample_01.fq"
	readsFilename="${outputDir}/sample_01.fq"
elif [ "$readsType" = "rlsim" ]; then
	echo "Simulating rlsim Reads"
	outputDir="rlsim_${outputDir}"
	mkdir -p "${outputDir}"
	rlsim/tools/sel "${transcriptomeLib}" > "${outputDir}/${transcriptomeLib}_expr.fa"
	rlsim/src/rlsim -n 1000 "${outputDir}/${transcriptomeLib}_expr.fa" > "${outputDir}/${transcriptomeLib}_frags.fa"
	cat "${outputDir}/${transcriptomeLib}_frags.fa"| simNGS/bin/simNGS rlsim/src/test/cov/s_4_0066.runfile > "${outputDir}/${transcriptomeLib}_reads.fq"
	readsFilename="${outputDir}/${transcriptomeLib}_reads.fq"
elif [ "$readsType" = "real" ]; then
	echo "Reading Real Reads"
	readsFilename=$4
else
	echo "Invalid reads type"
	exit
fi

readStats=$(python YanagiReadStats.py ${readsFilename})
meanReadLen=$(echo ${readStats} | cut -d ' ' -f1)
readLenSd=$(echo ${readStats} | cut -d ' ' -f2)

segmentsLib="seg_${transcriptomeLib}"

echo "Calculating Yanagi TPMs"
kallisto index -i "${segmentsLib}.idx" "${segmentsLib}"
kallisto quant -i "${segmentsLib}.idx" -o "${outputDir}/KallistoOutputYanagi" --single -l "$meanReadLen" -s "$readLenSd" --pseudobam "${readsFilename}" \
| python YanagiCount.py --single --stdin1 "${segmentsLib}" "${outputDir}/YanagiCountOutput" 100
python YanagiQuant.py -counts "${outputDir}/YanagiCountOutput/seg_counts.tsv" -readlen "$meanReadLen" -fraglen "$meanReadLen" -out "${outputDir}/YanagiTPM.tsv"

echo "Calculating Kallisto TPMs"
kallisto index -i "${transcriptomeLib}.idx" "${transcriptomeLib}"
kallisto quant -i "${transcriptomeLib}.idx" -o "${outputDir}/KallistoOutput" --single -l "$meanReadLen" -s "$readLenSd" "${readsFilename}"

if [ "$readsType" = "polyester" ] || [ "$readsType" = "rlsim" ]; then
	echo "Comparing to True TPMs"
	python YanagiTPMCalculator.py "$readsType" "${readsFilename}" "${transcriptomeLib}" "$meanReadLen" "${outputDir}/trueTPM.tsv"
	python YanagiQuantTest.py "${outputDir}/trueTPM.tsv" "${outputDir}/YanagiTPM.tsv" "${outputDir}/KallistoOutput/abundance.tsv"
fi