#!/bin/bash
#./YanagiPipeline.sh sim_reads ~/Penn_Internship Homo_sapiens_txome.fa YanagiSegmenter/homo_sapiens_segs_101.fa single polyester
#./YanagiPipeline.sh sim_reads ~/Penn_Internship Homo_sapiens_txome.fa create single rlsim

#./YanagiPipeline.sh sim_reads ~/Penn_Internship homo_sapiens_txome.fa_62.fa YanagiSegmenter/homo_sapiens/homo_sapiens_segs_151.fa_62.fa single polyester simulate

outputDir=$1
YanagiSuiteDir=$2
transcriptomeLib=$3
segmentsLib=$4 #can be create to create new segments
pairedEnd=$5
readsType=$6 #can be polyester, rlsim, or real
readsFilename1=$7 #can be simulate to simulate

if [ "$pairedEnd" != "single" ] && [ "$pairedEnd" != "paired" ]; then
	echo "Invalid reads type"
	exit
elif [ "$readsType" != "polyester" ] && [ "$readsType" != "rlsim" ] && [ "$readsType" != "real" ]; then
	echo "Invalid reads type"
	exit
fi

outputDir="${readsType}_${outputDir}"
mkdir -p "${outputDir}"

if [ "$readsFilename1" = "simulate" ]; then
	if [ "$pairedEnd" = "single" ]; then
		if [ "$readsType" = "polyester" ]; then
			echo "Simulating polyester Reads"
			Rscript PolyesterSimulateData.R "${transcriptomeLib}" "${outputDir}" single
			perl fasta_to_fastq.pl "${outputDir}/sample_01.fasta" > "${outputDir}/sample_01.fq"
			readsFilename1="${outputDir}/sample_01.fq"
		elif [ "$readsType" = "rlsim" ]; then
			echo "Simulating rlsim Reads"
			rlsim/tools/sel "${transcriptomeLib}" > "${outputDir}/${transcriptomeLib}_expr.fa"
			rlsim/src/rlsim -n 30000 "${outputDir}/${transcriptomeLib}_expr.fa" > "${outputDir}/${transcriptomeLib}_frags.fa"
			cat "${outputDir}/${transcriptomeLib}_frags.fa" \
			| simNGS/bin/simNGS -p single rlsim/src/test/cov/s_4_0066.runfile > "${outputDir}/${transcriptomeLib}_reads.fq"
			readsFilename1="${outputDir}/${transcriptomeLib}_reads.fq"
		fi
	elif [ "$pairedEnd" = "paired" ]; then
		if [ "$readsType" = "polyester" ]; then
			echo "Simulating polyester Reads"
			Rscript PolyesterSimulateData.R "${transcriptomeLib}" "${outputDir}" paired
			perl fasta_to_fastq.pl "${outputDir}/sample_01_1.fasta" > "${outputDir}/sample_01_1.fq"
			perl fasta_to_fastq.pl "${outputDir}/sample_01_2.fasta" > "${outputDir}/sample_01_2.fq"
			readsFilename1="${outputDir}/sample_01_1.fq"
			readsFilename2="${outputDir}/sample_01_2.fq"
		elif [ "$readsType" = "rlsim" ]; then
			echo "Simulating rlsim Reads"
			rlsim/tools/sel "${transcriptomeLib}" > "${outputDir}/${transcriptomeLib}_expr.fa"
			rlsim/src/rlsim -n 30000 "${outputDir}/${transcriptomeLib}_expr.fa" > "${outputDir}/${transcriptomeLib}_frags.fa"
			cat "${outputDir}/${transcriptomeLib}_frags.fa" \
			| simNGS/bin/simNGS -p paired rlsim/src/test/cov/s_4_0066.runfile > "${outputDir}/${transcriptomeLib}_reads.fq"
			python rlsimPairedReadSeparator.py "${outputDir}/${transcriptomeLib}_reads.fq"
			readsFilename1="${outputDir}/${transcriptomeLib}_reads_1.fq"
			readsFilename2="${outputDir}/${transcriptomeLib}_reads_2.fq"
		fi
	fi
else
	if [ "$pairedEnd" = "paired" ]; then
		readsFilename2=$8 #can be simulate to simulate
	fi
fi

if [ "$readsType" = "polyester" ]; then
	meanFragmentLength=250
	sdFragmentLength=25
elif [ "$readsType" = "rlsim" ]; then
	meanFragmentLength=174.46
	sdFragmentLength=19.33
fi

#Need to do fragment length
meanReadLen=$(python YanagiReadStats.py ${readsFilename1} fastq)
echo "$readsFilename1"
echo "$readsFilename2"

if [ "$segmentsLib" = "create" ]; then
	Rscript YanagiSegmenter/preprocessing_main.R "${YanagiSuiteDir}" hg38.fa Homo_sapiens.GRCh38.92.chr.gtf
	python YanagiSegmenter/yanagi_main.py "$meanReadLen"
	segmentsLib="homo_sapiens_segs_${meanReadLen}.fa"
fi

if [ "$pairedEnd" = "single" ]; then
	echo "Calculating Yanagi TPMs"
	kallisto index -i "${segmentsLib}.idx" "${segmentsLib}"
	kallisto quant -i "${segmentsLib}.idx" -o "${outputDir}/KallistoOutputYanagi" --single -l "$meanReadLen" -s 0.0001 --pseudobam "${readsFilename1}" \
	| python YanagiCount.py --single --stdin1 "${segmentsLib}" "${outputDir}/YanagiCountOutput" 1000
	python YanagiQuant.py -counts "${outputDir}/YanagiCountOutput/seg_counts.tsv" -readlen "$meanReadLen" -fraglen "$meanFragmentLength" -out "${outputDir}/YanagiTPM.tsv"

	echo "Calculating Kallisto TPMs"
	kallisto index -i "${transcriptomeLib}_kallisto.idx" "${transcriptomeLib}"
	kallisto quant -i "${transcriptomeLib}_kallisto.idx" -o "${outputDir}/KallistoOutput" --single -l "$meanFragmentLength" -s "$sdFragmentLength" "${readsFilename1}"

	echo "Calculating Salmon TPMs"
	salmon/bin/salmon index -t "${transcriptomeLib}" -i "${transcriptomeLib}_salmon"
	salmon/bin/salmon quant -i "${transcriptomeLib}_salmon" -l A -r "${readsFilename1}" -o "${outputDir}/SalmonOutput"

	if [ "$readsType" = "polyester" ] || [ "$readsType" = "rlsim" ]; then
		echo "Comparing to True TPMs"
		python YanagiTPMCalculator.py "$readsType" "${readsFilename1}" "${transcriptomeLib}" "$meanFragmentLength" "${outputDir}/trueTPM.tsv"
		python YanagiQuantTest.py "${outputDir}/trueTPM.tsv" "${outputDir}/YanagiTPM.tsv" "${outputDir}/KallistoOutput/abundance.tsv" "${outputDir}/SalmonOutput/quant.sf"
	fi
elif [ "$pairedEnd" = "paired" ]; then
	echo "Calculating Yanagi TPMs"
	kallisto index -i "${segmentsLib}.idx" "${segmentsLib}"
	kallisto quant -i "${segmentsLib}.idx" -o "${outputDir}/KallistoOutputYanagi" --single -l "$meanReadLen" -s 0.0001 --pseudobam "${readsFilename1}" \
	| python YanagiCount.py --paired --stdin1 "${segmentsLib}" "${outputDir}/YanagiCountOutput" 1000
	kallisto quant -i "${segmentsLib}.idx" -o "${outputDir}/KallistoOutputYanagi" --single -l "$meanReadLen" -s 0.0001 --pseudobam "${readsFilename2}" \
	| python YanagiCount.py --paired --stdin2 "${segmentsLib}" "${outputDir}/YanagiCountOutput" 1000
	python YanagiQuant.py -counts "${outputDir}/YanagiCountOutput/seg_counts.tsv" -readlen "$meanReadLen" -out "${outputDir}/YanagiTPM.tsv"

	echo "Calculating Kallisto TPMs"
	kallisto index -i "${transcriptomeLib}_kallisto.idx" "${transcriptomeLib}"
	kallisto quant -i "${transcriptomeLib}_kallisto.idx" -o "${outputDir}/KallistoOutput" "${readsFilename1}" "${readsFilename2}"

	echo "Calculating Salmon TPMs"
	salmon/bin/salmon index -t "${transcriptomeLib}" -i "${transcriptomeLib}_salmon"
	salmon/bin/salmon quant -i "${transcriptomeLib}_salmon" -l A -1 "${readsFilename1}" -2 "${readsFilename2}" -o "${outputDir}/SalmonOutput"

	if [ "$readsType" = "polyester" ] || [ "$readsType" = "rlsim" ]; then
		echo "Comparing to True TPMs"
		python YanagiTPMCalculator.py "$readsType" "${readsFilename1}" "${transcriptomeLib}" "$meanFragmentLength" "${outputDir}/trueTPM.tsv"
		python YanagiQuantTest.py "${outputDir}/trueTPM.tsv" "${outputDir}/YanagiTPM.tsv" "${outputDir}/KallistoOutput/abundance.tsv" "${outputDir}/SalmonOutput/quant.sf"
	fi
fi