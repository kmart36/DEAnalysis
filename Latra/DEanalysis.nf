nextflow.enable.dsl=2

workDir = '/home/kam071/DEAnalysis/Latra/work/'

process fastQC {
	conda 'fastqc'
	executor 'slurm'
	memory '8 GB'

	input:
	tuple val(sample_id), path(sample)

	script:
	"""
		mkdir -p /home/kam071/DEAnalysis/Latra/fastqc_files/${sample[0].simpleName}_fastqc
		fastqc --outdir /home/kam071/DEAnalysis/Latra/fastqc_files/${sample[0].simpleName}_fastqc ${sample[0]} --extract
	"""
}

//Do we want to save unpaired? Or just paired?

process trimmomatic {

	publishDir "trimmo_files", mode: 'copy'

	cpus 4
	memory '8 GB'
	executor 'slurm'

    input:
	tuple val(sample_id), path(sample)

	output:
	tuple val(sample_id), path('Pcorr*_{1,2}_paired.fq')

	script:
	"""
		echo ${sample_id}
		trimmomatic PE ${sample[0]} ${sample[1]} ${sample[0].simpleName}_paired.fq ${sample[0].simpleName}_unpaired.fq ${sample[1].simpleName}_paired.fq ${sample[1].simpleName}_unpaired.fq \
		ILLUMINACLIP:/software/apps/Trimmomatic/current/adapters/TruSeq3-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	"""
}

process kraken {
	
	publishDir "kraken_files", mode: 'copy'

	executor 'slurm'
	cpus 4
	memory '96 GB'
	queue 'medium'

	input:
	tuple val(sample_id), path(kraken)

	output:
	path("*.unclassified.out.fq")

	script:
	"""
		kraken2 --threads 16 --db /data/Kraken2/PlusPF/PlusPF_2023.db --paired \
        --classified-out ${sample_id}_R#.classified.out.fq \
        --unclassified-out ${sample_id}_R#.unclassified.out.fq \
        --confidence 0.5 \
        --output ${sample_id}.kraken.txt \
        --report ${sample_id}.kraken.report \
        ${kraken[0]} \
        ${kraken[1]}
	"""
}

workflow trimmo {
	take: 
		rawReads
	main:
		//fastQC(rawReads)
		trimmomatic(rawReads)
	emit:
		trimmomatic.out
}

workflow kraken2 {
	take:
		trim_files
	main:
		fastQC(trim_files)
		kraken(trim_files)
	emit:
		kraken.out
}

workflow {
	//rawReads = Channel.fromFilePairs('/home/sel025/LowerLab/2021_09_29_Svistunov_V_Latra_ant_BL_transcriptome_novogene/usftp21.novogene.com/raw_data/Pcorr*/Pcorr*_{1,2}.fq.gz', checkIfExists: true)
	trim = Channel.fromFilePairs('/home/kam071/DEAnalysis/Latra/trimmo_files/*_{1,2}_paired.fq', checkIfExists: true)
	kraken(trim)
	fastQC(kraken.out)
	//trimmo(rawReads)
	//kraken2(trimmo.out)
	//fastQC(kraken2.out)
}