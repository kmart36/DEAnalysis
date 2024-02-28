nextflow.enable.dsl=2

// Change this path to where you want your working directories to be generated. These will be cleared when the pipeline finishes, so
// this is not incredibly important. 
workDir = '/home/kam071/DEAnalysis/Latra/'

// Change this path to where you want your generated files to be stored (slash at the beginning, no slash at the end).
mainPath = '/home/kam071/DEAnalysis/Latra'

// Change this path to where your python scripts are located (find_pair.py and table_generator.py) 
// (slash at the beginning, no slash at the end).
pyPath = '/home/kam071/DEAnalysis/Latra'

// Change this path to where your raw reads are located. You may also have to change what's written in the rawReads line so it
// matches the format of your files (slash at the beginning, no slash at the end).
rawDataPath = '/home/sel025/LowerLab/2021_09_29_Svistunov_V_Latra_ant_BL_transcriptome_novogene/usftp21.novogene.com/raw_data/Pcorr*'
rawReads = Channel.fromFilePairs(rawDataPath + '/*_{1,2}.fq.gz', checkIfExists: true)
latra_query = '/home/kam071/DEAnalysis/Latra/query_ORs_V3.faa'

process fastQC {
	conda 'fastqc'
	executor 'slurm'
	memory '4 GB'

	input:
	tuple val(sample_id), path(sample)

	output:
	path("done.txt")

	script:
	"""
		mkdir -p $mainPath/fastqc_files/${sample[0].baseName}_fastqc
		fastqc --outdir $mainPath/fastqc_files/${sample[0].baseName}_fastqc ${sample[0]} --extract

		echo done > done.txt
	"""
}

process trimmomatic {

	publishDir "trimmo_files", mode: 'copy'

	cpus 4
	memory '4 GB'
	executor 'slurm'

    input:
	tuple val(sample_id), path(sample)

	output:
	tuple val(sample_id), path('*_{1,2}_paired.fq')

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

process table_gen {
	executor 'slurm'
	cpus 1
	memory '1 GB'
	queue 'short'

	input:
	path(done)

	output:
	path("kraken_done.txt")

	script:
	"""
		python $pyPath/table_generator.py
		echo done > done.txt
	"""
}

process find_pair {
	executor 'slurm'
	cpus 1
	memory '1 GB'
	queue 'short'

	input:
	path(done)

	output:
	path("done.txt")

	script:
	"""
		python $pyPath/find_pair.py
		echo done > done.txt
	"""
}

mode = 'norm'
left = file(mainPath + '/left_kraken.txt')
right = file(mainPath + '/right_kraken.txt')

process trinity {
	
	executor 'slurm'
	cpus 6
	memory '64 GB'
	queue 'medium'

	publishDir "trinity_files", mode: 'copy'

	input:
	path(done)

	output:
	path("*.Trinity.fasta")

	script:
	if (mode == 'infected')
		"""
			singularity exec /software/singularity-containers/2021-10-21-trinityrnaseq.v2.13.2.simg Trinity --seqType fq --left $in_left.text --right $in_right.text --max_memory 64G --CPU 8 --output trinity_files --full_cleanup
		"""
	else if (mode == 'noninfected')
		"""
			singularity exec /software/singularity-containers/2021-10-21-trinityrnaseq.v2.13.2.simg Trinity --seqType fq --left $nonin_left.text --right $nonin_right.text --max_memory 64G --CPU 8 --output trinity_files --full_cleanup
		"""
	else
		"""
			singularity exec /software/singularity-containers/2021-10-21-trinityrnaseq.v2.13.2.simg Trinity --seqType fq --left $left.text --right $right.text --max_memory 64G --CPU 8 --output trinity_files --full_cleanup

		"""
}

process busco {

	executor 'slurm'
	cpus 2
	memory '8 GB'

	publishDir "busco_files", mode: 'move'

	input:
	path(trin)

	script:
	"""
		cp -r /software/apps/augustus/current/config/ .
		busco -i $trin -o ${trin.simpleName}.busco -l endopterygota_odb10 -m transcriptome -f
	"""
}

process cdhit {

	executor 'slurm'
	cpus 8
	memory '64 GB'

	input:
	path(trin)

	output:
	path("trinity.cd-hit.results*")

	script:
	"""
		cd-hit-est -i $trin -o $mainPath/trinity.cd-hit.results -c 0.90 -n 8
	"""
}

process salmon {

	executor 'slurm'
	cpus 8
	memory '64 GB'

	script:
	"""
	"""
}

process blast_db {

	executor 'slurm'
	cpus 2
	memory '16 GB'

	input:
	path(cd_hit)

	output:
	path("Latr_db1*")

	script:
	"""
		makeblastdb -in $cd_hit -out Latr_db1 -dbtype nucl
	"""
}

process blast {

	input:
	path(db_1)

	script:
	"""
		blastp -query $latra_query -db $db_1 -out Latr_blast1_hits.faa -evalue 1 -outfmt 7
	"""
}


workflow trimmo {
	take: 
		rawReads
	main:
		fastQC(rawReads)
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

workflow trin {
	take:
		done
	main:
		trinity(done)
	emit:
		trinity.out
}


workflow {
	trimmo(rawReads) | kraken2 | fastQC | collectFile(name: 'kraken_done.txt', newLine: true) | table_gen | find_pair | trinity | cdhit | blast_db | blast
}