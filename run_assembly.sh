#!/bin/bash

#TODO: Ejecutar varios assemblies en paralelo en el caso que sea posible (no con trinity) en funcion de la cantidad de procesadores

# Default values
kmin=23
kmax=59
kstep=12
read_length=150
pair_end=0
fragment_size=0
fragment_size_sd=0
rlen_var=0
flen_var=0
n_cpus=1
assembly_soft="trinity"

# Ayuda y modo de empleo
mostrarAyuda () {
echo -e "USAGE:"
echo -e "\t$(echo $0 | awk -F'/' '{print $5}') [-a <string> | -l <number> | -M <number> | -m <number> | -p <number> | -d <number> | -s <number> | -c <number> | -L | -P | -h]"
echo -e "OPTIONS:"
echo -e "\t-a\tAssembly software (Available options: oases, trinity, abyss, soap, bridger, idba)."
echo -e "\t-p\tPair end fragment size (Default: Single end)."
echo -e "\t-d\tPair end fragment size standard deviation (Default: 0)."
echo -e "\t-l\tRead length (Default: 150 bps)."
echo -e "\t-L\tRead length is a variable (Default: False)."
echo -e "\t-P\tPair end fragment length is a variable (Default: False)."
echo -e "\t-m\tMinimum kmer value (Default: 23)."
echo -e "\t-M\tMaximum kmer value (Default: 59)."
echo -e "\t-s\tKmer step (Default: 12)."
echo -e "\t-c\tNumber of CPUs for each assembly (Default: 1)."
echo -e "\t-h\tShow this help."
}
# Procesamiento de parámetros
while getopts 'hLPm:M:s:l:p:a:c:d::' OPCION; do
  case $OPCION in
    l)
      read_length=$OPTARG;;
    a)
      assembly_soft=$OPTARG
      if [[ $assembly_soft == "oases" ]]; then
          n_as=1
      elif [[ $assembly_soft == "trinity" ]]; then
          n_as=2
      elif [[ $assembly_soft == "abyss" ]]; then
          n_as=3
      elif [[ $assembly_soft == "soap" ]]; then
          n_as=4
      elif [[ $assembly_soft == "bridger" ]]; then
          n_as=5
      elif [[ $assembly_soft == "idba" ]]; then
          n_as=6
      else
          echo "ERROR: The assembly software name is incorrect."
          mostrarAyuda
          exit 1
      fi
      log_file=$PWD/assembly_${assembly_soft}.log;;
    p)
      pair_end=1
      fragment_size=$OPTARG;;
    d)
      fragment_size_sd=$OPTARG;;
    m)
      kmin=$OPTARG;;
    M)
      kmax=$OPTARG;;
    s)
      kstep=$OPTARG;;
    c)
      n_cpus=$OPTARG;;
    L)
      rlen_var=1;;
    P)
      pair_end=1
      flen_var=1;;
    h)
      mostrarAyuda
      exit 0;;
    ?)
      echo "ERROR: Invalid arguments."
      mostrarAyuda
      exit 1;;
  esac
done

#fix kmax value only for trinity
if [[ $assembly_soft == "trinity" ]]; then
    if [[ $kmax -gt 32 ]]; then
        kmax=32
    fi
    if [[ $kmin -lt 23 ]]; then
        kmin=23
    fi
fi

# WRITE LOG FILE
echo "Assembly software: $assembly_soft" > $log_file
echo "kmer range: $kmin - $kmax" >> $log_file
echo "kmer step: $kstep" >> $log_file

# pair_end = 0 => single end assembly
# pair_end = 1 => pair end assembly
for dir in $(ls -d */); do
    if [[ $rlen_var -eq 1 ]]; then
        read_length=$(echo $dir | awk -F'_' '{print $2}' | sed 's/\///g')
    fi
    if [[ $flen_var -eq 1 ]]; then
        fragment_size=$(echo $dir | awk -F'_' '{print $2}' | sed 's/\///g')
    fi
    cd $dir
    base_dir=$PWD
    if [[ $pair_end -eq 0 ]]; then
        fastq_expr="*.fastq"
    else
        fastq_expr="*_1.fastq"
    fi
    for fq_file in $(ls -p -d $fastq_expr | fgrep -v '/'); do
        mkdir -p assembly_$fq_file
        cd assembly_$fq_file
            rm -rf $assembly_soft
            mkdir -p $assembly_soft
            cd $assembly_soft
                if [[ $pair_end -eq 0 ]]; then
                    #Single end assembly
                    case $n_as in
                        1) ## OASES ##
                            /usr/bin/time -o resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' oases_pipeline.py -g 27 -m $kmin -M $kmax -s $kstep -o k -d "-short -fastq $base_dir/$fq_file";
                            for k in $(seq $kmin $kstep $kmax); do
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
                            done;
                            fasta_size_filter.pl --min_size_cutoff=200 --input_file=kMerged/transcripts.fa --output_file=kMerged/transcripts_200.fa
                            rm k*/LastGraph k*/Graph2 k*/PreGraph k*/Roadmaps k*/Sequences k*/Log k*/contigs.fa k*/contig-ordering.txt k*/stats.txt;;
                        2) ## TRINITY ##
                            for k in $(seq $kmin $kstep $kmax); do
                                mkdir -p trinity_k_$k;
			        mkdir k_$k;
                                while [ ! -f trinity_k_$k/Trinity.fasta ]; do /usr/bin/time -o k_$k/resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' trinity-2.1.1 --seqType fq --max_memory 50G --single $base_dir/$fq_file --CPU $n_cpus --output trinity_k_$k/ --KMER_SIZE $k --min_contig_length $k --no_version_check; done
                                mv trinity_k_$k/Trinity.fasta k_$k/transcripts.fa;
                                rm -rf trinity_k_$k/;
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
                            done;;
                        3) ## ABYSS ##
                            abyss_path=''; abyss_prefix='';
                            for k in $(seq $kmin $kstep $kmax); do
                                abyss_path="$abyss_path aux_$k/transabyss-1.fa";
                                abyss_prefix="$abyss_prefix aux_$k/transabyss-1.";
                                mkdir -p aux_$k; mkdir -p k_$k;
                                /usr/bin/time -o k_$k/resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' transabyss --se $base_dir/$fq_file -k $k --outdir aux_$k --threads $n_cpus;
			        cp aux_$k/transabyss-1.fa k_$k/transcripts.fa;
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
                            done;
                            mkdir -p kMerged;
                            /usr/bin/time -o kMerged/resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' transabyss-merge --mink $kmin --maxk $kmax --prefix $abyss_prefix --threads $n_cpus $abyss_path ;
                            mv transabyss-merged.fa kMerged/transcripts.fa;
                            fasta_size_filter.pl --min_size_cutoff=200 --input_file=kMerged/transcripts.fa --output_file=kMerged/transcripts_200.fa;
                            rm -rf aux_*;;
                        4) ## SOAP ##
                            for k in $(seq $kmin $kstep $kmax); do
                                mkdir -p aux_$k; mkdir -p k_$k;
                                cd aux_$k;
                                echo -e "max_rd_len=$read_length\n[LIB]\nq=$base_dir/$fq_file" > config_file;
                                /usr/bin/time -o ../k_$k/resources.stats.1 --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' SOAPdenovo-Trans-127mer pregraph -s config_file -K $k -o transcripts -p $n_cpus; 
                                /usr/bin/time -o ../k_$k/resources.stats.2 --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' SOAPdenovo-Trans-127mer contig -g transcripts;
                                cd ..;
			        mv aux_$k/transcripts.contig k_$k/transcripts.fa;
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
			        rm -rf aux_$k;
                            done;;
                        5) ## BRIDGER ##
                            for k in $(seq $kmin $kstep $kmax); do
                                mkdir -p k_$k;
                                /usr/bin/time -o k_$k/resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' Bridger.pl --seqType fq --single $base_dir/$fq_file --kmer_length $k --CPU $n_cpus --output aux_$k;
			        mv aux_$k/Bridger.fasta k_$k/transcripts.fa;
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
			        rm -rf aux_$k;
                            done;;
                        6) ## IDBA ##
                            /usr/bin/time -o resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' idba_tran --long_read $base_dir/$fq_file --mink $kmin --maxk $kmax --step $kstep --num_threads $n_cpus --max_isoforms 1000 --out kaux;
                            for k in $(seq $kmin $kstep $kmax); do
                                mkdir -p k_$k;
                                mv kaux/contig-$k.fa k_$k/transcripts.fa;
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
                            done;
			    mkdir kMerged;
                            mv kaux/contig.fa kMerged/transcripts.fa;
                            fasta_size_filter.pl --min_size_cutoff=200 --input_file=kMerged/transcripts.fa --output_file=kMerged/transcripts_200.fa;
			    rm -rf kaux;;
                    esac
                else
                    #Pair end assembly
                    insert_size=$(( $fragment_size - 2*$read_length ))
                    fq_file2=$(echo $fq_file | sed 's/_1.fastq/_2.fastq/g')
                    case $n_as in
                        1) ## OASES ##
                            /usr/bin/time -o resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' oases_pipeline.py -g 27 -m $kmin -M $kmax -s $kstep -o k -d "-shortPaired -fastq -separate $base_dir/$fq_file $base_dir/$fq_file2" -p "-ins_length $fragment_size -ins_length_sd $fragment_size_sd"; 
                            for k in $(seq $kmin $kstep $kmax); do
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
                            done;
                            fasta_size_filter.pl --min_size_cutoff=201 --input_file=kMerged/transcripts.fa --output_file=kMerged/transcripts_200.fa;
                            rm k*/LastGraph k*/Graph2 k*/PreGraph k*/Roadmaps k*/Sequences k*/Log k*/contigs.fa k*/contig-ordering.txt k*/stats.txt;;
                        2) ## TRINITY ##
                            for k in $(seq $kmin $kstep $kmax); do
                                mkdir -p trinity_k_$k;
			        mkdir k_$k;
                                while [ ! -f trinity_k_$k.Trinity.fasta ]; do /usr/bin/time -o k_$k/resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' trinity-2.1.1 --seqType fq --max_memory 50G --left $base_dir/$fq_file --right $base_dir/$fq_file2 --CPU $n_cpus --output trinity_k_$k --full_cleanup --group_pairs_distance $(( $fragment_size + $fragment_size_sd )) --KMER_SIZE $k --no_version_check; done
                                mv trinity_k_$k.Trinity.fasta k_$k/transcripts.fa;
                                rm -rf trinity_k_$k/;
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
                            done;;
                        3) ## ABYSS ##
                            abyss_path=''; abyss_prefix='';
                            for k in $(seq $kmin $kstep $kmax); do
                                abyss_path="$abyss_path aux_$k/transabyss-1.fa";
                                abyss_prefix="$abyss_prefix aux_$k/transabyss-1.";
                                mkdir -p aux_$k; mkdir -p k_$k;
                                /usr/bin/time -o k_$k/resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' transabyss --pe $base_dir/$fq_file $base_dir/$fq_file2 -k $k --outdir aux_$k --threads $n_cpus; 
			        cp aux_$k/transabyss-1.fa k_$k/transcripts.fa;
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
                            done;
                            mkdir -p kMerged; 
                            /usr/bin/time -o kMerged/resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' transabyss-merge --mink $kmin --maxk $kmax --prefix $abyss_prefix --threads $n_cpus $abyss_path ;
                            mv transabyss-merged.fa kMerged/transcripts.fa;
                            fasta_size_filter.pl --min_size_cutoff=200 --input_file=kMerged/transcripts.fa --output_file=kMerged/transcripts_200.fa;
                            rm -rf aux_*;;
                        4) ## SOAP ##
                            for k in $(seq $kmin $kstep $kmax); do
                                mkdir -p aux_$k; mkdir -p k_$k;
                                cd aux_$k;
                                echo -e "max_rd_len=$read_length\n[LIB]\navg_ins=$fragment_size\nq1=$base_dir/$fq_file\nq2=$base_dir/$fq_file2" > $PWD/config_file; 
                                /usr/bin/time -o ../k_$k/resources.stats.1 --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' SOAPdenovo-Trans-127mer pregraph -s config_file -K $k -o transcripts -p $n_cpus; 
                                /usr/bin/time -o ../k_$k/resources.stats.2 --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' SOAPdenovo-Trans-127mer contig -g transcripts;
                                cd ..;
			        mv aux_$k/transcripts.contig k_$k/transcripts.fa;
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
			        rm -rf aux_$k;
                            done;;
                        5) ## BRIDGER ##
                            for k in $(seq $kmin $kstep $kmax); do
                                mkdir -p k_$k;
                                /usr/bin/time -o k_$k/resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' Bridger.pl --seqType fq --left $base_dir/$fq_file --right $base_dir/$fq_file2 --kmer_length $k --CPU $n_cpus --output aux_$k --pair_gap_length $insert_size; 
			        mv aux_$k/Bridger.fasta k_$k/transcripts.fa;
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
			        rm -rf aux_$k;
                            done;;
                        6) ## IDBA ##
                            fq2fa --merge $base_dir/$fq_file $base_dir/$fq_file2 reads.fasta
                            /usr/bin/time -o resources.stats --format='User time:%U\nCPU time:%S\nReal time:%E\nCPU Usage:%P\nMaximum memory usage:%M\nAverage memory usage:%k' idba_tran --long_read reads.fasta --mink $kmin --maxk $kmax --step $kstep --num_threads $n_cpus --max_isoforms 1000 --out kaux;
                            for k in $(seq $kmin $kstep $kmax); do
                                mkdir -p k_$k;
                                mv kaux/contig-$k.fa k_$k/transcripts.fa;
                                fasta_size_filter.pl --min_size_cutoff=200 --input_file=k_$k/transcripts.fa --output_file=k_$k/transcripts_200.fa;
                            done;
			    mkdir kMerged;
                            mv kaux/contig.fa kMerged/transcripts.fa;
                            fasta_size_filter.pl --min_size_cutoff=200 --input_file=kMerged/transcripts.fa --output_file=kMerged/transcripts_200.fa;
			    rm -rf kaux reads.fasta;;
                    esac  
                fi
            cd ..
        cd ..
    done
    cd ..    
done

