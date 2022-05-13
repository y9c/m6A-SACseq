#!/usr/bin/env sh
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-04-20 22:28

dbpath="https://raw.githubusercontent.com/chelab/db/main/reference_sequence"
species="human"
outdir="ref_test"
logfile="$outdir/indexing.log"
threads=16

if command -v axel >/dev/null 2>&1; then
  downloader="axel -n ${threads} -q -o"
else
  downloader="wget -q -O"
fi

usage_error() {
  echo >&2 "$(basename $0):  $1"
  exit 2
}
assert_argument() { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }

# One loop, nothing more.
if [ "$#" != 0 ]; then
  EOL=$(echo '\01\03\03\07')
  set -- "$@" "$EOL"
  while [ "$1" != "$EOL" ]; do
    opt="$1"
    shift
    case "$opt" in

    # Your options go here.
    -q | --quite) quite='true' ;;
    -s | --threads)
      assert_argument "$1" "$opt"
      threads="$1"
      shift
      ;;
    -s | --species)
      assert_argument "$1" "$opt"
      species="$1"
      shift
      ;;
    -i | --spikein)
      assert_argument "$1" "$opt"
      spikein="$1"
      shift
      ;;
    -o | --outdir)
      assert_argument "$1" "$opt"
      outdir="$1"
      shift
      ;;

    # Arguments processing. You may remove any unneeded line after the 1st.
    - | '' | [!-]*) set -- "$@" "$opt" ;;                             # positional argument, rotate to the end
    --*=*) set -- "${opt%%=*}" "${opt#*=}" "$@" ;;                    # convert '--name=arg' to '--name' 'arg'
    -[!-]?*) set -- $(echo "${opt#-}" | sed 's/\(.\)/ -\1/g') "$@" ;; # convert '-abc' to '-a' '-b' '-c'
    --) while [ "$1" != "$EOL" ]; do
      set -- "$@" "$1"
      shift
    done ;;                                             # process remaining arguments as positional
    -*) usage_error "unknown option: '$opt'" ;;         # catch misspelled options
    *) usage_error "this should NEVER happen ($opt)" ;; # sanity test for previous patterns

    esac
  done
  shift # $EOL
fi

download_db() {
  # Download the database (gz compressed file).
  local inurl="$1"
  local outfile="$2"
  if [ -f "${outfile}" ]; then
    echo "Reference file: ${outfile} exist. Do you want overwrite it? (y/N)"
    local yn="N"
    read yn
    if [ "$yn" = "${yn#[Nn]}" ]; then
      return
    fi
  fi
  echo "$(date -u)  Downloading db: ${outfile}..."
  ${downloader} ${outfile}.gz ${inurl} 2>&1 >>${logfile}
  gunzip -f ${outfile}.gz 2>&1 >>${logfile}
}

if [ -d "${outdir}" ]; then
  echo "Directory ${outdir} exist. Do you want overwrite it? (Y/n)"
  yn="Y"
  read yn
  if [ "$yn" != "${yn#[Nn]}" ]; then
    exit 0
  fi
else
  echo "The referece directory is not exist. Creating a new one..."
  mkdir -p "${outdir}"
fi

echo "$(date -u)  Start to build index..." >${logfile}

## prepare spike index
echo "$(date -u)  Preparing spike index..."
if [ -z ${spikein+x} ]; then
  cat <<EOF >${outdir}/spike_degenerate.fa
>probe_0
TATCTGTCTCGACGTNNANNGGCCTTTGCAACTAGAATTACACCATAATTGCT
>probe_25
TATCTGTCTCGACGTNNANNGGCATTCAAGCCTAGAATTACACCATAATTGCT
>probe_50
TATCTGTCTCGACGTNNANNGGCGAGGTGATCTAGAATTACACCATAATTGCT
>probe_75
TATCTGTCTCGACGTNNANNGGCTTCAACAACTAGAATTACACCATAATTGCT
>probe_100
TATCTGTCTCGACGTNNANNGGCGATGGTTTCTAGAATTACACCATAATTGCT
EOF
else
  cp $spike ${outdir}/spike_degenerate.fa
fi
makeblastdb -in ${outdir}/spike_degenerate.fa -dbtype nucl -out ${outdir}/spike_degenerate 2>&1 >>${logfile}
# expand ATGC
cat ${outdir}/spike_degenerate.fa |
  paste - - |
  awk 'BEGIN{Ns["A"]=1;Ns["T"]=2;Ns["G"]=3;Ns["C"]=4}{split($2,a,"NN");for(b1 in Ns)for(b2 in Ns)for(b3 in Ns)for(b4 in Ns)print $1"_"b1""b2"_"b3""b4"\n"a[1]""b1""b2""a[2]""b3""b4""a[3]}' >${outdir}/spike_expand.fa
bowtie2-build --threads {threads} ${outdir}/spike_expand.fa ${outdir}/spike_expand 1>>${logfile} 2>>${logfile}

# prepare contamination index
echo "$(date -u)  Preparing contamination index..."
download_db "${dbpath}/contamination.fa.gz" ${outdir}/contamination.fa
bowtie2-build --threads {threads} ${outdir}/contamination.fa ${outdir}/contamination 1>>${logfile} 2>>${logfile}

# prepare rRNA/ smallRNA/ genome index (base on different species)
if [ "$species" = "human" ]; then
  species_prefix="Homo_sapiens.GRCh38"
  url_gtf="http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz"
  url_fa="http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
elif [ "$species" = "mouse" ]; then
  species_prefix="Mus_musculus.GRCm38"
  url_gtf="http://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.gtf.gz"
  url_fa="http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
else
  echo "ERROR: Only support human/mouse in current version"
  exit 0
fi

# prepare rRNA index
echo "$(date -u)  Preparing rRNA (${species}) index..."
download_db "${dbpath}/${species_prefix}.rRNA.fa.gz" ${outdir}/rRNA_${species}.fa
bowtie2-build --threads {threads} ${outdir}/rRNA_${species}.fa ${outdir}/rRNA_${species} 1>>${logfile} 2>>${logfile}

# prepare smallRNA index
echo "$(date -u)  Preparing smallRNA (${species}) index..."

download_db "${dbpath}/${species_prefix}.smallRNA.fa.gz" ${outdir}/smallRNA_${species}.fa
bowtie2-build --threads {threads} ${outdir}/smallRNA_${species}.fa ${outdir}/smallRNA_${species} 1>>${logfile} 2>>${logfile}

# prepare genome index
echo "$(date -u)  Preparing genomne (${species}) index..."
download_db ${url_gtf} ${outdir}/genome_${species}.gtf
download_db ${url_fa} ${outdir}/genome_${species}.fa
# collpase gtf
#### python3 collapse_annotation.py gencode.v26.GRCh38.annotation.gtf gencode.v26.GRCh38.genes.gtf
# build star index
#### STAR (TODO: let docker to do this, since user might not hve STAR installed)

echo "\nCOPY THE CONFIGURE BELLOW"
echo "         ↓ ↓ ↓         \n"
echo '\033[0;32m'
echo "references:"
echo "  spike:"
echo "    fa: ${outdir}/spike_expand.fa"
echo "    bt2: ${outdir}/spike_expand"
echo "  spikeN:"
echo "    fa: ${outdir}/spike_degenerate.fa"
echo "    blast: ${outdir}/spike_degenerate"
echo "  contamination:"
echo "    fa: ${outdir}/contamination.fa"
echo "    bt2: ${outdir}/contamination"
echo "  rRNA:"
echo "    fa: ${outdir}/rRNA_${species}.fa"
echo "    bt2: ${outdir}/rRNA_${species}"
echo "  smallRNA:"
echo "    fa: ${outdir}/smallRNA_${species}.fa"
echo "    bt2: ${outdir}/smallRNA_${species}"
echo "  genome:"
echo "    fa: ${outdir}/genome_${species}.fa"
echo "    fai: ${outdir}/genome_${species}"
echo "    gtf: ${outdir}/genome_${species}.gtf"
echo "    gtf_collapse: ${outdir}/genome_collapse_${species}.gtf"
echo "    star: ${outdir}/genome_${species}"
echo '\033[0m'
