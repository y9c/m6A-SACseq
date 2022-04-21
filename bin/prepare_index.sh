#!/usr/bin/env sh
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-04-20 22:28

usage_error() {
  echo >&2 "$(basename $0):  $1"
  exit 2
}
assert_argument() { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }

species="human"
outdir="ref_test"
logfile="$outdir/indexing.log"

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

mkdir -p "${outdir}"

## prepare spike index
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
makeblastdb -in ${outdir}/spike_degenerate.fa -dbtype nucl -out ${outdir}/spike_degenerate 2>&1 >/dev/null
# expand ATGC
cat ${outdir}/spike_degenerate.fa |
  paste - - |
  awk 'BEGIN{Ns["A"]=1;Ns["T"]=2;Ns["G"]=3;Ns["C"]=4}{split($2,a,"NN");for(b1 in Ns)for(b2 in Ns)for(b3 in Ns)for(b4 in Ns)print $1"_"b1""b2"_"b3""b4"\n"a[1]""b1""b2""a[2]""b3""b4""a[3]}' >${outdir}/spike_expand.fa
bowtie2-build ${outdir}/spike_expand.fa ${outdir}/spike_expand 1>${logfile} 2>${logfile}

# prepare contamination index

# prepare rRNA/ smallRNA/ genome index (base on different species)
if [[ "$species" == "human" ]]; then
  # prepare rRNA index for human
  # prepare smallRNA index for human
  # prepare genome index for human
  echo ${species} > ${logfile}
elif [[ "$species" == "mouse" ]]; then
  # prepare rRNA index for mosue
  # prepare smallRNA index for mosue
  # prepare genome index for mosue
  echo ${species} > ${logfile}
else
  echo ${species} > ${logfile}
fi


echo "COPY THE CONFIGURE BELLOW"
echo "         ↓ ↓ ↓         \n"
echo '\033[0;32m'
echo "references:"
echo "  spike:"
echo "    fa: ${outdir}/spike_expand.fa"
echo "    bt2: ${outdir}/spike_expand"
echo "  spikeN:"
echo "    fa: ${outdir}/spike_degenerate.fa"
echo "    blast: ${outdir}/spike_degenerate"
echo "  rRNA:"
echo "    fa:"
echo "    bt2:"
echo "  smallRNA:"
echo "    fa:"
echo "    bt2:"
echo "  genome:"
echo "    fa:"
echo "    fai:"
echo "    gtf:"
echo "    gtf_collapse:"
echo "    star:"
echo "  contamination:"
echo "    fa:"
echo "    bt2:"
echo '\033[0m'
