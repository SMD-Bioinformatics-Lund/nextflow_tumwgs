#!/bin/bash
determine_shards_from_bam(){
    local bam step tag chr len pos end
    bam="$1"
    step="$2"
    pos=1
    samtools view -H $bam |\
    while read tag chr len; do
        [ $tag == '@SQ' ] || continue
        chr=$(expr "$chr" : 'SN:\(.*\)')
        len=$(expr "$len" : 'LN:\(.*\)')
        while [ $pos -le $len ]; do
            end=$(($pos + $step - 1))
            if [ $pos -lt 0 ]; then
                start=1
            else
                start=$pos
            fi
            if [ $end -gt $len ]; then
                echo -n "$chr:$start-$len "
                pos=$(($pos-$len))
                break
            else
                echo "$chr:$start-$end"
                pos=$(($end + 1))
            fi
        done
    done
    echo "NO_COOR"
}

determine_shards_from_dict(){
    local bam step tag chr len pos end
    dict="$1"
    step="$2"
    pos=1
    cat $dict |\
    while read tag chr len UR; do
        [ $tag == '@SQ' ] || continue
        chr=$(expr "$chr" : 'SN:\(.*\)')
        len=$(expr "$len" : 'LN:\(.*\)')
        while [ $pos -le $len ]; do
            end=$(($pos + $step - 1))
            if [ $pos -lt 0 ]; then
                start=1
            else
                start=$pos
            fi
            if [ $end -gt $len ];then
                echo -n "$chr:$start-$len "
                pos=$(($pos-$len))
                break
            else
                echo "$chr:$start-$end"
                pos=$(($end + 1))
            fi
        done
    done
    echo "NO_COOR"
}


determine_shards_from_fai(){
    local bam step tag chr len pos end
    fai="$1"
    step="$2"
    pos=1
    cat $fai |\
    while read chr len other; do
        while [ $pos -le $len ]; do
            end=$(($pos + $step - 1))
            if [ $pos -lt 0 ]; then
                start=1
            else
                start=$pos
            fi
            if [ $end -gt $len ]; then
                echo -n "$chr:$start-$len "
                pos=$(($pos-$len))
                break
            else
                echo "$chr:$start-$end"
                pos=$(($end + 1))
            fi
        done
    done
    echo "NO_COOR"
}

if [ $# -eq 2 ]; then
    filename=$(basename "$1")
    extension="${filename##*.}"
    if [ "$extension" = "fai" ]; then
        determine_shards_from_fai $1 $2
    elif [ "$extension" = "bam" ]; then
        determine_shards_from_bam $1 $2
    elif [ "$extension" = "dict" ]; then
        determine_shards_from_dict $1 $2
    fi
else
    echo "usage $0 file shard_size"
fi


##  ./generate_shards.sh GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai 100000000

## ./generate_shards.sh GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai 200000000