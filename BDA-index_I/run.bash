 #!/bin/bash

X=("english" "proteins" "sources")
Y=(16 32 64 128 256 512 1024)

for x in "${X[@]}"; do
    for y in "${Y[@]}"; do
        directory="./testFiles/${x}_${y}_files"
        mkdir -p "$directory"

        command="/usr/bin/time -v ./bda-index_I ./data/texts/${x}.200MB $y ./data/pattern/${x}_${y} ${directory}/${x}_${y}_out 960000 250000 ${directory}/${x}_${y}_idx"

        echo "Running command: $command >> ${x}_${y}_output.txt 2>&1"
        $command >> "${x}_${y}_output.txt" 2>&1

        echo "--------------------------------------"
    done
done

