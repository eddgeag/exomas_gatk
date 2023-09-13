#!/bin/bash

nohup Rscript wgs_flujo.R $1  > salida.out 2>&1 &

# Almacena la cadena de argumento en una variable
cadena="$1"

# Divide la cadena en variables usando el espacio como delimitador
IFS=" " read -ra variables <<< "$cadena"

for var in "${variables[@]}"; do
    # Construye y ejecuta el comando cp para cada variable
    nohup cp -r -v "/repositorio/exomas/pipeline/$var" /mnt/nas_mount/Respaldo_NAS/resultados > salida.out 2>&1 &
done

