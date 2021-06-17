#!/bin/bash

# Este script es para obtener el nombre de los bam que pertenecen a cada archivo (con más de 1000 reads). Crea un txt llamado réplicas con los nombres de las guías correspondientes. 

echo "Guia_ID	Bam_names"  > /mnt/c/Users/yolib/Documents/TFM/Linf_T/Datos/all_target_replicas.txt 

cd /mnt/c/Users/yolib/Documents/TFM/Linf_T/Datos/Datos

carpetas=$(ls)

#echo $carpetas
for file in $carpetas # voy a entrar en las carpetas con los datos una a una 
do 

  guide_name=$(echo $file | cut -d '-' -f2)
  
  cd $file # entro en la carpeta
  
  bam_files=$(ls | grep '.bam$') # obtengo también el nombre de los bam
  cadena_bam=$(echo '') # inicio la variable cadena_bam
  for bam_file in $bam_files  # recorro los bam files que hay en cada carpeta
	  do
	  	n_alignm=$(samtools view $bam_file | wc -l) # veo cuantos alineamientos tiene cada bam. No tienen header 
	  	if (("$n_alignm" >= "1000")) # si tiene más de 1000 alineamientos, lo guardo en cadena_bam. El nombre de cada archivo va separado por ";"
  		then
  			#echo $bam_file
  			cadena_bam=$(echo $cadena_bam";"$bam_file)
  			
  		fi
  		
	  done
   
   if ! [ -z $cadena_bam ] # si la cadena no está vacía
   then 
   	echo "$guide_name	$cadena_bam"  >> /mnt/c/Users/yolib/Documents/TFM/Linf_T/Datos/all_target_replicas.txt 
   	#echo $guide_name
   fi
   cd ..
   
done
