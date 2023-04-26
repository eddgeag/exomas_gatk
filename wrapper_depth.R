library(tools)
args <- commandArgs(trailingOnly = TRUE)

directorio <- args[1]
directorio_salida <- args[2]

archivos <- list.files(directorio,full.names = T,pattern="*.bam")
if(!dir.exists(directorio_salida)){
  dir.create(directorio_salida)
}
computar_cobertura <- function(archivos){
  
  for(archivo_n in 1:length(archivos)){
    archivo <- archivos[archivo_n]
    directorio_salida_2 <- file.path(directorio_salida,
                                     file_path_sans_ext(basename(archivo)))
    if(!dir.exists(directorio_salida_2)){
	dir.create(directorio_salida_2)
	}
    archivo <- archivos[archivo_n]
    
    command <- paste("./compute_depth.sh",archivo,
                     ">",
                     file.path(directorio_salida_2,paste0(file_path_sans_ext(basename(archivo))))
                     ,".txt")
    system(command = command)
       
  }
  
  
}

computar_cobertura(archivos)
