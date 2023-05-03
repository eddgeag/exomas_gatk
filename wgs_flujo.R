#### CARGAMOS LAS LIBRERIAS



load.libs <- c(
  "data.table",
  "limma",
  "grid",
  "egg",
  "Rbwa",
  "ggpubr",
  "tools",
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3",
  "RISmed",
  "ggplot2",
  "dplyr")

lapply(load.libs, require, character.only = TRUE)





control_calidad <- function(fastq_dir, output_dir) {
  ### Creamos directorio de QC si no existe
  if (!dir.exists(file.path(output_dir, "QC"))) {
    dir.create(file.path(output_dir, "QC"))
  }
  ### Si no existen archivos de salida
  if (length(list.files(file.path(output_dir, "QC"))) == 0) {
    command <-
      paste(
        "fastqc -t 4 ",
        paste0(fastq_dir, "/*.", unique(file_ext(
          list.files(fastq_dir)
        ))),
        "-o",
        file.path(output_dir, "QC")
      )
    system(command, intern = T)
  } else{
    message("Ya se ha hecho el control de Calidad")
  }
  
  
  
}

fn_exists_fasta <- function(folder_fasta) {
  extension = unlist(lapply(list.files(folder_fasta, pattern = "fa"), function(x)
    file_ext(x)))
  extension_fa <- extension[grep("^fa$", extension)]
  extension_fasta <- extension[grep("^fasta$", extension)]
  ## Es fa o fasta, y existe ?
  if (length(extension_fa) == 0 && length(extension_fasta) == 0) {
    stop("No existe archivo de referencia")
    
  } else if (length(extension_fa) != 0 ||
             length(extension_fasta) != 0) {
    if (length(extension_fa) != 0) {
      extension <- extension_fa
    } else if (length(extension_fasta) != 0) {
      extension <- extension_fasta
    }
    
    
  }
  ## el archivo fasta ?
  fasta_file <-
    list.files(folder_fasta,
               pattern = paste0(".", extension, "$"),
               full.names = T)
  return(fasta_file)
  
  
}



index_fasta_samtools <- function(folder_fasta = folder_fasta) {
  ## Vemos si existe el archivo con patron fa
  
  fasta_file <- fn_exists_fasta(folder_fasta)
  if (length(list.files(folder_fasta, pattern = "fai$")) == 0) {
    command <- paste("samtools faidx", fasta_file)
    print(command)
    system(command, intern = T)
  } else{
    message("Ya esta el index fai")
  }
  
  
  
  
}


index_bwa <- function(folder_fasta) {
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  extension <- c("amb", "ann", "bwt", "pac", "sa")
  
  if (length(file.exists(list.files(folder_fasta, extension[1], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[2], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[3], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[4], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[5], full.names = T))) ==
      0) {
    ## no existen bwa index
    print("Creando ficheros índices para bwa mem...")
    
    bwa_build_index(fasta_file, index_prefix = fasta_file)
    
  } else if (length(file.exists(list.files(folder_fasta, extension[1], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[2], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[3], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[4], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[5], full.names = T))) !=
             0) {
    message("Ya se han creado los ficheros para el alineamiento")
    
  }
  
  
}


bwamem <- function(fastq_dir = fastq_dir ,
                   folder_fasta = folder_fasta) {
  ### conseguimos el archivo fasta
  
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  ####fai index exist ?######
  
  if (!length(file.exists(list.files(dirname(fasta_file), "fai")))) {
    print("### generating fai index...")
    
    index_fasta_samtools(input_directory = reference_genome)
    
  }
  
  ## buscamos los archivos fastq
  fastq_files <- list.files(fastq_dir, full.names = T)
  ## Archivos fastq concatenados
  fastq_full_path_files <-
    c(fastq_1 = fastq_files[1], fastq_2 = fastq_files[2])
  ## Creamos directorio de mapeo
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  
  output_file_name <-
    file_path_sans_ext(output_file_name[length(output_file_name)])
  ### Creando el directorio de mapeo
  
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  if (!dir.exists(mapping_output_dir)) {
    dir.create(mapping_output_dir)
  }
  ### out file name sam file
  output_file_sam <-
    file.path(mapping_output_dir, paste0(output_file_name, ".sam"))
  ### out file name bam file
  
  
  output_file_bam <-
    file.path(mapping_output_dir, paste0(output_file_name, ".bam"))
  ### out file name bam sort file
  
  output_file_sorted_bam <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.bam"))
  ### Si no esta mapepeado, mapear
  
  if (length(list.files(mapping_output_dir)) == 0) {
    print("#### MAPPING...#####")
    bwa_mem(
      index_prefix = (fasta_file),
      fastq_files = fastq_full_path_files,
      sam_file = output_file_sam,
      type = "paired"
    )
    ### Si ya se ha mapeado pero no esta el bam, crearlo
    
    print("#### SAM TO BAM")
    
    command_sam_to_bam <-
      paste("samtools view -S -b",
            output_file_sam,
            "-o",
            output_file_bam)
    
    print(command_sam_to_bam)
    system(command_sam_to_bam, intern = T)
    ### Si ya esta el bam, sortearlo
    
    print("### BAM to sorted BAM")
    command_out_bam_sorted <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk SortSam -CREATE_INDEX true -INPUT",
        output_file_bam,
        "-OUTPUT",
        output_file_sorted_bam,
        "-SORT_ORDER coordinate -VALIDATION_STRINGENCY STRICT"
      )
    
    print(command_out_bam_sorted)
    system(command_out_bam_sorted, intern = T)
    ### AHORA SE PROCEDERIA A MERGE LOS BAMS EN CASO DE TRIO
    
  } else if (file.exists(output_file_sam) &&
             !file.exists(output_file_bam) &&
             !file.exists(output_file_sorted_bam)) {
    print("#### SAM TO BAM")
    
    command_sam_to_bam <-
      paste("samtools view -S -b",
            output_file_sam,
            "-o",
            output_file_bam)
    
    print(command_sam_to_bam)
    system(command_sam_to_bam, intern = T)
    ### Si ya esta el bam, sortearlo
    
    print("### BAM to sorted BAM")
    command_out_bam_sorted <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk SortSam -CREATE_INDEX true -INPUT",
        output_file_bam,
        "-OUTPUT",
        output_file_sorted_bam,
        "-SORT_ORDER coordinate -VALIDATION_STRINGENCY STRICT"
      )
    
    print(command_out_bam_sorted)
    system(command_out_bam_sorted, intern = T)
  }
  else if (file.exists(output_file_sam) &&
           file.exists(output_file_bam) &&
           !file.exists(output_file_sorted_bam)) {
    print("### BAM to sorted BAM")
    command_out_bam_sorted <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk SortSam -CREATE_INDEX true -INPUT",
        output_file_bam,
        "-OUTPUT",
        output_file_sorted_bam,
        "-SORT_ORDER coordinate -VALIDATION_STRINGENCY STRICT"
      )
    
  } else if (file.exists(output_file_sam) &&
             file.exists(output_file_bam) &&
             file.exists(output_file_sorted_bam)) {
    print("YA SE HA MAPEADO")
    
  }  else{
    stop(message("No se ha mapeado bien", call = T))
  }
  
  
}


markdups <- function(output_dir = output_dir,
                     fastq_dir = fastq_dir) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  
  ## Creamos directorio de mapeo, renombramos al archivo
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  ## quitamos la extension del archivo
  output_file_name <- file_path_sans_ext(output_file_name)
  ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  ## nombramos al arhivo bam sorteado
  bam_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.bam"))
  ## nombramos al archivo bam marcado con duplciados
  mark_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup.bam"))
  ## nombramos al archivo con las meteicas
  metrics_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup.txt"))
  ## Si no existe el archivo crearlo
  if (!file.exists(mark_file)) {
    command <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk  MarkDuplicates CREATE_INDEX=true INPUT=",
        bam_file,
        " VALIDATION_STRINGENCY=STRICT OUTPUT=",
        mark_file,
        "M=",
        metrics_file
      )
    print(command)
    system(command = command, intern = T)
    
  } else{
    message("Ya se han marcado los duplicados")
  }
  ## llamamos al comando
  
  
}



create_dict <- function(folder_fasta) {
  ## llamamos al archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  fasta_dict <- paste0(file_path_sans_ext(fasta_file), ".dict")
  if (!file.exists(fasta_dict)) {
    ## creamos el diccionario del fasta
    command <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk CreateSequenceDictionary -R",
        fasta_file
      )
    system(command, intern = T)
    
  } else{
    message("Ya se ha creado el diccionario fasta")
  }
  
  
}


creacion_readgroup <-
  function(output_dir = output_dir,
           fastq_dir = fastq_dir) {
    ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
    mapping_output_dir <- file.path(output_dir, "mapping_output")
    fastq_files <- list.files(fastq_dir, full.names = F)
    output_file_name <-
      unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
    ## quitamos la extension del archivo
    output_file_name <- file_path_sans_ext(output_file_name)
    ## obtenemos el archivo marcado con duplicados
    mark_file <-
      file.path(mapping_output_dir,
                paste0(output_file_name, ".sorted.mark_dup.bam"))
    ## quitamos la extension y renombramos al archivo de salida
    out_file <- paste0(file_path_sans_ext(mark_file), "_RG.bam")
    
    ## si ya se ha creado el archivo con grupo
    
    if (!file.exists(out_file)) {
      command <-
        paste(
          "~/tools/exomas_tools/gatk-4.3.0.0/gatk AddOrReplaceReadGroups I=",
          mark_file,
          "O=",
          out_file,
          "RGID=1 RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM=1"
        )
      system(command = command, intern = T)
      
    } else{
      message("Ya estan los grupos")
    }
    
    
    
    
  }




base_recalibrator <-
  function(folder_fasta,
           output_dir,
           folder_data_gatk,
           fastq_dir) {
    ## llamamos al knwon sites file vcf
    known_sites_file <-
      list.files(folder_data_gatk,
                 pattern = ".vcf$",
                 full.names = T)
    ## llamamos al archivo fasta
    fasta_file <- fn_exists_fasta(folder_fasta)
    ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
    mapping_output_dir <- file.path(output_dir, "mapping_output")
    fastq_files <- list.files(fastq_dir, full.names = F)
    output_file_name <-
      unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
    ## quitamos la extension del archivo
    output_file_name <- file_path_sans_ext(output_file_name)
    ## obtenemos el archivo marcado con duplicados
    RG_file <-
      file.path(mapping_output_dir,
                paste0(output_file_name, ".sorted.mark_dup_RG.bam"))
    ## quitamos la extension y renombramos al archivo de salida
    out_file <- file.path(dirname(RG_file), "recal_data.table")
    ## si no existe la tabla de recalibracion, la calculamos
    
    if (!file.exists(out_file)) {
      command <-
        paste(
          "~/tools/exomas_tools/gatk-4.3.0.0/gatk BaseRecalibrator -I",
          RG_file,
          " -R",
          fasta_file,
          " --known-sites",
          known_sites_file,
          " -O"  ,
          out_file
        )
      print(command)
      system(command = command, intern = T)
    } else{
      message("ya existe la tabla de base recalculator")
    }
    
    
  }


applybqsr <- function(folder_fasta, output_dir, fastq_dir) {
  ## llamamos al archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  ## quitamos la extension del archivo
  output_file_name <- file_path_sans_ext(output_file_name)
  ## obtenemos el archivo marcado con duplicados
  RG_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup_RG.bam"))
  
  recal_data.table <-
    file.path(dirname(RG_file), "recal_data.table")
  
  out_file <-
    file.path(
      dirname(recal_data.table),
      paste0(output_file_name, ".sorted.mark_dup_RG_bqsr.bam")
    )
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        " ~/tools/exomas_tools/gatk-4.3.0.0/gatk ApplyBQSR -I",
        RG_file ,
        "-R",
        fasta_file,
        " --bqsr-recal-file",
        recal_data.table,
        " -O",
        out_file
      )
    print(command)
    system(command = command, intern = T)
  } else{
    message("Ya se ha aplicado el bsqr")
  }
  
}


bam_statistics <- function(folder_fasta, fastq_dir, output_dir) {
  ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  ## bam file
  bam_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup_RG_bqsr.bam"))
  ## directorio de salida
  out_dir <-
    paste0(file_path_sans_ext(bam_file), ".CollectMultipleMetrics")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  out_file <- file.path(out_dir, "CollectMultipleMetrics")
  ## obtenemos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## verificacion de archivos en el directorio
  verificacion <- length(list.files(out_dir))
  ## comando estadistico
  
  if (verificacion < 2) {
    command <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk CollectMultipleMetrics R=",
        fasta_file,
        "I=",
        bam_file,
        "O=",
        out_file
      )
    
    system(command = command, intern = T)
    ## ahora juntamos las estadisticas
    command <- paste("multiqc", out_dir, "-o", out_dir)
    
    system(command = command, intern = T)
    
  } else{
    message("Ya se ha calculado la estadistica bam")
  }
  
  
}



haplotype_caller <- function(output_dir, folder_fasta, fastq_dir) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## recuperamos el ultimo archivo bam
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  ## bam file
  bam_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup_RG_bqsr.bam"))
  out_dir <- file.path(output_dir, "variantCalling")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  out_file <-
    file.path(out_dir, paste0(basename(output_file_name), ".g.vcf.gz"))
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk HaplotypeCaller -I",
        bam_file,
        "-R",
        fasta_file,
        "-ERC GVCF -O",
        out_file
      )
    system(command, intern = T)
  } else{
    message("Ya se han llamado a las variantes")
  }
  
  
}


genotypeGVCF <- function(folder_fasta, output_dir, fastq_dir) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## obtenemos el nombre del archivo
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  ## indicamos el archivo de entrada
  file_in <-
    file.path(output_dir,
              "variantCalling",
              paste0(output_file_name, ".g.vcf.gz"))
  ## indciamos el archivo de salida
  file_out <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz")
    )
  
  if (!file.exists(file_out)) {
    command <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk GenotypeGVCFs -R",
        fasta_file,
        "-V",
        file_in,
        "-O",
        file_out
      )
    system(command = command, intern = T)
  } else{
    message("Ya se ha calculado la probabilidad posterior de alelo no referente")
  }
  
  
}

variantRecallibrator <-
  function(fastq_dir,
           folder_fasta,
           folder_data_gatk,
           output_dir) {
    ## recuperamos el archivo fasta
    fasta_file <- fn_exists_fasta(folder_fasta)
    ## obtenemos el nombre del archivo
    fastq_files <- list.files(fastq_dir, full.names = F)
    output_file_name <-
      unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
    output_file_name <- file_path_sans_ext(output_file_name)
    
    in_file <-
      file.path(
        output_dir,
        "variantCalling",
        paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz")
      )
    
    snps_recal_file <-
      file.path(
        output_dir,
        "variantCalling",
        paste0(output_file_name, "_apply_genotypeGVCF.recal")
      )
    
    tranches_file <-
      file.path(
        output_dir,
        "variantCalling",
        paste0(output_file_name, "_apply_genotypeGVCF.g.tranches")
      )
    
    if (!file.exists(snps_recal_file) |
        !file.exists(tranches_file)) {
      command <-
        paste(
          "~/tools/exomas_tools/gatk-4.3.0.0/gatk VariantRecalibrator -V",
          in_file,
          " --resource:hapmap,known=false,training=true,truth=true,prior=15",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
          ),
          "--resource:omni,known=false,training=true,truth=true,prior=12",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
          ),
          "--resource:1000G,known=false,training=true,truth=false,prior=10",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
          ),
          "--resource:dbsnp,known=true,training=false,truth=false,prior=7",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
          ),
          "-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP",
          "-O",
          snps_recal_file,
          "--tranches-file",
          tranches_file
        )
      system(command = command, intern = T)
      
    } else{
      message("Ya se ha hecho el pre recalibrado de variantes")
    }
    
    
  }

applyVQSR <- function(folder_fasta, fastq_dir, output_dir) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## obtenemos el nombre del archivo
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  in_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz")
    )
  
  snps_recal_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.recal")
    )
  
  tranches_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.g.tranches")
    )
  
  out_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.vqsr.vcf")
    )
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk ApplyVQSR -R",
        fasta_file,
        "-V",
        in_file,
        "-O",
        out_file,
        "--truth-sensitivity-filter-level 99.0 --tranches-file",
        tranches_file,
        "--recal-file",
        snps_recal_file,
        "-mode SNP"
      )
    system(command = command, intern = T)
    
  } else{
    message("Ya se ha aplicado VQSR")
  }
  
}


variantFiltration <- function(folder_fasta, output_dir, fastq_dir) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## obtenemos el nombre del archivo
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  in_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.vqsr.vcf")
    )
  
  out_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(
        output_file_name,
        "_apply_genotypeGVCF.vqsr.varfilter.vcf"
      )
    )
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk VariantFiltration -R",
        fasta_file,
        "-V",
        in_file,
        "-O",
        out_file,
        "--filter-name LOW_depth10  --filter-expression 'DP< 10'"
      )
    system(command = command, intern = T)
    
    
  } else{
    "Ya se ha filtrado"
  }
  
  
  
}


analysisReady <- function(folder_fasta, output_dir, fastq_dir) {
  ## obtenemos el nombre del archivo
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  in_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(
        output_file_name,
        "_apply_genotypeGVCF.vqsr.varfilter.vcf"
      )
    )
  out_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(
        output_file_name,
        "_apply_genotypeGVCF.vqsr.varfilter.pass.vcf"
      )
    )
  
  if (!file.exists(out_file)) {
    command <-
      paste("bcftools view -f 'PASS,.' -O v -o", out_file, in_file)
    system(command = command, intern = T)
    
  } else{
    "ya se ha hecho el PASS filter"
  }
  
  
}



snpeff <- function(output_dir, folder_fasta) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  in_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.vqsr.vcf")
    )
  out_dir <- file.path(output_dir, "anotacion")
  out_file_1 <-
    file.path(out_dir, paste0(output_file_name, "anotacion_1.vcf"))
  out_file_2 <-
    file.path(out_dir, paste0(output_file_name, "anotacion_2.vcf"))
  out_file_3 <-
    file.path(out_dir, paste0(output_file_name, "anotacion_3.vcf"))
  
  campos_freq <-
    "~/tools/exomas_tools/snpEff/data/dbNSFP/hg38/dbNSFP4.1a.txt.gz"
  
  
  campos_clinvar <-
    "~/tools/exomas_tools/snpEff/data/clinvar/GRCh38/clinvar.vcf.gz"
  snpeff <- "~/tools/exomas_tools/snpEff/snpEff.jar"
  snpsift <- "~/tools/exomas_tools/snpEff/SnpSift.jar"
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  if (!file.exists(out_file_1)) {
    command <-
      paste("java -Xmx8g -jar",
            snpeff,
            "hg38",
            "-v",
            in_file,
            ">",
            out_file_1)
    system(command = command, intern = T)
  } else{
    message("Ya se a anotado 1")
  }
  if (file.exists(out_file_1) & !file.exists(out_file_2)) {
    command <- paste(
      "java -Xmx8g -jar",
      snpsift,
      "DbNsfp -db",
      campos_freq,
      "-f Uniprot_acc,Interpro_domain,1000Gp3_AC,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,SIFT_pred,SIFT_score,Polyphen2_HDIV_pred,Polyphen2_HDIV_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,MutationTaster_pred,MutationTaster_score,phastCons100way_vertebrate,phastCons100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP30way_mammalian_rankscore,phyloP17way_primate_rankscore,integrated_fitCons_score",
      out_file_1,
      "-v >",
      out_file_2
    )
    system(command = command, intern = T)
  } else{
    message("Ya se ha anotado nivel 1 y 2")
  }
  if (file.exists(out_file_1) &
      file.exists(out_file_2) & !file.exists(out_file_3)) {
    command <- paste(
      "java -Xmx8g -jar",
      snpsift,
      "annotate",
      campos_clinvar,
      out_file_2,
      "-v >",
      out_file_3
    )
    system(command = command, intern = T)
  } else{
    message("ya se ha anotado todo")
  }
  
  
  
}


export_tsv <- function(folder_fasta, output_dir, fastq_dir) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  
  in_file <-
    file.path(output_dir,
              "anotacion",
              paste0(output_file_name, "anotacion_3.vcf"))
  
  out_file <-
    file.path(output_dir,
              "anotacion",
              paste0(output_file_name, "anotacion_3.tsv"))
  if (!file.exists(out_file)) {
    command <-
      paste(
        "~/tools/exomas_tools/gatk-4.3.0.0/gatk VariantsToTable -R",
        fasta_file,
        "-V",
        in_file,
        "--show-filtered",
        "-O",
        out_file
      )
    system(command = command, intern = T)
    
  } else{
    message("Ya existe el archivo TSV")
  }
  
  
  
}

annovar <- function(output_dir, fastq_dir, folder_fasta) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  in_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(
        output_file_name,
        "_apply_genotypeGVCF.vqsr.varfilter.pass.vcf"
      )
    )
  out_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(
        output_file_name,
        "_apply_genotypeGVCF.vqsr.varfilter.pass_annotated"
      )
    )
  
  if (!file.exists(paste0(out_file, ".hg38_multianno.txt"))) {
    command <-
      paste(
        "~/tools/exomas_tools/annovar/table_annovar.pl",
        in_file,
        " ~/tools/exomas_tools/annovar/humandb/ -buildver hg38 -out myanno -remove -protocol refGene,clinvar_20220320,dbscsnv11,dbscsnv11,exac03,esp6500siv2_aa,gnomad211_exome,dbnsfp30a,avsnp150   -operation  g,f,f,f,f,f,f,f,f -nastring . -vcfinput -outfile",
        out_file
      )
    print(command)
    system(command)
    
    
  } else{
    print("YA ESTA ANOTADO")
  }
  
}

computo_cobertura <- function(output_dir, fastq_dir) {
  directorio_salida <- file.path(output_dir, "cobertura_y_stats")
  if (!dir.exists(directorio_salida)) {
    dir.create(directorio_salida)
  }
  fastq_files <- list.files(fastq_dir, full.names = F)
  
  bam_file <-
    list.files(
      file.path(output_dir, "mapping_output"),
      full.names = T,
      pattern = "sorted.mark_dup_RG_bqsr.bam"
    )
  print(bam_file)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "coverage", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  if (!file.exists(file.path(directorio_salida, paste0(output_file_name, ".txt")))) {
    command <- paste(
      "./compute_depth.sh",
      bam_file,
      ">",
      file.path(directorio_salida, paste0(output_file_name, ".txt"))
      ,
      ".txt"
    )
    
    system(command = command)
    
  } else{
    print("YA SE HA CALCULADO LA COBERTURA")
  }
  
  
}


computar_graficos <- function(output_dir, fastq_dir) {
  direcectorio_salida <- file.path(output_dir, "cobertura_y_stats")
  fastq_files <- list.files(fastq_dir, full.names = F)
  
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "coverage", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  in_file <- file_path_sans_ext(output_file_name)
  in_file <-
    file.path(output_dir, "cobertura_y_stats", paste0(in_file, ".txt"))
  
  archivo_salida <- file.path(
    output_dir,
    "cobertura_y_stats",
    paste0(output_file_name, "_graficos_cobertura.jpeg")
  )
  if (!file.exists(archivo_salida)) {
    cobertura <- read.delim(in_file, header = F)
    
    gcov <- cobertura[cobertura[, 1] == 'all',]
    ###
    longitud <- 300
    datos.pre <-
      data.frame(
        X = gcov[1:longitud, 2],
        Y = 100 * (1 - cumsum(gcov[1:longitud, 5])),
        Z = 100 * gcov[1:longitud, 5],
        relleno = gcov[1:longitud, 1]
      )
    datos.pre$relleno <- as.factor(datos.pre$relleno)
    p1 <-
      ggplot(data = datos.pre, aes(X, Y, fill = relleno)) + geom_line(color =
                                                                        "steelblue", linewidth = 2) + xlab("Profundidad de Cobertura") + ylab("Porcentaje de la región >= Profunidad") + theme(legend.position = "none") +
      theme_classic()
    p2 <-
      ggplot(data = datos.pre, aes(X, Z, fill = relleno)) + geom_col() + scale_fill_discrete(type =
                                                                                               "steelblue") + xlab("Profundidad de Cobertura") + ylab("Porcentaje de la región") +
      theme_article() + theme(legend.position = "none") + theme_classic()
    
    p3 <- ggarrange(p1, p2, ncol = 2)
    
    p4 <-
      annotate_figure(p3,
                      top = paste(
                        output_file_name,
                        "Profundidad media de cobertura:",
                        round(mean(cobertura[cobertura$V1 !=
                                               "all", "V7"]), 2),
                        "X"
                      ))
    
    
    
    ggsave(
      filename = file.path(
        output_dir,
        "cobertura_y_stats",
        paste0(output_file_name, "_graficos_cobertura.jpeg")
      ),
      plot = p4,
      bg = "white"
    )
  } else{
    print("YA SE REALIZARON LOS GRAFICOS")
  }
}

buscar_herencia <- function(df) {
  vector_hpo <- df$hpo  # Vector de la columna "hpo"
  vector_valores <-
    strsplit(vector_hpo, ";")  # Separar los valores por el delimitador ";"
  vector_resultado1 <-
    sapply(vector_valores, function(x)
      ifelse(grepl("dominant", x, ignore.case = T), "dominant", NA))
  vector_resultado1 <-
    unlist(lapply(vector_resultado1, function(x)
      ifelse(all(is.na(
        x
      )), NA, x[which(!is.na(x))])))
  vector_resultado2 <-
    sapply(vector_valores, function(x)
      ifelse(grepl("recessive", x, ignore.case = T), "recessive", NA))
  vector_resultado2 <-
    unlist(lapply(vector_resultado2, function(x)
      ifelse(all(is.na(
        x
      )), NA, x[which(!is.na(x))])))
  
  # data frame resultado
  df <- as.data.table(df)
  df <-
    bind_cols(dominante = vector_resultado1, recesivo = vector_resultado2, df)
  return(df)
}
pretratado <- function(output_dir, fastq_dir) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  in_file <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  in_file_tmp <- file_path_sans_ext(in_file)
  
  in_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(
        in_file_tmp,
        "_apply_genotypeGVCF.vqsr.varfilter.pass_annotated.hg38_multianno.txt"
      )
    )
  codigo <- unlist(strsplit(in_file_tmp, "_"))[1]
  
  directorio_salida <- file.path(output_dir, "post_process_results")
  directorio_salida <-
    file.path(directorio_salida,
              paste(codigo, "GATK_CODIGO_HPO", sep = "_"))
  out_file <-
    file.path(directorio_salida,
              paste(codigo, "GATK_CODIGO_HPO.csv", sep = "_"))
  if (!file.exists(out_file)) {
    hpo <-
      read.delim("./genes_to_phenotype.txt",
                 skip = 1,
                 header = F)[, c(2, 4)]
    X <-
      read.delim(
        in_file,
        header = T,
        na.strings = ".",
        row.names = NULL
      )
    
    QUAL <- X$Otherinfo2
    DP <- X$Otherinfo3
    
    
    
    m <-
      (apply(strsplit2(X[, "Otherinfo11"], ";"), 2, function(x)
        strsplit2(x, "=")[2]))
    m <- m[!m == ""] ## variables
    X.tmp <-
      strsplit2(X[, "Otherinfo11"], split = ";") # matriz spliteada
    variable <-
      matrix(NA, nrow = nrow(X), ncol = length(m)) ## matriz resultados
    colnames(variable) <- m
    indices <- vector("numeric", length = nrow(X))
    
    matriz_aux <-  matrix(NA, nrow = nrow(X), ncol = length(m))
    for (j in 1:length(m)) {
      m_i <- m[j]
      for (i in 1:nrow(X.tmp)) {
        pos_i <- grep(paste0("^", m_i, "="), X.tmp[i,])
        if (length(pos_i) == 0) {
          matriz_aux[i, j]  <- NA
        } else{
          matriz_aux[i, j] <- strsplit2(X.tmp[i, pos_i], "=")[2]
          
        }
      }
      
    }
    
    colnames(matriz_aux) <- m
    
    matriz_12 <- strsplit2(X$Otherinfo12, ":")
    matriz_13 <- strsplit2(X$Otherinfo13, ":")
    
    other_info <- matriz_13[, c(1, 2, 4)]
    colnames(other_info) <- unique(matriz_12[, c(1, 2, 4)])
    other_info <- as.data.table(other_info)
    midf <- X[,-grep("otherinfo", ignore.case = T, colnames(X))]
    w <- which(colnames(matriz_aux) %in% colnames(midf))
    
    colnames(midf)[grep("^AF$", colnames(midf))] <- "AF_gnomad"
    midf <-
      bind_cols(midf[, 1:5], as.data.table(matriz_aux), other_info, midf[, 6:ncol(midf)])
    
    campos_freq <- c(
      "ExAC_ALL",
      "ExAC_AFR",
      "ExAC_AMR",
      "ExAC_EAS",
      "ExAC_FIN",
      "ExAC_NFE",
      "ExAC_OTH",
      "ExAC_SAS",
      "esp6500siv2_aa",
      "AF_popmax",
      "AF_male",
      "AF_female",
      "AF_raw",
      "AF_afr",
      "AF_sas",
      "AF_amr",
      "AF_eas",
      "AF_nfe",
      "AF_fin",
      "AF_asj",
      "AF_oth",
      "non_topmed_AF_popmax",
      "non_neuro_AF_popmax",
      "non_cancer_AF_popmax",
      "controls_AF_popmax"
    )
    midf[, campos_freq] <- apply(midf[, campos_freq], 2, as.numeric)
    midf[, campos_freq] <-
      apply(midf[, campos_freq], 2, function(x)
        ifelse(is.na(x), 0, x))
    midf$freq <- rowMeans(midf[, campos_freq], na.rm = T)
    midf <- midf[,-which(colnames(midf) %in% campos_freq)]
    if (any(midf$GT == "0/0" | midf$GT == "0|0")) {
      stop(errorCondition(message = "Hay homocigotos regerefetes"))
    }
    cigosidad <-
      ifelse(midf$GT == "1/1" | midf$GT == "1|1", "HOMZ_ALT", "HETZ")
    
    midf <- as.data.table(midf)
    midf <- bind_cols(codigo = codigo, cigosidad = cigosidad, midf)
    
    
    hpo_ <- aggregate(V4 ~ V2, hpo, FUN = paste, collapse = ";")
    colnames(hpo_) <- c("Gene.refGene", "hpo")
    midf <- left_join(midf, hpo_, by = "Gene.refGene")

    if (!dir.exists(directorio_salida)) {
      dir.create(directorio_salida, recursive = T)
    }
    out_file <-
      file.path(directorio_salida,
                paste(codigo, "GATK_CODIGO_HPO.csv", sep = "_"))
    midf <- as.data.frame(buscar_herencia(midf))
    write.csv(midf, file = out_file, row.names = F)
  } else{
    print("YA SE HA PRETRATADO")
  }
  
}


filtrado_1 <- function(X, threshold) {
  X <-
    X[which(X$Func.refGene == "exonic;splicing" |
              X$Func.refGene == "exonic"),]
  X <-
    X[which(
      X$ExonicFunc.refGene == "nonsynonymous SNV" |
        X$ExonicFunc.refGene == "exonic" |
        X$ExonicFunc.refGene == "nonsynonymous SNV" |
        X$ExonicFunc.refGene == "frameshift deletion" |
        X$ExonicFunc.refGene == "frameshift insertion" |
        X$ExonicFunc.refGene == "startloss" |
        X$ExonicFunc.refGene == "startgain" |
        X$ExonicFunc.refGene == "stoploss"
    ),]
  
  X <- X[which(X[, grep("freq", colnames(X))] < threshold),]
  
  return(X)
  
}


filtrado_2 <- function(X) {
  X$CADD_phred <- as.numeric(X$CADD_phred)
  X <- X[which(X$CADD_phred >= 20),]
  X <-
    X[which((!is.na(X$CLNSIG) &
               !is.na(X$CLNDN) & !is.na(X$avsnp150)) == T),]
  
  return(X)
}

filtrado_3 <- function(X) {
  X <- X[which(grepl("benign", X$CLNSIG, ignore.case = T) == F),]
  return(X)
}


computo_frecuencias <- function(output_dir, fastq_dir) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  in_file <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  in_file_tmp <- file_path_sans_ext(in_file)
  codigo <- unlist(strsplit(in_file_tmp, "_"))[1]
  
  directorio <-
    file.path(output_dir,
              "post_process_results",
              paste0(codigo, "_GATK_CODIGO_HPO/"))
  
  archivo <- list.files(directorio, full.names = T,pattern = ".csv")
  
  
  X <- read.csv(archivo, header = T, row.names = NULL)
  
  
  todos <- readRDS("./todos.rds")
  
  
  
  if (!any(!names(todos) == codigo)) {
    todos[[codigo]] <- X
    
    
    todos <- lapply(todos, as.data.frame)
    common_cols <- Reduce(intersect, lapply(todos, colnames))
    
    X_new <- lapply(todos, function(df)
      df[, common_cols])
    
    
    
    df <- bind_rows(X_new) %>%
      select(codigo, cigosidad, AF, Chr, Start, End, Gene.refGene) %>% group_by(Chr, Start, End, Gene.refGene) %>%
      summarise(AF = (sum(cigosidad == "HETZ") + 2 * sum(cigosidad == "HOMZ_ALT")) /
                  (2 * n()))
    
    cromosomas <- c(paste0("chr", 1:22), "X", "Y")
    
    
    df <- bind_rows(X_new) %>%
      select(codigo, cigosidad, AF, Chr, Start, End, Gene.refGene) %>% group_by(Chr, Start, End, Gene.refGene) %>%
      summarise(
        paste_m = toString(codigo),
        n = n(),
        AF = (sum(cigosidad == "HETZ") + 2 * sum(cigosidad == "HOMZ_ALT")) / (2 *
                                                                                n())
      )
    
    df.final <- df[df$Chr %in% cromosomas,]
    
    
    write.csv(df.final, "./frecuencias_alelicas_lab.csv")
  } else{
    print("YA SE TIENE EN CUENTA LA FRECUENCIA ALELICA NUEVA")
  }
}

filtrado_general <- function(output_dir, fastq_dir) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  in_file <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  in_file_tmp <- file_path_sans_ext(in_file)
  codigo <- unlist(strsplit(in_file_tmp, "_"))[1]
  
  directorio_salida <- file.path(output_dir, "post_process_results")
  directorio_salida <-
    file.path(directorio_salida,
              paste(codigo, "GATK_CODIGO_HPO", sep = "_"))
  exoma.archivo <-
    file.path(directorio_salida,
              paste(codigo, "GATK_CODIGO_HPO.csv", sep = "_"))
  
  exoma <- read.csv(exoma.archivo, header = T, row.names = NULL)
  
  
  threshold <- 1
  
  retorno <- filtrado_1(exoma, threshold)
  retorno2 <- filtrado_2(retorno)
  retorno3 <- filtrado_3(retorno2)
  
  directorio_salida.crudo <- file.path(directorio_salida,"filtrado_sin_comparar")
  if(!dir.exists(directorio_salida.crudo)){
    dir.create(directorio_salida.crudo)
  }
  
  write.csv(retorno, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_1_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno2, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_2_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno3, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_3_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  
  threshold <- 0.02
  
  retorno <- filtrado_1(exoma, threshold)
  retorno2 <- filtrado_2(retorno)
  retorno3 <- filtrado_3(retorno2)
  
  write.csv(retorno, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_1_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno2, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_2_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno3, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_3_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  
  ### ahroa vamos con la comparativa
  
  freq_alelicas <- read.csv("./frecuencias_alelicas_lab.csv")
  freq_alelicas <-
    freq_alelicas[!duplicated(freq_alelicas[, c("Start", "End", "paste_m")]),]
  
  unicas <-
    freq_alelicas[grep(paste0("^", codigo, "$"), freq_alelicas$paste_m),]
  print(colnames(freq_alelicas))
  exoma.comparado <-
    left_join(
      exoma,
      unicas,
      by = c("Start", "End", "Chr", "Gene.refGene")
      ,
      relationship = "many-to-many"
    )
  
  threshold <- 1
  
  retorno <- filtrado_1(exoma.comparado, threshold)
  retorno2 <- filtrado_2(retorno)
  retorno3 <- filtrado_3(retorno2)
  directorio_salida.comparando <- file.path(directorio_salida,"comparando_bd_lab")
  if(!dir.exists(directorio_salida.comparando)){
    dir.create(directorio_salida.comparando)
  }
  write.csv(retorno, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_1_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno2, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_2_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno3, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_3_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  
  threshold <- 0.02
  
  retorno <- filtrado_1(exoma.comparado, threshold)
  retorno2 <- filtrado_2(retorno)
  retorno3 <- filtrado_3(retorno2)
  
  write.csv(retorno, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_1_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno2, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_2_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno3, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_3_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  
  
  
}


filtrado_vias <- function(output_dir,fastq_dir){
  fastq_files <- list.files(fastq_dir, full.names = F)
  in_file <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  in_file_tmp <- file_path_sans_ext(in_file)
  codigo <- unlist(strsplit(in_file_tmp, "_"))[1]
  
  directorio_salida <- file.path(output_dir, "post_process_results")
  directorio_salida <-
    file.path(directorio_salida,
              paste(codigo, "GATK_CODIGO_HPO", sep = "_"))
  
  exoma.archivo <-
    file.path(directorio_salida,
              paste(codigo, "GATK_CODIGO_HPO.csv", sep = "_"))
  
  exoma <- read.csv(exoma.archivo)
  freqs <- read.csv("./frecuencias_alelicas_lab.csv")
  
  genes_univ <- unique(freqs$Gene.refGene)
  exoma_goi <- unique(exoma$Gene.refGene)
  
  goi <- clusterProfiler::bitr(exoma_goi,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  universe <- clusterProfiler::bitr(genes_univ,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  
  
  
  ewp.up <- clusterProfiler::enrichWP(
    goi[,2],
    universe = universe[,2],
    organism = "Homo sapiens",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05, #p.adjust cutoff; relaxed for demo purposes
  )
  
  ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keyType = "ENTREZID")
  
  resultados.1 <- ewp.up@result
  resultados.1 <- resultados.1[order(resultados.1$pvalue),]
  
  p1 <- ggplot(resultados.1[1:20,], aes(x=Description, y=Count, fill=pvalue)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_continuous(low="blue", high="red") +
    labs(x = "", y = "", fill = "p.value") +
    theme(axis.text=element_text(size=11))  
  
  ggsave(filename=file.path(directorio_salida,"vias.jpeg"),plot=p1)
  
  ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keyType = "ENTREZID")
  resultados.1 <- ewp.up@result
  resultados.1 <- resultados.1[order(resultados.1$pvalue),]
  
  goi2 <- unlist(strsplit(resultados.1$geneID,"/"))
  goi2 <- clusterProfiler::bitr(goi2,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  
  dose <- DOSE::enrichDO(goi2[,2])
  dose <- DOSE::setReadable(dose, org.Hs.eg.db, keyType = "ENTREZID")
  
  resultado3 <- dose@result[order(dose@result$p.adjust),]
  p2 <- ggplot(resultado3[1:20,], aes(x=Description, y=Count, fill=-log(p.adjust))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_continuous(low="blue", high="red") +
    labs(x = "", y = "", fill = "p.adjust") +
    theme(axis.text=element_text(size=11))  
  ggsave(filename=file.path(directorio_salida,"patologias.jpeg"),plot=p2)
  
  
  pathos.list <- list()
  pathos <- resultado3
  for(l in 1:dim(pathos)[1]){
    repeticion <- length(unlist(strsplit(pathos$geneID[l],"/")))
    pathos.list[[l]] <- data.frame(patologia =rep(pathos$Description[l],repeticion),
                                   Gene.refGene = unlist(strsplit(pathos$geneID[l],"/")),
                                   p_value = rep(pathos$pvalue[l],repeticion),
                                   p_adjust = rep(pathos$p.adjust[l],repeticion))
    
    
  }
  
  pathos.df <- (Reduce(rbind,pathos.list))
  
  final <- left_join(exoma,pathos.df,by="Gene.refGene")
  
  
  threshold <- 1
  
  retorno <- filtrado_1(final, threshold)
  retorno2 <- filtrado_2(retorno)
  retorno3 <- filtrado_3(retorno2)
  
  directorio_salida.crudo <- file.path(directorio_salida,"filtrado_sin_comparar_vias")
  if(!dir.exists(directorio_salida.crudo)){
    dir.create(directorio_salida.crudo)
  }
  
  write.csv(retorno, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_1_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno2, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_2_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno3, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_3_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  
  threshold <- 0.02
  
  retorno <- filtrado_1(final, threshold)
  retorno2 <- filtrado_2(retorno)
  retorno3 <- filtrado_3(retorno2)
  
  write.csv(retorno, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_1_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno2, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_2_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno3, file = file.path(directorio_salida.crudo, paste(
    codigo,
    paste0("filtrado_3_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  
  ### ahroa vamos con la comparativa
  
  freq_alelicas <- read.csv("./frecuencias_alelicas_lab.csv")
  freq_alelicas <-
    freq_alelicas[!duplicated(freq_alelicas[, c("Start", "End", "paste_m")]),]
  
  unicas <-
    freq_alelicas[grep(paste0("^", codigo, "$"), freq_alelicas$paste_m),]
  print(colnames(freq_alelicas))
  exoma.comparado <-
    left_join(
      final,
      unicas,
      by = c("Start", "End", "Chr", "Gene.refGene")
      ,
      relationship = "many-to-many"
    )
  
  threshold <- 1
  
  retorno <- filtrado_1(exoma.comparado, threshold)
  retorno2 <- filtrado_2(retorno)
  retorno3 <- filtrado_3(retorno2)
  directorio_salida.comparando <- file.path(directorio_salida,"comparando_bd_lab_vias")
  if(!dir.exists(directorio_salida.comparando)){
    dir.create(directorio_salida.comparando)
  }
  write.csv(retorno, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_1_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno2, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_2_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno3, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_3_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  
  threshold <- 0.02
  
  retorno <- filtrado_1(exoma.comparado, threshold)
  retorno2 <- filtrado_2(retorno)
  retorno3 <- filtrado_3(retorno2)
  
  write.csv(retorno, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_1_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno2, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_2_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  write.csv(retorno3, file = file.path(directorio_salida.comparando, paste(
    codigo,
    paste0("filtrado_3_comparado_thresh_", threshold),
    ".csv",
    sep = "_"
  )))
  
  
  
  
}

muestras <- c("DX010-23")
for (muestra in muestras) {
  pipeline_dir <- "/repositorio/exomas/pipeline"
  fastq_dir <- file.path(pipeline_dir, muestra, "fastqfiles")
  output_dir <- file.path(pipeline_dir, muestra, "output_dirs")
  folder_fasta <- "/repositorio/exomas/datos/datos_gatk/referencia"
  folder_data_gatk <- "/repositorio/exomas/datos/datos_gatk"
  
  
  
  
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  
  
  
  
  control_calidad(fastq_dir, output_dir)
  ### Solo se corre una vez por el genoma
  index_fasta_samtools(folder_fasta = folder_fasta)
  ### Indice genoma corre una vez
  index_bwa(folder_fasta = folder_fasta)
  ### Mapeamos
  bwamem(fastq_dir = fastq_dir, folder_fasta = folder_fasta)
  ### marcamos duplicados
  markdups(output_dir = output_dir, fastq_dir = fastq_dir)
  ## creamos diccionario
  create_dict(folder_fasta)
  ## anadimos reaad group
  creacion_readgroup(output_dir, fastq_dir)
  ## Recalibramos
  base_recalibrator(folder_fasta, output_dir, folder_data_gatk, fastq_dir)
  ### aplicamos el recalibrado
  applybqsr(folder_fasta, output_dir, fastq_dir)
  ## estadisticas del pieline bam
  bam_statistics(folder_fasta, fastq_dir, output_dir)
  ## llamamos a las variantes
  haplotype_caller(output_dir, folder_fasta, fastq_dir)
  ## Calculamos la probabilidad posterior del alelo referente
  genotypeGVCF(folder_fasta, output_dir, fastq_dir)
  ## calculamos variant Recalibrator
  variantRecallibrator(fastq_dir, folder_fasta, folder_data_gatk, output_dir)
  ## apply VQSR
  applyVQSR(folder_fasta, fastq_dir, output_dir)
  ## primer filtraje
  variantFiltration(folder_fasta, output_dir, fastq_dir)
  ## preparamos el archivo listo para elanalisis
  analysisReady(folder_fasta, output_dir, fastq_dir)
  ##
  annovar(output_dir, fastq_dir, folder_fasta)
  
  computo_cobertura(output_dir, fastq_dir)
  
  computar_graficos(output_dir, fastq_dir)
  
  pretratado(output_dir, fastq_dir)
  
  computo_frecuencias(output_dir, fastq_dir)
  
  filtrado_general(output_dir, fastq_dir)
  
  filtrado_vias(output_dir,fastq_dir )
  
  print(paste("YA TERMINO LA MUESTRA ", muestra))
  
}
# ## anotamos
## snpeff(output_dir,folder_fasta)
# ## exportamos TSV
# export_tsv(folder_fasta,output_dir ,fastq_dir)
