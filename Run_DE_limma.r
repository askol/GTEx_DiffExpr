job.files = c()

ScriptDir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Scripts/"

if (!dir.exists(ScriptDir)){
    dir.create(ScriptDir)
}
setwd(ScriptDir)

tissues <- read.table(file = paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                          "data/support_files/all_v8_tissues_both_sexes.txt"), as.is=T)
tissues <- unlist(tissues)

for (tissue in tissues){

    job.file = paste0(tissue,".pbs")
    job.files = c(job.files, job.file)

    write(file = job.file, "#!/bin/bash")
    write(file = job.file, paste("#PBS -l nodes=1:ppn=1,mem=8gb"), append=TRUE)
    write(file=job.file, "#PBS -l walltime=36:00:00", append=TRUE)
    write(file = job.file,
          paste0("#PBS -o ", ScriptDir,
                 tissue,"_pbs.out"), append=TRUE)
    write(file = job.file, "#PBS -j oe", append=TRUE)
    write(file = job.file, paste0("#PBS -N Crt_dat_",tissue), append=TRUE)
    write(file = job.file, paste0("R CMD BATCH  '--args ",tissue,
              "' /gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Code/DE_limma.r ",
              ScriptDir, tissue,".out"),
              append=TRUE)    
    system(paste("chmod +x", job.file))
}

for (job in job.files){
    system(paste("qsub",job))
}
