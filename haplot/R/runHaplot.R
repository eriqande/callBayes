#' Run Haplotype Plot.
#'
#' Run Haplotype shiny plot
#' @export
#' @examples
#' runHaplotype()
runHaplotype <- function() {
  appDir <- system.file("shiny", "haPLOType", package = "haplot")
  if (appDir == "") {
    stop("Could not find shiny directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  runApp(appDir, display.mode = "normal")
}


#' Extract haplotype from alignment reads.
#'
#' The function \code{haplot} extracts haplotype from sequence alignment files through perl script 'hapture' and returns a summary table of the read depth and read quality associate with haplotype.
#'
#' @param run.label character vector. Run label to be used to display in haPLOType. Required
#' @param sam.path string. Directory path containings all sequence alignment files (SAM). Required
#' @param label.path string. Path to label file. This customized label file is a tab-separate file that contains entries of SAM file name, individual ID, and group label. Required
#' @param vcf.path string. Path to vcf file. Required
#' @param out.path string. Optional. If not specified, the intermediate files are created under \code{sam.path}, assuming that directory is granted for written permission.
#' @export
#' @examples
#' runHaplot("example 1", sam.path="data", label.path="data/label.txt", vcf.path="data/vcf.txt")
runHaplot <- function(run.label, sam.path, label.path, vcf.path,
  out.path=sam.path){


  run.label <- gsub(" +","_",run.label)

  haptureDir <- system.file("perl", "hapture", package = "haplot")


  # Need to check whether all path and files exist

  bash.cmd <- paste0("while IFS=$'\t' read -r -a line; do echo ",
    "\"perl ", haptureDir,
    " -v ", vcf.path, " ",
    " -s ", sam.path, "/${line[0]} ",
    " -i ${line[1]} ",
    " -g ${line[2]} ", " > ",
    out.path, "/", run.label, "_", "${line[1]}.summary &\"",
    "; done < ", label.path ,
    " > ", out.path, "/", "run_all.sh")


  system(bash.cmd)
  cat("Extracting Haplotyping ...")


   -v ${FOLDPATH}/chinook_noMNP_noComplex_noPriors.vcf -s ${FOLDPATH}/ch_flashed_${line[1]}_aln.sam  -i ${line[1]} -g '${line[0]}' > ${FOLDPATH}/chinook_panel1_${line[1]}.summary&";
  done < /home/biopipe/Genetics_Lab_Data/GTseq/chinook05172016/map/ch_microhaps_pops.txt > run_all.sh
  bash run_all.sh


  system("perl ${SRCPATH}/hap_recap.pl -v ${FOLDPATH}/chinook_noMNP_noComplex_noPriors.vcf -s ${FOLDPATH}/ch_flashed_${line[1]}_aln.sam  -i ${line[1]} -g '${line[0]}' > ${FOLDPATH}/chinook_panel1_${line[1]}.summary&")


  # paste0
  system("cat ...")


  summary.table.name<-"/home/biopipe/Genetics_Lab_Data/GTseq/chinook05172016/map/chinook_panel1.summary"
  run.label <- "chinook_05_17_2016"



  haplo.sum <- read.table(summary.table.name, stringsAsFactors = FALSE, sep="\t") %>%
    tbl_df

  colnames(haplo.sum) <- c("group", "id", "locus", "haplo", "depth", "logP.call", "logP.miscall")


  # clean the file: remove any haplotype with alignment character ("*" - unk, keep "_" - del)  grepl("[_|*]") &
  # remove any loci that does less than half of total individual
  # post remove any loci with highly variant haplotypes > 40

  num.id <- length(unique(haplo.sum$id))

  haplo.cleanup <- haplo.sum %>%
    #filter(grepl("[*]", haplo)) %>%
    group_by(locus, id) %>%
    mutate(n.haplo.per.indiv=n()) %>%
    ungroup() %>%
    group_by(locus) %>%
    mutate(n.indiv.per.locus = length(unique(id)), max.uniq.hapl=max(n.haplo.per.indiv)) %>%
    ungroup() %>%
    filter(n.indiv.per.locus > num.id/2, max.uniq.hapl < 40)  %>%
    select(group, id, locus, haplo, depth, logP.call, logP.miscall)


  haplo.add.balance <- haplo.cleanup %>%
    group_by(locus,id) %>%
    arrange(-depth) %>%
    mutate(allele.balance = depth/depth[1], rank=row_number() ) %>%
    ungroup()

  saveRDS(haplo.add.balance, paste0(run.label,".rds"))

  system("cp ...")


}
