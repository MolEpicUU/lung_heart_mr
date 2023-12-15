# https://github.com/MRCIEU/gwasvcf/blob/master/R/proxy.r
# https://github.com/explodecomputer/genetics.binaRies

library(genetics.binaRies)
library(magrittr)
library(rlang)

get_ld_proxies <- function(rsid, bfile, searchspace=NULL, tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1, out=tempfile())
{
  stopifnot(check_plink())
  searchspacename <- paste0(out, ".searchspace")
  targetsname <- paste0(out, ".targets")
  outname <- paste0(out, ".targets.ld.gz")
  utils::write.table(rsid, file=targetsname, row=FALSE, col=FALSE, qu=FALSE)
  if(!is.null(searchspace))
  {
    stopifnot(is.character(searchspace))
    utils::write.table(c(rsid, searchspace), file=searchspacename, row=F, col=F, qu=F)
    extract_param <- paste0(" --extract ", searchspacename)
  } else {
    extract_param <- " "
  }
  cmd <- paste0(options()[["tools_plink"]],
                " --bfile ", bfile, 
                extract_param,
                " --keep-allele-order ",
                " --r in-phase with-freqs gz",
                " --ld-snp-list ", targetsname,
                " --ld-window-kb ", tag_kb,
                " --ld-window ", tag_nsnp,
                " --out ", targetsname,
                " --threads ", threads,
                " 2>&1 > /dev/null"
  )
  message("Finding proxies...")
  system(cmd)
  
  ld <- data.table::fread(paste0("gunzip -c ", outname), header=TRUE) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(.data[["R"]]^2 > tag_r2) %>%
    dplyr::filter(.data[["SNP_A"]] != .data[["SNP_B"]]) %>%
    dplyr::mutate(PHASE=gsub("/", "", .data[["PHASE"]])) %>%
    dplyr::filter(nchar(.data[["PHASE"]]) == 4)
  unlink(searchspacename)
  unlink(targetsname)
  unlink(paste0(targetsname, c(".log", ".nosex")))
  unlink(outname)
  if(nrow(ld) == 0)
  {
    message("No proxies found")
    return(ld)
  }
  temp <- do.call(rbind, strsplit(ld[["PHASE"]], "")) %>% dplyr::as_tibble(.data, .name_repair="minimal")
  names(temp) <- c("A1", "B1", "A2", "B2")
  ld <- cbind(ld, temp) %>% dplyr::as_tibble(.data, .name_repair="minimal")
  # ld <- dplyr::arrange(ld, desc(abs(R))) %>%
  # 	dplyr::filter(!duplicated(SNP_A))
  ld <- dplyr::arrange(ld, dplyr::desc(abs(.data[["R"]])))
  message("Found ", nrow(ld), " proxies")
  return(ld)
}

check_plink <- function()
{
  if(is.null(options()[["tools_plink"]]))
  {
    message("'tools_plink' option is not set. See 'set_plink' for information.")
    return(FALSE)
  }
  filecheck <- file.exists(options()[["tools_plink"]])
  if(filecheck)
  {
    return(TRUE)
  }
  pathcheck <- any(sapply(strsplit(Sys.getenv("PATH"), split=":"), function(x) file.exists(file.path(x, options()[["tools_plink"]]))))
  if(pathcheck)
  {
    return(TRUE)
  }
  message("'tools_plink' option is not set. See 'set_plink' for information.")
  return(FALSE)
}

set_plink <- function(path="")
{
  if(is.null(path))
  {
    options(tools_plink = NULL)
  } else if(path == "")
  {
    a <- requireNamespace("genetics.binaRies")
    if(a)
    {
      message("Path not provided, using binaries in the explodecomputer/genetics.binaRies package")
      options(tools_plink = genetics.binaRies::get_plink_binary())
    } else {
      stop("Please provide a path to plink binary or run devtools::install_github('explodecomputer/genetics.binaRies')")
    }
  } else {
    options(tools_plink = path)
  }
}

set_bcftools <- function(path="")
{
  if(is.null(path))
  {
    options(tools_bcftools = NULL)
  } else if(path == "")
  {
    a <- requireNamespace("genetics.binaRies")
    if(a)
    {
      message("Path not provided, using binaries in the explodecomputer/genetics.binaRies package")
      options(tools_bcftools = genetics.binaRies::get_bcftools_binary())
    } else {
      stop("Please provide a path to bcftools binary or run devtools::install_github('explodecomputer/genetics.binaRies')")
    }
  } else {
    options(tools_bcftools = path)
  }
}

set_plink()