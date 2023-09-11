repo=BiocManager::repositories()
i = which(names(repo) == "CRAN")
pkgs.cran <- tools::package_dependencies('clusterProfiler', 
                            db=utils::available.packages(repo=repo[i]), 
                    which = c("Depends", "Imports"), reverse=TRUE)[[1]]

pkgs.bioc <- tools::package_dependencies('clusterProfiler', 
                            db=utils::available.packages(repo=repo[-i]), 
                            which = c("Depends", "Imports"), 
                            reverse=TRUE)[[1]]           

d1 <- data.frame(Package = pkgs.cran,
                Description = sapply(pkgs.cran, yulab.utils::packageTitle),
                repo = "CRAN"
            )
d2 <- data.frame(Package = pkgs.bioc,
                Description = sapply(pkgs.bioc, yulab.utils::packageTitle, repo="BioC"), 
                repo = "Bioconductor"
            )  
d <- rbind(d1, d2)   


gt::gt(d)

