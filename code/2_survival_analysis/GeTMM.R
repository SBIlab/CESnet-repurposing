##---------------------------------------------------
## GeTMM Normalization of RSEM expected counts
## TCGA patients
##---------------------------------------------------

## Import library
library(edgeR)
library(tibble)

## Initialize
#tcga_list = c('TCGA-BLCA')
tcga_list = c('TCGA-BRCA','TCGA-GBM','TCGA-OV','TCGA-LUAD','TCGA-UCEC','TCGA-KIRC',
              'TCGA-HNSC','TCGA-LGG','TCGA-THCA','TCGA-LUSC','TCGA-PRAD','TCGA-SKCM',
              'TCGA-COAD','TCGA-STAD','TCGA-BLCA','TCGA-LIHC','TCGA-CESC','TCGA-KIRP',
              'TCGA-SARC','TCGA-LAML','TCGA-ESCA','TCGA-PAAD','TCGA-PCPG','TCGA-READ',
              'TCGA-TGCT','TCGA-THYM','TCGA-KICH','TCGA-ACC','TCGA-MESO','TCGA-UVM',
              'TCGA-DLBC','TCGA-UCS','TCGA-CHOL')
fi_dir = '../../data/TCGA_patient/'
setwd(fi_dir)

## Import max gene length in kb
gene_length = read.table('biomart_ensembl_max_exon_length_in_kilo_bp.txt', header=TRUE, sep='\t', check.names=FALSE)


## GeTMM normalization for each cancer types
for (cancer in tcga_list){
  # set output directory
  tmp_fo_dir = paste0(fi_dir, cancer)
  setwd(tmp_fo_dir)
  
  print(paste0('testing ', cancer))
  
  if (file.exists('rna_seq_table.csv')==TRUE){
    if (file.exists('GeTMM_rna_seq.txt')==FALSE){
      # Import expression data (HTseq count)
      count = read.csv('rna_seq_table.csv', header=TRUE, check.names=FALSE)
      rownames(count) = NULL
      colnames(count) = c(names(count))
      colnames(count)[1] = 'Gene stable ID'
      
      # Calculate RPK
      common_column_name = intersect(names(gene_length), names(count))
      new_count = merge(gene_length, count, by=common_column_name, all.dataframe_1=TRUE)
      rpk <- (new_count[,3:ncol(new_count)]/new_count[,2]) # rpk
      #new_count = new_count[,-1] # remove gene length column in 'new_count'
      
      # ---------------
      # Compute GeTMM
      #tmp <- new_count
      #tmp$`Gene stable ID` <- NULL
      
      # Create DGEList object
      dgList <- DGEList(counts=rpk, genes=new_count$`Gene stable ID`)
      #dgList$samples
      #head(dgList$counts)
      #head(dgList$genes)
      
      # Filtering
      countsPerMillion <- cpm(dgList)
      #summary(countsPerMillion)
      countCheck <- countsPerMillion > 1
      
      keep<- which(rowSums(countCheck) >= 2)
      dgList <- dgList[keep,]
      #summary(cpm(dgList))
      
      # Normalization (TMM)
      dgList <- calcNormFactors(dgList, method='TMM')
      cps <- cpm(dgList) #, normalized.lib.sizes = TRUE)
      output <- data.frame(dgList$genes, cps, check.names=FALSE)
      
      cps2 <- cpm(dgList, log=TRUE, prior.count=1) #, normalized.lib.sizes = TRUE)
      output2 <- data.frame(dgList$genes, cps2, check.names=FALSE)
      
      # Return normalized count results
      write.table(output, file='GeTMM_rna_seq.txt', sep='\t',quote = FALSE, row.names=FALSE)  
      write.table(output2, file='GeTMM_log2_rna_seq.txt', sep='\t',quote = FALSE, row.names=FALSE)  
      
    }
    
  }
  
  # set directory to parent folder
  setwd(fi_dir)
}
