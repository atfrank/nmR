# (0) set working directory
setwd("~/GitSoftware/nmR/")

# (1) load library
library(nmR)

# (2) load chemical shift data
rnas <- unlist(strsplit("1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1SCL 1UUU 1XHP 1YSV 1Z2J 1ZC5 2FDT 2JWV 2K66 2KOC 2L3E 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LQZ 2LU0 2LUB 2LUN 2LV0 2M4W 2M5U 2M8K 2M12 2M21 2M22 2M24 2MEQ 2MFD 2MHI 2MIS 2MNC 2MXL 2N2O 2N2P 2N4L 2N6S 2N6T 2N6W 2N6X 2N7X 2NBY 2NBZ 2NC0 2NCI 2QH2 2QH4 2Y95 4A4S 4A4T 4A4U 5A17 5A18 5KQE", " "))
for (rna in rnas){
  try({
    predcs <- nmR::load_cs_data(csfile=paste("data/larmord_", rna, ".txt", sep = ""), accuracyFile = "data/larmord_accuracy_resname_nucleus.txt", atomBasedWeights = TRUE, names = c("conformation", "resid", "resname", "nucleus", "predCS", "ID"))
    expcs <-  nmR::load_cs_data(csfile=paste("data/observed_shifts_corrected_larmord_", rna, ".txt", sep = ""), names = c("resname", "resid", "nucleus", "expCS", "error"))
    
    # (3) merge chemical shift data
    cs <- merge(expcs, predcs, by = c("resname", "resid", "nucleus"))
    cs <- cs[order(cs$resid, cs$resname, cs$nucleus, cs$conformation),] # the be careful here; inspect the files written out below to ensure they correspond with your expectations
    
    # (4) write out vector, corresponding to observed shifts and matrix, corresponding to predict shifts
    tmp <- subset(cs, conformation==1)
    observed <- tmp$expCS
    weights <-  tmp$weight
    predicted <- matrix(cs$predCS, nrow = length(observed), byrow = TRUE)
    write.table(observed, file = paste("data/observed_vector_", rna, ".txt", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(weights, file = paste("data/weights_vector_", rna, ".txt", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(predicted, file = paste("data/predicted_matrix_", rna, ".txt", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE)
  })
}
