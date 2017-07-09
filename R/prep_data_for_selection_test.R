# (0) set working directory
setwd("~/GitSoftware/nmR/")

# (1) load library
library(nmR)

# (2) load chemical shift data
predcs <- nmR::load_cs_data(csfile="data/larmord_1SCL.txt", names = c("conformation", "resid", "resname", "nucleus", "predCS", "ID"))
expcs <-  nmR::load_cs_data(csfile="data/observed_shifts_corrected_larmord_1SCL.txt", names = c("resname", "resid", "nucleus", "expCS", "error"))

# (3) merge chemical shift data
cs <- merge(expcs, predcs, by = c("resname", "resid", "nucleus"))
cs <- cs[order(cs$conformation, cs$resid, cs$resname, cs$nucleus),]

# (4) write out vector, corresponding to observed shifts and matrix, corresponding to predict shifts
tmp <- subset(cs, conformation==1)
observed <- tmp$expCS
predicted <- matrix(cs$predCS, nrow = length(observed), byrow = TRUE)
write.table(observed, file = "data/observed_vector.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(predicted, file = "data/predicted_matrix.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
