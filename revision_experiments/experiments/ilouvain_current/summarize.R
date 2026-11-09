#!/usr/bin/env Rscript
out <- Sys.getenv("REVISION_DANCER_OUT", "results/ilouvain_current")
p <- read.csv(file.path(out, "proposed_results.csv"))
i <- read.csv(file.path(out, "ilouvain_results.csv"))
keys <- c("graph", "attribute_condition")
stopifnot(!anyDuplicated(p[c(keys, "repetition")]), !anyDuplicated(i[keys]))
summary <- do.call(rbind, lapply(split(p, interaction(p$graph, p$attribute_condition, drop = TRUE)), function(d) {
  data.frame(d[1L, keys], proposed_attempted = nrow(d), proposed_returned = sum(is.finite(d$ARI)),
    proposed_joint_stops = sum(d$stop_reason == "joint_stability"), proposed_caps = sum(d$stop_reason == "safety_cap"),
    proposed_ARI_mean = mean(d$ARI, na.rm = TRUE), proposed_ARI_SD = sd(d$ARI, na.rm = TRUE))
}))
summary <- merge(summary, i[c(keys, "ARI", "K_found", "external_code_sha256")], by = keys)
names(summary)[names(summary) == "ARI"] <- "ilouvain_ARI"
names(summary)[names(summary) == "K_found"] <- "ilouvain_K_found"
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
