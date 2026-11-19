#!/usr/bin/env Rscript
# Runs irace on a scenario file.   Rscript tuning/irace_run.R SCENARIO.txt [--check]
#   --check: only check the scenario (irace executes the target runner twice on the first configuration) and stop.
# Otherwise: irace(scenario); the final elite configurations are written to final_elites.tsv next to the log file.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript tuning/irace_run.R SCENARIO.txt [--check]")
suppressPackageStartupMessages(library(irace))
cat("irace", as.character(packageVersion("irace")), "|", R.version.string, "\n")
scenario <- readScenario(filename = args[1])
if ("--check" %in% args) {
  ok <- checkIraceScenario(scenario)
  quit(status = if (isTRUE(ok)) 0 else 1)
}
elites <- irace(scenario = scenario)
out <- file.path(dirname(scenario$logFile), "final_elites.tsv")
utils::write.table(elites, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("final elite configurations written to", out, "\n")
