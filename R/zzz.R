.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    crayon::green(
      "metid,
More information can be found at https://tidymass.github.io/metid/
If you use metid in you publication, please cite this publication:
Metabolic reaction network-based recursive metabolite annotation for untargeted metabolomics.
Authors: Xiaotao Shen (shenxt1990@163.com)
Maintainer: Xiaotao Shen."
    )
# cat(crayon::green(
#   c(
#     "                _    _____  ___ ",
#     " _ __ ___   ___| |_  \\_   \\/   \\",
#     "| '_ ` _ \\ / _ \\ __|  / /\\/ /\\ /",
#     "| | | | | |  __/ |_/\\/ /_/ /_// ",
#     "|_| |_| |_|\\___|\\__\\____/___,'  ",
#     "                                "
#   )
# ), sep = "\n")
  )
}

globalVariables(
  names = c(
    "MS1.peak.name",
    "MS2.spectra.name",
    "Candidate.number",
    "Compound.name",
    "CAS.ID",
    "HMDB.ID",
    "KEGG.ID",
    "Lab.ID",
    "Adduct",
    "mz.error",
    "mz.match.score",
    "RT.error",
    "RT.match.score",
    "CE",
    "SS",
    "Total.score",
    "Database",
    "Peak.name",
    "name",
    ".",
    "Cor.RT",
    "Ref.RT",
    "Degree.Span.MSE",
    "Poly.MSE",
    "adduct",
    "hilic.pos",
    "rp.pos",
    "hilic.neg",
    "rp.neg",
    "database",
    "ms1.match.ppm",
    "candidate.num",
    "threads",
    "Exp.mz",
    "Exp.intensity",
    "Lib.mz",
    "Lib.intensity",
    "extenstion",
    "Identification",
    "Name",
    "ExactMass",
    "DB#",
    "Formula",
    "mz",
    "RT",
    "mz.pos",
    "mz.neg",
    "Submitter",
    "Family",
    "Sub.pathway",
    "Note",
    "V1",
    "V2",
    "intensity"
  )
)