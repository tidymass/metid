

# ###HILIC Positive
# setwd("D:/study/database and library/inhouse/Metabolite database/RT correction/HILIC/POS")
# metabolite.info <- readr::read_csv("metabolite.info_HILIC.csv")
# metabolite.info$RT <- metabolite.info$RT/60
# colnames(metabolite.info)
# submitter <- pull(.data = metabolite.info, var = "Submitter")
# table(submitter)
# mix.b <- metabolite.info[which(submitter == "Mix_B"),c("Lab.ID", "Compound.name", "RT"),
#                          drop = FALSE]
#
# write.csv(mix.b, "mix_b/mix.b.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_B.xlsx",
#       cor.feature.table = "mix.b.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_b")
#
#
# mix.c <- metabolite.info[which(submitter == "Mix_C"),c("Lab.ID", "Compound.name", "RT"),
#                          drop = FALSE]
#
# write.csv(mix.c, "mix_c/mix.c.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_C.xlsx",
#       cor.feature.table = "mix.c.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_c")
#
#
# mix.d <- metabolite.info[which(submitter == "Mix_D"),c("Lab.ID", "Compound.name", "RT"),
#                          drop = FALSE]
#
# write.csv(mix.d, "mix_d/mix.d.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_D.xlsx",
#       cor.feature.table = "mix.d.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_d")
#
#
# mix.e <- metabolite.info[which(submitter == "Mix_E"),c("Lab.ID", "Compound.name", "RT"),
#                          drop = FALSE]
#
# write.csv(mix.e, "mix_e/mix.e.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_E.xlsx",
#       cor.feature.table = "mix.e.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_e")
#
#
# mix.f <- metabolite.info[which(submitter == "Mix_F"),c("Lab.ID", "Compound.name", "RT"),
#                          drop = FALSE]
#
# write.csv(mix.f, "mix_f/mix.f.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_F.xlsx",
#       cor.feature.table = "mix.f.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_f")
#
#
# mix.mm <- metabolite.info[which(submitter == "Mix_MM"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
#
# write.csv(mix.mm, "mix_mm/mix.mm.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_MM.xlsx",
#       cor.feature.table = "mix.mm.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_mm")
#
#
# pool.1 <- metabolite.info[which(submitter == "Pool_1"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
#
# write.csv(pool.1, "pool_1/pool.1.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_1.xlsx",
#       cor.feature.table = "pool.1.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_1")
#
#
# pool.2 <- metabolite.info[which(submitter == "Pool_2"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
#
# write.csv(pool.2, "pool_2/pool.2.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_2.xlsx",
#       cor.feature.table = "pool.2.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_2")
#
# pool.3 <- metabolite.info[which(submitter == "Pool_3"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(pool.3, "pool_3/pool.3.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_3.xlsx",
#       cor.feature.table = "pool.3.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_3")
#
#
# pool.4 <- metabolite.info[which(submitter == "Pool_4"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(pool.4, "pool_4/pool.4.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_4.xlsx",
#       cor.feature.table = "pool.4.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_4")
#
#
# pool.5 <- metabolite.info[which(submitter == "Pool_5"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(pool.5, "pool_5/pool.5.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_5.xlsx",
#       cor.feature.table = "pool.5.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_5")
#
#
# pool.6 <- metabolite.info[which(submitter == "Pool_6"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(pool.6, "pool_6/pool.6.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_6.xlsx",
#       cor.feature.table = "pool.6.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_6")
#
#
# pool.7 <- metabolite.info[which(submitter == "Pool_7"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(pool.7, "pool_7/pool.7.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_7.xlsx",
#       cor.feature.table = "pool.7.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_7")
#
#
# ipop <- metabolite.info[which(submitter == "iPOP"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(ipop, "ipop/ipop.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "iPOP.xlsx",
#       cor.feature.table = "ipop.csv",
#       return.predicted.rt = FALSE,
#       path = "ipop")
#
#
# exercise <- metabolite.info[which(submitter == "iPOP_Exercise"),c("Lab.ID", "Compound.name", "RT"),
#                         drop = FALSE]
# write.csv(ipop, "exercise/exercise.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Exercise.xlsx",
#       cor.feature.table = "exercise.csv",
#       return.predicted.rt = FALSE,
#       path = "exercise")
#
# ``
# file.name <- c("mix_b", "mix_c", "mix_d", "mix_e", "mix_f", "mix_mm",
#                "pool_1", "pool_2", "pool_3", "pool_4", "pool_5", "pool_6", "pool_7",
#                "ipop", "exercise")
#
# rt.table <- pbapply::pblapply(file.name, function(x){
#   readr::read_csv(file.path(x, "Feature.table.with.corrected.RT.csv"))
# })
#
# rt.table <- do.call(rbind, rt.table)
#
# plot(rt.table$RT - rt.table$RT.correction)
#
# plot(rt.table$RT * 60 - rt.table$RT.correction * 60)
#
# rt.correction <- rt.table$RT.correction
# names(rt.correction) <- rt.table$Lab.ID
#
# metabolite.info$RT[match(names(rt.correction), metabolite.info$Lab.ID)] <- rt.correction
# metabolite.info$RT <- metabolite.info$RT * 60
#
#
# write.csv(metabolite.info, "metabolite.info.pos.csv")
#
#
#
#
#
#
#
# ###HILIC Negative
# setwd("D:/study/database and library/inhouse/Metabolite database/RT correction/HILIC/NEG")
# metabolite.info <- readr::read_csv("metabolite.info_HILIC.csv")
# metabolite.info$RT <- metabolite.info$RT/60
# colnames(metabolite.info)
# submitter <- pull(.data = metabolite.info, var = "Submitter")
# table(submitter)
# mix.b <- metabolite.info[which(submitter == "Mix_B"),c("Lab.ID", "Compound.name", "RT"),
#                          drop = FALSE]
#
# write.csv(mix.b, "mix_b/mix.b.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_B.xlsx",
#       cor.feature.table = "mix.b.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_b")
#
#
# mix.c <- metabolite.info[which(submitter == "Mix_C"),c("Lab.ID", "Compound.name", "RT"),
#                          drop = FALSE]
#
# write.csv(mix.c, "mix_c/mix.c.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_C.xlsx",
#       cor.feature.table = "mix.c.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_c")
#
#
# mix.d <- metabolite.info[which(submitter == "Mix_D"),c("Lab.ID", "Compound.name", "RT"),
#                          drop = FALSE]
#
# write.csv(mix.d, "mix_d/mix.d.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_D.xlsx",
#       cor.feature.table = "mix.d.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_d")
#
#
# mix.e <- metabolite.info[which(submitter == "Mix_E"),c("Lab.ID", "Compound.name", "RT"),
#                          drop = FALSE]
#
# write.csv(mix.e, "mix_e/mix.e.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_E.xlsx",
#       cor.feature.table = "mix.e.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_e")
#
#
# mix.f <- metabolite.info[which(submitter == "Mix_F"),c("Lab.ID", "Compound.name", "RT"),
#                          drop = FALSE]
#
# write.csv(mix.f, "mix_f/mix.f.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_F.xlsx",
#       cor.feature.table = "mix.f.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_f")
#
#
# mix.mm <- metabolite.info[which(submitter == "Mix_MM"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
#
# write.csv(mix.mm, "mix_mm/mix.mm.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Mix_MM.xlsx",
#       cor.feature.table = "mix.mm.csv",
#       return.predicted.rt = FALSE,
#       path = "mix_mm")
#
#
# pool.1 <- metabolite.info[which(submitter == "Pool_1"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
#
# write.csv(pool.1, "pool_1/pool.1.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_1.xlsx",
#       cor.feature.table = "pool.1.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_1")
#
#
# pool.2 <- metabolite.info[which(submitter == "Pool_2"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
#
# write.csv(pool.2, "pool_2/pool.2.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_2.xlsx",
#       cor.feature.table = "pool.2.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_2")
#
# pool.3 <- metabolite.info[which(submitter == "Pool_3"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(pool.3, "pool_3/pool.3.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_3.xlsx",
#       cor.feature.table = "pool.3.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_3")
#
#
# pool.4 <- metabolite.info[which(submitter == "Pool_4"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(pool.4, "pool_4/pool.4.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_4.xlsx",
#       cor.feature.table = "pool.4.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_4")
#
#
# pool.5 <- metabolite.info[which(submitter == "Pool_5"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(pool.5, "pool_5/pool.5.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_5.xlsx",
#       cor.feature.table = "pool.5.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_5")
#
#
# pool.6 <- metabolite.info[which(submitter == "Pool_6"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(pool.6, "pool_6/pool.6.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_6.xlsx",
#       cor.feature.table = "pool.6.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_6")
#
#
# pool.7 <- metabolite.info[which(submitter == "Pool_7"),c("Lab.ID", "Compound.name", "RT"),
#                           drop = FALSE]
# write.csv(pool.7, "pool_7/pool.7.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Pool_7.xlsx",
#       cor.feature.table = "pool.7.csv",
#       return.predicted.rt = FALSE,
#       path = "pool_7")
#
# ipop <- metabolite.info[which(submitter == "iPOP"),c("Lab.ID", "Compound.name", "RT"),
#                         drop = FALSE]
# write.csv(ipop, "ipop/ipop.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "iPOP.xlsx",
#       cor.feature.table = "ipop.csv",
#       return.predicted.rt = FALSE,
#       path = "ipop")
#
#
# exercise <- metabolite.info[which(submitter == "iPOP_Exercise"),c("Lab.ID", "Compound.name", "RT"),
#                             drop = FALSE]
# write.csv(ipop, "exercise/exercise.csv", row.names = FALSE)
#
# correct_rt(ref.is.table = "ref.is.table.xlsx",
#       cor.is.table = "Exercise.xlsx",
#       cor.feature.table = "exercise.csv",
#       return.predicted.rt = FALSE,
#       path = "exercise")
#
#
# file.name <- c("mix_b", "mix_c", "mix_d", "mix_e", "mix_f", "mix_mm",
#                "pool_1", "pool_2", "pool_3", "pool_4", "pool_5", "pool_6", "pool_7",
#                "ipop", "exercise")
#
#
# rt.table <- pbapply::pblapply(file.name, function(x){
#   readr::read_csv(file.path(x, "Feature.table.with.corrected.RT.csv"))
# })
#
# rt.table <- do.call(rbind, rt.table)
#
# plot(rt.table$RT - rt.table$RT.correction)
#
# plot(rt.table$RT * 60 - rt.table$RT.correction * 60)
#
# rt.correction <- rt.table$RT.correction
# names(rt.correction) <- rt.table$Lab.ID
#
# metabolite.info$RT[match(names(rt.correction), metabolite.info$Lab.ID)] <- rt.correction
# metabolite.info$RT <- metabolite.info$RT * 60
#
#
# write.csv(metabolite.info, "metabolite.info.neg.csv")
#
#
# met.hilic.pos <-
#   readr::read_csv("D:/study/database and library/inhouse/Metabolite database/RT correction/HILIC/POS/metabolite.info.pos.csv")
#
# met.hilic.neg <-
#   readr::read_csv("D:/study/database and library/inhouse/Metabolite database/RT correction/HILIC/NEG/metabolite.info.neg.csv")
#
# dim(met.hilic.pos)
# dim(met.hilic.neg)
#
# plot(met.hilic.pos$RT, met.hilic.neg$RT)
# abline(0, 1)
#
#
# setwd("D:/study/database and library/inhouse/Metabolite database/database construction/RPLC")
# databaseConstruction(path = ".", version = "0.0.1",
#                      metabolite.info.name = "metabolite.info_RPLC.csv")
#
# load("msDatabase0.0.1")
# msDatabase_rplc0.0.1 <- msDatabase0.0.1
# msDatabase_rplc0.0.1
#
# ###replace RT
# rt.rplc <- apply(cbind(met.rplc.pos$RT, met.rplc.neg$RT), 1, mean)
# names(rt.rplc) <- met.rplc.pos$Lab.ID
# plot(msDatabase_rplc0.0.1@spectra.info$RT[match(names(rt.rplc), msDatabase_rplc0.0.1@spectra.info$Lab.ID)],
# rt.rplc)
#
# msDatabase_rplc0.0.1@spectra.info$RT[match(names(rt.rplc), msDatabase_rplc0.0.1@spectra.info$Lab.ID)] <-
#   rt.rplc
#
# save(msDatabase_rplc0.0.1, file = "msDatabase_rplc0.0.1")
#
#
#
#
#
# setwd("D:/study/database and library/inhouse/Metabolite database/database construction/HILIC")
# databaseConstruction(path = ".", version = "0.0.1",
#                      metabolite.info.name = "metabolite.info_HILIC.csv")
#
# load("msDatabase0.0.1")
# msDatabase_hilic0.0.1 <- msDatabase0.0.1
# msDatabase_hilic0.0.1
#
# ###replace RT
# rt.hilic <- apply(cbind(met.hilic.pos$RT, met.hilic.neg$RT), 1, mean)
# names(rt.hilic) <- met.hilic.pos$Lab.ID
# plot(msDatabase_hilic0.0.1@spectra.info$RT[match(names(rt.hilic), msDatabase_hilic0.0.1@spectra.info$Lab.ID)],
#      rt.hilic)
#
# msDatabase_hilic0.0.1@spectra.info$RT[match(names(rt.hilic), msDatabase_hilic0.0.1@spectra.info$Lab.ID)] <-
#   rt.hilic
#
# save(msDatabase_hilic0.0.1, file = "msDatabase_hilic0.0.1")


#' @title Correct RTs of metabolites in database using internal standard list
#' @description Correct RTs of metabolites in database using internal standard list.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param experiment.is.table Experiment internal standard table.
#' @param database.is.table Database internal standard table.
#' @param database Database, must be datasetClass object.
#' @param method polyline or loess.
#' @param poly A numeric vector.
#' @param degree A numeric vector.
#' @param path Work directory.
#' @return A metIdentifyClass object.
#' @export

correct_database_rt = function(experiment.is.table,
                               database.is.table,
                               database,
                               method = c("polyline", "loess"),
                               poly = c(1, 2, 3, 4, 5),
                               degree = c(1, 2),
                               path = "."){
  dir.create(path, showWarnings = FALSE)
  method <- match.arg(method)
  
  cor.feature.table <-
    as.data.frame(database@spectra.info[, c("Lab.ID", "Compound.name", "RT")])
  write.csv(
    cor.feature.table,
    file = file.path(path, "cor.feature.table.csv"),
    row.names = FALSE
  )
  
  predicted.rt <-
    correct_rt(
      ref.is.table = experiment.is.table,
      cor.is.table = database.is.table,
      cor.feature.table = "cor.feature.table.csv",
      method = method,
      poly = poly,
      degree = degree,
      path = path,
      return.predicted.rt = TRUE
    )
  unlink(x = file.path(path, "cor.feature.table.csv"))
  predicted.rt <- predicted.rt$RT.correction
  database@spectra.info$RT <- predicted.rt
  return(database)
}



correct_rt = function(ref.is.table,
                      cor.is.table,
                      cor.feature.table,
                      method = c("polyline", "loess"),
                      poly = c(1, 2, 3, 4, 5),
                      degree = c(1, 2),
                      path = ".",
                      return.predicted.rt = TRUE){
  # dir.create(path)
  method <- match.arg(method)
  
  ref.is.table <-
    readTable(file = file.path(path, ref.is.table))
  cor.is.table <-
    readTable(file = file.path(path, cor.is.table))
  cor.feature.table <-
    readTable(file = file.path(path, cor.feature.table))
  ref.is.table <- as.data.frame(ref.is.table)
  cor.is.table <- as.data.frame(cor.is.table)
  cor.feature.table <- as.data.frame(cor.feature.table)
  
  ref.is.table[, 3] <- as.numeric(ref.is.table[, 3])
  cor.is.table[, 3] <- as.numeric(cor.is.table[, 3])
  cor.feature.table[, 3] <-
    as.numeric(cor.feature.table[, 3])
  
  remain.idx <-
    which(!is.na(ref.is.table[, 3]) & !is.na(cor.is.table[, 3]))
  ref.is.table <-
    ref.is.table[remain.idx, , drop = FALSE]
  cor.is.table <-
    cor.is.table[remain.idx, , drop = FALSE]
  
  rt.table <- data.frame(ref = ref.is.table[, 3],
                         cor = cor.is.table[, 3],
                         stringsAsFactors = FALSE)
  
  rt1 <- rt.table$ref
  rt2 <- rt.table$cor
  
  if (method == "polyline") {
    poly <- poly[which(poly < nrow(rt.table) - 1)]
    mse <- bestpoly(
      x = rt.table$cor,
      y = rt.table$ref,
      poly = poly,
      path = path
    )
    idx <- poly[which.min(mse)]
    
    model <- lm(rt1 ~ poly(rt2, idx))
  } else{
    result <- bestloess(
      x = rt.table$cor,
      y = rt.table$ref,
      span.begin = 0.5,
      span.end = 1,
      span.step = 0.1,
      degree = degree,
      path = path
    )
    idx <- which.min(result[, 3])
    model <- loess(rt1 ~ rt2,
                   span = result[idx, 2],
                   degree = result[idx, 1])
    
  }
  
  raw.rt <- cor.feature.table$RT
  
  predict.rt <- predict(object = model,
                        newdata = data.frame(rt2 = raw.rt))
  predict.rt[raw.rt < min(rt2)] <-
    raw.rt[raw.rt < min(rt2)]
  predict.rt[raw.rt > max(rt2)] <-
    raw.rt[raw.rt > max(rt2)]
  
  temp.data <- data.frame(raw.rt, predict.rt,
                          stringsAsFactors = FALSE)
  
  # temp.data
  
  plot <-
    ggplot2::ggplot(data = temp.data, ggplot2::aes(x = raw.rt, y = predict.rt)) +
    ggplot2::geom_point(color = "black") +
    # ggplot2::geom_line() +
    ggplot2::geom_smooth(method = "loess",
                         color = "black",
                         fill = "gray") +
    ggplot2::geom_abline(intercept = 0,
                         slope = 1,
                         color = "red") +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    filename = file.path(path, "Raw.RT vs Predict.RT plot.pdf"),
    plot = plot,
    width = 8,
    height = 6
  )
  
  cor.feature.table <- data.frame(cor.feature.table,
                                  "RT.correction" = predict.rt,
                                  stringsAsFactors = FALSE)
  if (return.predicted.rt) {
    return(cor.feature.table)
  } else{
    write.csv(
      cor.feature.table,
      file = file.path(path, "Feature.table.with.corrected.RT.csv"),
      row.names = FALSE
    )
  }
}




###parameters optimization
bestpoly = function(x,
                    y,
                    poly = c(1, 2, 3, 4),
                    path = "."){
  pre.all <- list()
  for (i in seq_along(poly)) {
    # cat(i, " ")
    ##LOO
    pre <- NULL
    for (j in seq_along(x)) {
      # cat(j, " ")
      y1 <- y[-j]
      x1 <- x[-j]
      model <- lm(y1 ~ poly(x1, poly[i]))
      pre[j] <- predict(object = model,
                        newdata = data.frame(x1 = x[j]))
    }
    pre.all[[i]] <- pre
  }
  
  mse <-
    unlist(lapply(pre.all, function(x) {
      sum((y - x) ^ 2) / length(y)
    }))
  
  # y.lim1 <- 0.8*min(c(y, unlist(lapply(pre.all, min))))
  # y.lim2 <- 1.2*max(c(y, unlist(lapply(pre.all, max))))
  
  temp.data1 <- do.call(c, pre.all)
  temp.data1 <- data.frame(
    Cor.RT = rep(x, length(pre.all)),
    Ref.RT = temp.data1,
    Poly.MSE = as.character(rep(paste(
      poly, round(mse, 3), sep = ": "
    ), each = length(x))),
    stringsAsFactors = FALSE
  )
  
  temp.data2 <- data.frame(Cor.RT = x,
                           Ref.RT = y,
                           stringsAsFactors = FALSE)
  
  plot <-
    ggplot2::ggplot(data = temp.data2,
                    ggplot2::aes(x = Cor.RT, y = Ref.RT)) +
    ggplot2::geom_point(color = "black") +
    ggplot2::geom_line() +
    ggplot2::geom_smooth(method = "loess",
                         color = "black",
                         fill = "gray") +
    ggplot2::geom_point(data = temp.data1,
                        mapping = ggplot2::aes(x = Cor.RT, y = Ref.RT, color = Poly.MSE)) +
    ggplot2::geom_line(data = temp.data1,
                       mapping = ggplot2::aes(x = Cor.RT, y = Ref.RT, color = Poly.MSE)) +
    ggplot2::geom_smooth(
      data = temp.data1,
      mapping = ggplot2::aes(
        x = Cor.RT,
        y = Ref.RT,
        color = Poly.MSE,
        fill = Poly.MSE
      ),
      method = "loess"
    ) +
    # ggthemes::theme_pander() +
    ggplot2::theme_bw()
  ggplot2::ggsave(
    filename = file.path(path, "MSE.plot.pdf"),
    plot = plot,
    width = 8,
    height = 6
  )
  
  return(mse)
}



bestloess = function(x,
                     y,
                     span.begin = 0.5,
                     span.end = 1,
                     span.step = 0.1,
                     degree = c(1, 2),
                     path = "."){
  span <- seq(span.begin, span.end, span.step)
  para <- NULL
  for (i in seq_along(degree)) {
    para <- rbind(para, cbind(degree[i], span))
  }
  colnames(para) <- c("degree", "span")
  
  y <- y[order(x)]
  x <- sort(x)
  
  pre.all <- list()
  
  for (i in seq_len(nrow(para))) {
    temp.degree <- para[i, 1]
    temp.span <- para[i, 2]
    pre <- NULL
    # for(j in 2:(length(x) - 1)) {
    for (j in seq_along(x)) {
      y1 <- y
      x1 <- x
      # y1 <- y[-j]
      # x1 <- x[-j]
      model <-
        loess(y1 ~ x1, span = temp.span, degree = temp.degree)
      pre[j] <- predict(object = model,
                        newdata = data.frame(x1 = x[j]))
    }
    pre.all[[i]] <- pre[!is.na(pre)]
  }
  
  mse <-
    unlist(lapply(pre.all, function(x) {
      sum((y[2:(length(y) - 1)] - x) ^ 2) / length(y)
    }))
  result <- data.frame(para, mse, stringsAsFactors = FALSE)
  
  
  temp.data1 <- do.call(c, pre.all)
  para <- paste(para[, 1], para[, 2], sep = ";")
  
  temp.data1 <- data.frame(
    Cor.RT = rep(x, length(pre.all)),
    Ref.RT = temp.data1,
    Degree.Span.MSE = as.character(rep(paste(
      para, round(mse, 3), sep = ": "
    ), each = length(x))),
    stringsAsFactors = FALSE
  )
  
  
  temp.data2 <- data.frame(Cor.RT = x,
                           Ref.RT = y,
                           stringsAsFactors = FALSE)
  
  plot <-
    ggplot2::ggplot(data = temp.data2,
                    ggplot2::aes(x = Cor.RT, y = Ref.RT)) +
    ggplot2::geom_point(color = "black") +
    ggplot2::geom_line() +
    ggplot2::geom_smooth(method = "loess",
                         color = "black",
                         fill = "gray") +
    ggplot2::geom_point(
      data = temp.data1,
      mapping = ggplot2::aes(x = Cor.RT, y = Ref.RT, color = Degree.Span.MSE)
    ) +
    ggplot2::geom_line(
      data = temp.data1,
      mapping = ggplot2::aes(x = Cor.RT, y = Ref.RT, color = Degree.Span.MSE)
    ) +
    ggplot2::geom_smooth(
      data = temp.data1,
      mapping = ggplot2::aes(
        x = Cor.RT,
        y = Ref.RT,
        color = Degree.Span.MSE,
        fill = Degree.Span.MSE
      ),
      method = "loess"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    filename = file.path(path, "MSE.plot.pdf"),
    plot = plot,
    width = 8,
    height = 6
  )
  return(result)
} 
