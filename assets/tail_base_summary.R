options(warn=-1)
library(tidyverse)
library(openxlsx)
library(data.table)
library(reshape2)
library(lubridate)


#######
### args
args <- commandArgs(trailingOnly = T)

sample <- args[1]
gmc_file <- args[2]
profile_summary_file <- args[7]
doc_path <- paste0(args[8], "/")
plot_path <- paste0(args[9], "/")
bin_path <- paste0(args[10], "/")

dir.create(doc_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(args[3],sep = "\t", header = TRUE)

#######
### env
ref <- read.table(args[4])
ma <- read.table(args[5], sep = "\t")
ma1 <- read.table(args[6], sep = "\t")[-2]

ensure_base_columns <- function(df) {
  for (base in c("A", "C", "G", "T")) {
    if (!(base %in% names(df))) {
      df[[base]] <- 0
    }
  }
  df
}

#######
### function
tt_barplot <- function(mtfile) {
  df1 <- mtfile[, c(2, 13:16)]
  df1 <- gather(df1, base, percentage, "A_percentage":"G_percentage")
  df1$base <- gsub("A_percentage", "Tailing A", df1$base)
  df1$base <- gsub("T_percentage", "Tailing T", df1$base)
  df1$base <- gsub("C_percentage", "Tailing C", df1$base)
  df1$base <- gsub("G_percentage", "Tailing G", df1$base)
  df2 <- get(paste0(sample, "_trim_summary"))[, c(1, 4)]
  df2$base <- "Trimming"
  df2 <- df2[, c(1, 3, 2)]
  names(df2) <- c("tail_N", "base", "percentage")
  df2$percentage <- 0 - df2$percentage
  tmp <- rbind(df1, df2)
  tmp$tail_N <- as.factor(tmp$tail_N)

  P <- tmp %>%
    filter(!(tail_N %in% "0")) %>%
    ggplot() +
    aes(x = factor(tail_N, levels = unique(tail_N)), y = percentage, fill = base) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(`Tailing A` = "#00FFFF", `Tailing C` = "#FF0000", `Tailing G` = "#00FF00", `Tailing T` = "#0000FF", Trimming = "#FFA500")) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(labels = abs, breaks = seq(-100, 100, 20), limits = c(-100, 100)) +
    labs(x = "Tailing and trimming length(nt)", y = "Relative Percentage(%)", fill = NULL, title = sample)
}
extract_matrix <- function(df1) {
  # 将输入数据转换为data.table格式
  setDT(df1)

  # 使用lapply代替for循环，处理每个name
  result_list <- lapply(ref$V1, function(name) {
    len <- ref[which(ref$V1 == name), 4]

    PART1 <- df1[ID == name & aligned_length > len]
    if (nrow(PART1) > 0) {
      PART1[, miR_len := len]
      PART1[, tail_length := raw_len - len]
      PART1[, trim_length := 0]
      PART1[, tail_base := str_sub(raw_seq, start = raw_len - tail_length + 1)]
    }

    PART2 <- df1[ID == name & aligned_length <= len]
    if (nrow(PART2) > 0) {
      PART2[, miR_len := len]
      PART2[, tail_length := base]
      PART2[, trim_length := len - aligned_length]
      PART2[, tail_base := str_sub(raw_seq, start = raw_len - tail_length + 1)]
    }

    rbindlist(list(PART1, PART2), use.names = TRUE, fill = TRUE)
  })

  # 合并结果为一个数据表
  tail <- rbindlist(result_list, use.names = TRUE, fill = TRUE)
  return(tail)
}

#######
### pre-processing data
print(paste0("processing: ", sample))

meta$Sample <- gsub("-", "_", meta$Sample)
total <- meta[which(meta$Sample == sample), 3]
if (length(total) != 1) {
  stop(paste0("Expected one mapping-summary row for sample '", sample,
              "', found ", length(total)))
}

df <- fread(gmc_file, header = TRUE)
df <- df[, -c(3, 4)]
df1 <- separate(df, RNAME, into = c("sample", "rank", "raw_seq", "base"), sep = "-+", convert = TRUE,remove = FALSE)
df1$raw_len <- nchar(df1$raw_seq)
df1$aligned_length <- df1$raw_len - df1$base

# tailing and trimming summary file, there are two part you should know, detail can be found on notebook
tail <- extract_matrix(df1 = df1)

#######
### calculate tailing and trimming length

# template tail 1-10 summary
for (nt in 1:10) {

  if (nt == 1) {
    tmp_1 <- tail %>%
      filter(tail_length == nt, tail_base != "N") %>%
      select(6, 4, 7:12) %>%
      dplyr::arrange(ID, raw_seq) %>%
      left_join(ma, by = c("ID" = "V1"))

    tmp_1_summary <- tmp_1 %>%
      group_by(ID, tail_base) %>%
      count(name = "count") %>%
      ungroup() %>%
      mutate(percentage = round(((count / sum(count)) * 100), 4))

    tmp_2 <- tmp_1 %>%
      group_by_all() %>%
      count(name = "count") %>%
      select(ID, short_mechanism = V2, everything()) %>%
      ungroup() %>%
      mutate(percentage = round(((count / sum(count)) * 100), 4))

    tmp_1_summary1 <- spread(tmp_1_summary[1:3], key = tail_base, value = count, fill = 0) %>% ensure_base_columns()
    tmp_1_summary2 <- spread(tmp_1_summary[c(1:2, 4)], key = tail_base, value = percentage, fill = 0) %>% ensure_base_columns() %>%
      dplyr::rename("A_percentage"=A, "C_percentage"=C, "G_percentage"=G, "T_percentage"=T)

    tmp_summary <- fread(profile_summary_file, header = TRUE)
    tmp_summary <- tmp_summary %>%
      select("# ID", SUM, tail_1) %>%
      mutate("miRNA_tail_1/miRNA_total" = round(tail_1 / SUM * 100, 4))

    tail_1_summary <- tmp_1_summary1 %>%
      left_join(tmp_1_summary2, by = "ID") %>%
      left_join(tmp_summary, by = c("ID" = "# ID")) %>%
      left_join(ma, by = c("ID" = "V1")) %>%
      mutate(
        A_RPM = round((A / total * 10^6), 2),
        C_RPM = round((C / total * 10^6), 2),
        G_RPM = round((G / total * 10^6), 2),
        T_RPM = round((T / total * 10^6), 2)
      ) %>%
      select(ID, short_mechanism = V2, SUM, 2:5, ends_with("RPM"), ends_with("RPM"), `miRNA_tail_1/miRNA_total`)

    tail_list <- setNames(list(tmp_1_summary, tail_1_summary),
                          c(paste0("tail_", nt, "_total"), paste0("tail_", nt, "_summary")))
    write.xlsx(tail_list, paste0(doc_path, sample, "_tail_",nt,"_for_publish.xlsx"))
  } else {
    tmp_1 <- tail %>%
      filter(tail_length == nt, tail_base != "N") %>%
      select(6, 4, 7:12) %>%
      dplyr::arrange(ID, raw_seq) %>%
      left_join(ma, by = c("ID" = "V1"))

    tmp_2 <- tmp_1 %>%
      group_by_all() %>%
      count(name = "count") %>%
      ungroup() %>%
      mutate(percentage = round(((count / sum(count)) * 100), 4))

    tail_1_summary <- tmp_1 %>%
      group_by(ID, tail_base) %>%
      summarize(count = n()) %>%
      mutate(part_percentage = round(((count / sum(count)) * 100), 4) ) %>%
      ungroup() %>%
      mutate(global_percentage = round(((count / sum(count)) * 100), 4))

    tmp_summary <- fread(profile_summary_file, header = TRUE)
    tmp_summary <- tmp_summary %>%
      left_join( ma, by = c("# ID" = "V1")) %>%
      select("# ID", short_mechanism = V2, SUM, tail = paste0("tail_", nt)) %>% mutate(!!paste0("miRNA_tail_", nt, "/miRNA_total") := round( tail / SUM * 100, 4))

  }

  tail_list <- setNames(
    list(tmp_summary, tail_1_summary, tmp_2),
    c(paste("tail", nt, "summary",sep = "_"), paste("tail", nt, "detail",sep = "_"), paste("tail", nt, "raw",sep = "_"))
  )
  assign(paste(sample, "tail",nt,"summary",sep = "_"), tmp_2)
  write.xlsx(tail_list, paste0(doc_path, sample, "_tail_",nt,"_summary.xlsx"))
}

# non-template tail 1-10 summary, from here group miRNA, because non-template tail don't effect hairpin sequence
mirna_count <- tail %>%
  left_join(ma1, by = c("ID" = "V1")) %>%
  group_by(V3) %>%
  tally()

for (nt in 1:10) {

  if (nt == 1) {

    tmp_1 <- tail %>%
      filter(base != 0, tail_length == nt, tail_base != "N") %>%
      dplyr::arrange(ID, raw_seq) %>%
      left_join(ma1, by = c("ID" = "V1")) %>%
      select(ID = V3, raw_seq, short_mechanism = V4, 7:12)

    tmp_1_summary <- tmp_1 %>%
      group_by(ID, tail_base) %>%
      count(name = "count") %>%
      ungroup() %>%
      mutate(percentage = round(((count / sum(count)) * 100), 4))

    tmp_2 <- tmp_1 %>%
      group_by_all() %>%
      count(name = "count") %>%
      ungroup() %>%
      mutate(percentage = round(((count / sum(count)) * 100), 4))

    tmp_1_summary1 <- spread(tmp_1_summary[1:3], key = tail_base, value = count, fill = 0) %>% ensure_base_columns()
    tmp_1_summary2 <- spread(tmp_1_summary[c(1:2, 4)], key = tail_base, value = percentage, fill = 0) %>% ensure_base_columns() %>% dplyr::rename("A_percentage"=A, "C_percentage"=C, "G_percentage"=G, "T_percentage"=T)
    tmp_1_summary3 <- aggregate(count ~ ID, data = tmp_1_summary, FUN = sum)


    tmp_summary <- tmp_1_summary3 %>%
      left_join(mirna_count, by = c("ID" = "V3")) %>%
      select("ID", SUM = n, tail_1 = count) %>%
      mutate("miRNA_tail_1/miRNA_total" = round(tail_1 / SUM * 100, 4))

    tail_1_summary <- tmp_1_summary1 %>%
      left_join(tmp_1_summary2, by = "ID") %>%
      left_join(tmp_summary, by = "ID") %>%
      left_join(ma1[!duplicated(ma1$V3),], by = c("ID" = "V3")) %>%
      mutate(
        A_RPM = round((A / total * 10^6), 2),
        C_RPM = round((C / total * 10^6), 2),
        G_RPM = round((G / total * 10^6), 2),
        T_RPM = round((T / total * 10^6), 2)
      ) %>%
      select(ID, short_mechanism = V4, SUM, A,C,G,T, ends_with("RPM"), ends_with("RPM"), `miRNA_tail_1/miRNA_total`)

    tail_list <- setNames(list(tmp_1_summary, tail_1_summary),
                          c(paste0("tail_", nt, "_total"), paste0("tail_", nt, "_summary")))
    write.xlsx(tail_list, paste0(doc_path, sample, "_tail_",nt,"_nt_for_publish.xlsx"))
  } else {
    tmp_1 <- tail %>%
      filter(base != 0, tail_length == nt, tail_base != "N") %>%
      dplyr::arrange(ID, raw_seq) %>%
      left_join(ma1, by = c("ID" = "V1")) %>%
      select(ID = V3, raw_seq, short_mechanism = V4, 7:12)

    tmp_2 <- tmp_1 %>%
      group_by_all() %>%
      count(name = "count") %>%
      ungroup() %>%
      mutate(percentage = round(((count / sum(count)) * 100), 4))

    tail_1_summary <- tmp_1 %>%
      group_by(ID, tail_base) %>%
      summarize(count = n()) %>%
      mutate(part_percentage = round(((count / sum(count)) * 100), 4)) %>%
      ungroup() %>%
      mutate(global_percentage = round(((count / sum(count)) * 100), 4))

    tmp_1_summary <- tmp_1 %>%
      group_by(ID, tail_base) %>%
      count(name = "count") %>%
      aggregate(count ~ ID, FUN = sum)

    tmp_summary <- tmp_1_summary %>%
      left_join(mirna_count, by = c("ID" = "V3")) %>%
      select("ID", SUM = n, tail = count) %>%
      mutate(!!paste0("miRNA_tail_", nt, "/miRNA_total") := round( tail / SUM * 100, 4))
  }
  tail_list <- setNames(
    list(tmp_summary, tail_1_summary, tmp_2),
    c(paste("tail", nt, "summary",sep = "_"), paste("tail", nt, "detail",sep = "_"), paste("tail", nt, "raw",sep = "_"))
  )
  assign(paste(sample, "tail",nt,"summary",sep = "_"), tmp_2)
  write.xlsx(tail_list, paste0(doc_path,sample, "_tail_",nt,"_summary.xlsx"))
}

#####
# summary file, here summary results is unique id.
tail <- tail[!duplicated(tail$RNAME),]

mirna_count <- tail %>%
  left_join(ma1, by = c("ID" = "V1")) %>%
  group_by(V3) %>%
  tally()

# non-template tail 1-10 summary, duplicated.
for (nt in 1:10) {

  if (nt == 1) {

    tmp_1 <- tail %>%
      filter(base != 0, tail_length == nt, tail_base != "N") %>%
      dplyr::arrange(ID, raw_seq) %>%
      left_join(ma1, by = c("ID" = "V1")) %>%
      select(ID = V3, raw_seq, short_mechanism = V4, 7:12)

    tmp_1_summary <- tmp_1 %>%
      group_by(ID, tail_base) %>%
      count(name = "count") %>%
      ungroup() %>%
      mutate(percentage = round(((count / sum(count)) * 100), 4))

    tmp_2 <- tmp_1 %>%
      group_by_all() %>%
      count(name = "count") %>%
      ungroup() %>%
      mutate(percentage = round(((count / sum(count)) * 100), 4))

    tmp_1_summary1 <- spread(tmp_1_summary[1:3], key = tail_base, value = count, fill = 0) %>% ensure_base_columns()
    tmp_1_summary2 <- spread(tmp_1_summary[c(1:2, 4)], key = tail_base, value = percentage, fill = 0) %>% ensure_base_columns() %>% dplyr::rename("A_percentage"=A, "C_percentage"=C, "G_percentage"=G, "T_percentage"=T)
    tmp_1_summary3 <- aggregate(count ~ ID, data = tmp_1_summary, FUN = sum)

    tmp_summary <- tmp_1_summary3 %>%
      left_join(mirna_count, by = c("ID" = "V3")) %>%
      select("ID", SUM = n, tail_1 = count) %>%
      mutate("miRNA_tail_1/miRNA_total" = round(tail_1 / SUM * 100, 4))

    tail_1_summary <- tmp_1_summary1 %>%
      left_join(tmp_1_summary2, by = "ID") %>%
      left_join(tmp_summary, by = "ID") %>%
      left_join(ma1[!duplicated(ma1$V3),], by = c("ID" = "V3")) %>%
      mutate(
        A_RPM = round((A / total * 10^6), 2),
        C_RPM = round((C / total * 10^6), 2),
        G_RPM = round((G / total * 10^6), 2),
        T_RPM = round((T / total * 10^6), 2)
      ) %>%
      select(ID, short_mechanism = V4, SUM, A,C,G,T, ends_with("RPM"), ends_with("percentage"), `miRNA_tail_1/miRNA_total`)

    tail_list <- setNames(list(tmp_1_summary, tail_1_summary),
                          c(paste0("tail_", nt, "_total"), paste0("tail_", nt, "_summary")))
    write.xlsx(tail_list, paste0(doc_path,sample,  "_tail_",nt,"_nt_duplicated_for_publish.xlsx"))
  } else {
    tmp_1 <- tail %>%
      filter(base != 0, tail_length == nt, tail_base != "N") %>%
      dplyr::arrange(ID, raw_seq) %>%
      left_join(ma1, by = c("ID" = "V1")) %>%
      select(ID = V3, raw_seq, short_mechanism = V4, 7:12)

    tmp_2 <- tmp_1 %>%
      group_by_all() %>%
      count(name = "count") %>%
      ungroup() %>%
      mutate(percentage = round(((count / sum(count)) * 100), 4))

    tail_1_summary <- tmp_1 %>%
      group_by(ID, tail_base) %>%
      summarize(count = n()) %>%
      mutate(part_percentage = round(((count / sum(count)) * 100), 4)) %>%
      ungroup() %>%
      mutate(global_percentage = round(((count / sum(count)) * 100), 4))

    tmp_1_summary <- tmp_1 %>%
      group_by(ID, tail_base) %>%
      count(name = "count") %>%
      aggregate(count ~ ID, FUN = sum)

    tmp_summary <- tmp_1_summary %>%
      left_join(mirna_count, by = c("ID" = "V3")) %>%
      select("ID", SUM = n, tail = count) %>%
      mutate(!!paste0("miRNA_tail_", nt, "/miRNA_total") := round( tail / SUM * 100, 4))


  }
  tail_list <- setNames(
    list(tmp_summary, tail_1_summary, tmp_2),
    c(paste("tail", nt, "summary",sep = "_"), paste("tail", nt, "detail",sep = "_"), paste("tail", nt, "raw",sep = "_"))
  )
  assign(paste(sample, "tail",nt,"summary",sep = "_"), tmp_2)
  write.xlsx(tail_list, paste0(doc_path, sample,  "_tail_",nt,"_duplicated_summary.xlsx"))
}

# trimming 1-10 summary
for (nt in 1:10) {
  tmp_1 <- tail %>%
    filter(trim_length == nt, tail_base != "N") %>%
    dplyr::arrange(ID, raw_seq) %>%
    left_join(ma1, by = c("ID" = "V1"))%>%
    select(ID = V3, raw_seq, short_mechanism = V4, 7:12)

  tmp_2 <- tmp_1 %>%
    group_by_all() %>%
    count(name = "count") %>%
    ungroup() %>%
    mutate(percentage = round(((count / sum(count)) * 100), 4))

  tmp_3 <- tmp_1 %>%
    group_by(ID) %>%
    tally(name = "count") %>%
    left_join(mirna_count, by=c("ID"="V3")) %>%
    dplyr::rename("SUM" = n) %>%
    mutate(part_percentage = round(((count / sum(count)) * 100), 4),
           global_percentage = round(((count / SUM) * 100), 4))

  trimming_list <- setNames(
    list(tmp_3,tmp_2),
    c(paste("trimming",nt,"summary",sep = "_"), paste("trimming",nt,"detail",sep = "_"))
  )

  assign(paste(sample, "trimming", nt, "summary",sep = "_"), tmp_1)
  write.xlsx(trimming_list, paste0(doc_path, sample, "_trimming_",nt,"_duplicated_summary.xlsx"))
}

tail_summary <- data.frame()
for (i in 0:10) {
  if (i == 0) {
    tmp_0 <- tail %>%
      filter(tail_length == 0) %>%
      group_by(base) %>%
      count(name = "count")
    tmp_0$A <- 0
    tmp_0$T <- 0
    tmp_0$C <- 0
    tmp_0$G <- 0
    tmp_0$tail_N <- 0
    tmp_0$library <- sample
    tmp_0 <- tmp_0[, c(8, 7, 2, 3:6)]
    tail_summary <- rbind(tail_summary, tmp_0)
  } else {
    tmp <- tail %>%
      filter(tail_length == i) %>%
      select(tail_base)
    base <- c("A", "T", "C", "G")
    num_base <- vector(mode = "integer", length = 4)
    for (j in seq_along(base)) {
      num_aux <- str_count(tmp$tail_base, base[j])
      num_base[j] <- sum(num_aux)
    }
    names(num_base) <- base
    tmp_num <- as.data.frame(num_base / i)
    tmp_num <- as.data.frame(t(tmp_num))
    tmp_num$count <- rowSums(tmp_num[, 1:4])
    tmp_num$tail_N <- i
    tmp_num$library <- sample
    tmp_num <- tmp_num[, c(7:5, 1:4)]
    tail_summary <- rbind(tail_summary, tmp_num)
  }
}

# tail summary
tail_summary <- tail_summary %>%
  mutate(
    A_RPM = round((A / total * 10^6), 2),
    T_RPM = round((T / total * 10^6), 2),
    C_RPM = round((C / total * 10^6), 2),
    G_RPM = round((G / total * 10^6), 2),
    percentage = round(((count / sum(count)) * 100), 2),
    A_percentage = round(((A / sum(count)) * 100), 2),
    T_percentage = round(((T / sum(count)) * 100), 2),
    C_percentage = round(((C / sum(count)) * 100), 2),
    G_percentage = round(((G / sum(count)) * 100), 2)
  )

# no-template tail summary
tail_nt_summary <- data.frame()
for (i in 0:10) {
  if (i == 0) {
    tmp_0 <- tail %>%
      filter(tail_length == 0) %>%
      group_by(base) %>%
      count(name = "count")
    tmp_0$A <- 0
    tmp_0$T <- 0
    tmp_0$C <- 0
    tmp_0$G <- 0
    tmp_0$tail_N <- 0
    tmp_0$library <- sample
    tmp_0 <- tmp_0[, c(8, 7, 2, 3:6)]
    tail_nt_summary <- rbind(tail_nt_summary, tmp_0)
  } else {
    tmp <- tail %>%
      filter(base != 0 & tail_length == i) %>%  # That's the difference from up here, when base = 0, tail_length != 0, means this miRNA has been tailing but this base same as genome sequence
      select(tail_base)
    base <- c("A", "T", "C", "G")
    num_base <- vector(mode = "integer", length = 4)
    for (j in seq_along(base)) {
      num_aux <- str_count(tmp$tail_base, base[j])
      num_base[j] <- sum(num_aux)
    }
    names(num_base) <- base
    tmp_num <- as.data.frame(num_base / i)
    tmp_num <- as.data.frame(t(tmp_num))
    tmp_num$count <- rowSums(tmp_num[, 1:4])
    tmp_num$tail_N <- i
    tmp_num$library <- sample
    tmp_num <- tmp_num[, c(7:5, 1:4)]
    tail_nt_summary <- rbind(tail_nt_summary, tmp_num)
  }
}

tail_nt_summary <- tail_nt_summary %>%
  mutate(
    A_RPM = round((A / total * 10^6), 2),
    T_RPM = round((T / total * 10^6), 2),
    C_RPM = round((C / total * 10^6), 2),
    G_RPM = round((G / total * 10^6), 2),
    percentage = round(((count / sum(count)) * 100), 2),
    A_percentage = round(((A / sum(count)) * 100), 2),
    T_percentage = round(((T / sum(count)) * 100), 2),
    C_percentage = round(((C / sum(count)) * 100), 2),
    G_percentage = round(((G / sum(count)) * 100), 2)
  )

# trim summary
trim_summary <- tail %>%
  group_by(trim_length) %>%
  count(name = "count") %>%
  rename(trim_N = trim_length)
trim_summary$RPM <- round((trim_summary$count / total * 10^6), 2)
trim_summary$percentage <- round(((trim_summary$count / sum(trim_summary$count)) * 100), 2)
trim_summary$library <- sample

# output summary file
assign(paste0(sample, "_tail_nt_summary"), tail_nt_summary)
write.xlsx(get(paste0(sample, "_tail_nt_summary")), paste0(doc_path, sample,  "_tail_nt_summary.xlsx"))
assign(paste0(sample, "_tail_summary"), tail_summary)
write.xlsx(get(paste0(sample, "_tail_summary")), paste0(doc_path, sample, "_tail_summary.xlsx"))
assign(paste0(sample, "_trim_summary"), trim_summary)
write.xlsx(get(paste0(sample, "_trim_summary")), paste0(doc_path, sample, "_trim_summary.xlsx"))

# plot
P <- tt_barplot(tail_summary)
ggsave(paste0(plot_path, sample, "_tailing_trimming_length.pdf"), P, width = 8, height = 6)
P_nt <- tt_barplot(tail_nt_summary)
ggsave(paste0(plot_path, sample, "_tailing_trimming_nt_length.pdf"), P_nt, width = 8, height = 6)

#####

save.image(paste0(bin_path, sample, "_", lubridate::today(), ".RData"))
