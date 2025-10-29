# quality check after fastp

pacman::p_load(tidyverse, patchwork, pacman)

summary_to_df = function(x) {
  path = x
  tmp = jsonlite::read_json(path = path) |> pluck("summary")
  list_rbind(tmp[3:4] |> map(as.data.frame), names_to = "step") |>
    mutate(sample_lane = {{x}} |>
             str_remove("^json/") |>
             str_remove("_R1_001.fastq.json"))
}

fastp_data =
  fs::dir_map("json/",
              function(x) summary_to_df(x)) |>
  list_rbind() |>
  separate_wider_delim(sample_lane, names = c("sample", "lane"), delim = "_L00") |>
  mutate(sample = str_remove(sample, "_S[0-9]*")) |>
  group_by(sample, step) |>
  summarise(across(total_reads:total_bases, sum), across(q20_bases:gc_content, mean))

bind_rows(fastp_data, fastp_data, .id = "fastp_run") |>
  filter(step %in% "after_filtering") |>
  ggplot() +
  geom_col(aes(y = sample, x = total_reads, fill = fastp_run),
           position = "dodge")

fastp_data2 |>
  ggplot() +
  geom_col(aes(y = sample, x = total_reads, fill = step),
           position = "dodge")

# quality check after kallisto and kmer size comparison: 17 wins!

# k21 =
#   fs::dir_ls("kall21_2510/", glob = "*.json") |>
#   map(~jsonlite::read_json(.x) |> as.data.frame()) |>
#   list_rbind(names_to = "path") |>
#   mutate(sample = str_remove_all(
#     path, "kall[0-9]*_[0-9]*/|_S[0-9]*.json"))
#
# k23 =
#   fs::dir_ls("kall23_2510/", glob = "*.json") |>
#   map(~jsonlite::read_json(.x) |> as.data.frame()) |>
#   list_rbind(names_to = "path") |>
#   mutate(sample = str_remove_all(
#     path, "kall[0-9]*_[0-9]*/|_S[0-9]*.json"))
#
# k27 =
#   fs::dir_ls("kall27_2510/", glob = "*.json") |>
#   map(~jsonlite::read_json(.x) |> as.data.frame()) |>
#   list_rbind(names_to = "path") |>
#   mutate(sample = str_remove_all(
#     path, "kall[0-9]*_[0-9]*/|_S[0-9]*.json"))

k31 =
  fs::dir_ls("kall31_0811/", glob = "*.json") |>
  map(~jsonlite::read_json(.x) |> as.data.frame()) |>
  list_rbind(names_to = "path") |>
  mutate(sample = str_remove_all(
    path, "kall[0-9]*_[0-9]*/|_S[0-9]*.json"))

old_k31 =
  fs::dir_ls("kall31_2510/", glob = "*.json") |>
  map(~jsonlite::read_json(.x) |> as.data.frame()) |>
  list_rbind(names_to = "path") |>
  mutate(sample = str_remove_all(
    path, "kall[0-9]*_[0-9]*/|_S[0-9]*.json"))

kallisto_data =
  list("k31" = k31, "old_k31" = old_k31) |>
  list_rbind(names_to = "kmer_size")
 #pivot_wider(id_cols = sample, names_from = kmer_size, values_from = where(is.numeric))
  #left_join(fastp_data |>
              # filter(step %in% "after_filtering") |>
              # transmute(sample, q20_rate, q30_rate,
              #        read_length = mean(read1_mean_length, read2_mean_length)))

ggplot(kallisto_data,
       aes(y = p_pseudoaligned,
           x = kmer_size)) +
  geom_point() +
  geom_line(aes(group = sample)) +
  scale_color_gradient(low = "steelblue", high = "coral")

lm(p_pseudoaligned ~ kmer_size, data = kallisto_data) |>
  anova()

rm(fastp_data, k15, k17, k19, k21, kallisto_data, summary_to_df)
