## rnaseq analysis Stefanie Taqhi

# load essentials

source("scripts/help_fun.R", echo = TRUE)

# LOADING -----------------------------------------------------------------

# load B.vel BAKTA annotation

bvel =
  read_annotation("bvel") |>
  left_join(rbh |> select(-blic), by = join_by(target_id == bvel))
  # left_join(read_tsv("bakta/bvel_clusters.tsv", col_select = 2:3) |>
  #             mutate(target_id = str_remove(protein_id, "STRG0A36PTC.")))


# load B.lic BAKTA annotation

blic =
  read_annotation("blic") |>
  left_join(rbh |> select(-bvel), by =  join_by(target_id == blic))
  # left_join(read_tsv("bakta/blic_clusters.tsv", col_select = 2:3) |>
  #             mutate(target_id = str_remove(protein_id, "STRG0A68UKT.")))


# combine annotations (mmseqs rbh) and fill/forge gene names

background =
  bind_rows(blic, bvel) |>
  mutate(
    Product = na_if(Product,
                    "hypothetical protein") |>
      coalesce(paste0("hp.", ifelse(is.na(IDX), row_number(),
                                    IDX))),
    Gene =
      coalesce(Gene,
               case_when(
                 !is.na(gene_extract(Product)) ~
                   gene_extract(Product),
                 str_detect(Product, "hp.") ~ Product,
                 .default = str_remove_all(Product, "[:punct:]") |>
                   str_to_lower() |>
                   abbreviate())))

# load kallisto transcript counts data and merge everything in one big og' table

# we also compute raw transcriptional probability as a function of transcript
# length and total transcribed area length (we use it later)

fastp_data =
  fs::dir_map("json/",
              function(x) summary_to_df(x)) |>
  list_rbind() |>
  separate_wider_delim(sample_lane, names = c("sample", "lane"), delim = "_L00") |>
  mutate(sample = str_remove(sample, "_S[0-9]*")) |>
  summarise(across(total_reads:total_bases, sum), across(q20_bases:gc_content, mean), .by = c(sample, step))

fastp_data |>
  filter(!sample %in%  c("61", "62", "63", "Undetermined")) |>
  ungroup() |>
  summarise(avg_reads = median(total_reads), .by = step)

fastp_data |>
  ggplot() +
  geom_col(aes(y = sample, x = total_reads, fill = step),
           position = "dodge")

kallisto_data =
  fs::dir_ls("kall31_0811/", glob = "*.json") |>
  map(~jsonlite::read_json(.x) |> as.data.frame()) |>
  list_rbind(names_to = "path") |>
  mutate(id = str_c("id.", str_remove_all(path,
                                 "kall[0-9]*_[0-9]*/|_S[0-9]*.json")),
         path = str_replace(path, ".json", ".tsv"))

raw_data =
  kallisto_data |>
  pull(path) |>
  as.list() |>
  set_names() |>
  map(read_tsv) |>
  list_rbind(names_to = "id") |>
  mutate(id = str_c("id.", str_remove_all(id,
                             "kall.*/|_S[0-9]*.tsv")),
         counts = ceiling(est_counts)) |>
  complete(target_id, id, fill = list(counts = 0)) |>
  left_join(background, multiple = "any") |>
  filter(!str_detect(Product,
                     "ibosom|tRNA|rRNA|transfer-messenger|elongation|initiation") &
           sum(counts) > 1, .by = Gene) |>
  summarise(counts = sum(counts),
            .by = c(Gene, id))

metadata =
  read_tsv("metadata_ST.csv") |>
  mutate(id = str_c("id.", id),
         sample = factor(sample, c("Plim", "BMM", "BC+")),
         timepoint = str_remove(timepoint, " h") |> as.numeric(),
         strain = factor(strain, c("B.velezensis",
                                   "B.licheniformis")),
         media = factor(media, c("Control", "Treatment")),
         batch = str_c("b.", batch)) |>
  left_join(kallisto_data |> select(id, p_pseudoaligned)) |>
  add_row(id = str_c("id.", 61:63)) |>
  mutate(group = str_c(timepoint, strain, sample, sep = "; "))


metadata |>
  count(group, timepoint)

metadata |>
  left_join(
    raw_data |>
      filter(sum(counts > 5) > 5, .by = Gene) |>
      summarise(shannon = quick_shannon(counts),
                lib_size = sum(counts),
                med = median(counts),
                IQR = IQR(counts),
                sparsity = sparsity(counts),
                valid = sum(counts >= 5)/n(),
                .by = id)) |>
  filter(lib_size > 1e4) |>
  # ggplot(aes(x = shannon, y = log(lib_size), color = str_c(timepoint, sample))) +
  # geom_smooth(alpha = .5, method  = "lm", se = F) +
  # geom_point() +
  # ggpp::stat_quadrant_counts(xintercept = 6, yintercept = 10) +
  # ggpp::geom_quadrant_lines(xintercept = 6, yintercept = 10)
  count(strain, sample, timepoint)

# implementation of robust-clr.
# Here we have a zero-offset based on percentage of library space
# This "penalize" 0s from high-intensity samples as it is more likely
# they are true zeros

# with replacement

clean_data =
  raw_data |>
  left_join(metadata) |>
  filter(sum(counts > 5) > 5, .by = Gene) |>
  filter(sum(counts) > 10000, .by = id) |>
  mutate(
    pseudo = case_match(counts,
                        0 ~ NA,
                        .default = counts),
    mean_nz = mean(log(pseudo), na.rm = T),
    pseudo = replace_na(pseudo, 1/exp(1)),
    LR = log(pseudo) - mean_nz,
    .by = id)

# sdv1 = clean_data |>
#   pivot_wider(id_cols = id, names_from = Gene, values_from = pLR) |>
#   column_to_rownames("id") |>
#   prcomp()
#
# sdv2 = clean_data |>
#   pivot_wider(id_cols = id, names_from = Gene, values_from = LR) |>
#   column_to_rownames("id") |>
#   prcomp()
#
# vegan::protest(sdv1, sdv2)

n_distinct(clean_data$Gene)

summarise(clean_data,
          M = max(LR),
          m = min(LR),
          mu = mean(LR),
          .by = id)

clean_data |>
  summarise(spread = sd(LR),
            location = mean(LR),
            .by = c(Gene, strain)) |>
  ggplot(aes(x = location, y = spread)) +
  geom_point(aes(color = strain)) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
  scale_color_manual(values = pal_bac) +
  geom_hline(yintercept = 0, color = "hotpink2", linewidth = .3) +
  geom_vline(xintercept = 0, color = "hotpink2", linewidth = .3) +
  geom_text(data = clean_data |> summarise(sparsity = sum(counts == 0)/n()),
            aes(x = 0, y = 0, label = round(sparsity, 2))) +
  ggpp::stat_dens2d_labels(aes(label = Gene), keep.fraction = 0.01) +
  facet_wrap(~ strain)

meta =
  clean_data |>
  summarise(
    lib_size = sum(counts),
    sparsity = sparsity(counts),
    .by = c(id, sample, timepoint,
            strain, group, media, batch,
            p_pseudoaligned)) |>
  mutate(
    time_char = str_c(timepoint, "hpo") |>
      factor(levels = c("6hpo", "48hpo")),
    bias = lib_size/sum(lib_size),
    sample_time = str_c(timepoint, sample, sep = "; ") |>
      factor(levels = c("6; BMM", "6; Plim", "6; BC+",
                        "48; BMM", "48; Plim", "48; BC+")),
    group = str_c(strain, sample, sep = "; ") |>
      factor(levels = c("B.velezensis; BMM", "B.velezensis; Plim", "B.velezensis; BC+",
                        "B.licheniformis; BMM", "B.licheniformis; Plim", "B.licheniformis; BC+")),
    group_time = str_c(timepoint, strain, sample, sep = "; ") |>
      factor(levels = c("6; B.velezensis; BMM", "6; B.velezensis; Plim",
                        "6; B.velezensis; BC+", "48; B.velezensis; BMM",
                        "48; B.velezensis; Plim", "48; B.velezensis; BC+",
                        "6; B.licheniformis; BMM", "6; B.licheniformis; Plim",
                        "6; B.licheniformis; BC+", "48; B.licheniformis; BMM",
                        "48; B.licheniformis; Plim", "48; B.licheniformis; BC+")))

meta |>
  count(strain, sample, timepoint)


background |>
  summarise(.by = Gene,
            n = n(),
            across(where(is.character),
                   ~ str_c(.x, collapse = "; "))) |> View()

my_data =
  left_join(clean_data |>
              select(-group), meta) |>
  left_join(
    background |>
      summarise(
        gene_copies = n(),
        across(where(is.character),
               ~ str_replace_na(.x) |>
                 str_c(collapse = "; ")),
        .by = Gene),
    relationship = "many-to-many")

ggplot(my_data |>
         mutate(rn = row_number(),
                .by = id) |>
         filter(counts > 0) |>
         mutate(VLR = var(LR), .by = rn) |>
         group_by(rn, timepoint, strain, sample) |>
         summarise(VLR = mean(VLR), LR = mean(LR))) +
  geom_path(aes(x = rn,
                y = LR,
                color = strain,
                alpha = abs(VLR),
                linewidth = abs(VLR))) +
  facet_grid(rows = vars(group)) +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "LR") +
  scale_color_manual(values = pal_bac) +
  scale_linewidth(range = c(0.05, 0.95)) +
  ggh4x::facet_nested(timepoint + sample ~ strain, scales = "free_x")

ggsave("horizon.png", width = 11, height = 5.5)

# clean environment

gc()
