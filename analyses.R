### RNAseq script for Stefanie
## Roberto Siani
# 11.11.23

# SETUP -------------------------------------------------------------------

source("scripts/load_results.R", echo = TRUE)

rm(background, blic, bvel, clean_data, kallisto_data, metadata, raw_data, rbh)

pal_extended <-
  c(
    monochromeR::generate_palette("cornflowerblue", "go_lighter", 6),
    monochromeR::generate_palette("orange", "go_lighter", 6)
  )

# FULL - ANALYSES --------------------------------------------------------

# singular value decomposition

res_svd <-
  my_data |>
  pivot_wider(
    id_cols = id,
    names_from = Gene,
    values_from = LR
  ) |>
  column_to_rownames("id") |>
  as.matrix() |>
  prcomp()

svd_df <-
  broom::tidy(res_svd) |>
  left_join(broom::tidy(res_svd, "d")) |>
  mutate(name = str_c(PC, " ", round(percent, 2))) |>
  pivot_wider(
    id_cols = row,
    names_prefix = "PC"
  ) |>
  left_join(meta, join_by(row == id))

# visualization

a <- plot_svd_nice2(svd_df, 1, 2)
# theme(legend.position = "none")

vegan::adonis2(dist(res_svd$x) ~ sample * strain * timepoint,
  data = meta[match(
    rownames(res_svd$x),
    meta$id
  ), ],
  by = "terms",
  permutations = 999
) |>
  broom::tidy() |>
  gt::gt(caption = "PERMANOVA") |>
  gt::opt_stylize(color = "blue")

res_svd_vel <-
  my_data |>
  filter(strain %in% "B.velezensis" & str_detect(target_id, "bv")) |>
  pivot_wider(
    id_cols = id,
    names_from = Gene,
    values_from = LR
  ) |>
  column_to_rownames("id") |>
  as.matrix() |>
  prcomp()

svd_df_vel <-
  broom::tidy(res_svd_vel) |>
  left_join(broom::tidy(res_svd_vel, "d")) |>
  mutate(name = str_c(PC, " ", round(percent, 2))) |>
  pivot_wider(
    id_cols = row,
    names_prefix = "PC"
  ) |>
  left_join(meta, join_by(row == id)) |>
  mutate(across(where(is.factor), droplevels))


vegan::adonis2(dist(res_svd_vel$x) ~ sample * timepoint,
  data = meta[match(
    rownames(res_svd_vel$x),
    meta$id
  ), ],
  by = "terms",
  permutations = 999
) |>
  broom::tidy() |>
  gt::gt(caption = "B.velezensis") |>
  gt::opt_stylize(color = "blue")


# visualization

b <- plot_svd_nice2(svd_df_vel, 1, 2) +
  theme(legend.position = "none")

res_svd_blic <-
  my_data |>
  filter(strain %in% "B.licheniformis" & str_detect(target_id, "bl")) |>
  pivot_wider(
    id_cols = id,
    names_from = Gene,
    values_from = LR
  ) |>
  column_to_rownames("id") |>
  as.matrix() |>
  prcomp()

svd_df_blic <-
  broom::tidy(res_svd_blic) |>
  left_join(broom::tidy(res_svd_blic, "d")) |>
  mutate(name = str_c(PC, " ", round(percent, 2))) |>
  pivot_wider(
    id_cols = row,
    names_prefix = "PC"
  ) |>
  left_join(meta, join_by(row == id)) |>
  mutate(across(where(is.factor), droplevels))


vegan::adonis2(dist(res_svd_blic$x) ~ sample * timepoint,
  data = meta[match(
    rownames(res_svd_blic$x),
    meta$id
  ), ],
  by = "terms",
  permutations = 999
) |>
  broom::tidy() |>
  gt::gt(caption = "B.licheniformis") |>
  gt::opt_stylize(color = "blue")

# visualization

c <-
  plot_svd_nice2(svd_df_blic, 1, 2) +
  scale_color_manual(values = pal_extended[7:12]) +
  scale_fill_manual(values = pal_extended[7:12]) +
  theme(legend.position = "none")
#
# svd_df_blic |>
#   summarise(across(starts_with("PC"), ~ coin::kruskal_test(.x ~ group_time) |>
#                      coin::pvalue())) |>
#   unlist() |> sort()
#
# c = plot_svd_nice(svd_df_blic, 1, 2) +
#   scale_color_manual(values = pal_extended[4:6]) +
#   scale_fill_manual(values = pal_extended[4:6]) +
#   theme(legend.position = "none")

a + (b / c) +
  plot_layout() +
  plot_annotation(tag_levels = "a")

# export manually cuz patchwork is dumb

# model --------------------------------------------------------------------

# mixed model with flagella as predictor and mediator of response to different
# media type. Strain are used as grouping variable and we account for
# unequal variances and differences in group sizes

res_glm <-
  my_data |>
  nest(.by = c(
    Gene,
    target_id,
    Product,
    KO,
    Description,
    strain,
    pgp,
    starts_with("level")
  )) |>
  mutate(
    mod = map(data,
      .progress = T,
      \(.x) glm(LR ~ sample * time_char,
        contrasts = list(sample = "contr.sum"),
        data = .x
      )
    ),
    augmented = map(mod, broom::augment),
    tidied = map(mod, broom::tidy)
  )

res_glm |>
  pluck("mod") |>
  map(residuals) |>
  unlist() |>
  plot()

# conditional predictions and comparisons for different focal variables:
# effects due to flagella
# differences in effects due to flagella in the differents media
# differences in effects due to media between flagellated and unflagellated

comp_by_media <-
  res_glm |>
  mutate(
    a_comp =
      map(mod,
        .progress = T,
        \(.x) marginaleffects::avg_comparisons(.x,
          by = "time_char",
          # by = c("strain", "time_char"),
          variables = list(sample = "pairwise")
        )
      )
  )

# comp_by_strain =
#   res_glm |>
#   mutate(
#     a_comp =
#       map(mod,
#           .progress = T,
#           \(.x) marginaleffects::avg_comparisons(.x,
#                                                  variables = "strain",
#                                                  by = c("sample", "time_char"))))

comp_by_time <-
  res_glm |>
  mutate(
    a_comp =
      map(mod,
        .progress = T,
        \(.x) marginaleffects::avg_comparisons(.x,
          variables = "time_char",
          # by = c("strain", "sample"),
          by = "sample"
        )
      )
  )

all_res <- list(
  "general" = res_glm,
  "media" = comp_by_media,
  # "strain" = comp_by_strain,
  "time" = comp_by_time
)

saveRDS(all_res, "all_results_separated_by_strain.RDS")

all_res <- readRDS("all_results_separated_by_strain.RDS")

# effects -----------------------------------------------------------------

tidied_results <-
  all_res |>
  pluck("general") |>
  filter((strain %in% "B.velezensis" & str_detect(target_id, "bv")) |
    (strain %in% "B.licheniformis" & str_detect(target_id, "bl"))) |>
  unnest(tidied) |>
  mutate(FDR = p.adjust(p.value, "fdr"), .by = strain)

with(tidied_results, table(term, FDR < 0.05))

# calculate adjusted qvalues, local false discovery rates and significance

res_comp_by_time <-
  all_res |>
  pluck("time") |>
  filter((str_detect(target_id, "bl") & strain %in% "B.licheniformis") |
    (str_detect(target_id, "bv") & strain %in% "B.velezensis")) |>
  unnest(a_comp) |>
  qvalue_adj()

res_comp_by_time |>
  filter(significant) |>
  count(strain, sample, side)

res_comp_by_media <-
  all_res |>
  pluck("media") |>
  unnest(a_comp) |>
  filter((str_detect(target_id, "bl") & strain %in% "B.licheniformis") |
    (str_detect(target_id, "bv") & strain %in% "B.velezensis")) |>
  qvalue_adj() |>
  mutate(contrast = str_c(minuend, subtrahend, sep = "_"))

res_comp_by_media |>
  filter(significant) |>
  count(time_char, strain, side)

res_comp_by_media |>
  filter(significant) |>
  mutate(group = str_c(time_char, strain, side)) |>
  summarise(group = list(group), .by = Gene) |>
  ggplot(aes(x = group)) +
  geom_bar() +
  ggupset::scale_x_upset(n_intersections = 15)

# res_comp_by_strain =
#   all_res |>
#   pluck("strain") |>
#   unnest(a_comp) |>
#   qvalue_adj()
#
# res_comp_by_strain |>
#   filter(significant) |>
#   count(time_char, sample, side) |>
#   gt::gt(caption = "FDR < 0.05") |>
#   gt::opt_stylize(color = "pink")

pal_sample <- c(
  "ns" = "#c5c5c5",
  "Plim" = "#c7e9b4",
  "BC+" = "#41b6c4",
  "BMM" = "#225ea8"
)

int_genes <-
  c(
    "tatB", "arsR", "arsC", "tagD", "phoH", "GntR", "ArsR", "batp",
    "ldh", "lldP", "ugpQ", "pstA", "pstB", "pstC", "phoA", "nptA",
    "arK", "narI", "narJ", "narH", "narK", "spoVB", "flgI", "flgK",
    "flgM", "flgC", "fliW", "fliJ", "fliM", "flgP", "ugpE", "caiA",
    "aceB", "cydB", "cydC", "cydD", "nirD", "ugpB", "citT", "gabT",
    "gabD", "putP", "this", "thiF", "pit", "buk", "sulP", "spxH",
    "gbsB"
  )


res_comp_by_media |>
  ggplot(aes(x = estimate, y = s.value, color = side)) +
  geom_point() +
  ggrepel::geom_text_repel(
    data = ~ filter(.x, Gene %in% int_genes, significant, estimate > 0),
    aes(label = Gene), color = "#555555", size = 5,
    min.segment.length = 0,
    fontface = "italic", family = "Fira Sans",
    max.overalps = Inf, box.padding = 0.3,
    max.time = 1, max.iter = 1e4,
    point.padding = 0,
    force_pull = 0,
    hjust = 0,
    nudge_x = 5 - (filter(res_comp_by_media, Gene %in% int_genes, significant, estimate > 0) |>
      pluck("estimate")),
    direction = "y"
  ) +
  ggrepel::geom_text_repel(
    data = ~ filter(.x, Gene %in% int_genes, significant, estimate < 0),
    aes(label = Gene), color = "#555555", size = 5,
    min.segment.length = 0,
    hjust = 1,
    fontface = "italic", family = "Fira Sans",
    max.overalps = Inf, box.padding = 0.3,
    max.time = 1, max.iter = 1e4,
    point.padding = 0,
    force_pull = 0,
    nudge_x = -5 - (filter(res_comp_by_media, Gene %in% int_genes, significant, estimate < 0) |>
      pluck("estimate")),
    direction = "y"
  ) +
  xlim(-6, 6) +
  ggh4x::facet_nested_wrap(~ strain + contrast + time_char, nrow = 2) +
  scale_color_manual(values = pal_sample) +
  ggpp::stat_group_counts()

my_ggsave("new_res/volcano2", width = 36, height = 24)

ggsave("figures/new_res/volcano2.svg", width = 36, height = 24)

write_csv(res_comp_by_media, "figures/new_res/res_comp_by_media.csv")
# write_csv(res_comp_by_strain, "res_comp_by_strain.csv")
write_csv(res_comp_by_time, "figures/new_res/res_comp_by_time.csv")

# res_comp_by_strain |>
#   ggplot(aes(x = estimate, y = s.value, color = side)) +
#   geom_point(shape = 19) +
#   ggh4x::facet_nested_wrap(~ sample + minuend + time_char, nrow = 1) +
#   scale_color_manual(values = c(pal_bac, "#757575")) +
#   ggpp::stat_group_counts()
#
# my_ggsave("volcano_strain_vs")

stack_venn <-
  map2(
    c(
      "Plim", "Plim",
      "BC+", "BC+",
      "BMM", "BMM"
    ),
    c(
      "6hpo", "48hpo",
      "6hpo", "48hpo",
      "6hpo", "48hpo"
    ),
    ~ ggvenn::ggvenn(
      list(
        B.velezensis = res_comp_by_media |>
          filter(side %in% {{ .x }} &
            time_char %in% {{ .y }} &
            strain %in% c("B.velezensis")) %>% pull(Gene),
        B.licheniformis = res_comp_by_media |>
          filter(side %in% {{ .x }} &
            time_char %in% {{ .y }} &
            strain %in% c("B.licheniformis")) |> pull(Gene)
      ),
      fill_color = c("cornflowerblue", "orange"),
      fill_alpha = .3,
      text_size = 10,
      set_name_size = 15,
      show_percentage = F,
      stroke_color = "#333333",
      set_name_color = "#333333",
      text_color = "#333333"
    ) + ggtitle(str_c({{ .x }}, {{ .y }}, sep = ": ")) +
      theme(title = element_text(size = 20))
  )

(stack_venn[[1]] + stack_venn[[2]]) /
  (stack_venn[[3]] + stack_venn[[4]]) /
  (stack_venn[[5]] + stack_venn[[6]])

# stack_venn2 =
#   map2(c("Plim", "Plim",
#          "BC+", "BC+" ,
#          "BMM", "BMM"),
#        c("6hpo", "48hpo",
#          "6hpo", "48hpo",
#          "6hpo", "48hpo"),
#        ~ ggvenn::ggvenn(
#          list(
#            B.velezensis = res_comp_by_time |>
#              filter(sample %in% {{ .x }} &
#                       side %in% {{ .y}} &
#                       strain %in% c("B.velezensis")) %>% pull(Gene),
#            B.licheniformis = res_comp_by_time |>
#              filter(sample %in% {{ .x }} &
#                       side %in% {{ .y}} &
#                       strain %in% c("B.licheniformis")) |> pull(Gene)),
#          fill_color = c("cornflowerblue", "orange"),
#          fill_alpha = .3,
#          text_size = 10,
#          set_name_size = 15,
#          show_percentage = F,
#          stroke_color = "#333333",
#          set_name_color = "#333333",
#          text_color = "#333333"
#        ) + ggtitle(str_c({{ .y }}, {{ .x }}, sep = ": ")) +
#          theme(title = element_text(size = 20)))
#
# ((stack_venn2[[1]] + stack_venn2[[2]]) /
#   (stack_venn2[[3]] + stack_venn2[[4]]) /
#   (stack_venn2[[5]] + stack_venn2[[6]]) )|
#   ((stack_venn[[1]] + stack_venn[[2]]) /
#   (stack_venn[[3]] + stack_venn[[4]]) /
#   (stack_venn[[5]] + stack_venn[[6]]))

my_ggsave("new_res/venn")

# enrichment --------------------------------------------------------------

# nothing to see here

enricher <- clusterProfiler::enricher

data(kegg_category, package = "clusterProfiler")

T2N <-
  kegg_category |>
  mutate(KO = str_c("K", id)) |>
  select(KO, subcategory)

T2G <-
  kegg_category |>
  mutate(KO = str_c("K", id)) |>
  select(category, KO)

longer_media <-
  res_comp_by_media |>
  separate_longer_delim(
    cols = c(pgp, starts_with("level")),
    delim = "; "
  ) |>
  filter(pgp != "NA") |>
  pivot_longer(c(level2:level5),
    names_to = "level",
    values_to = "category"
  )



longer_time |>
  pull(level) |>
  unique() |>
  set_names() |>
  map(
    ~ clusterProfiler::compareCluster(
      Gene ~ side + strain + time_char,
      data = res_comp_by_media |>
        filter(significant),
      fun = "enricher",
      TERM2GENE = longer_media |>
        filter(level %in% .x) |>
        select(category, Gene) |>
        distinct(),
      pvalueCutoff = 1
    ) |> pluck("compareClusterResult")
  ) |>
  list_rbind(names_to = "level") |>
  View()

ce <-
  clusterProfiler::compareCluster(
    Gene ~ side + strain + sample,
    data = longer_time |>
      filter(significant),
    fun = "enricher",
    universe = my_data |> filter(str_detect(pgp, "P")) |> pull(Gene) |> unique(),
    TERM2GENE = longer_time |> select(pgp, Gene),
    TERM2NAME = longer_time |> select(pgp, level),
    pvalueCutoff = 1,
    qvalueCutoff = 0.05
  )

ce |>
  pluck("compareClusterResult") |>
  View()

## creating output for string

dir.create("enrichement_test_string")

# by media

dir.create("enrichement_test_string/by_media")

res_comp_by_media |>
  filter(significant) |>
  select(strain, side, time_char, contrast) |>
  count(strain, side, time_char, contrast)

pwalk(
  res_comp_by_media |>
    filter(significant) |>
    select(strain, side, time_char, contrast) |>
    distinct(),
  \(side, strain, time_char, contrast) res_comp_by_media |>
    filter(
      strain %in% {{ strain }},
      side %in% {{ side }},
      time_char %in% {{ time_char }},
      contrast %in% {{ contrast }},
      significant
    ) |>
    select(target_id) |>
    separate_rows(target_id, sep = "; ") |>
    filter(str_detect(
      target_id,
      if_else(strain %in% "B.licheniformis",
        "bl",
        "bv"
      )
    )) |>
    distinct() |>
    write_csv(
      str_c(
        "enrichement_test_string/by_media/",
        strain, "_",
        side, "_",
        contrast, "_",
        time_char, ".csv"
      ),
      col_names = F
    )
)

# by strain

dir.create("enrichement_test_string/by_strain")

res_comp_by_strain |>
  filter(significant) |>
  select(sample, side, time_char) |>
  count(side, sample, time_char)

pwalk(
  res_comp_by_strain |>
    filter(significant) |>
    select(sample, side, time_char) |>
    distinct(),
  \(sample, side, time_char) res_comp_by_strain |>
    filter(
      sample %in% {{ sample }},
      side %in% {{ side }},
      time_char %in% {{ time_char }},
      significant
    ) |>
    select(target_id) |>
    separate_rows(target_id, sep = "; ") |>
    filter(str_detect(
      target_id,
      if_else(side %in% "B.licheniformis",
        "bl",
        "bv"
      )
    )) |>
    distinct() |>
    write_csv(
      str_c(
        "enrichement_test_string/by_strain/",
        side, "_", sample, "_", time_char, ".csv"
      ),
      col_names = F
    )
)

# by time

dir.create("enrichement_test_string/by_time")

res_comp_by_time |>
  filter(significant) |>
  select(strain, sample, side) |>
  count(strain, sample, side)

pwalk(
  res_comp_by_time |>
    filter(significant) |>
    select(strain, sample, side) |>
    distinct(),
  \(strain, sample, side) res_comp_by_time |>
    filter(
      sample %in% {{ sample }},
      side %in% {{ side }},
      strain %in% {{ strain }},
      significant
    ) |>
    select(target_id) |>
    separate_rows(target_id, sep = "; ") |>
    filter(str_detect(
      target_id,
      if_else(strain %in% "B.licheniformis",
        "bl",
        "bv"
      )
    )) |>
    distinct() |>
    write_csv(
      str_c(
        "enrichement_test_string/by_time/",
        strain, "_", side, "_", sample, ".csv"
      ),
      col_names = F
    )
)

# viz2 --------------------------------------------------------------------



my_data |>
  separate_longer_delim(
    cols = c(pgp, starts_with("level")),
    delim = "; "
  ) |>
  filter(
    Gene %in% c(
      pull(
        res_comp_by_media |>
          filter(significant),
        Gene
      ),
      pull(
        res_comp_by_time |>
          filter(significant),
        Gene
      )
    ),
    str_detect(level3, "PHOSPHATE_SOLUBILIZATION")
  ) |>
  ggplot() +
  geom_raster(aes(x = id, y = Gene, fill = LR)) +
  ggh4x::facet_nested(level4 ~ strain + sample + time_char, scales = "free", space = "free") +
  rcartocolor::scale_fill_carto_c(palette = "Tropic") +
  theme(
    axis.text.x = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.ticks = element_blank(),
    strip.clip = "off"
  )

ggsave("heatmap.png", width = 24, height = 30)

target_list_by_media <-
  readxl::read_xlsx("Interesting_Genes.xlsx", 3) |>
  pull(Gene) |>
  unique()
target_list_by_time <-
  readxl::read_xlsx("Interesting_Genes.xlsx", 4) |>
  pull(Gene) |>
  unique()

target_list <-
  readxl::read_xlsx("Interesting_Genes.xlsx", 2) |>
  separate_longer_delim(cols = Gene, delim = ", ") |>
  right_join(my_data) |>
  drop_na(category)

target_list |>
  ggplot() +
  geom_raster(aes(x = id, y = Gene, fill = LR)) +
  ggh4x::facet_nested(category ~ strain + sample + time_char, scales = "free", space = "free") +
  rcartocolor::scale_fill_carto_c(palette = "Tropic") +
  theme(
    axis.text = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.ticks = element_blank(),
    strip.clip = "off"
  )

ggsave("target_genes.png", width = 24, height = 14)

target <-
  list(
    "pho" = c("phoD", "phoE", "phoA", "phoH"),
    "gln" = c("glnK", "glnK", "glnA", "glnQ", "glnR", "glnA"),
    "amt" = c("amtB"),
    "ure" = c("ureC", "5_ureB_sRNA", "ureB", "ureA"),
    "pts" = c("pstB", "pstA", "pstC", "pstS"),
    "nir" = c("nirB", "nirD"),
    "proton" = c(
      "fadH", "fmnn", "azoR1", "azoR2", "ndhd", "nuoL", "yjlD",
      "ndh", "ndfo", "ptnd", "yugJ", "yugK", "ptno", "ywnB",
      "wrbA", "ndpfo", "ndphdo", "yfmJ", "qor", "yjgC", "afma",
      "malK", "ugpC", "tadA"
    )
  )

plot_stack <-
  target |>
  imap(~ my_data |>
    filter(Gene %in% .x) |>
    ggplot(aes(
      y = group_time,
      x = LR,
      color = group_time
    )) +
    geom_jitter(width = 0, alpha = 0.5) +
    geom_boxplot(fill = NA) +
    scale_color_manual(values = pal_extended) +
    ggtitle(label = .y) +
    theme(legend.position = "none"))

plot_stack |>
  imap(~ ggsave(str_c(.y, ".png"), plot = .x, ))

res_comp_by_strain |>
  filter(!is.na(KO)) |>
  left_join(kegg_category |> mutate(KO = str_c("K", id))) |>
  mutate(
    category =
      category |>
        as.factor() |>
        fct_lump_n(5, ties.method = "average") |>
        fct_na_value_to_level("Other") |>
        fct_infreq(), .by = sample
  ) |>
  ggplot() +
  geom_bar(
    aes(
      y = category,
      fill = side
    ),
    position = "dodge"
  ) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    panel.spacing.y = unit(0, "lines"),
    panel.spacing.x = unit(0, "lines"),
    strip.text.y = element_text(angle = 0, hjust = 0)
  ) +
  scale_fill_manual(values = pal_bac) +
  scale_y_discrete(labels = \(.x) str_wrap(.x, width = 15)) +
  scale_x_log10() +
  ggh4x::facet_nested_wrap(~ sample + time_char, scales = "free", nrow = 3)

my_ggsave("enrichment_strain")

res_comp_by_media |> # enrichment of drugs and metabolism in
  filter(significant) |>
  left_join(kegg_category |> mutate(KO = str_c("K", id))) |>
  mutate(
    subcategory =
      subcategory |>
        as.factor() |>
        fct_lump_n(5) |>
        fct_na_value_to_level("Other") |>
        fct_infreq(), .by = c(strain, contrast)
  ) |>
  ggplot() +
  geom_bar(
    aes(
      y = subcategory,
      fill = side
    ),
    position = "dodge"
  ) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    panel.spacing.y = unit(0, "lines"),
    panel.spacing.x = unit(0, "lines"),
    strip.text.y = element_text(angle = 0, hjust = 0)
  ) +
  scale_fill_manual(values = pal_sample) +
  scale_y_discrete(labels = \(.x) str_wrap(.x, width = 15)) +
  ggh4x::facet_nested_wrap(~ strain + contrast + time_char,
    scales = "free",
    nrow = 2
  )

my_ggsave("enrichment_media")


# network -----------------------------------------------------------------

p_load(propr, tidygraph, ggraph)

meta |>
  count(strain, time_char)

# nested_propr =
#   my_data |>
#   filter((str_detect(target_id, "bl") & strain %in% "B.licheniformis") |
#            (str_detect(target_id, "bv") & strain %in% "B.velezensis")) |>
#   nest(.by  = c(strain, sample)) |>
#   mutate(
#     wide_data =
#       map(data,
#           ~ pivot_wider(.x,
#                         id_cols = id,
#                         values_from = counts,
#                         names_from = Gene) |>
#             column_to_rownames("id")),
#     res_propd =
#       map(wide_data,
#           ~ propd(
#             .x,
#             group = meta[match(rownames(.x),
#                                meta$id),]$time_char,
#             alpha = 0.01,
#             p = 99) |>
#             setActive(what = "theta_d") |>
#             updateCutoffs(cutoff = seq(0.05, 0.95, 0.05))))


nested_propr <-
  my_data |>
  filter((str_detect(target_id, "bl") & strain %in% "B.licheniformis") |
    (str_detect(target_id, "bv") & strain %in% "B.velezensis")) |>
  nest(.by = c(strain, time_char)) |>
  mutate(
    wide_data =
      map(
        data,
        ~ filter(.x, all(counts > 0), .by = Gene) |>
          pivot_wider(
            id_cols = id,
            values_from = LR,
            names_from = Gene
          ) |>
          column_to_rownames("id")
      ),
    res_propr =
      map(
        wide_data,
        ~ propr(
          .x,
          ivar = NA,
          metric = "rho",
          p = 99
        ) |>
          updateCutoffs(
            cutoff = seq(0.6, 0.95, 0.05),
            ncores = 6
          )
      )
  )

nested_propr[, 1]

nested_propr |>
  pluck("res_propr", 2, "fdr")

net_licheniformis_first <-
  nested_propr |>
  pluck("res_propr", 1) |>
  getResults() |>
  filter(propr > .90) |>
  select(
    from = Pair,
    to = Partner,
    weights = propr
  ) |>
  as_tbl_graph(directed = F) |>
  mutate(
    fast_greedy = group_fast_greedy(weights = weights),
    color = group_color(),
    degree = centrality_degree(weights = weights),
    hub = centrality_hub(weights = weights)
  ) |>
  activate("nodes") |>
  filter(degree > 1)

ggraph(net_licheniformis_first, layout = "fr") +
  geom_edge_link(aes(alpha = weights)) +
  geom_node_point(
    aes(
      size = hub,
      color = as.factor(fast_greedy)
    ),
    alpha = 1 / 3
  ) +
  theme_graph(plot_margin = margin(0, 0, 0, 0)) +
  theme(legend.position = "none") +
  scale_alpha(range = c(0.1, 0.7))

net_velezensis <-
  nested_propr |>
  pluck("res_propr", 2) |>
  getResults() |>
  select(
    from = Pair,
    to = Partner,
    weights = propr
  ) |>
  as_tbl_graph(directed = F) |>
  mutate(group = group_fast_greedy(weights = weights)) |>
  activate("nodes") |>
  as_tibble()


my_data |>
  filter(strain %in% "B.velezensis") |>
  left_join("")


wide_data <-
  pivot_wider(my_data,
    id_cols = Gene,
    values_from = LR,
    names_from = id
  ) |>
  column_to_rownames("Gene")

res_tsne <-
  Rtsne::Rtsne(
    wide_data,
    perplexity = 90,
    max_iter = 1e3,
    check_duplicates = F,
    initial_dims = 30
  )

data.frame(
  Gene = rownames(wide_data),
  D = res_tsne$Y
) |>
  left_join(res_comp_by_media) |>
  ggplot() +
  geom_point(aes(D.1, D.2,
    alpha = abs(estimate),
    color = side
  )) +
  ggh4x::facet_nested_wrap(~ time_char + strain) +
  scale_color_manual(values = pal_sample)

data.frame(
  Gene = rownames(wide_data),
  D = res_tsne$Y
) |>
  left_join(res_comp_by_strain) |>
  ggplot() +
  geom_point(aes(D.1, D.2,
    alpha = abs(estimate),
    color = side
  )) +
  ggh4x::facet_nested_wrap(~ time_char + sample) +
  scale_color_manual(values = pal_bac)

data.frame(
  Gene = rownames(wide_data),
  D = res_tsne$Y
) |>
  left_join(res_comp_by_time) |>
  ggplot() +
  geom_point(aes(D.1, D.2,
    alpha = abs(estimate),
    color = side
  )) +
  ggh4x::facet_nested_wrap(~ strain + sample)


wide_bv <-
  my_data |>
  filter(str_detect(target_id, "bv") & strain %in% "B.velezensis") |>
  pivot_wider(
    id_cols = Gene,
    values_from = LR,
    names_from = id
  ) |>
  column_to_rownames("Gene")

res_tsne_vel <-
  Rtsne::Rtsne(
    wide_bv,
    perplexity = 30,
    max_iter = 1e3,
    check_duplicates = F,
    initial_dims = 10
  )

data.frame(
  Gene = rownames(wide_bv),
  D = res_tsne_vel$Y
) |>
  left_join(res_comp_by_time |>
    filter(strain %in% "B.velezensis")) |>
  ggplot() +
  geom_point(aes(D.1, D.2,
    alpha = abs(estimate),
    color = side
  )) +
  ggh4x::facet_nested_wrap(~ strain + sample)

data.frame(
  Gene = rownames(wide_bv),
  D = res_tsne_vel$Y
) |>
  left_join(res_comp_by_media |>
    filter(strain %in% "B.velezensis")) |>
  ggplot() +
  geom_point(aes(D.1, D.2,
    alpha = abs(estimate),
    color = side
  )) +
  ggh4x::facet_nested_wrap(~ strain + time_char) +
  scale_color_manual(values = pal_sample)

wide_bl <-
  my_data |>
  filter(str_detect(target_id, "bl") & strain %in% "B.licheniformis") |>
  pivot_wider(
    id_cols = Gene,
    values_from = LR,
    names_from = id
  ) |>
  column_to_rownames("Gene")

res_tsne_lic <-
  Rtsne::Rtsne(
    wide_bl,
    perplexity = 30,
    max_iter = 1e3,
    check_duplicates = F,
    initial_dims = 10
  )

data.frame(
  Gene = rownames(wide_bl),
  D = res_tsne_lic$Y
) |>
  left_join(res_comp_by_time |>
    filter(strain %in% "B.licheniformis")) |>
  ggplot() +
  geom_point(aes(D.1, D.2,
    alpha = abs(estimate),
    color = side
  )) +
  ggh4x::facet_nested_wrap(~ strain + sample)

data.frame(
  Gene = rownames(wide_bl),
  D = res_tsne_lic$Y
) |>
  left_join(res_comp_by_media |>
    filter(strain %in% "B.licheniformis")) |>
  ggplot() +
  geom_point(aes(D.1, D.2,
    alpha = abs(estimate),
    color = side
  )) +
  ggh4x::facet_nested_wrap(~ strain + time_char) +
  scale_color_manual(values = pal_sample)
