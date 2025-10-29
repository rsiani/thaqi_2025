pacman::p_load(tidyverse, patchwork, pacman)

theme_set(
  theme_bw() +
    theme(
      plot.margin = margin(3, 5, 0, 3),
      text = element_text(size = 30,
                          family = "Fira Sans",
                          color = "#555555"),
      axis.line = element_line(color = "#555555", linewidth = .5),
      axis.ticks = element_line(color = "#555555", linewidth = .5),
      axis.title.x = element_text(hjust = 1),
      axis.title.y = element_text(hjust = 1),
      legend.position = "bottom",
      legend.justification = c(1, 0),
      legend.margin = margin(0, 5, 3, 5),
      legend.title = element_blank(),
      strip.background = element_blank(),
      panel.spacing.x = unit(.5, "lines"),
      panel.spacing.y = unit(.5, "lines"),
    )
)

pal_bac = c("B.velezensis" = "cornflowerblue",
            "B.licheniformis" = "orange")


my_ggsave =
  function(.x, ...) {
    ggsave(filename = str_c("figures/", .x, ".svg", sep = ""),
           width = 11,
           height = 11,
           device = "svg",
           units = "in"
    )
  }



summary_to_df = function(x) {
  path = x
  tmp = jsonlite::read_json(path = path) |> pluck("summary")
  list_rbind(tmp[3:4] |> map(as.data.frame), names_to = "step") |>
    mutate(sample_lane = {{x}} |>
             str_remove("^json/") |>
             str_remove("_R1_001.fastq.json"))
}

quick_shannon =
  function(x) {
    freq = x/sum(x)
    freq = freq[freq > 0]
    y = -sum(freq * log2(freq))
    return(y)
  }

gene_extract =
  \(.x) str_extract(.x, "\\p{alpha}\\p{lower}{2}\\d?\\p{upper}\\b")

plot_svd_nice =
  function(my.data, PCx, PCy){
    ggplot(
      data =
        my.data,
      aes(
        x = .data[[names(my.data)[PCx + 1]]],
        y = .data[[names(my.data)[PCy + 1]]]
      )) +
      geom_point(
        aes(
          shape = time_char,
          color = strain,
        ),
        size = 5,
        stroke = 1,
        alpha = 1/3,
        show.legend = F
      ) +
      scale_shape_manual(values = c(1, 2)) +
      scale_color_manual(
        values = pal_bac,
        aesthetics = c("color")) +
      ggforce::geom_mark_hull(
        aes(fill = group_time,
            color = after_scale(colorspace::darken(fill)),
            # label = group_time
            ),
        linewidth = 0.2,
        alpha = .1,
        show.legend = F,
        concavity = 10) +
      scale_fill_manual(
        values = pal_extended) +
      ggside::geom_xsidedensity(aes(fill = group_time),
                                alpha = 1/2,
                                color = NA,
                                show.legend = T) +
      ggside::geom_ysidedensity(aes(fill = group_time),
                                alpha = 1/2,
                                color = NA,
                                show.legend = T) +
      ggside::scale_xsidey_continuous(labels = NULL) +
      ggside::scale_ysidex_continuous(labels = NULL)
  }


plot_svd_nice2 =
  function(my.data, PCx, PCy){
    ggplot(
      data =
        my.data,
      aes(
        x = .data[[names(my.data)[PCx + 1]]],
        y = .data[[names(my.data)[PCy + 1]]]
      )) +
      scale_color_manual(values = pal_bac) +
      scale_fill_manual(
        values = pal_extended) +
      geom_point(
        aes(
          shape = sample_time,
          color = strain,
        ),
        size = 7,
        stroke = 1,
        alpha = .9,
        show.legend = T
      ) +
      # ggforce::geom_mark_hull(
      #   aes(fill = group_time,
      #       # color = after_scale(colorspace::darken(fill)),
      #       # label = group_time
      #       ),
      #   linewidth = 0,
      #   alpha = .25,
      #   label.margin = margin(.1, .1, .1, .1, "lines"),
      #   label.fontsize = 10,
      #   label.buffer = unit(1, "lines"),
      #   expand = unit(1, "lines"),
      #   show.legend = F) +
      scale_shape_manual(values = c(0, 1, 2,
                                    15, 16, 17))
      # ggside::geom_xsidedensity(aes(fill = group_time),
      #                           alpha = 1/2,
      #                           color = "white",
      #                           show.legend = F) +
      # ggside::geom_ysidedensity(aes(fill = group_time),
      #                           alpha = 1/2,
      #                           color = "white",
      #                           show.legend = F) +
      # ggside::scale_xsidey_continuous(labels = NULL) +
      # ggside::scale_ysidex_continuous(labels = NULL)
  }

qvalue_adj =
  function(.x) {
    separate_wider_delim(.x,
                         contrast,
                         names = c("minuend", "subtrahend"),
                         delim = " - ") |>
      mutate(
        # minuend = str_extract(minuend, "(?<=\\().*(?=\\))"),
        # subtrahend = str_extract(subtrahend, "(?<=\\().*(?=\\))"),
        FDR = p.adjust(p.value, "fdr"),
        significant = if_else(FDR < 0.05, T, F),
        side = case_when(
          statistic > 0 & significant ~ minuend,
          statistic < 0 & significant ~ subtrahend,
          .default = "ns"))
  }

read_interpro =
  function(filename) {
    read_tsv(filename,
             col_names = c("LR", "MD5", "Length", "Analysis",
                           "Signature_accession", "Signature_description",
                           "Start", "End", "score",
                           "Status", "Date", "IP_accession", "IP_description"),
             na = "-") |>
      group_by(LR) |>
      summarise(
        InterPro = str_c(str_unique(paste(IP_accession,
                                          IP_description,
                                          sep = ": ")),
                         collapse = "; ") |>
          str_remove_all("NA: NA(; )?") |>
          na_if(""),
        Full = str_c(str_unique(paste(Signature_accession,
                                      Signature_description,
                                      sep = ": ")),
                     collapse = "; ") |>
          str_remove_all("NA: NA(; )?") |>
          na_if(""))
  }

read_pfam =
  function(filename) {
    read_tsv(filename,
             col_names = c("target_id", "MD5", "Length", "Analysis",
                           "Signature_accession", "Signature_description",
                           "Start", "End", "score",
                           "Status", "Date", "IP_accession", "IP_description"),
             na = "-") |>
      group_by(target_id) |>
      summarise(
        InterPro = str_c(str_unique(IP_accession),
                         collapse = "; ") |>
          str_remove_all("NA(; )?") |>
          na_if(""),
        Pfam = str_c(str_unique(Signature_accession),
                     collapse = "; ") |>
          str_remove_all("NA(; )?") |>
          na_if(""),
        Description = str_c(str_unique(Signature_description, ignore_case = T),
                            collapse = "; ") |>
          str_remove_all("NA(; )?") |>
          na_if(""))
  }

rbh =
  read_tsv("bakta/rbh.tsv",
           col_names = c("bvel", "blic"),
           col_select = 1:2) |>
  full_join(read_tsv("bakta/blic.tsv",
                     col_select = 6,
                     skip = 5),
            by = c("blic" = "Locus Tag"),
            na_matches = "never") |>
  full_join(read_tsv("bakta/bvel.tsv",
                     col_select = 6,
                     skip = 5),
            by = c("bvel" = "Locus Tag"),
            na_matches = "never") |>
  filter(!(is.na(bvel) & is.na(blic))) |>
  mutate(IDX = row_number())

read_fwf("bakta/bvel_annotation.txt") |>
  set_names(c("pgp", "target_id")) |>
  fill(pgp) |>
  left_join(
    read_tsv(str_c("bakta/", "bvel", "_pgp.txt")) |>
      separate_wider_delim(level6,
                           names = c("level6", "pgp"),
                           delim = "->"))

read_pgp =
  function(ID) {
    read_fwf(str_c("bakta/", {{ID}}, "_annotation.txt")) |>
      set_names(c("pgp", "target_id")) |>
      fill(pgp) |>
      left_join(
        read_tsv(str_c("bakta/", {{ID}}, "_pgp.txt")) |>
          separate_wider_delim(level6,
                               names = c("level6", "pgp"),
                               delim = "->"))
  }

read_annotation =
  function(ID) {
    read_tsv(str_c("bakta/", {{ID}}, ".tsv"), skip = 5) |>
      select("target_id" = `Locus Tag`, Gene, Product) |>
      filter(!is.na(target_id)) |>
      left_join(read_tsv(str_c("bakta/", {{ID}}, "_ko_definition.txt"),
                         col_names = c("target_id", "KO", "Description"),
                         col_select = 1:3)) |>
      left_join(read_pgp({{ID}}))
  }


sparsity = \(.x) sum(.x == 0)/n()

# entropy_estimate =
#   rep(1, 300) |>
#   as.list() |>
#   map(.progress = T,
#       \(.x) tidy_LR |>
#         filter(counts > 0) |>
#         uncount(counts) |>
#         slice_sample(n = 100, replace = T, by = sample) |>
#         count(Gene, sample) |>
#         summarise(shannon = quick_shannon(n), .by = sample)) |>
#   list_rbind()

correct_batches =
  function(.x) {
    mat_data =
      pivot_wider(.x,
                  id_cols = Gene,
                  names_from = id,
                  values_from = counts) |>
      column_to_rownames("Gene") |>
      as.matrix()

    metadata_ordered =
      metadata[match(colnames(mat_data), metadata$id),]

    batch_correction =
      sva::ComBat_seq(
        counts = mat_data,
        batch = metadata_ordered$batch,
        group = str_c(metadata_ordered$strain,
                      metadata_ordered$sample)
      )

    tidy_bc =
      batch_correction |>
      as.data.frame() |>
      rownames_to_column("Gene") |>
      pivot_longer(-Gene,
                   values_to = "corrected",
                   names_to = "id")

    return(tidy_bc)

  }





