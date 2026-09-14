library(data.table)
library(grid)

options(bitmapType = "cairo")

# ============================================================
# 1. Parse command-line arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- which(args == flag)

  if (length(idx) == 0 || idx == length(args)) {
    return(NA_character_)
  }

  args[idx + 1]
}

mode       <- get_arg("--mode")
base_dir   <- get_arg("--base")
trait      <- get_arg("--trait")
out_dir    <- get_arg("--outdir")
trait_list <- get_arg("--trait-list")


# ============================================================
# 2. Shared functions
# ============================================================

get_element_name <- function(trait_name) {
  sub(
    "_(log|no_outlier)$",
    "",
    trait_name
  )
}


find_column <- function(dat, possible_names) {

  matched <- possible_names[
    possible_names %in% names(dat)
  ]

  if (length(matched) == 0) {
    stop(
      "Required column not found. Tried: ",
      paste(possible_names, collapse = ", "),
      "\nAvailable columns: ",
      paste(names(dat), collapse = ", ")
    )
  }

  matched[1]
}


read_mlma <- function(mlma_file) {

  if (!file.exists(mlma_file)) {
    stop("MLMA file not found: ", mlma_file)
  }

  cat("Reading:\n")
  cat(mlma_file, "\n")

  df <- fread(
    mlma_file,
    showProgress = TRUE
  )

  chr_col <- find_column(df, c("Chr", "CHR", "chr"))
  snp_col <- find_column(df, c("SNP", "snp"))
  bp_col  <- find_column(df, c("bp", "BP", "pos", "POS"))
  p_col   <- find_column(df, c("p", "P"))

  gwas <- data.table(
    SNP = as.character(df[[snp_col]]),
    CHR = suppressWarnings(as.integer(df[[chr_col]])),
    BP  = suppressWarnings(as.numeric(df[[bp_col]])),
    P   = suppressWarnings(as.numeric(df[[p_col]]))
  )

  rm(df)
  gc()

  gwas <- gwas[
    !is.na(CHR) &
      !is.na(BP) &
      !is.na(P) &
      CHR >= 1 &
      CHR <= 22 &
      BP >= 0 &
      P > 0 &
      P <= 1
  ]

  if (nrow(gwas) == 0) {
    stop("No valid GWAS records retained from: ", mlma_file)
  }

  setorder(gwas, CHR, BP)

  cat("Variants retained:", nrow(gwas), "\n")
  cat(
    "Chromosomes:",
    paste(sort(unique(gwas$CHR)), collapse = ", "),
    "\n"
  )

  gwas
}


prepare_manhattan <- function(gwas) {

  chr_info <- gwas[
    ,
    .(
      chr_length = max(BP, na.rm = TRUE)
    ),
    by = CHR
  ]

  setorder(chr_info, CHR)

  chr_info[
    ,
    offset := shift(
      cumsum(chr_length),
      fill = 0
    )
  ]

  plot_data <- merge(
    gwas,
    chr_info[, .(CHR, offset)],
    by = "CHR",
    all.x = TRUE,
    sort = FALSE
  )

  plot_data[
    ,
    BP_cumulative := BP + offset
  ]

  plot_data[
    ,
    LOGP := -log10(
      pmax(P, .Machine$double.xmin)
    )
  ]

  setorder(plot_data, CHR, BP)

  axis_data <- plot_data[
    ,
    .(
      axis_position = (
        min(BP_cumulative, na.rm = TRUE) +
          max(BP_cumulative, na.rm = TRUE)
      ) / 2
    ),
    by = CHR
  ]

  setorder(axis_data, CHR)

  list(
    data = plot_data,
    axis = axis_data
  )
}


draw_manhattan <- function(gwas, title_text) {

  prepared <- prepare_manhattan(gwas)

  plot_data <- prepared$data
  axis_data <- prepared$axis

  point_colours <- ifelse(
    plot_data$CHR %% 2 == 0,
    "grey65",
    "grey25"
  )

  point_colours[
    plot_data$P < 5e-8
  ] <- "red"

  ymax <- max(
    plot_data$LOGP,
    na.rm = TRUE
  )

  if (!is.finite(ymax)) {
    ymax <- 10
  }

  ymax <- max(
    9,
    ceiling(ymax * 1.05)
  )

  plot(
    x = plot_data$BP_cumulative,
    y = plot_data$LOGP,
    type = "n",
    xaxt = "n",
    xlab = "Chromosome",
    ylab = expression(-log[10](italic(P))),
    main = title_text,
    ylim = c(0, ymax),
    bty = "l",
    las = 1,
    cex.main = 1.5,
    cex.lab = 1.2,
    cex.axis = 1
  )

  points(
    x = plot_data$BP_cumulative,
    y = plot_data$LOGP,
    pch = 20,
    cex = 0.28,
    col = point_colours
  )

  # Genome-wide significance threshold
  abline(
    h = -log10(5e-8),
    col = "red",
    lwd = 1.5,
    lty = 2
  )

  # Suggestive significance threshold
  abline(
    h = -log10(1e-6),
    col = "blue",
    lwd = 1.3,
    lty = 2
  )

  axis(
    side = 1,
    at = axis_data$axis_position,
    labels = axis_data$CHR,
    tick = FALSE,
    line = -0.5,
    cex.axis = 0.9
  )
}


draw_qq <- function(p_values, title_text) {

  p_values <- p_values[
    is.finite(p_values) &
      p_values > 0 &
      p_values <= 1
  ]

  p_values <- sort(p_values)

  n <- length(p_values)

  if (n == 0) {
    stop("No valid P values for QQ plot.")
  }

  expected <- -log10(ppoints(n))

  observed <- -log10(
    pmax(p_values, .Machine$double.xmin)
  )

  plot_limit <- max(
    expected,
    observed,
    na.rm = TRUE
  )

  plot(
    expected,
    observed,
    pch = 20,
    cex = 0.35,
    col = "black",
    xlab = expression(Expected~~-log[10](italic(P))),
    ylab = expression(Observed~~-log[10](italic(P))),
    main = title_text,
    xlim = c(0, plot_limit),
    ylim = c(0, plot_limit),
    asp = 1,
    bty = "l",
    las = 1,
    cex.main = 1.5,
    cex.lab = 1.2,
    cex.axis = 1
  )

  abline(
    a = 0,
    b = 1,
    col = "red",
    lwd = 1.5,
    lty = 2
  )
}


draw_png_pdf_page <- function(
  manhattan_file,
  qq_file
) {

  if (!requireNamespace("png", quietly = TRUE)) {
    stop(
      "The R package 'png' is required for PDF assembly."
    )
  }

  if (!file.exists(manhattan_file)) {
    stop("Manhattan PNG missing: ", manhattan_file)
  }

  if (!file.exists(qq_file)) {
    stop("QQ PNG missing: ", qq_file)
  }

  manhattan_image <- png::readPNG(
    manhattan_file
  )

  qq_image <- png::readPNG(
    qq_file
  )

  grid.newpage()

  # Manhattan plot on upper part
  pushViewport(
    viewport(
      x = 0.5,
      y = 0.68,
      width = 0.97,
      height = 0.60
    )
  )

  grid.raster(
    manhattan_image,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = TRUE
  )

  popViewport()

  # QQ plot on lower part
  pushViewport(
    viewport(
      x = 0.5,
      y = 0.23,
      width = 0.42,
      height = 0.42
    )
  )

  grid.raster(
    qq_image,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = TRUE
  )

  popViewport()
}


# ============================================================
# 3. Single-trait plotting mode
# ============================================================

if (mode == "single") {

  if (any(is.na(c(base_dir, trait, out_dir)))) {
    stop(
      "Single mode usage:\n",
      "Rscript script.R ",
      "--mode single ",
      "--base <GWAS_BASE> ",
      "--trait <TRAIT> ",
      "--outdir <OUTPUT_DIR>"
    )
  }

  trait_dir <- file.path(
    base_dir,
    trait
  )

  mlma_file <- file.path(
    trait_dir,
    paste0(trait, "_MLMA.mlma")
  )

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  element_name <- get_element_name(trait)

  manhattan_png <- file.path(
    out_dir,
    paste0(trait, "_Manhattan.png")
  )

  qq_png <- file.path(
    out_dir,
    paste0(trait, "_QQ.png")
  )

  cat("========================================\n")
  cat("Mode: single\n")
  cat("Trait:", trait, "\n")
  cat("Title:", element_name, "\n")
  cat("MLMA:", mlma_file, "\n")
  cat("Output:", out_dir, "\n")
  cat("========================================\n")

  gwas <- read_mlma(mlma_file)

  # Manhattan PNG
  png(
    filename = manhattan_png,
    width = 14,
    height = 6,
    units = "in",
    res = 300,
    type = "cairo"
  )

  par(
    mar = c(5, 5, 3.5, 1.5)
  )

  draw_manhattan(
    gwas = gwas,
    title_text = element_name
  )

  dev.off()

  # QQ PNG
  png(
    filename = qq_png,
    width = 7,
    height = 7,
    units = "in",
    res = 300,
    type = "cairo"
  )

  par(
    mar = c(5, 5, 3.5, 1.5)
  )

  draw_qq(
    p_values = gwas$P,
    title_text = element_name
  )

  dev.off()

  cat("Manhattan plot:\n")
  cat(manhattan_png, "\n")

  cat("QQ plot:\n")
  cat(qq_png, "\n")

  cat("Finished:", trait, "\n")

  quit(
    save = "no",
    status = 0
  )
}


# ============================================================
# 4. Combined-PDF mode
# ============================================================

if (mode == "combine") {

  if (any(is.na(c(base_dir, trait_list, out_dir)))) {
    stop(
      "Combine mode usage:\n",
      "Rscript script.R ",
      "--mode combine ",
      "--base <GWAS_BASE> ",
      "--trait-list <TRAIT_LIST> ",
      "--outdir <OUTPUT_DIR>"
    )
  }

  if (!file.exists(trait_list)) {
    stop("Trait list not found: ", trait_list)
  }

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  traits <- readLines(trait_list)

  traits <- trimws(traits)
  traits <- traits[nzchar(traits)]

  if (length(traits) == 0) {
    stop("Trait list is empty.")
  }

  trait_table <- data.table(
    trait = traits
  )

  trait_table[
    ,
    element := get_element_name(trait)
  ]

  trait_table[
    ,
    analysis_type := fifelse(
      grepl("_no_outlier$", trait),
      "no_outlier",
      fifelse(
        grepl("_log$", trait),
        "log",
        "other"
      )
    )
  ]

  # no_outlier has priority over log
  trait_table[
    ,
    priority := fifelse(
      analysis_type == "no_outlier",
      1L,
      fifelse(
        analysis_type == "log",
        2L,
        3L
      )
    )
  ]

  trait_table[
    ,
    manhattan_png := file.path(
      base_dir,
      trait,
      "plots",
      paste0(trait, "_Manhattan.png")
    )
  ]

  trait_table[
    ,
    qq_png := file.path(
      base_dir,
      trait,
      "plots",
      paste0(trait, "_QQ.png")
    )
  ]

  missing_manhattan <- trait_table[
    !file.exists(manhattan_png)
  ]

  missing_qq <- trait_table[
    !file.exists(qq_png)
  ]

  if (nrow(missing_manhattan) > 0) {
    print(
      missing_manhattan[
        ,
        .(trait, manhattan_png)
      ]
    )

    stop("One or more Manhattan PNG files are missing.")
  }

  if (nrow(missing_qq) > 0) {
    print(
      missing_qq[
        ,
        .(trait, qq_png)
      ]
    )

    stop("One or more QQ PNG files are missing.")
  }

  # ----------------------------------------------------------
  # PDF 1: all 40 analyses
  # ----------------------------------------------------------

  setorder(
    trait_table,
    element,
    priority,
    trait
  )

  all_pdf <- file.path(
    out_dir,
    paste0(
      "All_",
      nrow(trait_table),
      "_analyses_Manhattan_QQ.pdf"
    )
  )

  cat("Creating all-analysis PDF:\n")
  cat(all_pdf, "\n")

  pdf(
    file = all_pdf,
    width = 14,
    height = 10,
    onefile = TRUE,
    useDingbats = FALSE
  )

  for (i in seq_len(nrow(trait_table))) {

    cat(
      "All:",
      i,
      "/",
      nrow(trait_table),
      trait_table$trait[i],
      "\n"
    )

    draw_png_pdf_page(
      manhattan_file = trait_table$manhattan_png[i],
      qq_file = trait_table$qq_png[i]
    )
  }

  dev.off()

  # ----------------------------------------------------------
  # PDF 2: one result per element
  #
  # Prefer no_outlier over log
  # ----------------------------------------------------------

  preferred_table <- trait_table[
    order(element, priority, trait),
    .SD[1],
    by = element
  ]

  setorder(
    preferred_table,
    element
  )

  selection_file <- file.path(
    out_dir,
    "Preferred_analysis_selected_for_each_element.tsv"
  )

  fwrite(
    preferred_table[
      ,
      .(
        Element = element,
        Selected_analysis = trait,
        Analysis_type = analysis_type,
        Manhattan_file = manhattan_png,
        QQ_file = qq_png
      )
    ],
    selection_file,
    sep = "\t"
  )

  preferred_pdf <- file.path(
    out_dir,
    paste0(
      "Preferred_",
      nrow(preferred_table),
      "_elements_Manhattan_QQ.pdf"
    )
  )

  cat("Creating preferred-analysis PDF:\n")
  cat(preferred_pdf, "\n")

  pdf(
    file = preferred_pdf,
    width = 14,
    height = 10,
    onefile = TRUE,
    useDingbats = FALSE
  )

  for (i in seq_len(nrow(preferred_table))) {

    cat(
      "Preferred:",
      i,
      "/",
      nrow(preferred_table),
      preferred_table$trait[i],
      "\n"
    )

    draw_png_pdf_page(
      manhattan_file = preferred_table$manhattan_png[i],
      qq_file = preferred_table$qq_png[i]
    )
  }

  dev.off()

  cat("========================================\n")
  cat("PDF generation completed\n")
  cat("All analyses:", nrow(trait_table), "\n")
  cat("Unique elements:", nrow(preferred_table), "\n")
  cat("\nAll-analysis PDF:\n")
  cat(all_pdf, "\n")
  cat("\nPreferred-analysis PDF:\n")
  cat(preferred_pdf, "\n")
  cat("\nSelection table:\n")
  cat(selection_file, "\n")
  cat("========================================\n")

  quit(
    save = "no",
    status = 0
  )
}


# ============================================================
# 5. Invalid mode
# ============================================================

stop(
  "Invalid or missing --mode.\n",
  "Valid modes are:\n",
  "  --mode single\n",
  "  --mode combine\n"
)