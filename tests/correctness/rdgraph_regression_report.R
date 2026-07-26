#!/usr/bin/env Rscript

find.repo.root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(cur, "DESCRIPTION")) &&
        dir.exists(file.path(cur, "R"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not locate repository root from ", start, call. = FALSE)
    }
    cur <- parent
  }
}

html.escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

format.metric.value <- function(x) {
  if (is.numeric(x)) {
    ifelse(is.na(x), "", formatC(x, digits = 4L, format = "fg"))
  } else {
    html.escape(x)
  }
}

html.table <- function(df, columns = names(df)) {
  if (is.null(df) || nrow(df) == 0L) return("<p>No rows.</p>")
  df <- df[, columns, drop = FALSE]
  header <- paste0("<th>", html.escape(names(df)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(df)), function(i) {
    cells <- vapply(df, function(col) format.metric.value(col[[i]]), character(1L))
    paste0("<tr>", paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
  }, character(1L))
  paste0(
    "<table><thead><tr>", header, "</tr></thead><tbody>",
    paste(rows, collapse = "\n"),
    "</tbody></table>"
  )
}

relative.path <- function(path, base.dir) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  base.dir <- normalizePath(base.dir, winslash = "/", mustWork = TRUE)
  if (startsWith(path, paste0(base.dir, "/"))) {
    return(substring(path, nchar(base.dir) + 2L))
  }
  path
}

html.relative.href <- function(path, base.dir) {
  html.escape(relative.path(path, base.dir))
}

source.correctness.file <- function(repo.root, file) {
  source(file.path(repo.root, "tests", "correctness", file), chdir = TRUE)
}

fit.rdgraph.variant <- function(data, y, fit.args, graph.model = "iknn", graph.spec = NULL) {
  elapsed <- system.time({
    fit.or.error <- tryCatch({
      args <- c(list(X = data$X, y = y), fit.args)
      if (!is.null(graph.spec)) {
        args$k <- 2L
        args$adj.list <- graph.spec$adj.list
        args$weight.list <- graph.spec$weight.list
        args$use.counting.measure <- TRUE
        args$apply.geometric.pruning <- FALSE
        args$max.ratio.threshold <- 0
        args$threshold.percentile <- 0
      }
      do.call(gflowx::fit.rdgraph.regression, args)
    }, error = function(e) e)
  })[["elapsed"]]

  if (inherits(fit.or.error, "error")) {
    list(y = y, fit = NULL, error = fit.or.error, elapsed = elapsed, graph.model = graph.model)
  } else {
    list(y = y, fit = fit.or.error, error = NULL, elapsed = elapsed, graph.model = graph.model)
  }
}

run.rdgraph.case <- function(case) {
  set.seed(case$seed)
  data <- case$generate()

  variants <- list()
  for (response.name in names(data$responses)) {
    y <- as.double(data$responses[[response.name]])
    variants[[paste(response.name, "iknn", sep = "::")]] <-
      fit.rdgraph.variant(data, y, case$fit.args, graph.model = "iknn")

    if (identical(response.name, names(data$responses)[[1L]])) {
      oracle.graph <- rdgraph.case.oracle.graph(case, data)
      if (!is.null(oracle.graph)) {
        variants[[paste(response.name, oracle.graph$model, sep = "::")]] <-
          fit.rdgraph.variant(
            data,
            y,
            case$fit.args,
            graph.model = oracle.graph$model,
            graph.spec = oracle.graph
          )
        case$oracle.graph.command <- oracle.graph$command
      }
    }
  }

  result <- list(case = case, data = data, variants = variants)
  result$metrics <- rdgraph.case.metrics(result)
  result$k.sweep <- rdgraph.run.k.sweep(case, data)
  result
}

render.case.section <- function(case.result, report.dir) {
  plot.tags <- vapply(case.result$plots, function(plot.info) {
    rel <- relative.path(plot.info$path, report.dir)
    paste0(
      "<figure><img src=\"", html.escape(rel), "\" alt=\"", html.escape(plot.info$title),
      "\"><figcaption>", html.escape(plot.info$title), "</figcaption></figure>"
    )
  }, character(1L))

  metric.columns <- c(
    "response_variant",
    "graph_model",
    "fit_status",
    "rmse_raw",
    "rmse_fit",
    "rmse_improvement",
    "rmse_improvement_fraction",
    "cor_fit_truth",
    "var_ratio_fit_truth",
    "boundary_rmse_fit",
    "fit_elapsed_sec",
    "n_edges"
  )

  paste0(
    "<section class=\"case\">",
    "<h2>", html.escape(case.result$case$title), "</h2>",
    "<p><code>", html.escape(case.result$case$id), "</code> | group: ",
    html.escape(case.result$case$group), "</p>",
    "<h3>Reproducible Commands</h3>",
    "<p>Data generation:</p><pre><code>", html.escape(case.result$case$data.command %||% ""), "</code></pre>",
    "<p>Conditional expectation fit:</p><pre><code>", html.escape(case.result$case$fit.command %||% ""), "</code></pre>",
    if (!is.null(case.result$case$oracle.graph.command)) {
      paste0(
        "<p>Oracle graph fit:</p><pre><code>",
        html.escape(case.result$case$oracle.graph.command),
        "</code></pre>"
      )
    } else "",
    html.table(case.result$metrics, columns = metric.columns),
    paste(plot.tags, collapse = "\n"),
    render.k.sweep.section(case.result, report.dir),
    "</section>"
  )
}

render.k.sweep.section <- function(case.result, report.dir) {
  ks <- case.result$k.sweep
  if (is.null(ks)) return("")

  plot.tags <- if (length(ks$plots)) {
    vapply(ks$plots, function(plot.info) {
      rel <- relative.path(plot.info$path, report.dir)
      paste0(
        "<figure><img src=\"", html.escape(rel), "\" alt=\"", html.escape(plot.info$title),
        "\"><figcaption>", html.escape(plot.info$title), "</figcaption></figure>"
      )
    }, character(1L))
  } else character(0)

  paste0(
    "<h3>k Selection Measures</h3>",
    "<p><strong>Oracle RMSE</strong>: RMSE between fitted values and the known synthetic truth. Lower is better and is available here because this report uses synthetic data.</p>",
    "<p><strong>GCV</strong>: generalized cross-validation reported by <code>gflowx::fit.rdgraph.regression()</code> at the selected iteration. Lower is better.</p>",
    "<p><strong>Min-normalized GCV</strong>: <code>GCV(k) / min_k GCV(k)</code>. Lower is better; the best possible value is 1.</p>",
    "<p><strong>Robust-normalized GCV</strong>: <code>(GCV(k) - median_k GCV(k)) / MAD_k GCV(k)</code>, with IQR/SD/constant fallbacks. Lower is better.</p>",
    "<p><strong>Graph structural stability</strong>: adjacent-k graphs are compared by graph edit distance and Jensen-Shannon divergence of degree distributions. Smaller values indicate more stable adjacent graph structure.</p>",
    "<p><strong>Largest connected component fit</strong>: when a k graph is disconnected, the report refits <code>gflowx::fit.rdgraph.regression()</code> on the largest connected component only, using the induced weighted graph. LCC metrics are computed only on vertices inside that component and are reported separately from full-graph metrics.</p>",
    "<p><strong>Component-wise fit</strong>: when multiple components are scientifically meaningful, the report also refits each connected component with at least 10 vertices and computes aggregate component-wise metrics over all fitted components.</p>",
    "<h4>k Selection Diagnostics</h4>",
    html.table(ks$selection),
    render.lcc.sweep.section(ks),
    render.component.sweep.section(ks),
    "<h4>k Sweep Summary</h4>",
    html.table(ks$summary),
    "<h4>Adjacent-k Stability</h4>",
    html.table(ks$stability),
    paste(plot.tags, collapse = "\n")
  )
}

render.lcc.sweep.section <- function(ks) {
  s <- ks$summary
  if (is.null(s) || !nrow(s) || !"lcc_fit_status" %in% names(s)) return("")
  lcc <- s[!is.na(s$lcc_fit_status) & s$lcc_fit_status != "connected_graph", , drop = FALSE]
  if (!nrow(lcc)) return("")

  columns <- intersect(c(
    "k",
    "component_count",
    "largest_component_size",
    "largest_component_fraction",
    "second_largest_component_size",
    "lcc_fit_status",
    "lcc_n",
    "lcc_fraction",
    "lcc_rmse_fit",
    "lcc_cor_fit_truth",
    "lcc_gcv",
    "lcc_elapsed_sec",
    "lcc_edge_count",
    "lcc_error_message"
  ), names(lcc))

  paste0(
    "<h4>Largest Connected Component Fits</h4>",
    "<p>Rows appear only for disconnected k graphs. Points outside the largest component are excluded from these fits and from the LCC-only metrics.</p>",
    html.table(lcc, columns = columns)
  )
}

render.component.sweep.section <- function(ks) {
  s <- ks$summary
  if (is.null(s) || !nrow(s) || !"component_fit_status" %in% names(s)) return("")
  component <- s[!is.na(s$component_fit_status) & s$component_fit_status != "connected_graph", , drop = FALSE]
  if (!nrow(component)) return("")

  columns <- intersect(c(
    "k",
    "component_count",
    "largest_component_size",
    "largest_component_fraction",
    "second_largest_component_size",
    "component_fit_status",
    "component_fit_components",
    "component_fit_n",
    "component_fit_fraction",
    "component_fit_rmse",
    "component_fit_cor_truth",
    "component_fit_elapsed_sec",
    "component_fit_error_message"
  ), names(component))

  paste0(
    "<h4>Component-Wise Fits</h4>",
    "<p>Rows appear for disconnected k graphs. Each connected component with at least 10 vertices is fit separately, and aggregate metrics are computed over the fitted vertices.</p>",
    html.table(component, columns = columns)
  )
}

rdgraph.report.style <- function() {
  paste0(
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:32px;line-height:1.45;color:#1f2933;}",
    "h1,h2{line-height:1.15;} code{background:#f2f4f7;padding:2px 4px;border-radius:4px;}",
    "table{border-collapse:collapse;margin:16px 0;width:100%;font-size:13px;}",
    "th,td{border:1px solid #d9dee7;padding:6px 8px;text-align:left;vertical-align:top;}",
    "th{background:#f6f8fb;} figure{margin:22px 0;} img{max-width:100%;border:1px solid #d9dee7;}",
    "figcaption{font-size:13px;color:#52616f;margin-top:6px;} .case{border-top:2px solid #d9dee7;padding-top:18px;margin-top:28px;}",
    "pre{background:#f6f8fb;padding:12px;overflow:auto;font-size:12px;}",
    ".report-nav{margin:0 0 18px 0;} .report-nav a{margin-right:14px;}"
  )
}

write.rdgraph.report <- function(case.results,
                                 metrics,
                                 report.path,
                                 report.dir,
                                 report.title = "Rdgraph Regression Correctness Report",
                                 report.description = NULL) {
  session.text <- paste(utils::capture.output(utils::sessionInfo()), collapse = "\n")
  generated.at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  sections <- vapply(case.results, render.case.section, character(1L), report.dir = report.dir)

  metric.columns <- c(
    "round_id",
    "case_id",
    "response_variant",
    "graph_model",
    "fit_status",
    "rmse_raw",
    "rmse_fit",
    "rmse_improvement",
    "cor_fit_truth",
    "fit_elapsed_sec"
  )
  metric.columns <- intersect(metric.columns, names(metrics))

  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\">",
    "<title>", html.escape(report.title), "</title>",
    "<style>", rdgraph.report.style(), "</style></head><body>",
    "<nav class=\"report-nav\"><a href=\"index.html\">Report index</a>",
    "<a href=\"rdgraph_regression_correctness.html\">Combined report</a></nav>",
    "<h1>", html.escape(report.title), "</h1>",
    "<p>Generated at ", html.escape(generated.at), ".</p>",
    if (!is.null(report.description)) paste0("<p>", html.escape(report.description), "</p>") else "",
    "<p>This report is visual-first: plots are intended for human review, while the metric table records candidate values for later automation.</p>",
    "<h2>Metric Summary</h2>",
    html.table(metrics, columns = metric.columns),
    paste(sections, collapse = "\n"),
    "<h2>Session Info</h2><pre>", html.escape(session.text), "</pre>",
    "</body></html>"
  )

  writeLines(html, report.path)
}

write.rdgraph.report.index <- function(case.results, metrics, report.records, index.path, report.dir) {
  generated.at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  rows <- do.call(rbind, lapply(report.records, function(rec) {
    data.frame(
      report = paste0("<a href=\"", html.relative.href(rec$path, report.dir), "\">", html.escape(rec$title), "</a>"),
      cases = rec$n.cases,
      description = html.escape(rec$description %||% ""),
      stringsAsFactors = FALSE
    )
  }))
  html.table.raw <- function(df) {
    header <- paste0("<th>", names(df), "</th>", collapse = "")
    body <- vapply(seq_len(nrow(df)), function(i) {
      cells <- vapply(df, function(col) as.character(col[[i]]), character(1L))
      paste0("<tr>", paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
    }, character(1L))
    paste0("<table><thead><tr>", header, "</tr></thead><tbody>", paste(body, collapse = "\n"), "</tbody></table>")
  }

  case.rows <- do.call(rbind, lapply(case.results, function(res) {
    data.frame(
      round = res$case$round$title %||% "",
      case_id = res$case$id,
      title = res$case$title,
      group = res$case$group,
      stringsAsFactors = FALSE
    )
  }))

  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\">",
    "<title>Rdgraph Regression Correctness Report Index</title>",
    "<style>", rdgraph.report.style(), "</style></head><body>",
    "<h1>Rdgraph Regression Correctness Report Index</h1>",
    "<p>Generated at ", html.escape(generated.at), ".</p>",
    "<p>Open a round-specific report for focused review, or the combined report for a full-suite pass.</p>",
    "<h2>Reports</h2>",
    html.table.raw(rows),
    "<h2>Cases</h2>",
    html.table(case.rows),
    "</body></html>"
  )

  writeLines(html, index.path)
}

write.rdgraph.round.reports <- function(case.results, metrics, report.dir) {
  round.ids <- unique(vapply(case.results, function(res) res$case$round$id %||% "unassigned", character(1L)))
  records <- list()

  for (round.id in round.ids) {
    keep <- vapply(case.results, function(res) identical(res$case$round$id, round.id), logical(1L))
    round.results <- case.results[keep]
    round.meta <- round.results[[1L]]$case$round
    round.metrics <- metrics[metrics$round_id == round.id, , drop = FALSE]
    report.path <- file.path(report.dir, round.meta$file)
    write.rdgraph.report(
      round.results,
      round.metrics,
      report.path,
      report.dir,
      report.title = round.meta$title,
      report.description = round.meta$description
    )
    metrics.path <- file.path(report.dir, paste0(tools::file_path_sans_ext(round.meta$file), "_metrics.csv"))
    utils::write.csv(round.metrics, metrics.path, row.names = FALSE)
    records[[length(records) + 1L]] <- list(
      title = round.meta$title,
      description = round.meta$description,
      path = report.path,
      metrics.path = metrics.path,
      n.cases = length(round.results)
    )
  }

  records
}

main <- function() {
  repo.root <- find.repo.root()
  setwd(repo.root)

  source.correctness.file(repo.root, "rdgraph_metrics.R")
  source.correctness.file(repo.root, "rdgraph_graphs.R")
  source.correctness.file(repo.root, "rdgraph_k_sweep.R")
  source.correctness.file(repo.root, "rdgraph_cases.R")
  source.correctness.file(repo.root, "rdgraph_plots.R")

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The pkgload package is required to run the correctness report.", call. = FALSE)
  }
  pkgload::load_all(repo.root, quiet = TRUE)

  report.dir <- file.path(repo.root, "tests", "manual", "reports", "rdgraph-regression-correctness")
  plot.dir <- file.path(report.dir, "plots")
  if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

  cases <- rdgraph.regression.correctness.cases()
  case.results <- vector("list", length(cases))

  for (i in seq_along(cases)) {
    message(sprintf("[%d/%d] Running %s", i, length(cases), cases[[i]]$id))
    case.results[[i]] <- run.rdgraph.case(cases[[i]])
    case.results[[i]]$plots <- rdgraph.write.case.plots(case.results[[i]], plot.dir)
    if (!is.null(case.results[[i]]$k.sweep)) {
      case.results[[i]]$k.sweep$plots <- rdgraph.write.k.sweep.plots(case.results[[i]], plot.dir)
    }
  }

  metrics <- do.call(rbind, lapply(case.results, `[[`, "metrics"))
  metrics.path <- file.path(report.dir, "metrics.csv")
  utils::write.csv(metrics, metrics.path, row.names = FALSE)

  report.path <- file.path(report.dir, "rdgraph_regression_correctness.html")
  write.rdgraph.report(
    case.results,
    metrics,
    report.path,
    report.dir,
    report.title = "Rdgraph Regression Correctness Report: Combined",
    report.description = "Combined report containing all currently implemented correctness cases."
  )
  report.records <- write.rdgraph.round.reports(case.results, metrics, report.dir)
  report.records <- c(list(list(
    title = "Combined Report",
    description = "Full-suite report containing all rounds.",
    path = report.path,
    metrics.path = metrics.path,
    n.cases = length(case.results)
  )), report.records)
  index.path <- file.path(report.dir, "index.html")
  write.rdgraph.report.index(case.results, metrics, report.records, index.path, report.dir)

  message("Wrote index: ", index.path)
  message("Wrote report: ", report.path)
  for (rec in report.records[-1L]) {
    message("Wrote round report: ", rec$path)
  }
  message("Wrote metrics: ", metrics.path)
}

if (identical(environment(), globalenv())) {
  main()
}
