## ----setup, include=FALSE-----------------------------------------------------
set.seed(0)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 8,
  fig.height = 8,
  out.extra = 'style="display:block; margin:auto"',
  dev = "CairoPNG"
)

## ----gap_incl, echo=FALSE, warning=FALSE--------------------------------------
htmltools::tags$iframe(
  src = "https://jinghuazhao.github.io/R/vignettes/gap.html",
  width = "100%",
  height = "600px",
  style = "border:none;"
)

