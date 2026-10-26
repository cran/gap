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
htmltools::tagList(
  htmltools::tags$iframe(
    id = "gapframe",
    src = "https://jinghuazhao.github.io/R/vignettes/gap.html",
    style = "
      width: 100%;
      border: none;
      display: block;
      min-height: 200px;
    "
  ),
  htmltools::tags$script(htmltools::HTML("
window.addEventListener('message', function(event) {
  if (!event.data || event.data.type !== 'setHeight') return;
  var iframe = document.getElementById('gapframe');
  if (!iframe) return;
  iframe.style.height = Math.max(event.data.height, 200) + 'px';
});
  "))
)

