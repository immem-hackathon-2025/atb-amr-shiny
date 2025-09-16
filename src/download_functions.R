IconButton <- function(outputId, type, ...) {
  if (type == 'data_dl') {
    s.class <- 'shiny-download-link'
    icon <- icon('table')
  } else if (type == 'graph_modal') {
    s.class <- 'action-button'
    icon <- icon('chart-area')
  } else {
    stop('Got bad IconButton type')
  }
  aTag <- tags$a(
    id=outputId,
    class=paste('btn btn-default', s.class),
    style='padding: 2px 4px;',
    href='',
    target='_blank',
    download=NA,
    icon,
    '',
    ...
  )
}