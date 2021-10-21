#' Create a balance table
#'
#' This function creates a standard balance tables as used in economics.
#' It takes a data set (coerced to a data.table), identifying variables,
#' a group (or treatment) variable, and variables to balance across groups.
#' Optionally, it may save the output to a LaTeX file and change the names of
#' the balanced variables.
#'
#' @param data Data set to be used. Can be a panel (one row per observation-period), with
#' variables to be balanced as columns if that is the case, set `panel` to `TRUE`.
#' Alternatively, `data` may be in long format, with each row a unique combination of
#' observation, period and variable. A column `value` should be specified, which will hold
#' the unique value of the row.
#' @param id.vars The identifying variables, uniquely identifying an observation and year.
#' @param bal.vars The variables to balance across groups.
#' @param group.var The variable identifying the two groups, treatment and control. Should be
#' binary, with 0 identifying the control group and 1 the treatment.
#' @param panel By default `F`. Switch to `T` is data set is a panel with variables to be
#' balanced as columns.
#' @param new.names.bal.vars A named vector used to rename the variables to be balanced.
#' @param save.path By default `NULL`. Indicate a path where the balance table should be saved.
#' @return A `gt` table object.
#' @export
btable <- function(data, id.vars, bal.vars, group.var,
                   panel = F, new.names.bal.vars = NULL,
                   save.path = NULL) {

  # Data should already be circumscribed to the desired period!

  ### Check if data is data.table and convert ---------------------------
  if(!is.data.table(data)) data = as.data.table(data)

  ### If the data in panel form, transform to long ---------------------------
  if (panel == T){

    # Select columns to keep
    data.clean = data[, c(id.vars, bal.vars, group.var)]

    # Melt the panel to long
    data.l = melt(data.clean,
                  id.vars = c(id.vars, group.var),
                  measure.vars = bal.vars,
                  variable.name = 'var')
  } else {
    data.l = data
  }

  ### Group computation ---------------------------
  ## N
  btab_gr.N = data.l[, .(N = .N),
                     by = c(group.var, 'var')]
  # Cast to wide
  btab.N.w = dcast(btab_gr.N, as.formula(str_glue("var ~ {group.var}")), value.var = c('N'))
  setkey(btab.N.w, var)

  ## Mean and SE
  btab_gr.mse = data.l[, .(mean = mean(value),
                           se = stderr(value)),
                       by = c(group.var, 'var')]
  # Reshapes
  btab.mse.l = melt(btab_gr.mse, id.vars = c(group.var, 'var'), variable.name = 'stat')

  btab.mse.w = dcast(btab.mse.l, as.formula(str_glue("var + stat ~ {group.var}")),
                     value.var = c('value'))
  setkey(btab.mse.w, var)

  ## Merge and compute mean difference
  btab.gr = btab.mse.w[btab.N.w]
  setnames(btab.gr, c('0', '1', 'i.0', 'i.1'), c('g_0', 'g_1', 'N_0', 'N_1'))
  # Mean differences
  btab.gr[stat == 'mean', diff := g_1-g_0]
  # Prepare for merge
  setkey(btab.gr, var, stat)

  ### Overall computations ---------------------------
  ## N
  btab_all.N = data[, .(N_all = .N),
                    by = 'var']
  setkey(btab_all.N, var)

  ## Mean and SE
  btab_all.mse = data[, .(mean = mean(value),
                          se = stderr(value)),
                      by = c('var')]
  # Reshapes
  btab_all.mse.l = melt(btab_all.mse, id.vars = c('var'), variable.name = 'stat', value.name = 'g_all')
  setkey(btab_all.mse.l, var, stat)



  ### Merge Group and Overall ---------------------------
  btab = btab.gr[btab_all.mse.l]
  btab = btab[btab_all.N]
  # Erase some cells
  btab[stat == 'se', `:=` (N_0 = NA, N_1 = NA, N_all = NA)]

  ### Compute ttests ---------------------------
  ## t-test loop and storage
  for (v in bal.vars) {
    if (match(v, bal.vars) == 1){ # Create a list if first item of variable list
      ttests = vector(mode = 'list', length = length(bal.vars))
      names(ttests) = bal.vars # Name the list
    }
    # Store t-test in list
    ttests[v] = t.test(as.formula(str_glue("value ~ {group.var}")), data = data[var == v])$p.value
  }

  ## Transform in data.table and clean up
  ttests_dt = as.data.table(ttests)
  ttests_dt_t = data.table::transpose(ttests_dt, keep.names = 'var')
  ttests_dt_t[, stat := 'se']
  setnames(ttests_dt_t, 'V1', 'p_val')
  setkey(ttests_dt_t, var, stat)

  ### Merge t-tests ---------------------------
  btab.fin = merge(btab, ttests_dt_t, by = c('var', 'stat'), all.x = T)
  btab.fin[stat == 'se', diff := p_val][, p_val := NULL] # Erase some cells
  btab.fin[, stat := NULL] # Remove stat column
  btab.fin[(1:nrow(btab.fin)%%2 == 0), var := ''] # Remove variable name for second lineline

  ### Build the table ---------------------------
  ## Change variable names if provided
  if(!is.null(new.names.bal.vars)){
    for(n in names(new.names.bal.vars)){
      btab.fin[var == n, var := new.names.bal.vars[n]]
    }
  }

  gt_tb =
    gt(btab.fin) %>%
    tab_spanner(
      label = 'Control',
      columns = c(N_0, g_0)
    ) %>%
    tab_spanner(
      label = 'Treatment',
      columns = c(N_1, g_1)
    ) %>%
    tab_spanner(
      label = 'Overall',
      columns = c(N_all, g_all)
    ) %>%
    tab_spanner(
      label = 'T -- C',
      columns = c(diff)
    ) %>%
    cols_label(
      var = md('*Variable*'),
      N_0 = 'N',
      g_0 = 'mean/se',
      N_1 = 'N',
      g_1 = 'mean/se',
      N_all = 'N',
      g_all = 'mean/se',
      diff = 'diff./p-value'
    ) %>% cols_move_to_end(diff) %>%
    fmt_missing(everything(), missing_text = " ") %>%
    fmt_number(
      columns = matches('g_|diff'),
      rows = !is.na(N_0),
      decimals = 3
    ) %>%
    fmt_number(
      columns = matches('g_|diff'),
      rows = is.na(N_0),
      decimals = 3,
      pattern = '({x})'
    ) %>%
    cols_align(
      align = 'center'
    ) %>%
    cols_align(
      align = 'left',
      columns = 'var'
    )

  ### Display balance table in viewer  ---------------------------
  return(gt_tb)

  ### If file path defined, save  ---------------------------
  if (!is.null(save.path)) {
    gtsave(gt_tb, filename = save.path)
  }
}
