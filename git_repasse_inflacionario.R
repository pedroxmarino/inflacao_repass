library(rbcb)
library(dplyr)
library(tsibble)
library(purrr)
library(fabletools)
library(feasts)
library(rio)
library(tidyr)
library(restriktor)
library(forecast)
library(broom)
library(ggplot2)

exp_ipca <- rbcb::get_twelve_months_inflation_expectations(
  indic = "IPCA",
  end_date = Sys.Date()
)

View(exp_ipca)

#Trimestralização das expectativas pela média das medianas
exp_ipca_aux <- exp_ipca |> 
  dplyr::filter(base==0, smoothed =="S") |> 
  dplyr::group_by(data = tsibble::yearquarter(date)) |> 
  dplyr::summarise(ipca_exp_12m = mean(median), .groups = "drop")

# Inflação total e de preços livres
dados_ipca <- rbcb::get_series(
  code = c("ipca_total" = 433, "ipca_livres"= 11428),
  start_date="1999-12-01"
) |> 
  purrr::reduce(.f=dplyr::full_join, by="date") |> # reduce aplica função a cada elemento do conjunto
  dplyr::group_by(data=tsibble::yearquarter(date)) |> 
  dplyr::summarise(
    dplyr::across(
      .cols = c("ipca_total", "ipca_livres"),
      .fns = ~(prod((.x/100)+1)-1)*100),
    .groups ="drop"
  )

# ajustar sazonalmente IPCA Livres
ipca_livres_sa <- dados_ipca |> 
  dplyr::filter(data>=tsibble::yearquarter("1999 Q4"), !is.na(ipca_livres)) |> 
  tsibble::as_tsibble(index="data") |> 
  fabletools::model(x11 = feasts::X_13ARIMA_SEATS(ipca_livres ~x11())) |> 
  fabletools::components() |> 
  dplyr::select("data", "ipca_livres_sa" = "season_adjust")


# juntar com os dados do IPCA
dados_ipca <- dplyr::left_join(
  x= dados_ipca,
  y = ipca_livres_sa,
  by = "data"
)

# Coletando IC-BR
dados_ic <- rbcb::get_series(code = c("ic_br"=27574), start_date="1998-01-01") |> 
  dplyr::group_by(data= tsibble::yearquarter(date)) |> 
  dplyr::summarise(
    ic_br = (log(ic_br / dplyr::lag(ic_br, n = 2)) * 100) |> 
      dplyr::last(),
    .groups = "drop"
  )

# Coletando hiato do Porduto do IFI
hiato_ifi <- rio::import(
  file = paste0("https://www12.senado.leg.br/ifi/dados/arquivos/estimativas-do-hiato-do",
                "-produto-ifi/@@download/file/Hiato%20do%20Produto%20IFI.xlsx"),
  formate = "xlsx",
  setclass = "tibble",
  sheet = "Hiato do Produto",
  skip = 1
) |> 
  dplyr::mutate(
    data= tsibble::yearquarter(`Trim-Ano`),
    hiato = `Hiato` * 100,
    .keep = "none"
  )

# Coletando hiato do produto Bacen
hiato_bcb <- rio::import(
  file = paste0(
    "https://www.bcb.gov.br/content/ri/relatorioinflacao/",
    "202303/ri202303anp.xlsx"
  ),
  format = "xlsx",
  setclass = "tibble",
  sheet = "Graf 2.2.4",
  skip=8
) |> 
  dplyr::mutate(
    data=tsibble::yearquarter(`Trimestre`),
    hiato =`Hiato`,
    .keep = "none"
  )

# Reunindo os dados
dados_reg <- purrr::reduce(
  .x= list(dados_ipca, dados_ic, hiato_ifi, exp_ipca_aux),
  .f = dplyr::full_join,
  by= "data"
) |> 
  dplyr::arrange(data) |> 
  dplyr::mutate(
    ipca_lag1 = dplyr::lag(ipca_total,1),
    ipca_lag2 = dplyr::lag(ipca_total,2),
    ic_lag1 = dplyr::lag(hiato,1),
    hiato_lag3 = dplyr::lag(hiato, 3)
  ) |> 
  dplyr::filter(data>=tsibble::yearquarter("2002 Q1")) |> 
  tidyr::drop_na()

# Estimando os modelos
modelo_irrestrito <- lm(
  formula = ipca_livres ~ -1 + ipca_lag1 + ipca_lag2 + ipca_exp_12m + hiato_lag3 + ic_lag1,
  data = dados_reg
)

modelo_restrito <- restriktor::restriktor(
  object = modelo_irrestrito,
  constraints = "ipca_lag1 + ipca_lag2 + ipca_exp_12m + ic_lag1 == 1"
)

summary(modelo_restrito)

# Olhando resíduos
forecast::checkresiduals(modelo_restrito)


irr_tidy <- broom::tidy(modelo_irrestrito)
res_tidy <- dplyr::tibble(
  term = irr_tidy$term,
  estimate = modelo_restrito$b.restr,
  std.error = summary(modelo_restrito)$coefficients[, 2],
  low = estimate - 1.64*std.error,
  high = estimate + 1.64*std.error
)

res_tidy |>
  dplyr::filter(!term %in% c("quarter1", "quarter2", "quarter3", "quarter4")) |>
  ggplot2::ggplot() + 
  ggplot2::aes(x = estimate, y = term, xmin = low, xmax = high, height = 0) +
  ggplot2::geom_point() +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::geom_errorbarh() +
  ggplot2::scale_x_continuous(breaks = seq(-0.5, 1, 0.25)) +
  ggplot2::labs(
    title = "Coeficientes estimados da Curva de Phillips",
    subtitle = "Modelo com restrição de verticalidade",
    x = "",
    y = ""
  )

res_plot <- dados_reg |>
  dplyr::transmute(
    data = data,
    resid = modelo_restrito$residuals[, 1],
    ipca_total_lag1 = coef(modelo_restrito)["ipca_lag1"] * ipca_lag1,
    ipca_total_lag2 = coef(modelo_restrito)["ipca_lag2"] * ipca_lag1,
    ipca_exp_12m = coef(modelo_restrito)["ipca_exp_12m"] * ipca_exp_12m,
    hiato_lag3 = coef(modelo_restrito)["hiato_lag3"] * hiato_lag3,
    ic_br_lag1 = coef(modelo_restrito)["ic_lag1"] * ic_lag1
  ) |>
  tidyr::pivot_longer(cols = -"data", names_to = "term", values_to = "contribution")

ggplot2::ggplot(res_plot) +
  ggplot2::aes(x = as.Date(data), y = contribution, fill = term) +
  ggplot2::geom_col() +
  ggplot2::labs(
    title = "Contribuição de cada variável no IPCA livres",
    x = NULL,
    y = NULL,
    fill = NULL
  )


