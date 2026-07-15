# Volatility Forecasting with GARCH Models under Robust Loss Functions

Sebastiano Pinotti — July 2026

This project compares twelve GARCH-family specifications: GARCH(1,1), GARCH-in-mean, EGARCH, and APARCH, each estimated under Gaussian, Student-t, and GED innovations, on two assets with deliberately different risk profiles: Goldman Sachs (GS) and the iShares Silver Trust (SLV). Models are estimated on 2020–2024 daily returns and evaluated on one-step-ahead rolling volatility forecasts over 2025–2026, a test period that includes an extreme SLV return near −30%.

## Headline finding

The evaluation loss function is not a neutral implementation detail. Under the non-robust MAE loss, asymmetric GED models appear significantly superior for SLV out-of-sample, contradicting the in-sample evidence. Under the proxy-robust QLIKE and MSE losses of Patton (2011), the effect vanishes entirely: no pairwise Diebold–Mariano comparison is significant at the 5% level. The apparent result was an artifact of the evaluation criterion: MAE's optimal forecast is the conditional *median* of the squared return, far below the conditional variance for fat-tailed returns, so it systematically rewards downward-biased volatility forecasts. For GS, by contrast, the asymmetric Student-t advantage survives robust evaluation: EGARCH+t significantly outperforms every symmetric specification under QLIKE.

Full methodology, results, and discussion: [`garch_volatility_report.pdf`](garch_volatility_report.pdf).

## Reproducing the results

Requirements: R ≥ 4.0 with `xts`, `readxl`, `fBasics`, `rugarch`, `FinTS`, `sandwich`.

The script runs given two Excel files of daily closing prices (`data/PriceHistory_GS.xlsx`, `data/PriceHistory_SLV.xlsx`). The rolling forecast step (12 models × 2 assets × 58 re-estimations) is the expensive part and is cached to `output/cache/rolls.rds` on first run. Delete the cache after changing specifications, refit frequency, or data. 

## Data

Daily closing prices for GS and SLV, January 2, 2020 – February 27, 2026, obtained from LSEG Refinitiv. The raw files are not redistributed due to licensing terms. All results are reproducible from any daily close series for the same tickers and period (e.g. Yahoo Finance): daily closes for liquid large caps and ETFs are identical across vendors up to rounding.

## Key reference

Patton, A.J. (2011). Volatility forecast comparison using imperfect volatility proxies. *Journal of Econometrics*, 160(1), 246–256.

## License

MIT, see [LICENSE](LICENSE).

