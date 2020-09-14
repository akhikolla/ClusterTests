## ---- eval = FALSE------------------------------------------------------------
#  data(psid, package = "bife")
#  head(psid)

## ---- eval = FALSE------------------------------------------------------------
#  ##    ID LFP KID1 KID2 KID3     INCH AGE TIME
#  ## 1:  1   1    1    1    1 58807.81  26    1
#  ## 2:  1   1    1    0    2 41741.87  27    2
#  ## 3:  1   1    0    1    2 51320.73  28    3
#  ## 4:  1   1    0    1    2 48958.58  29    4
#  ## 5:  1   1    0    1    2 53634.62  30    5
#  ## 6:  1   1    0    0    3 50983.13  31    6

## ---- eval = FALSE------------------------------------------------------------
#  library(bife)
#  stat <- bife(LFP ~ KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID, psid, "probit")
#  summary(stat)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID
#  ##
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.7144874  0.0562414 -12.704   < 2e-16 ***
#  ## KID2      -0.4114824  0.0515524  -7.982  1.44e-15 ***
#  ## KID3      -0.1298719  0.0415477  -3.126   0.00177 **
#  ## log(INCH) -0.2417756  0.0541720  -4.463  8.08e-06 ***
#  ## AGE        0.2319779  0.0375351   6.180  6.40e-10 ***
#  ## I(AGE^2)  -0.0028846  0.0004989  -5.781  7.41e-09 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 6058.88,
#  ## null deviance= 8152.05,
#  ## nT= 5976, N= 664
#  ##
#  ## ( 7173 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6
#  ##
#  ## Average individual fixed effect= -1.121

## ---- eval = FALSE------------------------------------------------------------
#  apes_stat <- get_APEs(stat)
#  summary(apes_stat)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.0927845  0.0176177  -5.267  1.39e-07 ***
#  ## KID2      -0.0534358  0.0158732  -3.366  0.000762 ***
#  ## KID3      -0.0168654  0.0132361  -1.274  0.202595
#  ## log(INCH) -0.0313974  0.0165434  -1.898  0.057713 .
#  ## AGE        0.0301250  0.0117372   2.567  0.010269 *
#  ## I(AGE^2)  -0.0003746  0.0001563  -2.397  0.016518 *
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## ---- eval = FALSE------------------------------------------------------------
#  stat_bc <- bias_corr(stat)
#  summary(stat_bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID
#  ##
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.6309003  0.0555075 -11.366   < 2e-16 ***
#  ## KID2      -0.3635502  0.0511328  -7.110  1.16e-12 ***
#  ## KID3      -0.1149806  0.0413488  -2.781   0.00542 **
#  ## log(INCH) -0.2139642  0.0536616  -3.987  6.68e-05 ***
#  ## AGE        0.2052755  0.0373054   5.503  3.74e-08 ***
#  ## I(AGE^2)  -0.0025520  0.0004962  -5.144  2.70e-07 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 6058.88,
#  ## null deviance= 8152.05,
#  ## nT= 5976, N= 664
#  ##
#  ## ( 7173 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6
#  ##
#  ## Average individual fixed effect= -0.969

## ---- eval = FALSE------------------------------------------------------------
#  apes_stat_bc <- get_APEs(stat_bc)
#  summary(apes_stat_bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.1015555  0.0171068  -5.937  2.91e-09 ***
#  ## KID2      -0.0585204  0.0156695  -3.735  0.000188 ***
#  ## KID3      -0.0185083  0.0131282  -1.410  0.158593
#  ## log(INCH) -0.0344416  0.0162926  -2.114  0.034520 *
#  ## AGE        0.0330430  0.0116438   2.838  0.004542 **
#  ## I(AGE^2)  -0.0004108  0.0001551  -2.648  0.008101 **
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## ---- eval = FALSE------------------------------------------------------------
#  library(data.table)
#  setDT(psid)
#  setkey(psid, ID, TIME)
#  psid[, LLFP := shift(LFP), by = ID]

## ---- eval = FALSE------------------------------------------------------------
#  dyn <- bife(LFP ~ LLFP + KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID, psid, "probit")
#  dyn_bc <- bias_corr(dyn, L = 1L)
#  summary(dyn_bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ LLFP + KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) |
#  ##     ID
#  ##
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## LLFP       1.0025730  0.0473072  21.193   < 2e-16 ***
#  ## KID1      -0.4741452  0.0679084  -6.982  2.91e-12 ***
#  ## KID2      -0.1958494  0.0625931  -3.129  0.001754 **
#  ## KID3      -0.0754092  0.0505115  -1.493  0.135461
#  ## log(INCH) -0.1947087  0.0621154  -3.135  0.001721 **
#  ## AGE        0.2009782  0.0477735   4.207  2.59e-05 ***
#  ## I(AGE^2)  -0.0024145  0.0006293  -3.837  0.000125 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 4774.57,
#  ## null deviance= 6549.14,
#  ## nT= 4792, N= 599
#  ##
#  ## ( 1461 observation(s) deleted due to missingness )
#  ## ( 6896 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6
#  ##
#  ## Average individual fixed effect= -1.939

## ---- eval = FALSE------------------------------------------------------------
#  apes_dyn_bc <- get_APEs(dyn_bc)
#  summary(apes_dyn_bc)

## ---- eval = FALSE------------------------------------------------------------
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## LLFP       0.1826348  0.0179590  10.170   < 2e-16 ***
#  ## KID1      -0.0752484  0.0191149  -3.937  8.26e-05 ***
#  ## KID2      -0.0310820  0.0176830  -1.758    0.0788 .
#  ## KID3      -0.0119677  0.0143610  -0.833    0.4047
#  ## log(INCH) -0.0309009  0.0170945  -1.808    0.0707 .
#  ## AGE        0.0318959  0.0132689   2.404    0.0162 *
#  ## I(AGE^2)  -0.0003832  0.0001743  -2.198    0.0279 *
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

