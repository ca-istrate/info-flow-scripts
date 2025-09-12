import logging
import os

import pandas as pd
import numpy as np
from scipy.stats import norm

from common.docquier import docquier
from common.pairing import Pairing


logger = logging.getLogger("script")


def default_output_file(args) -> None:
    x, _ = os.path.splitext(os.path.basename(args.x))
    y, _ = os.path.splitext(os.path.basename(args.y))
    vars(args)["output_file"] = f"{x}-{y}-{args.start}-{args.end}-{args.shift}-analysis.out"
    return vars(args)["output_file"]


def calculate_p_value(T, x):
    if T < 0:
        return norm.cdf(x)
    else:
        return 1 - norm.cdf(x)



def compute_kalpha(alpha):
    if 1 - alpha == 0.95:
        k_alpha = 1.96
    elif 1 - alpha == 0.90:
        k_alpha = 1.65
    elif 1 - alpha == 0.85:
        k_alpha = 1.44
    elif 1 - alpha == 0.80:
        k_alpha = 1.28
    elif 1 - alpha == 0.75:
        k_alpha = 1.15
    else:
        logger.error("Eroare: te rog introdu alta valoare pentru alpha!")
        logger.error("Valori posibile: 0.05, 0.10, 0.15, 0.20, 0.25.")
        exit(1)

    return k_alpha


def hadamard_mean(N, A, B):
    return sum(np.multiply(np.array(A), np.array(B))) / N


def variation(X):
    dX = []
    for i in range(len(X) - 1):
          dX.append(X[i + 1] - X[i])
    return dX


def theoretical(data, args):
    logger.debug("Running theoretical calculation.")

    alpha = args.alpha 
    logger.debug(f"\n \u03B1 \t(alpha) \t- nivelul de semnificatie:   {alpha}")
    logger.debug(f" 1 - \u03B1 \t(1 - alpha) \t- coeficientul de incredere: {1 - alpha}\n")

    k_alpha = compute_kalpha(alpha)

    N, X1, X2 = len(data.x), data.x, data.y

    dX1 = variation(X1)
    dX2 = variation(X2)

    # Compute expected values (or means):
    E_X1 = sum(X1[0 : N - 1]) / (N - 1)
    E_X2 = sum(X2[0 : N - 1]) / (N - 1)

    E_dX1 = sum(dX1) / (N - 1)
    E_dX2 = sum(dX2) / (N - 1)

    E_X1X1 = hadamard_mean(N - 1, X1[0 : N - 1], X1[0 : N - 1])
    E_X2X2 = hadamard_mean(N - 1, X2[0 : N - 1], X2[0 : N - 1])
    E_X1X2 = hadamard_mean(N - 1, X1[0 : N - 1], X2[0 : N - 1])

    E_X1dX1 = hadamard_mean(N - 1, X1[0 : N - 1], dX1)
    E_X1dX2 = hadamard_mean(N - 1, X1[0 : N - 1], dX2)
    E_X2dX1 = hadamard_mean(N - 1, X2[0 : N - 1], dX1)
    E_X2dX2 = hadamard_mean(N - 1, X2[0 : N - 1], dX2)

    # Compute the covariance coefficients:
    C11 = E_X1X1 - E_X1 * E_X1
    C12 = E_X1X2 - E_X1 * E_X2
    C22 = E_X2X2 - E_X2 * E_X2

    C1d1 = E_X1dX1 - E_X1 * E_dX1
    C1d2 = E_X1dX2 - E_X1 * E_dX2
    C2d1 = E_X2dX1 - E_X2 * E_dX1
    C2d2 = E_X2dX2 - E_X2 * E_dX2

    # Compute information flow indicators:
    T1_2 = (C22 * C12 * C1d2 - C12 * C12 * C2d2) / (C22 * C22 * C11 - C22 * C12 * C12)
    T2_1 = (C11 * C12 * C2d1 - C12 * C12 * C1d1) / (C11 * C11 * C22 - C11 * C12 * C12)

    # Compute trusted intervals:
    a11 = (C22 * C1d1 - C12 * C2d1) / (C11 * C22 - C12 * C12)
    a12 = (C11 * C2d1 - C12 * C1d1) / (C11 * C22 - C12 * C12)
    f1  = E_dX1 - a11 * E_X1 - a12 * E_X2
    Q1  = 0
    for i in range(0, N - 1):
        tmp = dX1[i] - (f1 + a11 * X1[i] + a12 * X2[i])
        Q1 += tmp * tmp
    b1 = np.sqrt(Q1 / (N - 1))

    i11 = N / (b1 * b1)
    i12 = i11 * E_X1
    i13 = i11 * E_X2
    i14 = 0

    i21 = i12
    i22 = i11 * E_X1X1
    i23 = i11 * E_X1X2
    i24 = 0

    i31 = i13
    i32 = i23
    i33 = i11 * E_X2X2
    i34 = 0

    i41 = i14
    i42 = i24
    i43 = i34
    i44 = 2 * i11

    I_Fischer = np.array([[i11, i12, i13, i14], \
                  [i21, i22, i23, i24], \
                  [i31, i32, i33, i34], \
                  [i41, i42, i43, i44]])
    I_Fischer_inv = np.linalg.inv(I_Fischer)

    sigma_a12_patrat = I_Fischer_inv[2][2]

    dist_interval_2_1 = np.abs(k_alpha * C12 / C11 * np.sqrt(sigma_a12_patrat))

    low_2_1   = T2_1 - dist_interval_2_1
    upper_2_1 = T2_1 + dist_interval_2_1

    ##################################################################################

    a21 = (C22 * C1d2 - C12 * C2d2) / (C11 * C22 - C12 * C12)
    a22 = (C11 * C2d2 - C12 * C1d2) / (C11 * C22 - C12 * C12)
    f2  = E_dX2 - a21 * E_X1 - a22 * E_X2
    Q2  = 0
    for i in range(0, N - 1):
        tmp = dX2[i] - (f2 + a21 * X1[i] + a22 * X2[i])
        Q2 += tmp * tmp
    b2 = np.sqrt(Q2 / (N - 1))


    i11 = N / (b2 * b2)
    i12 = i11 * E_X1
    i13 = i11 * E_X2
    i14 = 0

    i21 = i12
    i22 = i11 * E_X1X1
    i23 = i11 * E_X1X2
    i24 = 0

    i31 = i13
    i32 = i23
    i33 = i11 * E_X2X2
    i34 = 0

    i41 = i14
    i42 = i24
    i43 = i34
    i44 = 2 * i11

    I_Fischer = np.array([[i11, i12, i13, i14], \
                  [i21, i22, i23, i24], \
                  [i31, i32, i33, i34], \
                  [i41, i42, i43, i44]])
    I_Fischer_inv = np.linalg.inv(I_Fischer)

    sigma_a21_patrat = I_Fischer_inv[1][1]

    dist_interval_1_2 = np.abs(k_alpha * C12 / C22 * np.sqrt(sigma_a21_patrat))

    low_1_2   = T1_2 - dist_interval_1_2
    upper_1_2 = T1_2 + dist_interval_1_2

    # Analyze data:
    T1_2_x = T1_2 / (C12 / C22 * np.sqrt(sigma_a21_patrat))
    T1_2_pval = calculate_p_value(T1_2, T1_2_x)

    T2_1_x = T2_1 / (C12 / C11 * np.sqrt(sigma_a12_patrat))
    T2_1_pval = calculate_p_value(T2_1, T2_1_x)

    return pd.DataFrame({
        "T21_theo": T2_1,
        "T21_int_dist": dist_interval_2_1,
        "T21_int_lower": low_2_1,
        "T21_int_upper": upper_2_1,
        "T21_x": T2_1_x,
        "T21_pval": T2_1_pval,
        "T12_theo": T1_2,
        "T12_int_dist": dist_interval_1_2,
        "T12_int_lower": low_1_2,
        "T12_int_upper": upper_1_2,
        "T12_x": T1_2_x,
        "T12_pval": T1_2_pval,
    }, index=[0])

def single_row_output(args, data, row, matrices):
    def significance_text(typ, T, x, p_value, alpha):
        logger.info(f"Test semnificatie T {typ}")
        logger.info(f"Valoarea lui x = {x}")
        if T < 0:
            if alpha >= p_value:
                logger.info(f"T {typ} este semnificativ < 0 corespunzator nivelului de " 
                f"semnificatie \u03B1 (alpha) = {alpha} si valorii p_value = " 
                f"{p_value:.20f} .")
            else:
                logger.info(f" T {typ} nu difera semnificativ fata de 0. p_value = "
                f"{p_value:.20f}  > {alpha} (alpha)")
        else: # T >= 0
            p_value = 1 - norm.cdf(x)
    
            if alpha >= p_value:
                logger.info(f" T {typ} este semnificativ > 0 corespunzator nivelului de " 
                f"semnificatie \u03B1 (alpha) = {alpha} si valorii p_value = "
                f"{p_value:.20f} .")
            else:
                logger.info(f" T {typ} nu difera semnificativ fata de 0. p_value = "
                f"{p_value:.20f}  > {alpha} (alpha)")
        logger.info("")

        _1st = "1"
        _2nd = "2"

        if typ == "2->1":
            _1st = "2"
            _2nd = "1"
    
        if 0.10 <= p_value < 0.15:
            logger.info(f"\tExista o cauzalitate slaba intre factorul X{_1st} si factorul X{_2nd}. (0.10 <= p_value = {p_value:.20f} < 0.15)")
        elif 0.05 <= p_value < 0.10:
            logger.info(f"\tExista o cauzalitate normala intre factorul X{_1st} si factorul X{_2nd}. (0.05 <= p_value = {p_value:.20f} < 0.10)")
        elif p_value < 0.05:
            logger.info(f"\tExista o cauzalitate puternica intre factorul X{_1st} si factorul X{_2nd}. (p_value = {p_value:.20f} < 0.05)")
        logger.info("")

    logger.info(f"Rezultate teoretice\n")
    logger.info(f" T 1->2 : {row.loc[0, 'T12_theo']} \u00b1 {row.loc[0, 'T12_int_dist']}")
    logger.info(f"\tInterval de incredere: [{row.loc[0, 'T12_int_lower']}; {row.loc[0, 'T12_int_upper']}]\n")

    significance_text("1->2", row.loc[0, "T12_theo"], row.loc[0, "T12_x"], row.loc[0, "T12_pval"], args.alpha)

    logger.info(f" T 2->1 : {row.loc[0, 'T21_theo']} \u00b1 {row.loc[0, 'T21_int_dist']}")
    logger.info(f"\tInterval de incredere: [{row.loc[0, 'T21_int_lower']}; {row.loc[0, 'T21_int_upper']}]\n")

    significance_text("2->1", row.loc[0, "T21_theo"], row.loc[0, "T21_x"], row.loc[0, "T21_pval"], args.alpha)

    logger.info(f"Rezultate bootstraping\n")

    # (T, tau, R, error_T, error_tau, error_R, sig_T, sig_tau, sig_R, pval_T, pval_tau, pval_R)
    logger.info("Matricea T")
    logger.info(matrices[0])
    logger.info("\n")
    logger.info("Matricea tau")
    logger.info(matrices[1])
    logger.info("\n")
    logger.info("Matricea R")
    logger.info(matrices[2])
    logger.info("\n")
    logger.info("Matricea erorii lui T")
    logger.info(matrices[3])
    logger.info("\n")
    logger.info("Matricea erorii lui tau")
    logger.info(matrices[4])
    logger.info("\n")
    logger.info("Matricea erorii lui R")
    logger.info(matrices[5])
    logger.info("\n")
    logger.info("Matricea significance lui T")
    logger.info(matrices[6])
    logger.info("\n")
    logger.info("Matricea significance lui tau")
    logger.info(matrices[7])
    logger.info("\n")
    logger.info("Matricea significance lui R")
    logger.info(matrices[8])
    logger.info("\n")
    logger.info("Matricea pvalue lui T")
    logger.info(matrices[9])
    logger.info("\n")
    logger.info("Matricea pvalue lui tau")
    logger.info(matrices[10])
    logger.info("\n")
    logger.info("Matricea pvalue lui R")
    logger.info(matrices[11])
    logger.info("\n")


def analyze(args, data):
    row = theoretical(data, args)
    doc_row, matrices = docquier(data, args)

    header = pd.DataFrame({
        "x_time_start": data.x_time[0],
        "x_time_end": data.x_time[-1],
        "y_time_start": data.y_time[0],
        "y_time_end": data.y_time[-1]
    }, index=[0])

    final_row = pd.concat([header, row, doc_row], axis=1)

    return final_row, matrices


def main(args) -> None:
    logger.info("Running info-flow calculation.")
    
    logger.debug(f"Input arguments: {args}")
    
    data = Pairing.from_files(args.x, args.y, start=args.start, end=args.end, shift=args.shift)
    single_row, matrices = analyze(args, data)

    single_row_output(args, data, single_row, matrices)
