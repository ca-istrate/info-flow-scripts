# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function to compute the transfer of information from variable xj to variable xi (T) and the corresponding normalization (tau)
Multivariate time series
Liang (2021, 'Normalized multivariate time series causality analysis and causal graph reconstruction')

Compute Pearson correlation coefficient R

Compute the error in T21, tau21 and R based on the normal bootstrap with replacement (resampling of variables)

Created: 20/04/2021
Last updated: 22/04/2021

@author: David Docquier
"""

import logging

import pandas as pd
import numpy as np
from sklearn.utils import resample


logger = logging.getLogger("script")


def compute_liang_nvar(x, dt, n_iter, conf):
    # Function to compute absolute transfer of information from xj to xi (T)
    def compute_liang_index(detC, Deltajk, Ckdi, Cij, Cii):
        T = (1. / detC) * np.sum(Deltajk * Ckdi) * (
                    Cij / Cii)  # absolute rate of information flowing from xj to xi (nats per unit time) (equation (14))
        return T

    # Function to compute relative transfer of information from xj to xi (tau)
    def compute_liang_index_norm(detC, Deltaik, Ckdi, T_all, Tii, gii, Cii, Tji):
        selfcontrib = (1. / detC) * np.sum(Deltaik * Ckdi)  # contribution from itself (equation (15))
        transfer = np.sum(np.abs(T_all)) - np.abs(Tii)  # transfer contribution (equation (20))
        noise = 0.5 * gii / Cii  # noise contribution
        Z = np.abs(selfcontrib) + transfer + np.abs(noise)  # normalizer (equation (20))
        tau = 100. * Tji / Z  # relative rate of information flowing from xj to xi (%) (equation (19))
        return tau

    # Function to test significance (based on the confidence interval)
    def compute_sig(var, error, conf):
        if (var - conf * error < 0. and var + conf * error < 0.) or (
                var - conf * error > 0. and var + conf * error > 0.):
            sig = 1
        else:
            sig = 0
        z = abs(var / error)  # z score
        pval = np.exp(
            -0.717 * z - 0.416 * z ** 2.)  # p value for 95% confidence interval (https://www.bmj.com/content/343/bmj.d2304)
        return sig, pval

    # Dimensions
    nvar = np.size(x, 0)  # number of variables
    N = np.size(x, 1)  # length of the time series (number of observations)

    # Compute tendency dx
    k = 1  # k = 1 (or 2 for highly chaotic and densely sampled systems)
    dx = np.zeros((nvar, N))  # initialization of dx (to have the same number of time steps as x)
    for i in np.arange(nvar):
        dx[i, 0:N - k] = (x[i, k:N] - x[i, 0:N - k]) / (k * dt)  # Euler forward finite difference of x (equation (7))

    # Compute covariances and matrix determinant
    C = np.cov(x)  # covariance matrix
    dC = np.empty_like(C) * 0.
    for i in np.arange(nvar):
        for j in np.arange(nvar):
            dC[j, i] = (np.sum((x[j, :] - np.nanmean(x[j, :])) * (dx[i, :] - np.nanmean(dx[i, :])))) / (
                        N - 1.)  # covariance between x and dx
    detC = np.linalg.det(C)  # matrix determinant

    # Compute cofactors
    Delta = np.linalg.inv(C).T * detC  # cofactor matrix (https://en.wikipedia.org/wiki/Minor_(linear_algebra))

    # Compute absolute transfer of information (T) and correlation coefficient
    T = np.zeros((nvar, nvar))
    R = np.zeros((nvar, nvar))
    for i in np.arange(nvar):
        for j in np.arange(nvar):
            T[j, i] = compute_liang_index(detC, Delta[j, :], dC[:, i], C[i, j], C[
                i, i])  # compute T (transfer of information from xj to xi) and create matrix
            R[j, i] = C[i, j] / np.sqrt(
                C[i, i] * C[j, j])  # compute correlation coefficient and create correlation matrix

    # Compute noise terms
    g = np.zeros(nvar)
    for i in np.arange(nvar):
        a1k = np.dot(np.linalg.inv(C), dC[:,
                                       i])  # compute a1k coefficients based on matrix-vector product (see beginning of page 4 in Liang (2014))
        f1 = np.nanmean(dx[i, :])
        for k in np.arange(nvar):
            f1 = f1 - a1k[k] * np.nanmean(x[k, :])
        R1 = dx[i, :] - f1
        for k in np.arange(nvar):
            R1 = R1 - a1k[k] * x[k, :]
        Q1 = np.sum(R1 ** 2.)
        g[i] = Q1 * dt / N  # equation (10)

    # Compute relative transfer of information (tau)
    tau = np.zeros((nvar, nvar))
    for i in np.arange(nvar):
        for j in np.arange(nvar):
            tau[j, i] = compute_liang_index_norm(detC, Delta[i, :], dC[:, i], T[:, i], T[i, i], g[i], C[i, i],
                                                 T[j, i])  # compute tau and create matrix

    # Compute error in Tji and tauji using bootstrap with replacement
    boot_T = np.zeros((n_iter, nvar, nvar))
    boot_tau = np.zeros((n_iter, nvar, nvar))
    boot_R = np.zeros((n_iter, nvar, nvar))

    for it in np.arange(n_iter):  # loop over realizations

        # Resample x and dx
        index = np.arange(N)
        boot_index = resample(index, replace=True)
        boot_x = np.zeros((nvar, N))
        boot_dx = np.zeros((nvar, N))
        for t in np.arange(N):
            boot_x[:, t] = x[:, boot_index[t]]
            boot_dx[:, t] = dx[:, boot_index[t]]

        # Compute covariances and matrix determinant based on resampled variables
        boot_C = np.cov(boot_x)
        boot_dC = np.empty_like(boot_C) * 0.
        for i in np.arange(nvar):
            for j in np.arange(nvar):
                boot_dC[j, i] = (np.sum(
                    (boot_x[j, :] - np.nanmean(boot_x[j, :])) * (boot_dx[i, :] - np.nanmean(boot_dx[i, :])))) / (N - 1.)
        boot_detC = np.linalg.det(boot_C)

        # Compute cofactors based on resampled variables
        boot_Delta = np.linalg.inv(boot_C).T * boot_detC

        # Compute absolute transfer of information (T) and correlation coefficient based on resampled variables
        for i in np.arange(nvar):
            for j in np.arange(nvar):
                boot_T[it, j, i] = compute_liang_index(boot_detC, boot_Delta[j, :], boot_dC[:, i], boot_C[i, j],
                                                       boot_C[i, i])
                boot_R[it, j, i] = boot_C[i, j] / np.sqrt(boot_C[i, i] * boot_C[j, j])

        # Compute noise terms based on resampled variables
        boot_g = np.zeros(nvar)
        for i in np.arange(nvar):
            a1k = np.dot(np.linalg.inv(boot_C), boot_dC[:, i])
            f1 = np.nanmean(boot_dx[i, :])
            for k in np.arange(nvar):
                f1 = f1 - a1k[k] * np.nanmean(boot_x[k, :])
            R1 = boot_dx[i, :] - f1
            for k in np.arange(nvar):
                R1 = R1 - a1k[k] * boot_x[k, :]
            Q1 = np.sum(R1 ** 2.)
            boot_g[i] = Q1 * dt / N  # equation (10)

        # Compute relative transfer of information (tau) based on resampled variables
        for i in np.arange(nvar):
            for j in np.arange(nvar):
                boot_tau[it, j, i] = compute_liang_index_norm(boot_detC, boot_Delta[i, :], boot_dC[:, i],
                                                              boot_T[it, :, i], boot_T[it, i, i], boot_g[i],
                                                              boot_C[i, i], boot_T[it, j, i])

    # Compute error in T21, tau21 and R (standard deviation of boostraped values)
    error_T = np.nanstd(boot_T, axis=0)
    error_tau = np.nanstd(boot_tau, axis=0)
    error_R = np.nanstd(boot_R, axis=0)

    # Calculate the matrices for significance and pvalue
    sig_T = np.zeros((nvar, nvar))
    sig_tau = np.zeros((nvar, nvar))
    sig_R = np.zeros((nvar, nvar))
    pval_T = np.zeros((nvar, nvar))
    pval_tau = np.zeros((nvar, nvar))
    pval_R = np.zeros((nvar, nvar))

    for i in np.arange(nvar):
        for j in np.arange(nvar):
            sig_T[i, j], pval_T[i, j] = compute_sig(T[i, j], error_T[i, j], conf)
            sig_tau[i, j], pval_tau[i, j] = compute_sig(tau[i, j], error_tau[i, j], conf)
            sig_R[i, j], pval_R[i, j] = compute_sig(R[i, j], error_R[i, j], conf)

    return pd.DataFrame({
        "T11_doc": T[0,0],
        "T12_doc": T[0,1],
        "T21_doc": T[1,0],
        "T22_doc": T[1,1],
        "tau11": tau[0,0],
        "tau12": tau[0,1],
        "tau21": tau[1,0],
        "tau22": tau[1,1],
        "R11": R[0,0],
        "R12": R[0,1],
        "R21": R[1,0],
        "R22": R[1,1],
        "error_T11_doc": error_T[0,0],
        "error_T12_doc": error_T[0,1],
        "error_T21_doc": error_T[1,0],
        "error_T22_doc": error_T[1,1],
        "error_tau11": error_tau[0,0],
        "error_tau12": error_tau[0,1],
        "error_tau21": error_tau[1,0],
        "error_tau22": error_tau[1,1],
        "error_R11": error_R[0,0],
        "error_R12": error_R[0,1],
        "error_R21": error_R[1,0],
        "error_R22": error_R[1,1],
        "sig_T11_doc": sig_T[0,0],
        "sig_T12_doc": sig_T[0,1],
        "sig_T21_doc": sig_T[1,0],
        "sig_T22_doc": sig_T[1,1],
        "sig_tau11": sig_tau[0,0],
        "sig_tau12": sig_tau[0,1],
        "sig_tau21": sig_tau[1,0],
        "sig_tau22": sig_tau[1,1],
        "sig_R11": sig_R[0,0],
        "sig_R12": sig_R[0,1],
        "sig_R21": sig_R[1,0],
        "sig_R22": sig_R[1,1],
        "pval_T11_doc": pval_T[0,0],
        "pval_T12_doc": pval_T[0,1],
        "pval_T21_doc": pval_T[1,0],
        "pval_T22_doc": pval_T[1,1],
        "pval_tau11": pval_tau[0,0],
        "pval_tau12": pval_tau[0,1],
        "pval_tau21": pval_tau[1,0],
        "pval_tau22": pval_tau[1,1],
        "pval_R11": pval_R[0,0],
        "pval_R12": pval_R[0,1],
        "pval_R21": pval_R[1,0],
        "pval_R22": pval_R[1,1]
    }, index=[0]), (T, tau, R, error_T, error_tau, error_R, sig_T, sig_tau, sig_R, pval_T, pval_tau, pval_R)


def docquier(data, args):
    logger.debug("Running bootstraping method")
    x = np.array([data.x, data.y])
    return compute_liang_nvar(x, 1, args.bootstrap_iter, args.bootstrap_conf_val)
