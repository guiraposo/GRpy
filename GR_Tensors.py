{
  "nbformat": 4,
  "nbformat_minor": 0,
  "metadata": {
    "colab": {
      "provenance": [],
      "authorship_tag": "ABX9TyObTDDIpq3V1i5BLQ0H2WX6",
      "include_colab_link": true
    },
    "kernelspec": {
      "name": "python3",
      "display_name": "Python 3"
    },
    "language_info": {
      "name": "python"
    }
  },
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "view-in-github",
        "colab_type": "text"
      },
      "source": [
        "<a href=\"https://colab.research.google.com/github/guiraposo/GRpy/blob/main/GR_Tensors.py\" target=\"_parent\"><img src=\"https://colab.research.google.com/assets/colab-badge.svg\" alt=\"Open In Colab\"/></a>"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "CmmhRtrHzMJU"
      },
      "outputs": [],
      "source": [
        "def christoffel_symbols(metric, coords):\n",
        "    dims = metric.shape[0]\n",
        "    christoffel = [[[0 for _ in range(dims)] for _ in range(dims)] for _ in range(dims)]\n",
        "    inv_metric = metric.inv()\n",
        "    for k in range(dims):\n",
        "        for i in range(dims):\n",
        "            for j in range(dims):\n",
        "                sum_term = 0\n",
        "                for l in range(dims):\n",
        "                    partial_1 = sp.simplify(sp.diff(metric[j, l], coords[i]))\n",
        "                    partial_2 = sp.simplify(sp.diff(metric[i, l], coords[j]))\n",
        "                    partial_3 = sp.simplify(sp.diff(metric[i, j], coords[l]))\n",
        "                    sum_term += inv_metric[k, l] * (partial_1 + partial_2 - partial_3)\n",
        "                christoffel[k][i][j] = sp.simplify(sum_term / 2)\n",
        "    return christoffel\n",
        "def riemann_tensor(christoffel, coords):\n",
        "    dims = len(coords)\n",
        "    R = [[[sp.zeros(dims, dims) for _ in range(dims)] for _ in range(dims)] for _ in range(dims)]\n",
        "    for rho in range(dims):\n",
        "        for sigma in range(dims):\n",
        "            for mu in range(dims):\n",
        "                for nu in range(dims):\n",
        "                    term1 = sp.simplify(sp.diff(christoffel[rho][nu][sigma], coords[mu]))\n",
        "                    term2 = sp.simplify(sp.diff(christoffel[rho][mu][sigma], coords[nu]))\n",
        "                    term3 = sp.simplify(sum(christoffel[rho][mu][lambda1] * christoffel[lambda1][nu][sigma] for lambda1 in range(dims)))\n",
        "                    term4 = sp.simplify(sum(christoffel[rho][nu][lambda1] * christoffel[lambda1][mu][sigma] for lambda1 in range(dims)))\n",
        "                    R[rho][sigma][mu][nu] = sp.simplify(term1 - term2 + term3 - term4)\n",
        "    return R\n",
        "def ricci_tensor(riemann, coords):\n",
        "    dims = len(coords)\n",
        "    ricci = sp.zeros(dims, dims)\n",
        "    for mu in range(dims):\n",
        "        for nu in range(dims):\n",
        "            ricci[mu, nu] = sp.simplify(sum(riemann[lambda1][mu][lambda1][nu] for lambda1 in range(dims)))\n",
        "    return ricci\n",
        "def ricci_scalar(ricci, metric):\n",
        "    inv_metric = metric.inv()\n",
        "    dims = metric.shape[0]\n",
        "    scalar = 0\n",
        "    for mu in range(dims):\n",
        "        for nu in range(dims):\n",
        "            scalar += sp.simplify(inv_metric[mu, nu] * ricci[mu, nu])\n",
        "    return sp.simplify(scalar)\n",
        "def einstein_tensor(ricci_tensor, ricci_scalar, metric):\n",
        "    dims = metric.shape[0]\n",
        "    einstein = sp.zeros(dims, dims)\n",
        "    for mu in range(dims):\n",
        "        for nu in range(dims):\n",
        "            einstein[mu, nu] = ricci_tensor[mu, nu] - 1/2 * metric[mu, nu] * ricci_scalar\n",
        "    return sp.simplify(einstein)\n",
        "def stress_energy_tensor_perfect_fluid(rho, p, u, metric):\n",
        "    dims = metric.shape[0]\n",
        "    stress_energy = sp.zeros(dims, dims)\n",
        "    for mu in range(dims):\n",
        "        for nu in range(dims):\n",
        "            stress_energy[mu, nu] = (rho + p) * u[mu] * u[nu] + p * metric[mu, nu]\n",
        "    return sp.simplify(stress_energy)\n",
        "def einsteins_equations(einstein_tensor, stress_energy_tensor_perfect_fluid, metric):\n",
        "    dims = metric.shape[0]\n",
        "    einstein_eqs = sp.zeros(dims, dims)\n",
        "    for mu in range(dims):\n",
        "        for nu in range(dims):\n",
        "            einstein_eqs[mu, nu] = einstein_tensor[mu, nu] - 8 * sp.pi * stress_energy_tensor_perfect_fluid[mu,nu]\n",
        "    # Einstein's Equations: G_{\\mu\\nu} - 8πT_{\\mu\\nu} = 0 (In geometric units)\n",
        "    # We return a tensor that should be all zeros if the equations are satisfied\n",
        "    return sp.simplify(einstein_eqs)\n"
      ]
    }
  ]
}