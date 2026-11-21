# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Published ephemerides and flare peak times of ASASSN-14ko (fixed before any new data were read).

Cycle number n counts outbursts from the first ASAS-SN peak (JD 2456861.2, July 2014),
the indexing used by Payne et al. 2021 (ApJ 910, 125; arXiv:2009.03321) and, as checked
in test_ko.py, by Payne et al. 2023 (ApJ 951, 134; arXiv:2206.11278, posted June 2022),
whose timing fit and forecasts are called "PAYNE22" here after the arXiv date.
The September 2020 peak is from Payne et al. 2022 (ApJ 926, 142; arXiv:2104.06414).
"""

TARGET = {
    "name": "ASASSN-14ko",
    "host": "ESO 253-G003",
    "ra_deg": 81.32533,     # 05h25m18.08s
    "dec_deg": -46.00583,   # -46d00m21.0s
}

# Payne et al. 2023 (arXiv:2206.11278), eq. 1 and the fit quoted in the abstract.
PAYNE22 = {
    "t0_mjd": 2456852.0 - 2400000.5,   # JD 2456852 (+10/-11)
    "t0_err": 10.5,
    "P0": 115.2,                       # days (+1.3/-1.2)
    "P0_err": 1.25,
    "Pdot": -0.0026,                   # dimensionless dP/dt
    "Pdot_err": 0.0006,
}

# Payne et al. 2023 (arXiv:2206.11278) forecasts of the next five optical peaks (MJD), n = 25..29.
PAYNE22_FORECASTS = {25: 59639.7, 26: 59747.5, 27: 59855.9, 28: 59962.1, 29: 60069.0}


def payne_time(n, t0=None, P0=None, Pdot=None):
    """Peak time of cycle n for a constant period derivative (Payne et al. 2023, eq. 1)."""
    t0 = PAYNE22["t0_mjd"] if t0 is None else t0
    P0 = PAYNE22["P0"] if P0 is None else P0
    Pdot = PAYNE22["Pdot"] if Pdot is None else Pdot
    return t0 + n * P0 + 0.5 * n * n * P0 * Pdot + (n ** 3) * P0 * Pdot * Pdot / 6.0


def _jd(jd):
    return round(jd - 2400000.5, 3)


# Published optical peak times: (n, MJD, 1-sigma error in days, instrument, reference).
# Payne+2021 Table 3 errors are asymmetric; the mean of the two sides is used.
PUBLISHED_PEAKS = [
    (0, _jd(2456861.2), 6.6, "ASAS-SN V", "Payne+2021"),
    (1, _jd(2456971.1), 9.1, "ASAS-SN V", "Payne+2021"),
    (2, _jd(2457083.8), 3.6, "ASAS-SN V", "Payne+2021"),
    (4, _jd(2457311.8), 3.8, "ASAS-SN V", "Payne+2021"),
    (5, _jd(2457419.4), 7.7, "ASAS-SN V", "Payne+2021"),
    (7, _jd(2457651.9), 3.0, "ASAS-SN V", "Payne+2021"),
    (8, _jd(2457761.4), 6.4, "ASAS-SN V", "Payne+2021"),
    (9, _jd(2457879.4), 5.15, "ASAS-SN V/g", "Payne+2021"),
    (10, _jd(2457988.4), 4.1, "ASAS-SN V/g", "Payne+2021"),
    (11, _jd(2458096.2), 0.9, "ASAS-SN V/g", "Payne+2021"),
    (12, _jd(2458210.4), 0.6, "ASAS-SN V/g", "Payne+2021"),
    (14, _jd(2458435.7), 1.2, "ASAS-SN g", "Payne+2021"),
    (15, _jd(2458547.1), 0.25, "ASAS-SN g", "Payne+2021"),
    (17, _jd(2458768.5), 1.0, "ASAS-SN g", "Payne+2021"),
    (18, _jd(2458878.4), 0.45, "ASAS-SN g", "Payne+2021"),
    (19, _jd(2458987.5), 1.6, "ASAS-SN g", "Payne+2021"),
    (20, 59097.76, 1.0, "ASAS-SN g", "Payne+2022 ApJ 926 142 (Sept 2020; error assumed 1.0 d)"),
    (21, 59205.7, 0.19, "ASAS-SN g", "Payne+2023 arXiv:2206.11278"),
    (22, 59315.2, 0.55, "ASAS-SN g", "Payne+2023 arXiv:2206.11278"),
    (23, 59433.3, 4.5, "ASAS-SN g", "Payne+2023 (weather gap; Swift B gives 59423.9)"),
    (24, 59532.2, 0.9, "ASAS-SN g", "Payne+2023 arXiv:2206.11278"),
]

# TESS peak times from Payne+2023 (arXiv:2206.11278, Table 2), used only to validate our TESS pipeline.
PUBLISHED_TESS_PEAKS = [
    (14, _jd(2458435.50), 0.12, 4),    # sector 4
    (21, _jd(2459206.13), 0.10, 33),   # sector 33
]

# Gate G1: cycles with a published ASAS-SN g-band peak to better than 1 day.
GATE_CYCLES = (15, 17, 18, 21, 22, 24)

FIRST_NEW_CYCLE = 25   # first cycle after the last published peak
