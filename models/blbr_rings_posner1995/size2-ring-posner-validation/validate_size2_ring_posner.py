#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import shutil
import struct
import subprocess
import tempfile
import zlib
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


ROOT = Path(__file__).resolve().parent
BNGL_PATH = ROOT / "size2_ring_posner_validation.bngl"
PNG_PATH = ROOT / "size2_ring_posner_validation.png"
JSON_PATH = ROOT / "size2_ring_posner_validation_summary.json"


BNG2_CANDIDATES = [
    Path("/Users/wish/Simulations/BioNetGen-2.9.3/BNG2.pl"),
]

NFSIM_CANDIDATES = [
    Path("/Users/wish/Code/nfsim/build/NFsim"),
    Path("/Users/wish/Simulations/BioNetGen-2.9.3/bin/NFsim"),
]


RT = 30
ST = 2 * RT
KF = 1.6605390671738464e-05
KPX = 0.0016666666666666668
KOFF = 0.01
JP = 1.0
JM = KOFF
J2 = JP / JM

T_END = 600.0
N_STEPS = 240
OUTPUT_DT = T_END / N_STEPS
T_GRID = [i * OUTPUT_DT for i in range(N_STEPS + 1)]


DOSES = (
    {"name": "low", "title": "LOW DOSE LT 150", "lt": 150, "seed": 1101},
    {"name": "peak", "title": "PEAK DOSE LT 331", "lt": 331, "seed": 1201},
    {"name": "high", "title": "HIGH DOSE LT 450", "lt": 450, "seed": 1301},
)


FONT_5X7 = {
    " ": ("00000", "00000", "00000", "00000", "00000", "00000", "00000"),
    "0": ("01110", "10001", "10011", "10101", "11001", "10001", "01110"),
    "1": ("00100", "01100", "00100", "00100", "00100", "00100", "01110"),
    "2": ("01110", "10001", "00001", "00010", "00100", "01000", "11111"),
    "3": ("11110", "00001", "00001", "01110", "00001", "00001", "11110"),
    "4": ("00010", "00110", "01010", "10010", "11111", "00010", "00010"),
    "5": ("11111", "10000", "10000", "11110", "00001", "00001", "11110"),
    "6": ("00110", "01000", "10000", "11110", "10001", "10001", "01110"),
    "7": ("11111", "00001", "00010", "00100", "01000", "01000", "01000"),
    "8": ("01110", "10001", "10001", "01110", "10001", "10001", "01110"),
    "9": ("01110", "10001", "10001", "01111", "00001", "00010", "11100"),
    "A": ("01110", "10001", "10001", "11111", "10001", "10001", "10001"),
    "B": ("11110", "10001", "10001", "11110", "10001", "10001", "11110"),
    "C": ("01110", "10001", "10000", "10000", "10000", "10001", "01110"),
    "D": ("11100", "10010", "10001", "10001", "10001", "10010", "11100"),
    "E": ("11111", "10000", "10000", "11110", "10000", "10000", "11111"),
    "F": ("11111", "10000", "10000", "11110", "10000", "10000", "10000"),
    "G": ("01110", "10001", "10000", "10111", "10001", "10001", "01110"),
    "H": ("10001", "10001", "10001", "11111", "10001", "10001", "10001"),
    "I": ("01110", "00100", "00100", "00100", "00100", "00100", "01110"),
    "J": ("00001", "00001", "00001", "00001", "10001", "10001", "01110"),
    "K": ("10001", "10010", "10100", "11000", "10100", "10010", "10001"),
    "L": ("10000", "10000", "10000", "10000", "10000", "10000", "11111"),
    "M": ("10001", "11011", "10101", "10101", "10001", "10001", "10001"),
    "N": ("10001", "10001", "11001", "10101", "10011", "10001", "10001"),
    "O": ("01110", "10001", "10001", "10001", "10001", "10001", "01110"),
    "P": ("11110", "10001", "10001", "11110", "10000", "10000", "10000"),
    "Q": ("01110", "10001", "10001", "10001", "10101", "10010", "01101"),
    "R": ("11110", "10001", "10001", "11110", "10100", "10010", "10001"),
    "S": ("01111", "10000", "10000", "01110", "00001", "00001", "11110"),
    "T": ("11111", "00100", "00100", "00100", "00100", "00100", "00100"),
    "U": ("10001", "10001", "10001", "10001", "10001", "10001", "01110"),
    "V": ("10001", "10001", "10001", "10001", "10001", "01010", "00100"),
    "W": ("10001", "10001", "10001", "10101", "10101", "10101", "01010"),
    "X": ("10001", "10001", "01010", "00100", "01010", "10001", "10001"),
    "Y": ("10001", "10001", "01010", "00100", "00100", "00100", "00100"),
    "Z": ("11111", "00001", "00010", "00100", "01000", "10000", "11111"),
}


def locate_executable(candidates: Sequence[Path], fallback_name: str) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    found = shutil.which(fallback_name)
    if found:
        return Path(found)
    raise FileNotFoundError(f"Could not find executable: {fallback_name}")


def run_checked(cmd: Sequence[str], cwd: Path, quiet: bool = True) -> None:
    stdout = subprocess.DEVNULL if quiet else None
    stderr = subprocess.STDOUT if quiet else None
    subprocess.run(list(cmd), cwd=cwd, check=True, stdout=stdout, stderr=stderr)


def parse_gdat(path: Path) -> Dict[str, List[float]]:
    with path.open() as handle:
        header = handle.readline().strip().lstrip("#").split()
        data = {name: [] for name in header}
        for line in handle:
            if not line.strip():
                continue
            parts = line.split()
            for name, value in zip(header, parts):
                data[name].append(float(value))
    return data


def normalize_series(raw: Dict[str, List[float]]) -> Dict[str, List[float]]:
    def pick(*names: str) -> List[float]:
        for name in names:
            if name in raw:
                return raw[name]
        raise KeyError(f"None of {names} found in gdat columns: {sorted(raw)}")

    return {
        "time": pick("time"),
        "Y1": pick("Y1"),
        "Y2": pick("Y2()", "Y2"),
        "Z2": pick("Z2()", "Z2"),
    }


def mean_trajectories(runs: Sequence[Dict[str, List[float]]]) -> Dict[str, List[float]]:
    keys = runs[0].keys()
    mean = {}
    for key in keys:
        cols = zip(*(run[key] for run in runs))
        mean[key] = [sum(values) / len(runs) for values in cols]
    return mean


def cjn_equilibrium(j: int, n: int, s_free: float, y1: float, y2: float) -> float:
    sc = s_free + y1 + 2.0 * y2
    if sc <= 0.0:
        return 0.0
    return (
        0.5
        * sc
        * math.comb(2, j)
        * (y1 / sc) ** j
        * (s_free / sc) ** (2 - j)
        * (2.0 * y2 / sc) ** (n - 1)
    )


def equilibrium_residual(s_free: float, l_free: float, lt_total: int) -> Tuple[float, float]:
    if s_free <= 0.0 or l_free <= 0.0:
        return (1.0e6, 1.0e6)
    k1_eq = KF / KOFF
    k2_eq = KPX / KOFF
    y1 = 2.0 * k1_eq * l_free * s_free
    y2 = 0.5 * k2_eq * s_free * y1
    c12 = cjn_equilibrium(1, 2, s_free, y1, y2)
    z2 = 0.5 * J2 * c12
    return (
        ST - (s_free + y1 + 2.0 * y2 + 2.0 * z2),
        lt_total - (l_free + y1 + y2 + z2),
    )


def solve_equilibrium(lt_total: int) -> Dict[str, float]:
    x = [0.6 * ST, 0.4 * lt_total]
    for _ in range(40):
        f0, f1 = equilibrium_residual(x[0], x[1], lt_total)
        norm = max(abs(f0), abs(f1))
        if norm < 1.0e-12:
            break

        h_s = 1.0e-6 * max(1.0, abs(x[0]))
        h_l = 1.0e-6 * max(1.0, abs(x[1]))
        fs0, fs1 = equilibrium_residual(x[0] + h_s, x[1], lt_total)
        fl0, fl1 = equilibrium_residual(x[0], x[1] + h_l, lt_total)

        a = (fs0 - f0) / h_s
        b = (fl0 - f0) / h_l
        c = (fs1 - f1) / h_s
        d = (fl1 - f1) / h_l
        det = a * d - b * c
        if abs(det) < 1.0e-20:
            raise RuntimeError("Equilibrium solver encountered a singular Jacobian")

        dx0 = (-d * f0 + b * f1) / det
        dx1 = (c * f0 - a * f1) / det

        step = 1.0
        improved = False
        while step > 1.0e-6:
            trial = [x[0] + step * dx0, x[1] + step * dx1]
            if trial[0] > 0.0 and trial[1] > 0.0:
                t0, t1 = equilibrium_residual(trial[0], trial[1], lt_total)
                if max(abs(t0), abs(t1)) < norm:
                    x = trial
                    improved = True
                    break
            step *= 0.5
        if not improved:
            x = [x[0] + dx0, x[1] + dx1]

    s_free, l_free = x
    k1_eq = KF / KOFF
    k2_eq = KPX / KOFF
    y1 = 2.0 * k1_eq * l_free * s_free
    y2 = 0.5 * k2_eq * s_free * y1
    c01 = cjn_equilibrium(0, 1, s_free, y1, y2)
    c11 = cjn_equilibrium(1, 1, s_free, y1, y2)
    c21 = cjn_equilibrium(2, 1, s_free, y1, y2)
    c02 = cjn_equilibrium(0, 2, s_free, y1, y2)
    c12 = cjn_equilibrium(1, 2, s_free, y1, y2)
    c22 = cjn_equilibrium(2, 2, s_free, y1, y2)
    z2 = 0.5 * J2 * c12

    return {
        "S": s_free,
        "L": l_free,
        "Y1": y1,
        "Y2": y2,
        "Z2": z2,
        "R2": z2 / 2.0,
        "C01": c01,
        "C11": c11,
        "C21": c21,
        "C02": c02,
        "C12": c12,
        "C22": c22,
    }


def posner_rhs(state: Sequence[float], lt_total: int) -> List[float]:
    y1, y2, z2, c01, c11, c21, c02, c12, c22 = state

    s_free = max(ST - y1 - 2.0 * y2 - 2.0 * z2, 0.0)
    l_free = max(lt_total - y1 - y2 - z2, 0.0)

    dy1 = 2.0 * KF * l_free * s_free - KOFF * y1 - KPX * s_free * y1 + 2.0 * KOFF * y2 - JP * c12 + 2.0 * JM * z2
    dy2 = KPX * s_free * y1 - 2.0 * KOFF * y2 - JP * c12 + 2.0 * JM * z2
    dz2 = 2.0 * JP * c12 - 4.0 * JM * z2

    dc01 = -4.0 * KF * l_free * c01 + KOFF * c11 - 2.0 * KPX * c01 * y1 + KOFF * (s_free - 2.0 * c01 - c11)
    dc11 = (
        4.0 * KF * l_free * c01
        - KOFF * c11
        - 2.0 * KF * l_free * c11
        + 2.0 * KOFF * c21
        - KPX * c11 * (y1 + s_free)
        + KOFF * (s_free + y1 - 2.0 * c01 - 2.0 * c11 - 2.0 * c21)
    )
    dc21 = (
        2.0 * KF * l_free * c11
        - 2.0 * KOFF * c21
        - 2.0 * KPX * c21 * s_free
        + KOFF * (y1 - c11 - 2.0 * c21)
    )
    dc02 = (
        -4.0 * KF * l_free * c02
        + KOFF * c12
        - 2.0 * KPX * c02 * y1
        + 2.0 * KPX * c01 * c11
        - 2.0 * KOFF * c02
        + KOFF * (s_free - 2.0 * c01 - c11 - 2.0 * c02 - c12)
    )
    dc12 = (
        4.0 * KF * l_free * c02
        - KOFF * c12
        - 2.0 * KF * l_free * c12
        + 2.0 * KOFF * c22
        - KPX * c12 * (y1 + s_free)
        + 4.0 * KPX * c01 * c21
        + KPX * (c11 ** 2)
        - 2.0 * KOFF * c12
        + KOFF * (s_free + y1 - 2.0 * c01 - 2.0 * c11 - 2.0 * c21 - 2.0 * c02 - 2.0 * c12 - 2.0 * c22)
        - JP * c12
        + 2.0 * JM * z2
    )
    dc22 = (
        2.0 * KF * l_free * c12
        - 2.0 * KOFF * c22
        - 2.0 * KPX * c22 * s_free
        + 2.0 * KPX * c11 * c21
        - 2.0 * KOFF * c22
        + KOFF * (y1 - c11 - 2.0 * c21 - c12 - 2.0 * c22)
    )

    return [dy1, dy2, dz2, dc01, dc11, dc21, dc02, dc12, dc22]


def rk4_step(state: Sequence[float], dt: float, lt_total: int) -> List[float]:
    k1 = posner_rhs(state, lt_total)
    s2 = [a + 0.5 * dt * b for a, b in zip(state, k1)]
    k2 = posner_rhs(s2, lt_total)
    s3 = [a + 0.5 * dt * b for a, b in zip(state, k2)]
    k3 = posner_rhs(s3, lt_total)
    s4 = [a + dt * b for a, b in zip(state, k3)]
    k4 = posner_rhs(s4, lt_total)
    return [a + dt * (b + 2.0 * c + 2.0 * d + e) / 6.0 for a, b, c, d, e in zip(state, k1, k2, k3, k4)]


def integrate_posner(lt_total: int, dt: float = 0.05) -> Dict[str, List[float]]:
    sample_every = int(round(OUTPUT_DT / dt))
    if abs(sample_every * dt - OUTPUT_DT) > 1.0e-12:
        raise ValueError("dt must divide the output interval exactly")

    state = [0.0] * 9
    state[3] = float(RT)

    traj = {"time": [0.0], "Y1": [0.0], "Y2": [0.0], "Z2": [0.0]}
    for sample_idx in range(1, N_STEPS + 1):
        for _ in range(sample_every):
            state = rk4_step(state, dt, lt_total)
        traj["time"].append(sample_idx * OUTPUT_DT)
        traj["Y1"].append(state[0])
        traj["Y2"].append(state[1])
        traj["Z2"].append(state[2])
    return traj


class Canvas:
    def __init__(self, width: int, height: int, background: Tuple[int, int, int] = (255, 255, 255)) -> None:
        self.width = width
        self.height = height
        self.pixels = bytearray(background * (width * height))

    def set_pixel(self, x: int, y: int, color: Tuple[int, int, int]) -> None:
        if 0 <= x < self.width and 0 <= y < self.height:
            idx = 3 * (y * self.width + x)
            self.pixels[idx : idx + 3] = bytes(color)

    def fill_rect(self, x: int, y: int, w: int, h: int, color: Tuple[int, int, int]) -> None:
        for yy in range(y, y + h):
            if yy < 0 or yy >= self.height:
                continue
            start = 3 * (yy * self.width + max(x, 0))
            end = 3 * (yy * self.width + min(x + w, self.width))
            if start < end:
                row = bytearray()
                for _ in range((end - start) // 3):
                    row.extend(color)
                self.pixels[start:end] = row

    def draw_line(self, x0: float, y0: float, x1: float, y1: float, color: Tuple[int, int, int], width: int = 1) -> None:
        dx = x1 - x0
        dy = y1 - y0
        steps = int(max(abs(dx), abs(dy), 1))
        for i in range(steps + 1):
            x = int(round(x0 + dx * i / steps))
            y = int(round(y0 + dy * i / steps))
            half = max(0, width // 2)
            for xx in range(x - half, x + half + 1):
                for yy in range(y - half, y + half + 1):
                    self.set_pixel(xx, yy, color)

    def draw_polyline(self, points: Sequence[Tuple[float, float]], color: Tuple[int, int, int], width: int = 1) -> None:
        for p0, p1 in zip(points, points[1:]):
            self.draw_line(p0[0], p0[1], p1[0], p1[1], color, width=width)

    def draw_dashed_horizontal(self, x0: int, x1: int, y: int, color: Tuple[int, int, int], dash: int = 7, gap: int = 5) -> None:
        x = x0
        while x < x1:
            self.draw_line(x, y, min(x + dash, x1), y, color)
            x += dash + gap

    def text_width(self, text: str, scale: int = 2) -> int:
        return len(text) * 6 * scale

    def draw_text(self, x: int, y: int, text: str, color: Tuple[int, int, int], scale: int = 2) -> None:
        cursor = x
        for char in text.upper():
            glyph = FONT_5X7.get(char, FONT_5X7[" "])
            for row, row_bits in enumerate(glyph):
                for col, bit in enumerate(row_bits):
                    if bit == "1":
                        self.fill_rect(cursor + col * scale, y + row * scale, scale, scale, color)
            cursor += 6 * scale


def png_chunk(kind: bytes, data: bytes) -> bytes:
    return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)


def write_png(path: Path, canvas: Canvas) -> None:
    raw = bytearray()
    for y in range(canvas.height):
        raw.append(0)
        start = 3 * y * canvas.width
        stop = 3 * (y + 1) * canvas.width
        raw.extend(canvas.pixels[start:stop])

    ihdr = struct.pack(">IIBBBBB", canvas.width, canvas.height, 8, 2, 0, 0, 0)
    idat = zlib.compress(bytes(raw), level=9)
    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + png_chunk(b"IHDR", ihdr)
        + png_chunk(b"IDAT", idat)
        + png_chunk(b"IEND", b"")
    )


def map_points(xs: Sequence[float], ys: Sequence[float], left: int, top: int, width: int, height: int, x_max: float, y_max: float) -> List[Tuple[float, float]]:
    mapped = []
    for x_val, y_val in zip(xs, ys):
        x = left + (x_val / x_max) * width
        y = top + height - (y_val / y_max) * height
        mapped.append((x, y))
    return mapped


def nice_top(value: float) -> float:
    if value <= 0.0:
        return 1.0
    exponent = math.floor(math.log10(value))
    base = value / (10 ** exponent)
    if base <= 1.0:
        nice = 1.0
    elif base <= 2.0:
        nice = 2.0
    elif base <= 5.0:
        nice = 5.0
    else:
        nice = 10.0
    return nice * (10 ** exponent)


def plot_validation(results: Dict[str, Dict[str, object]], out_path: Path) -> None:
    width = 1800
    height = 1280
    canvas = Canvas(width, height)

    black = (10, 10, 10)
    gray = (130, 130, 130)
    light_gray = (225, 225, 225)
    blue = (33, 99, 171)
    red = (184, 42, 42)
    pale = (245, 245, 245)

    canvas.draw_text(360, 26, "POSNER 1995 SIZE 2 RING VALIDATION", black, scale=4)
    canvas.draw_text(430, 68, "BNGL NFSIM VS KINETIC THEORY AND EQUILIBRIUM THEORY", black, scale=2)

    legend_y = 108
    canvas.draw_line(260, legend_y, 340, legend_y, blue, width=3)
    canvas.draw_text(352, legend_y - 9, "POSNER ODE", black, scale=2)
    canvas.draw_line(620, legend_y, 700, legend_y, red, width=3)
    canvas.draw_text(712, legend_y - 9, "NFSIM MEAN", black, scale=2)
    canvas.draw_dashed_horizontal(1010, 1090, legend_y, gray, dash=10, gap=6)
    canvas.draw_text(1102, legend_y - 9, "POSNER EQ", black, scale=2)

    left_margin = 88
    right_margin = 44
    top_margin = 160
    bottom_margin = 70
    col_gap = 34
    row_gap = 38
    panel_w = int((width - left_margin - right_margin - 2 * col_gap) / 3)
    panel_h = int((height - top_margin - bottom_margin - 2 * row_gap) / 3)

    obs_order = ("Y1", "Y2", "Z2")
    y_limits = {}
    for obs in obs_order:
        peak = 0.0
        for dose in results.values():
            theory = dose["theory"]
            nf_mean = dose["nf_mean"]
            equilibrium = dose["equilibrium"]
            peak = max(peak, max(theory[obs]), max(nf_mean[obs]), equilibrium[obs])
        y_limits[obs] = nice_top(1.1 * peak + 0.25)

    for col_idx, obs in enumerate(obs_order):
        title_x = left_margin + col_idx * (panel_w + col_gap) + 12
        canvas.draw_text(title_x, 132, obs, black, scale=3)

    for row_idx, dose in enumerate(DOSES):
        block = results[dose["name"]]
        theory = block["theory"]
        nf_mean = block["nf_mean"]
        equilibrium = block["equilibrium"]

        for col_idx, obs in enumerate(obs_order):
            x0 = left_margin + col_idx * (panel_w + col_gap)
            y0 = top_margin + row_idx * (panel_h + row_gap)
            canvas.fill_rect(x0, y0, panel_w, panel_h, pale)
            canvas.draw_line(x0, y0, x0 + panel_w, y0, black)
            canvas.draw_line(x0, y0 + panel_h, x0 + panel_w, y0 + panel_h, black)
            canvas.draw_line(x0, y0, x0, y0 + panel_h, black)
            canvas.draw_line(x0 + panel_w, y0, x0 + panel_w, y0 + panel_h, black)

            y_max = y_limits[obs]
            for frac in (0.25, 0.5, 0.75):
                yy = int(round(y0 + panel_h * (1.0 - frac)))
                canvas.draw_dashed_horizontal(x0 + 1, x0 + panel_w - 1, yy, light_gray, dash=4, gap=5)

            for tick in (0, 150, 300, 450, 600):
                xx = int(round(x0 + panel_w * (tick / T_END)))
                canvas.draw_line(xx, y0 + panel_h, xx, y0 + panel_h + 6, black)
                tick_label = str(tick)
                canvas.draw_text(xx - canvas.text_width(tick_label, scale=1) // 2, y0 + panel_h + 10, tick_label, black, scale=1)

            y_ticks = (0.0, y_max / 2.0, y_max)
            for tick_val in y_ticks:
                yy = int(round(y0 + panel_h - panel_h * (tick_val / y_max)))
                canvas.draw_line(x0 - 6, yy, x0, yy, black)
                tick_label = str(int(round(tick_val)))
                canvas.draw_text(x0 - 8 - canvas.text_width(tick_label, scale=1), yy - 4, tick_label, black, scale=1)

            eq_y = int(round(y0 + panel_h - panel_h * (equilibrium[obs] / y_max)))
            canvas.draw_dashed_horizontal(x0 + 1, x0 + panel_w - 1, eq_y, gray, dash=8, gap=5)

            theory_pts = map_points(theory["time"], theory[obs], x0, y0, panel_w, panel_h, T_END, y_max)
            nf_pts = map_points(nf_mean["time"], nf_mean[obs], x0, y0, panel_w, panel_h, T_END, y_max)
            canvas.draw_polyline(theory_pts, blue, width=2)
            canvas.draw_polyline(nf_pts, red, width=2)

            if col_idx == 0:
                canvas.draw_text(x0 + 12, y0 + 12, dose["title"], black, scale=2)
            if row_idx == 2:
                canvas.draw_text(x0 + panel_w - 70, y0 + panel_h + 34, "TIME", black, scale=2)

    write_png(out_path, canvas)


def bngl_output_prefix(suffix: str) -> str:
    return f"{BNGL_PATH.stem}_{suffix}"


def run_bngl_actions() -> None:
    bng2 = locate_executable(BNG2_CANDIDATES, "BNG2.pl")
    run_checked([str(bng2), str(BNGL_PATH.name)], cwd=ROOT)


def run_nfsim_ensemble(xml_path: Path, n_runs: int, seed0: int) -> Dict[str, List[float]]:
    nfsim = locate_executable(NFSIM_CANDIDATES, "NFsim")
    runs = []
    with tempfile.TemporaryDirectory(prefix="size2_nfsim_", dir=str(ROOT)) as tempdir:
        tmp = Path(tempdir)
        for offset in range(n_runs):
            seed = seed0 + offset
            gdat_path = tmp / f"{xml_path.stem}_seed{seed}.gdat"
            cmd = [
                str(nfsim),
                "-xml",
                str(xml_path),
                "-sim",
                str(T_END),
                "-oSteps",
                str(N_STEPS),
                "-seed",
                str(seed),
                "-ogf",
                "-bscb",
                "-utl",
                "5",
                "-o",
                str(gdat_path),
            ]
            run_checked(cmd, cwd=ROOT)
            runs.append(normalize_series(parse_gdat(gdat_path)))
    return mean_trajectories(runs)


def summarize_block(theory: Dict[str, List[float]], nf_mean: Dict[str, List[float]], equilibrium: Dict[str, float], lt_total: int, n_runs: int) -> Dict[str, object]:
    summary = {
        "LT": lt_total,
        "n_runs": n_runs,
        "equilibrium": {obs: equilibrium[obs] for obs in ("Y1", "Y2", "Z2")},
        "final_ode": {obs: theory[obs][-1] for obs in ("Y1", "Y2", "Z2")},
        "final_nfsim_mean": {obs: nf_mean[obs][-1] for obs in ("Y1", "Y2", "Z2")},
        "rms_nfsim_vs_ode": {},
    }
    for obs in ("Y1", "Y2", "Z2"):
        diffs = [a - b for a, b in zip(nf_mean[obs], theory[obs])]
        summary["rms_nfsim_vs_ode"][obs] = math.sqrt(sum(d * d for d in diffs) / len(diffs))
    return summary


def run_validation_suite(n_runs: int = 96) -> Dict[str, object]:
    run_bngl_actions()

    results = {}
    summary = {}
    for dose_idx, dose in enumerate(DOSES):
        prefix = bngl_output_prefix(dose["name"])
        xml_path = ROOT / f"{prefix}.xml"
        gdat_path = ROOT / f"{prefix}.gdat"
        if not xml_path.exists():
            raise FileNotFoundError(f"Expected BNGL-generated XML not found: {xml_path}")
        if not gdat_path.exists():
            raise FileNotFoundError(f"Expected BNGL-generated gdat not found: {gdat_path}")

        equilibrium = solve_equilibrium(dose["lt"])
        theory = integrate_posner(dose["lt"])
        nf_mean = run_nfsim_ensemble(xml_path, n_runs=n_runs, seed0=5000 + 1000 * dose_idx)
        single_run = normalize_series(parse_gdat(gdat_path))

        results[dose["name"]] = {
            "lt": dose["lt"],
            "title": dose["title"],
            "theory": theory,
            "equilibrium": equilibrium,
            "nf_mean": nf_mean,
            "single_run": single_run,
            "xml_path": str(xml_path),
            "gdat_path": str(gdat_path),
        }
        summary[dose["name"]] = summarize_block(theory, nf_mean, equilibrium, dose["lt"], n_runs)

    plot_validation(results, PNG_PATH)
    payload = {
        "bundle_dir": str(ROOT),
        "bngl_path": str(BNGL_PATH),
        "png_path": str(PNG_PATH),
        "summary": summary,
    }
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    payload = run_validation_suite()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
