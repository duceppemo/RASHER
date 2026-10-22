"""RASHER logo and the RASHER/BACoN pipeline diagrams, light and dark variants (run from docs/images/)."""

LIGHT = {"text": "#1d2433", "muted": "#5d6678", "box": "#f6f7f9", "edge": "#d7dbe2", "arrow": "#8a93a3"}
DARK = {"text": "#e6e8ec", "muted": "#9aa3b2", "box": "#1d2129", "edge": "#363c47", "arrow": "#737d8e"}
THEMES = {
    "light": LIGHT,
    "dark": DARK,
    # One adaptive file: colors are CSS variables, flipped by a media query inside the SVG. Works wherever
    # the SVG is shown in an <img>, GitHub wikis included (the #gh-*-mode-only trick does not work there).
    "auto": {k: f"var(--{k})" for k in LIGHT},
}


def style(theme):
    if theme != "auto":
        return ""
    inner = "; ".join(f"--{k}: {v}" for k, v in LIGHT.items())
    dark = "; ".join(f"--{k}: {v}" for k, v in DARK.items())
    return (f"\n  <style>:root {{ {inner} }} "
            f"@media (prefers-color-scheme: dark) {{ :root {{ {dark} }} }}</style>")
MEAT, MEAT_EDGE, FAT = "#d94f3d", "#9c2f21", "#faeedd"
FONT = "system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"


def wave_points(y, x0, x1, amp, periods=2.5, step_px=6):
    """Points along a sine from x0 to x1 at baseline y."""
    import math
    n = max(2, int((x1 - x0) / step_px))
    return [(x0 + (x1 - x0) * i / n, y + amp * math.sin(2 * math.pi * periods * i / n)) for i in range(n + 1)]


def polyline(points):
    return " ".join(f"L{x:.1f},{y:.1f}" for x, y in points)


def ribbon(x0=110, x1=490, top=40, thickness=52, amp=15):
    """A closed wavy strip of constant vertical thickness: a rasher seen flat, with cut ends."""
    upper = wave_points(top, x0, x1, amp)
    lower = wave_points(top + thickness, x1, x0, amp)[::1]
    lower = [(x, y) for x, y in reversed(wave_points(top + thickness, x0, x1, amp))]
    d = f"M{upper[0][0]:.1f},{upper[0][1]:.1f} " + polyline(upper[1:])
    d += f" L{lower[0][0]:.1f},{lower[0][1]:.1f} " + polyline(lower[1:]) + " Z"
    return d


def streak(y, x0=124, x1=476, amp=15):
    pts = wave_points(y, x0, x1, amp)
    return f"M{pts[0][0]:.1f},{pts[0][1]:.1f} " + polyline(pts[1:])


def logo(theme):
    t = THEMES[theme]
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 600 220" width="600" height="220"
  font-family="{FONT}">
  <title>RASHER: Bait, Assemble and Compare Illumina short reads</title>{style(theme)}
  <path d="{ribbon()}" fill="{MEAT}" stroke="{MEAT_EDGE}" stroke-width="3" stroke-linejoin="round"/>
  <path d="{streak(56)}" fill="none" stroke="{FAT}" stroke-width="10" stroke-linecap="butt"/>
  <path d="{streak(76)}" fill="none" stroke="{FAT}" stroke-width="10" stroke-linecap="butt"/>
  <text x="300" y="165" text-anchor="middle" font-size="54" font-weight="800" letter-spacing="10"
    fill="{t['text']}">RASHER</text>
  <text x="300" y="196" text-anchor="middle" font-size="15" letter-spacing="1"
    fill="{t['muted']}">Bait, Assemble and Compare Illumina short reads</text>
</svg>
"""


def box(x, w, n, title, lines, t, h=110, y=40):
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="{t["box"]}" stroke="{t["edge"]}"/>']
    out.append(f'<circle cx="{x + 24}" cy="{y + 26}" r="12" fill="{MEAT}"/>')
    out.append(f'<text x="{x + 24}" y="{y + 31}" text-anchor="middle" font-size="14" font-weight="700" '
               f'fill="#ffffff">{n}</text>')
    out.append(f'<text x="{x + 44}" y="{y + 31}" font-size="17" font-weight="700" '
               f'fill="{t["text"]}">{title}</text>')
    for i, line in enumerate(lines):
        out.append(f'<text x="{x + 16}" y="{y + 56 + i * 19}" font-size="12.5" fill="{t["muted"]}">{line}</text>')
    return "\n  ".join(out)


def chip(x, y, w, label, t):
    return (f'<rect x="{x}" y="{y}" width="{w}" height="36" rx="18" fill="none" stroke="{MEAT}" '
            f'stroke-width="1.6"/>\n  <text x="{x + w / 2}" y="{y + 23}" text-anchor="middle" font-size="14" '
            f'font-weight="600" fill="{t["text"]}">{label}</text>')


def arrow(x0, x1, t, y=95):
    return (f'<line x1="{x0}" y1="{y}" x2="{x1 - 9}" y2="{y}" stroke="{t["arrow"]}" stroke-width="2"/>'
            f'<path d="M{x1},{y} l-10,-5 v10 z" fill="{t["arrow"]}"/>')


def pipeline(theme, assemble_lines, filter_lines, tool):
    t = THEMES[theme]
    parts = [chip(12, 77, 74, "reads", t), arrow(86, 118, t)]
    parts.append(box(118, 168, 1, "Bait", ["minimap2 (-x sr)" if tool == "rasher" else "minimap2 (map-ont)",
                                           "or BBDuk"], t))
    parts.append(arrow(286, 318, t))
    parts.append(box(318, 172, 2, "Filter", filter_lines, t))
    parts.append(arrow(490, 522, t))
    parts.append(box(522, 230, 3, "Assemble", assemble_lines, t))
    parts.append(arrow(752, 784, t))
    parts.append(box(784, 230, 4, "Compare", ["SNPs: minimap2, SKA2 or Parsnp", "tree: FastTree or IQ-TREE"], t))
    parts.append(arrow(1014, 1046, t))
    parts.append(chip(1046, 77, 82, "report", t))
    body = "\n  ".join(parts)
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1142 190" width="1142" height="190"
  font-family="{FONT}">
  <title>{'RASHER' if tool == 'rasher' else 'BACoN'} pipeline: bait, filter, assemble, \
compare</title>{style(theme)}
  {body}
</svg>
"""


def write(name, text):
    with open(name, "w") as fh:
        fh.write(text)


for theme in ("light", "dark", "auto"):
    write("rasher_logo.svg" if theme == "auto" else f"rasher_logo_{theme}.svg", logo(theme))
    write("rasher_pipeline.svg" if theme == "auto" else f"rasher_pipeline_{theme}.svg", pipeline(
        theme, ["templated: samtools consensus", "de novo: SPAdes or SKESA"],
        ["fastp (quality, adapters),", "then even subsampling"], "rasher"))
    write("bacon_pipeline.svg" if theme == "auto" else f"bacon_pipeline_{theme}.svg", pipeline(
        theme, ["templated: samtools consensus", "de novo: Flye or myloasm"],
        ["Filtlong: the best reads,", "up to the target depth"], "bacon"))
print("written")
