"""The RASHER logo (a circular genome with paired reads) and both pipeline diagrams (run from docs/images/)."""

LIGHT = {"text": "#1d2433", "muted": "#5d6678", "box": "#f6f7f9", "edge": "#d7dbe2", "arrow": "#8a93a3",
         "reads": "#8a93a3", "ref": "#5d6678"}
DARK = {"text": "#e6e8ec", "muted": "#9aa3b2", "box": "#1d2129", "edge": "#363c47", "arrow": "#737d8e",
        "reads": "#737d8e", "ref": "#9aa3b2"}
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


def arc(cx, cy, r, a0, a1):
    """An SVG arc from angle a0 to a1 (degrees)."""
    import math
    x0, y0 = cx + r * math.cos(math.radians(a0)), cy + r * math.sin(math.radians(a0))
    x1, y1 = cx + r * math.cos(math.radians(a1)), cy + r * math.sin(math.radians(a1))
    return f"M{x0:.1f},{y0:.1f} A{r},{r} 0 0 1 {x1:.1f},{y1:.1f}"


# The icon: a circular reference genome with paired short-read arcs aligned on two tracks; SNPs as dots
READ_PAIRS = [(52, (-80, -58), (-52, -30)), (52, (-5, 17), (23, 45)), (52, (95, 117), (123, 145)),
              (64, (-130, -108), (-102, -80)), (64, (-40, -18), (-12, 10)), (64, (60, 82), (88, 110)),
              (64, (150, 172), (178, 200))]
SNPS = [(52, 34), (64, -91), (52, 106)]


def logo(theme):
    import math
    t = THEMES[theme]
    cx, cy = 85, 82
    parts = [f'<circle cx="{cx}" cy="{cy}" r="40" fill="none" stroke="{t["ref"]}" stroke-width="5"/>']
    for r, mate1, mate2 in READ_PAIRS:
        for a0, a1 in (mate1, mate2):
            parts.append(f'<path d="{arc(cx, cy, r, a0, a1)}" fill="none" stroke="{t["reads"]}" '
                         f'stroke-width="8" stroke-linecap="round"/>')
    for r, a in SNPS:
        x, y = cx + r * math.cos(math.radians(a)), cy + r * math.sin(math.radians(a))
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5.5" fill="{MEAT}"/>')
    icon = "\n  ".join(parts)
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 620 168" width="620" height="168"
  font-family="{FONT}">
  <title>RASHER: Bait, Assemble and Compare Illumina short reads</title>{style(theme)}
  {icon}
  <text x="185" y="92" font-size="56" font-weight="800" letter-spacing="7" fill="{t['text']}">RASHER</text>
  <text x="185" y="122" font-size="14.5" letter-spacing="0.8" fill="{t['muted']}">Bait, Assemble and Compare \
Illumina short reads</text>
</svg>
"""


def box(x, w, n, title, lines, t, h=112, y=40):
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="{t["box"]}" stroke="{t["edge"]}"/>']
    out.append(f'<circle cx="{x + 25}" cy="{y + 27}" r="13" fill="{MEAT}"/>')
    out.append(f'<text x="{x + 25}" y="{y + 32}" text-anchor="middle" font-size="15" font-weight="700" '
               f'fill="#ffffff">{n}</text>')
    out.append(f'<text x="{x + 46}" y="{y + 33}" font-size="19" font-weight="700" '
               f'fill="{t["text"]}">{title}</text>')
    for i, line in enumerate(lines):
        out.append(f'<text x="{x + 16}" y="{y + 60 + i * 21}" font-size="14" fill="{t["muted"]}">{line}</text>')
    return "\n  ".join(out)


def chip(x, y, w, label, t):
    return (f'<rect x="{x}" y="{y}" width="{w}" height="38" rx="19" fill="none" stroke="{MEAT}" '
            f'stroke-width="1.6"/>\n  <text x="{x + w / 2}" y="{y + 25}" text-anchor="middle" font-size="15.5" '
            f'font-weight="600" fill="{t["text"]}">{label}</text>')


def arrow(x0, x1, t, y=95):
    return (f'<line x1="{x0}" y1="{y}" x2="{x1 - 9}" y2="{y}" stroke="{t["arrow"]}" stroke-width="2"/>'
            f'<path d="M{x1},{y} l-10,-5 v10 z" fill="{t["arrow"]}"/>')


def pipeline(theme, assemble_lines, filter_lines, tool):
    t = THEMES[theme]
    parts = [chip(12, 76, 80, "reads", t), arrow(92, 124, t)]
    parts.append(box(124, 182, 1, "Bait", ["minimap2 (-x sr)" if tool == "rasher" else "minimap2 (map-ont)",
                                           "or BBDuk"], t))
    parts.append(arrow(306, 338, t))
    parts.append(box(338, 206, 2, "Filter", filter_lines, t))
    parts.append(arrow(544, 576, t))
    parts.append(box(576, 252, 3, "Assemble", assemble_lines, t))
    parts.append(arrow(828, 860, t))
    parts.append(box(860, 252, 4, "Compare", ["SNPs: minimap2, SKA2 or Parsnp", "tree: FastTree or IQ-TREE"], t))
    parts.append(arrow(1112, 1144, t))
    parts.append(chip(1144, 76, 90, "report", t))
    body = "\n  ".join(parts)
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1246 192" width="1246" height="192"
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
