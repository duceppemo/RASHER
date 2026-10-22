"""RASHER logo and the RASHER/BACoN pipeline diagrams, light and dark variants (run from docs/images/)."""

THEMES = {
    "light": {"text": "#1d2433", "muted": "#5d6678", "box": "#f6f7f9", "edge": "#d7dbe2", "arrow": "#8a93a3"},
    "dark": {"text": "#e6e8ec", "muted": "#9aa3b2", "box": "#1d2129", "edge": "#363c47", "arrow": "#737d8e"},
}
MEAT, MEAT_EDGE, FAT = "#c2402f", "#8c2b1f", "#f4dfc8"
FONT = "system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"


def wave(y, x0=90, x1=510, amp=13):
    """A smooth 2.5-period wave from x0 to x1 at baseline y."""
    step = (x1 - x0) / 5
    d = f"M{x0},{y} "
    up = True
    for i in range(5):
        a, b = x0 + i * step, x0 + (i + 1) * step
        cy = y - amp if up else y + amp
        d += f"C{a + step * 0.35:.0f},{cy} {b - step * 0.35:.0f},{cy} {b:.0f},{y} "
        up = not up
    return d


def logo(theme):
    t = THEMES[theme]
    strip = wave(62)
    fat1 = wave(50, 110, 490, 10)
    fat2 = wave(74, 110, 490, 10)
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 600 210" width="600" height="210"
  font-family="{FONT}">
  <title>RASHER: Bait, Assemble and Compare Illumina short reads</title>
  <path d="{strip}" fill="none" stroke="{MEAT_EDGE}" stroke-width="54" stroke-linecap="round"/>
  <path d="{strip}" fill="none" stroke="{MEAT}" stroke-width="46" stroke-linecap="round"/>
  <path d="{fat1}" fill="none" stroke="{FAT}" stroke-width="7" stroke-linecap="round" opacity="0.95"/>
  <path d="{fat2}" fill="none" stroke="{FAT}" stroke-width="7" stroke-linecap="round" opacity="0.95"/>
  <g fill="{FAT}">
    <circle cx="216" cy="62" r="5"/>
    <circle cx="300" cy="62" r="5"/>
    <circle cx="384" cy="62" r="5"/>
  </g>
  <text x="300" y="155" text-anchor="middle" font-size="54" font-weight="800" letter-spacing="10"
    fill="{t['text']}">RASHER</text>
  <text x="300" y="186" text-anchor="middle" font-size="15" letter-spacing="1"
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
  <title>{'RASHER' if tool == 'rasher' else 'BACoN'} pipeline: bait, filter, assemble, compare</title>
  {body}
</svg>
"""


def write(name, text):
    with open(name, "w") as fh:
        fh.write(text)


for theme in ("light", "dark"):
    write(f"rasher_logo_{theme}.svg", logo(theme))
    write(f"rasher_pipeline_{theme}.svg", pipeline(
        theme, ["templated: samtools consensus", "de novo: SPAdes or SKESA"],
        ["fastp (quality, adapters),", "then even subsampling"], "rasher"))
    write(f"bacon_pipeline_{theme}.svg", pipeline(
        theme, ["templated: samtools consensus", "de novo: Flye or myloasm"],
        ["Filtlong: the best reads,", "up to the target depth"], "bacon"))
print("written")
