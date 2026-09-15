#!/usr/bin/env python3
"""
Methods schematic for the multi-scale habitat association analysis.

Draws the data-to-inference workflow for the multi-scale habitat association
analysis: two preparation lanes (response variable / temporal aggregation, and
landscape predictors / spatial extraction / predictor screening) converging on
model fitting, then the estimated quantity.

The result panel shows one real scale-of-effect profile twice: as a coefficient
curve across the 16 spatial extents, and as the corresponding row of the domain
bar chart (Figure 3), so a reader can map one representation onto the other.

REPRODUCIBILITY
    The example profile is read from the analysis output when that file is
    available, so the figure tracks the fitted models rather than a transcription:

        output/scale_analysis/domain_classification_full.csv

    If the file is absent the script falls back to values embedded below, which
    were taken from the 2026-08-31 knit. When the file IS present the script
    compares the two and stops if they disagree, so a silent drift between the
    figure and the analysis is not possible.

    Radii in the nested-buffer glyph are mapped to log(extent), matching the log
    x-axis used elsewhere in the figure. The 16 extents are geometric (each
    ~1.414x the previous), so on a log scale they are evenly spaced. This is
    deliberate: drawn true to scale, the eight finest radii would all fall
    within about 4 px of the centre.

USAGE
    From repo/Python/:
        python methods_schematic.py ../output/methods_schematic/methods_schematic

    Writes OUTPUT_STEM.svg (default: methods_schematic.svg). If the optional
    dependency cairosvg is installed, also writes OUTPUT_STEM.png at 1700 px wide.
    SVG is the deliverable; the PNG is for on-screen checking only.

DEPENDENCIES
    Python 3.8+ standard library. cairosvg is optional (PNG export only).
"""

import base64
import csv
import math
import os
import struct
import sys

# ---------------------------------------------------------------------------
# Data: the example profile shown in the result panel
# ---------------------------------------------------------------------------

# Which profile to illustrate. Horse / October / human use intensity is used
# because it is the paper's flagship reversal: positive at 160 m, negative at
# 7.2 km, so the extent chosen decides the sign of the inference.
SPECIES, PERIOD, COVARIATE = "Horse", "Oct", "hu"

# Path to the analysis output, relative to this script (repo/Python/).
CSV_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "output", "scale_analysis", "domain_classification_full.csv",
)

# Fallback values, transcribed from the 2026-08-31 knit. Each tuple is
# (buffer radius in m, standardised coefficient, domain sign) and `CI_FALLBACK`
# holds the 95% confidence limits keyed by radius. Domain sign is '+' or '-'
# where the interval excludes zero in a consistent direction, '0' otherwise.
PROF_FALLBACK = [
    (40, 0.311, "0"), (56, 0.311, "0"), (80, 0.313, "0"), (113, 0.328, "0"),
    (160, 0.339, "+"), (226, 0.320, "0"), (319, 0.286, "0"), (451, 0.284, "0"),
    (638, 0.266, "0"), (903, 0.240, "0"), (1277, 0.129, "0"), (1806, -0.005, "0"),
    (2554, -0.167, "0"), (3612, -0.350, "-"), (5108, -0.586, "-"), (7223, -0.631, "-"),
]
CI_FALLBACK = {
    40: (-0.017, 0.638), 56: (-0.017, 0.639), 80: (-0.015, 0.642),
    113: (-0.003, 0.659), 160: (0.001, 0.676), 226: (-0.019, 0.660),
    319: (-0.054, 0.626), 451: (-0.056, 0.625), 638: (-0.071, 0.603),
    903: (-0.095, 0.574), 1277: (-0.153, 0.412), 1806: (-0.297, 0.287),
    2554: (-0.460, 0.127), 3612: (-0.672, -0.028), 5108: (-0.968, -0.204),
    7223: (-1.074, -0.189),
}


def load_profile(path=CSV_PATH):
    """Return (prof, ci) for the illustrated profile.

    Reads `domain_classification_full.csv` if present and verifies it against
    the embedded fallback; raises SystemExit on disagreement so the figure can
    never silently drift from the analysis. Falls back to the embedded values
    (with a warning) when the file is not available, e.g. when the script is
    run outside the project tree.
    """
    if not os.path.exists(path):
        sys.stderr.write(
            "NOTE: %s not found; using the embedded 2026-08-31 values.\n" % path
        )
        return PROF_FALLBACK, CI_FALLBACK

    rows = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if (r["species"] == SPECIES and r["period"] == PERIOD
                    and r["covariate"] == COVARIATE):
                rows.append(r)
    if len(rows) != 16:
        raise SystemExit(
            "Expected 16 extents for %s / %s / %s, found %d in %s"
            % (SPECIES, PERIOD, COVARIATE, len(rows), path)
        )

    rows.sort(key=lambda r: float(r["extent_radius_m"]))
    prof, ci = [], {}
    for r in rows:
        m = int(round(float(r["extent_radius_m"])))
        prof.append((m, round(float(r["estimate"]), 3), r["domain_sign"]))
        ci[m] = (round(float(r["conf.low"]), 3), round(float(r["conf.high"]), 3))

    # Guard against silent drift between the figure and the fitted models.
    if [p[0] for p in prof] != [p[0] for p in PROF_FALLBACK]:
        raise SystemExit("Buffer radii in the CSV do not match the embedded set.")
    for (m, est, sign), (fm, fest, fsign) in zip(prof, PROF_FALLBACK):
        if sign != fsign or abs(est - fest) > 0.0015:
            raise SystemExit(
                "Profile has changed since 2026-08-31 (extent %d m: %s %.3f vs "
                "embedded %s %.3f). Update PROF_FALLBACK/CI_FALLBACK after "
                "confirming the new knit is the one you want to publish." % (m, sign, est, fsign, fest)
            )
    return prof, ci


prof, ci = load_profile()

# ---------------------------------------------------------------------------
# Palette
#   Box fills form a luminance sequence (greyscale 249 / 235 / 221 / 191) so the
#   four categories separate under any colour vision deficiency and in greyscale.
#   Domain colours match the manuscript's Figure 3 exactly: they were sampled
#   from the rendered domain_bar_chart.png.
# ---------------------------------------------------------------------------
BLUE,RED,GREY='#0072B2','#B22222','#b0b0b0'
A_DATA,A_SPEC,A_RES='#333333','#333333','#333333'
S_SPEC='#7396bd'
F_MF,S_MF='#9ec5e8','#1f4e79'
F_DATA,F_SPEC,F_RES='#fafafa','#cfe0f3','#e2f0e7'
INK,MUTE,EDGE='#1a1a1a','#5a5a5a','#a2aab2'
W,H=1123,794

# --- Canvas crop --------------------------------------------------------------
# The drawing is laid out on the W x H grid above, then the visible canvas is
# cropped to the drawn content via the viewBox, with an equal MARGIN on all four
# sides. Nothing inside the figure moves; only the visible canvas changes.
#
# CONTENT_* is the ink bounding box in user units, measured from a 1:1 render
# (cairosvg at output_width=1123) on 2026-09-03: x 13.35 to 1102.65 (Data
# sources box edge to the right edge of the result panel), y 75.35 to 690.5
# (top of the upper box row to the last footnote line). Stroke widths included.
# Y1 tracks gy below: it is the second footnote baseline, gy + 46, plus descender.
# RE-MEASURE THESE IF ANY BOX MOVES OR THE FOOTNOTES CHANGE, or the crop will
# clip content or leave uneven margins.
CONTENT_X0,CONTENT_X1=13.35,1102.65
CONTENT_Y0,CONTENT_Y1=75.35,690.5
MARGIN=7.0                      # equal on all four sides
CROP_X=math.floor(CONTENT_X0-MARGIN)
CROP_Y=math.floor(CONTENT_Y0-MARGIN)
CROP_W=math.ceil(CONTENT_X1+MARGIN-CROP_X)
CROP_H=math.ceil(CONTENT_Y1+MARGIN-CROP_Y)

o=[]
# --- Species icons -----------------------------------------------------------
# The same three silhouettes used in the domain bar chart (Figure 3), tinted to
# the same colours, so a reader carries one visual key between the two figures.
# Embedded as base64 data URIs so the SVG stays self-contained. The source PNGs
# are greyscale-plus-alpha silhouettes and are never modified on disk; the tint
# is applied to a copy in memory. An SVG flood filter would be tidier but
# cairosvg ignores filters, so the PNG export would lose the colour.
ICON_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "input")
SPECIES_ICONS = [("horse.png", "#5C2E0A"),   # dark chestnut brown
                 ("cow.png",   "#000000"),   # black
                 ("elk.png",   "#708090")]   # slate blue-grey
# Alpha caps from S7_Main_analysis.Rmd (ensure_transparent alpha_value), so the
# icons composite to the same shade as in Figure 3.
ICON_ALPHA = {"horse.png": 0.60, "cow.png": 0.95, "elk.png": 0.60}

def _png_size(path):
    """Width and height from the PNG IHDR chunk (stdlib only)."""
    with open(path, "rb") as fh:
        return struct.unpack(">II", fh.read(24)[16:24])

def _icon_uri(path, colour):
    """Base64 data URI for one icon, prepared exactly as the domain bar chart
    (Figure 3) prepares it in S7_Main_analysis.Rmd (ensure_transparent, then
    colorize_icon):

      1. near-white pixels (grey > 0.9) become fully transparent;
      2. every other pixel's alpha is capped at ICON_ALPHA (0.60 horse and
         elk, 0.95 cattle);
      3. the silhouette is filled with the species colour.

    Step 1 matters for the elk: its PNG has opaque near-white pixels between
    the antler tines. Tinting by alpha alone (the approach used until
    2026-09-15) filled that gap in as a blotch.

    Requires Pillow. Without it the icon is embedded unchanged and renders as a
    black silhouette (with the same antler blotch): the figure is still
    correct, only not matched to Figure 3.
    """
    try:
        import io
        from PIL import Image
        la = Image.open(path).convert("LA")
        grey, alpha = la.split()
        cap = int(round(255 * ICON_ALPHA.get(os.path.basename(path), 0.60)))
        white = grey.point(lambda v: 255 if v > 0.9 * 255 else 0)
        alpha = Image.composite(Image.new("L", la.size, 0),
                                alpha.point(lambda a: min(a, cap)), white)
        rgb = tuple(int(colour[i:i + 2], 16) for i in (1, 3, 5))
        tinted = Image.new("RGBA", la.size, rgb + (255,))
        tinted.putalpha(alpha)
        buf = io.BytesIO()
        tinted.save(buf, format="PNG")
        data = buf.getvalue()
    except ImportError:
        with open(path, "rb") as fh:
            data = fh.read()
    return "data:image/png;base64," + base64.b64encode(data).decode("ascii")

def esc(s): return s.replace('&','&amp;').replace('<','&lt;')
def t(x,y,s,size=10.5,fill=INK,anchor='start',weight='normal'):
    o.append(f'<text x="{x}" y="{y}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" '
             f'font-weight="{weight}" font-family="Helvetica,Arial,sans-serif">{esc(s)}</text>')
LH=15.0
def bh(lines,gh=0): return 30+LH*len(lines)+12+gh
def box(x,y,w,title,lines,fill=F_DATA,stroke=EDGE,glyph=None,gh=0):
    h=bh(lines,gh)
    o.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" fill="{fill}" stroke="{stroke}" stroke-width="1.3"/>')
    t(x+13,y+23,title,12.5,INK,'start','bold')
    yy=y+42
    for ln in lines:
        if ln.startswith('~'):
            t(x+13,yy,ln[1:],10.3,MUTE)
        elif ln.startswith('  '):
            t(x+25,yy,ln.strip(),10.3,MUTE)
        else:
            o.append(f'<circle cx="{x+17}" cy="{yy-3.5}" r="1.9" fill="{MUTE}"/>')
            t(x+25,yy,ln,10.3,MUTE)
        yy+=LH
    if glyph is not None: glyph(x,w,yy-6)
    return h
def arrow(pts,col,wdt=2.3):
    d='M'+' L'.join(f'{a},{b}' for a,b in pts)
    o.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{wdt}" marker-end="url(#a_{col[1:]})"/>')
o.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{CROP_W}" height="{CROP_H}" '
         f'viewBox="{CROP_X} {CROP_Y} {CROP_W} {CROP_H}">')
o.append('<defs>')
for c in dict.fromkeys((A_DATA,A_SPEC,A_RES)):
    o.append(f'<marker id="a_{c[1:]}" viewBox="0 0 10 10" refX="9.2" refY="5" markerWidth="3.9" '
             f'markerHeight="3.9" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="{c}"/></marker>')
o.append('</defs>')
# White ground plus the outer frame. Inset by half the stroke so the whole
# stroke falls inside the viewBox instead of being clipped at the canvas edge.
# Stroke width and corner radius match the boxes (1.3 and rx=6).
BORDER_W=1.3
_bi=BORDER_W/2.0
o.append(f'<rect x="{CROP_X+_bi}" y="{CROP_Y+_bi}" width="{CROP_W-BORDER_W}" '
         f'height="{CROP_H-BORDER_W}" rx="6" fill="#ffffff" stroke="{INK}" '
         f'stroke-width="{BORDER_W}"/>')
HORSE,CATTLE,ELK='#1b9e77','#d95f02','#7570b3'
def glyph_species(bx,bw,gy):
    """Row of the three species silhouettes, matched to Figure 3."""
    ih = 26.0
    specs = []
    for fn, col in SPECIES_ICONS:
        path = os.path.join(ICON_DIR, fn)
        w, h = _png_size(path)
        specs.append((path, col, ih * w / h))
    gap = 14.0
    total = sum(sp[2] for sp in specs) + gap * (len(specs) - 1)
    x = bx + (bw - total) / 2.0
    for path, col, iw in specs:
        o.append(f'<image x="{x:.1f}" y="{gy:.1f}" width="{iw:.1f}" height="{ih:.1f}" '
                 f'href="{_icon_uri(path, col)}" preserveAspectRatio="xMidYMid meet"/>')
        x += iw + gap

def glyph_ladder(bx,bw,gy):
    months=['May','Jun','Jul','Aug','Sep','Oct']
    gx0=bx+40; cw=(bw-52)/6.0
    for i,m in enumerate(months):
        t(gx0+i*cw+cw/2,gy+8,m,6.4,MUTE,'middle')
    rows=[('Horse',HORSE,[(0,1),(1,1),(2,1),(3,1),(4,1),(5,1)]),
          ('Cattle',CATTLE,[(2,2),(4,2)]),
          ('Elk',ELK,[(0,3),(3,3)])]
    for j,(nm,col,bins) in enumerate(rows):
        yy=gy+13+j*15
        t(bx+37,yy+8.5,nm,7.4,MUTE,'end')
        if nm=='Cattle':
            o.append(f'<rect x="{gx0:.1f}" y="{yy}" width="{2*cw-2:.1f}" height="11" rx="2" fill="#efefef" stroke="#e2e2e2" stroke-width="0.7"/>')
            t(gx0+cw-1,yy+7.8,'no data',6.2,MUTE,'middle')
        for st,ln in bins:
            o.append(f'<rect x="{gx0+st*cw:.1f}" y="{yy}" width="{ln*cw-2:.1f}" height="11" rx="2" '
                     f'fill="{col}" fill-opacity="0.30" stroke="{col}" stroke-width="0.9"/>')
def glyph_rings(bx,bw,gy):
    cx=bx+bw/2; cy=gy+74
    ext=[40,56,80,113,160,226,319,451,638,903,1277,1806,2554,3612,5108,7223]
    for i,m in enumerate(ext):
        r=9+(64-9)*(math.log10(m/40)/math.log10(7223/40))
        o.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="none" stroke="#4a6fa5" '
                 f'stroke-opacity="{min(0.25+0.042*i,0.85):.2f}" stroke-width="0.9"/>')
    o.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="3.0" fill="{INK}"/>')
    o.append(f'<line x1="{cx:.1f}" y1="{cy:.1f}" x2="{cx+64:.1f}" y2="{cy:.1f}" stroke="{MUTE}" stroke-width="0.9" stroke-dasharray="2 2"/>')
    t(cx+7,cy-5,'40 m',7.4,MUTE)
    t(cx+62,cy-5,'7.2 km',7.4,MUTE,'end')
    t(cx,cy+78,'16 log-spaced radii',7.6,MUTE,'middle')

CW=178
X1,X2,X3,X4=196,400,604,None
TOPY,BOTY=76,414

box(14,236,150,'Data sources',
    ['Camera traps, 276','  locations, 2022 to 2025','PlanetScope NDVI, 3 m','AVI vegetation strata','Harvest layer, DEM'],F_DATA)
box(X1,TOPY,CW,'Response variable',
    ['Independent detections of','  horses, cattle, and elk,','  30 minute threshold',
     'The camera-month is the','  unit of analysis','Effort offset: log camera-days'],
    F_DATA,glyph=glyph_species,gh=34)
box(X2,TOPY,CW,'Temporal aggregation',
    ['~Set by detection rate','Monthly for horses','Bimonthly for cattle','Quarterly for elk'],F_SPEC,S_SPEC,glyph_ladder,62)
h_mf=box(X3,TOPY,CW,'Model fitting',
    ['Negative binomial GLM','Spatial autocovariate','Stratum fixed effect','Effort offset','Every species and period','  refitted at all 16 extents','160 base model fits'],F_MF,S_MF)
box(X1,BOTY,CW,'Landscape predictors',
    ['Herbaceous NDVI','Human use intensity','Proportion herbaceous,','  shrub, wetland and','  harvest land cover','Terrain ruggedness'],F_DATA)
box(X2,BOTY-74,CW,'Spatial extraction',
    ['16 nested buffers, 40 m','  to 7,223 m radius','Matched-scale design:','  all predictors in a model','  at the same extent'],F_SPEC,S_SPEC,glyph_rings,150)
box(X3,BOTY,CW,'Predictor screening',
    ['Events-per-predictor','  budget','AIC screening at fixed','  model size','Substitution for','  excluded predictors'],F_SPEC,S_SPEC)

# ---------------- result panel ----------------
RX,RY,RW,RH=800,76,302,516
F2POS,F2NEG,F2ZERO='#0072B2','#660303','#D9D9D9'
o.append(f'<rect x="{RX}" y="{RY}" width="{RW}" height="{RH}" rx="6" fill="{F_RES}" stroke="{A_RES}" stroke-width="1.3"/>')
t(RX+14,RY+23,'Scale of effect',12.5,INK,'start','bold')
for i,ln in enumerate(['Coefficient profile across the 16 extents; non-zero',
                       'domains; the peak extent is the scale of effect']):
    t(RX+14,RY+42+i*14,ln,10.3,MUTE)

labs=['40 m','56 m','80 m','113 m','160 m','226 m','319 m','451 m','638 m','903 m',
      '1.3 km','1.8 km','2.6 km','3.6 km','5.1 km','7.2 km']
px0,px1=RX+46,RX+286
py0,py1=RY+122,RY+228
lo,hi=math.log10(40),math.log10(7223)
Xm=lambda m: px0+(math.log10(m)-lo)/(hi-lo)*(px1-px0)
ymin,ymax=-1.25,0.85
Ym=lambda v: py1-(v-ymin)/(ymax-ymin)*(py1-py0)
step=(px1-px0)/15.0
def diamond(cx,cy,r,fill):
    o.append(f'<polygon points="{cx:.1f},{cy-r} {cx+r:.1f},{cy} {cx:.1f},{cy+r} {cx-r:.1f},{cy}" '
             f'fill="{fill}" stroke="{INK}" stroke-width="1.1" stroke-linejoin="round"/>')
def star(cx,cy,r,fill,sw=1.1):
    pts=[]
    for k in range(10):
        ang=-math.pi/2+k*math.pi/5; rr=r if k%2==0 else r*0.42
        pts.append(f'{cx+rr*math.cos(ang):.1f},{cy+rr*math.sin(ang):.1f}')
    o.append(f'<polygon points="{" ".join(pts)}" fill="{fill}" stroke="{INK}" stroke-width="{sw}" stroke-linejoin="round"/>')

t(px0,py0-26,'Horse, human use intensity, October',9.6,MUTE)
t(px0,py0-12,'each point is a separate model fitted at that extent',8.8,MUTE)
o.append(f'<rect x="{Xm(160)-step/2:.1f}" y="{py0}" width="{step:.1f}" height="{py1-py0}" fill="{F2POS}" fill-opacity="0.13"/>')
o.append(f'<rect x="{Xm(3612)-step/2:.1f}" y="{py0}" width="{3*step:.1f}" height="{py1-py0}" fill="{F2NEG}" fill-opacity="0.17"/>')
band=' '.join(f'{Xm(m):.1f},{Ym(ci[m][1]):.1f}' for m,e,sg in prof)+' '+' '.join(f'{Xm(m):.1f},{Ym(ci[m][0]):.1f}' for m,e,sg in reversed(prof))
o.append(f'<polygon points="{band}" fill="#808080" fill-opacity="0.16"/>')
o.append(f'<line x1="{px0}" y1="{Ym(0):.1f}" x2="{px1}" y2="{Ym(0):.1f}" stroke="{MUTE}" stroke-width="0.9" stroke-dasharray="3 3"/>')
o.append('<path d="'+' '.join(('M' if i==0 else 'L')+f'{Xm(m):.1f},{Ym(e):.1f}' for i,(m,e,sg) in enumerate(prof))+f'" fill="none" stroke="{INK}" stroke-width="1.6"/>')
for m,e,sg in prof:
    o.append(f'<circle cx="{Xm(m):.1f}" cy="{Ym(e):.1f}" r="2.8" fill="{F2POS if sg=="+" else F2NEG if sg=="-" else GREY}"/>')
diamond(Xm(160),Ym(0.339),5.2,F2POS)
star(Xm(7223),Ym(-0.631),7.0,F2NEG,1.0)
o.append(f'<line x1="{px0}" y1="{py0}" x2="{px0}" y2="{py1}" stroke="{INK}" stroke-width="0.9"/>')
o.append(f'<line x1="{px0}" y1="{py1}" x2="{px1}" y2="{py1}" stroke="{INK}" stroke-width="0.9"/>')
t(px0-6,Ym(0.5),'+',10,MUTE,'end'); t(px0-6,Ym(-0.7),'−',10,MUTE,'end')
o.append(f'<text x="{px0-24}" y="{(py0+py1)/2:.0f}" font-size="9.5" fill="{MUTE}" text-anchor="middle" '
         f'font-family="Helvetica,Arial,sans-serif" transform="rotate(-90 {px0-24} {(py0+py1)/2:.0f})">coefficient</text>')
for (m,e,sg),lab in zip(prof,labs):
    x=Xm(m); y=py1+12
    o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="6.6" fill="{MUTE}" text-anchor="end" '
             f'font-family="Helvetica,Arial,sans-serif" transform="rotate(-45 {x:.1f} {y:.1f})">{lab}</text>')

sy=py1+108
o.append(f'<line x1="{(px0+px1)/2:.0f}" y1="{py1+50:.0f}" x2="{(px0+px1)/2:.0f}" y2="{sy-16:.0f}" '
         f'stroke="{A_RES}" stroke-width="2.6" marker-end="url(#a_{A_RES[1:]})"/>')
for m,val in ((160,0.339),(7223,-0.631)):
    o.append(f'<line x1="{Xm(m):.1f}" y1="{Ym(val)+9:.1f}" x2="{Xm(m):.1f}" y2="{py1:.1f}" '
             f'stroke="{MUTE}" stroke-width="0.9" stroke-dasharray="2 3"/>')
    o.append(f'<line x1="{Xm(m):.1f}" y1="{py1+46:.1f}" x2="{Xm(m):.1f}" y2="{sy-3:.1f}" '
             f'stroke="{MUTE}" stroke-width="0.9" stroke-dasharray="2 3"/>')
BARH=13.0
for m,e,sg in prof:
    o.append(f'<line x1="{Xm(m):.1f}" y1="{sy-5:.1f}" x2="{Xm(m):.1f}" y2="{sy+BARH+5:.1f}" stroke="#ececec" stroke-width="0.8"/>')
runs=[]
for i,(m,e,sg) in enumerate(prof):
    if runs and runs[-1][0]==sg and runs[-1][2]==i-1: runs[-1][2]=i
    else: runs.append([sg,i,i])
for sg,i0,i1 in runs:
    x0=Xm(prof[i0][0])-step/2+1.1; x1=Xm(prof[i1][0])+step/2-1.1
    col={'+':F2POS,'-':F2NEG,'0':F2ZERO}[sg]
    o.append(f'<rect x="{x0:.1f}" y="{sy:.1f}" width="{x1-x0:.1f}" height="{BARH}" rx="{BARH/2:.1f}" fill="{col}"/>')
diamond(Xm(160),sy+BARH/2,7.0,F2POS)
star(Xm(7223),sy+BARH/2,9.5,F2NEG)
t(px0-step/2-4,sy+BARH/2+3,'October',8.8,INK,'end')
for (m,e,sg),lab in zip(prof,labs):
    x=Xm(m); y=sy+BARH+13
    o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="6.6" fill="{MUTE}" text-anchor="end" '
             f'font-family="Helvetica,Arial,sans-serif" transform="rotate(-45 {x:.1f} {y:.1f})">{lab}</text>')
t((px0+px1)/2,sy+BARH+52,'Spatial extent (buffer radius)',8.6,MUTE,'middle')

t((px0+px1)/2,sy+BARH+72,'Example profile as in Figure 3',9.8,INK,'middle')
ly=sy+BARH+94
for i,(c,lab) in enumerate([(F2POS,'positive domain'),(F2NEG,'negative domain'),(F2ZERO,'zero (CI spans zero)')]):
    yy=ly+i*15
    o.append(f'<rect x="{RX+16}" y="{yy-8}" width="16" height="8" rx="4" fill="{c}"/>')
    t(RX+38,yy,lab,9.4,MUTE)
diamond(RX+24,ly+45-3.5,7.0,F2POS); t(RX+38,ly+45,'peak within a domain',9.4,MUTE)
star(RX+24,ly+61-3.5,9.5,F2NEG);   t(RX+38,ly+61,'global peak: the scale of effect',9.4,MUTE)

# ---------------- arrows ----------------
arrow([(164,266),(172,266),(172,TOPY+55),(X1-3,TOPY+55)],A_DATA)
arrow([(164,300),(172,300),(172,BOTY+55),(X1-3,BOTY+55)],A_DATA)
arrow([(X1+CW,TOPY+55),(X2-4,TOPY+55)],A_DATA)
arrow([(X1+CW,BOTY+55),(X2-4,BOTY+55)],A_DATA)
arrow([(X2+CW,TOPY+55),(X3-4,TOPY+55)],A_SPEC)
arrow([(X2+CW,BOTY+55),(X3-4,BOTY+55)],A_SPEC)
arrow([(X3+CW/2,BOTY-4),(X3+CW/2,TOPY+h_mf+8)],A_SPEC)
arrow([(X2+CW-46,TOPY+bh(['a']*4,62)),(X2+CW-46,332),(X3+46,332),(X3+46,BOTY-4)],A_SPEC)
t(X2+CW-40,314,'the pooling period decides',8.8,MUTE)
t(X2+CW-40,325,'the predictor budget',8.8,MUTE)
arrow([(X3+CW,TOPY+55),(RX-6,TOPY+55)],A_RES)

# ---------------- legend + caption ----------------
# gy is the baseline of the legend labels; the swatches sit at gy-9.
# The deepest box (Spatial extraction) ends at y = 607, so this leaves a 25 px
# clear band. It was 76 px until 2026-09-03, when two-thirds of that white space
# was removed to shorten the figure for page fit. If a box row gets deeper,
# raise this or the legend will collide with it.
gy=643
LEGEND=[(F_DATA,EDGE,'Data'),(F_SPEC,S_SPEC,'Analytical step'),
        (F_MF,S_MF,'Analysis'),(F_RES,A_RES,'Result')]
GAP=42
def lw(lab): return 0.52*10.4*len(lab)
total=sum(20+7+lw(l) for _,_,l in LEGEND)+GAP*(len(LEGEND)-1)
x=(W-total)/2.0
for f,st,lab in LEGEND:
    o.append(f'<rect x="{x:.1f}" y="{gy-9}" width="20" height="14" rx="3" fill="{f}" stroke="{st}" stroke-width="1.3"/>')
    t(x+27,gy+3,lab,10.4,INK)
    x+=20+7+lw(lab)+GAP
t(W/2,gy+30,'Each species-period is modelled separately at all 16 spatial extents, giving 160 base model fits.',10,MUTE,'middle')
t(W/2,gy+46,'The scale of effect is the extent at which a predictor’s coefficient is strongest, not the best-fitting model.',10,MUTE,'middle')
o.append('</svg>')

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def main():
    stem = sys.argv[1] if len(sys.argv) > 1 else "methods_schematic"
    svg_path = stem + ".svg"
    with open(svg_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(o))
    print("wrote", svg_path)

    try:
        import cairosvg
    except ImportError:
        print("cairosvg not installed; SVG only (this is fine for submission).")
        return
    png_path = stem + ".png"
    cairosvg.svg2png(url=svg_path, write_to=png_path, output_width=1700)
    print("wrote", png_path)


if __name__ == "__main__":
    main()
