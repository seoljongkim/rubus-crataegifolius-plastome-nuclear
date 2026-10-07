#!/usr/bin/env python3
"""compose_fig2.py — Fig 2 합성 (a,b 상단 / c,d,e 하단). 패널 라벨은 MDPI 규약대로 소문자."""
import os
from PIL import Image, ImageDraw, ImageFont
BASE = os.environ.get("RUBUS_BASE", "/home/spenkim/rubus_cp")
P = f"{BASE}/figures_revision"
FONT = "/mnt/c/Windows/Fonts/timesbd.ttf"

panels = [("a", "Fig2A_SSR_type.png"), ("b", "Fig2B_SSR_proportion.png"),
          ("c", "Fig2C_tandem_length.png"), ("d", "Fig2D_dispersed_type.png"),
          ("e", "Fig2E_dispersed_length.png")]
ims = {k: Image.open(f"{P}/{f}").convert("RGB") for k, f in panels}

PAD, LAB = 40, 200                      # 패널 간격 / 라벨 영역 높이
def scale(im, h):
    w = int(im.width * h / im.height)
    return im.resize((w, h), Image.LANCZOS)

TOP_H = 1500                           # 상단 행 높이
BOT_H = 1350                           # 하단 행 높이
top = [scale(ims["a"], TOP_H), scale(ims["b"], TOP_H)]
bot = [scale(ims[k], BOT_H) for k in "cde"]

tw = sum(i.width for i in top) + PAD
bw = sum(i.width for i in bot) + PAD * 2
W = max(tw, bw)
H = LAB + TOP_H + PAD + LAB + BOT_H

canvas = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(canvas)
font = ImageFont.truetype(FONT, 150)

def row(imgs, keys, y0, x0):
    x = x0
    for im, k in zip(imgs, keys):
        d.text((x, y0 - LAB + 20), f"({k})", font=font, fill=(0, 0, 0))
        canvas.paste(im, (x, y0))
        x += im.width + PAD

row(top, "ab", LAB, (W - tw) // 2)
row(bot, "cde", LAB + TOP_H + PAD + LAB, (W - bw) // 2)

canvas.save(f"{P}/Fig2_repeat_elements.png", dpi=(300, 300))
canvas.save(f"{P}/Fig2_repeat_elements.tif", dpi=(300, 300), compression="tiff_lzw")
print(f"합성 완료: {canvas.size}  → {P}/Fig2_repeat_elements.png|tif")
