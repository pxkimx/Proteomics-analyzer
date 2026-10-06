"""Draw the app icon (a mass spectrum on a dark teal tile) and write AppIcon.icns."""
import os
import subprocess
import sys
import tempfile

from PIL import Image, ImageDraw, ImageFilter


def draw(size: int) -> Image.Image:
    s = size
    img = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    tile = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    d = ImageDraw.Draw(tile)
    m = int(s * 0.06)
    d.rounded_rectangle((m, m, s - m, s - m), radius=int(s * 0.22), fill=(26, 38, 56, 255))
    glow = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    gd = ImageDraw.Draw(glow)
    gd.ellipse((s * .15, s * .1, s * .85, s * .8), fill=(217, 69, 43, 90))
    glow = glow.filter(ImageFilter.GaussianBlur(s * 0.12))
    tile = Image.alpha_composite(tile, glow)
    mask = Image.new("L", (s, s), 0)
    ImageDraw.Draw(mask).rounded_rectangle((m, m, s - m, s - m), radius=int(s * 0.22), fill=255)
    img.paste(tile, (0, 0), mask)
    d = ImageDraw.Draw(img)
    base = s * 0.76
    d.line((s * .16, base, s * .84, base), fill=(240, 195, 106, 140), width=max(2, int(s * .008)))
    bars = [(.22, .20), (.31, .42), (.40, .26), (.50, .62), (.60, .34), (.69, .48), (.78, .22)]
    for i, (x, h) in enumerate(bars):
        top = base - h * s * 0.75
        col = (240, 195, 106, 255) if i == 3 else ((217, 69, 43, 255) if i % 2 else (232, 135, 58, 255))
        d.line((s * x, base, s * x, top), fill=col, width=int(s * 0.045))
        r = s * 0.0225
        d.ellipse((s * x - r, top - r, s * x + r, top + r), fill=col)
    return img


def main(out: str) -> None:
    with tempfile.TemporaryDirectory() as t:
        iconset = os.path.join(t, "AppIcon.iconset")
        os.makedirs(iconset)
        master = draw(1024)
        for base in (16, 32, 128, 256, 512):
            master.resize((base, base), Image.LANCZOS).save(f"{iconset}/icon_{base}x{base}.png")
            master.resize((base * 2, base * 2), Image.LANCZOS).save(f"{iconset}/icon_{base}x{base}@2x.png")
        subprocess.run(["iconutil", "-c", "icns", iconset, "-o", out], check=True)


if __name__ == "__main__":
    main(sys.argv[1])
