from pathlib import Path

here = Path(__file__).resolve().parent

VERSION = (here / "VERSION").read_text(encoding="utf-8").strip()
