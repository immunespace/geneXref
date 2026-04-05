"""Download a pre-built geneXref database from GitHub Releases."""

from __future__ import annotations

import json
import shutil
import urllib.request
from pathlib import Path

GITHUB_REPO = "immunespace/geneXref"
DB_FILENAME = "geneXref.db"
DEFAULT_DB_DIR = Path.home() / ".geneXref"


def _latest_release_asset_url() -> tuple[str, str]:
    """Return (download_url, tag) for the geneXref.db asset in the latest release."""
    api_url = f"https://api.github.com/repos/{GITHUB_REPO}/releases/latest"
    req = urllib.request.Request(api_url, headers={"Accept": "application/vnd.github+json"})
    with urllib.request.urlopen(req) as resp:
        release = json.loads(resp.read())

    tag = release["tag_name"]
    for asset in release.get("assets", []):
        if asset["name"] == DB_FILENAME:
            return asset["browser_download_url"], tag

    raise RuntimeError(
        f"No {DB_FILENAME} asset found in release {tag}. "
        f"See https://github.com/{GITHUB_REPO}/releases"
    )


def download_db(dest_dir: str | Path | None = None) -> Path:
    """Download the latest pre-built database from GitHub Releases.

    Parameters
    ----------
    dest_dir : str or Path, optional
        Directory to save the database to. Defaults to ``~/.geneXref/``.

    Returns
    -------
    Path
        Path to the downloaded database file.
    """
    dest_dir = Path(dest_dir) if dest_dir else DEFAULT_DB_DIR
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest_path = dest_dir / DB_FILENAME

    url, tag = _latest_release_asset_url()
    print(f"Downloading {DB_FILENAME} from release {tag}...")
    print(f"  {url}")

    with urllib.request.urlopen(url) as resp:
        with open(dest_path, "wb") as f:
            shutil.copyfileobj(resp, f)

    size_mb = dest_path.stat().st_size / (1024 * 1024)
    print(f"Saved to {dest_path} ({size_mb:.1f} MB)")
    return dest_path


def main():
    download_db()


if __name__ == "__main__":
    main()
