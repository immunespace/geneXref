def __getattr__(name):
    if name == "GeneMapper":
        from .mapper import GeneMapper
        return GeneMapper
    if name == "download_db":
        from .download import download_db
        return download_db
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = ["GeneMapper", "download_db"]
