import nswave.version as version


def test_nswave_version_exists():
    assert hasattr(version, "__version__")
    assert isinstance(version.__version__, str)
    assert version.__version__.strip() != ""
