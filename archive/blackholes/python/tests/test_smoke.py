import bhwave.version as version


def test_version_string():
    assert isinstance(version.__version__, str)
    assert version.__version__
