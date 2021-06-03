import pytest

@pytest.fixture
def create_metadata_fixed(tmp_path):
    d = tmp_path / "sub"
    d.makedir()
    
