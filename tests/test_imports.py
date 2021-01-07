import pytest


class TestImports:
    @pytest.mark.client
    def test_client(self):
        from chemserve._rdkit_imports import Chem, cairosvg, Image

        assert Chem is None
        assert cairosvg is None
        assert Image is None

    @pytest.mark.server
    def test_server(self):
        from chemserve._rdkit_imports import Chem, cairosvg, Image

        assert Chem is not None
        assert cairosvg is not None
        assert Image is not None


if __name__ == "__main__":
    pytest.main()
