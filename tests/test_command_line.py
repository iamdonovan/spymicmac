import sys
import unittest
import spymicmac.tools.register_ortho
import spymicmac.tools.join_hexagon


class TestCommandLine(unittest.TestCase):

    def test_register_ortho(self):
        with self.assertRaises(SystemExit):
            sys.argv.append('-h')
            spymicmac.tools.register_ortho.main()

    def test_join_hexagon(self):
        with self.assertRaises(SystemExit):
            sys.argv.append('-h')
            spymicmac.tools.join_hexagon.main()


if __name__ == '__main__':
    unittest.main()
