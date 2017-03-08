import unittest

from . import test_suite


unittest.TextTestRunner(verbosity=1).run(test_suite)
