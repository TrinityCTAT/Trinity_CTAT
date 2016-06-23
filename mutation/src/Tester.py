
__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2015"
__credits__ = ["Timothy Tickle", "Brian Haas"]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"


import FunctionalTester
import ScriptTester
import unittest


# Calls all unit tests as a regression suite.
suite = unittest.TestSuite()
# To run FunctionalTester. Please update paths at the top of the script.
#suite.addTest(FunctionalTester.suite())
suite.addTest(ScriptTester.suite())

runner = unittest.TextTestRunner()
runner.run(suite)
