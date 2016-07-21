""" Test Myrheim-Meyer dimension"""
import dagology as dag
import unittest
import os

class TestMyrheimMeyer(unittest.TestCase):
    def test_MMD_function():
        f_dict = dag.get_mmd_data(2)
        self.assertEqual(f_dict[0.25], 2.0)
    

	def test_MMD_write():
	    test_fp = './tmp.pkl'
	    f_dict = dag.get_mmd_data(2, filepath=test_fp)
	    f_dict_read = dag.get_mmd_data(2, filepath=test_fp)
	    self.assertTrue(len(f_dict_read) > 1)
	    # tidy up
	    os.remove(test_fp)
    
    
if __name__ == "__main__":
    unittest.main()