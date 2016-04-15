""" Test Myrheim-Meyer dimension"""
import dagology as dag

def test_MMD_function():
    f_dict = dag.get_mmd_data(2)
    assert f_dict[0.25] == 2.0
    
def test_MMD_write():
    test_fp = './tmp.pkl'
    f_dict = dag.get_mmd_data(2, filepath=test_fp)
    f_dict_read = dag.get_mmd_data(2, filepath=test_fp)
    print len(f_dict_read)
    
def main():
    test_MMD_write()
    
if __name__ == "__main__":
    main()
