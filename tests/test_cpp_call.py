from bisr_dpw_cpp_routines import add, strongly_connected_components

def test_add_call():

    assert(add(1,2)==3)

def test_dict_print():

    d = {}

    d[0] = [1,3,5]
    d[1] = [4,7,9]

    strongly_connected_components(d, d)

if __name__=='__main__':
    test_dict_print()
