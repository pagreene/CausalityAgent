from causality_agent import _resource_dir
import database_handler

db = database_handler.DatabaseHandler(_resource_dir)

def test_find_causality_targets_akt():
    def check(res):
        assert res == []
    source = {'id': 'AKT1',
              'rel': 'is-phosphorylated-by'}
    db.find_causality_targets(source, check)


def test_find_causality_targets_braf():
    def check(res):
        res1 = res[0]
        assert res1.get('rel') == 'is-phosphorylated-by'
        assert res1.get('id2') == 'MAPK1'
        assert res1.get('res2') == 'T'
        assert res1.get('pos2') == '185T'
    source = {'id': 'BRAF',
              'pSite': 'S365S',
              'rel': 'is-phosphorylated-by'}
    db.find_causality_targets(source, check)


def test_find_next_correlation_akt():
    # TODO: Implement checking the content of the results here
    res = db.find_next_unexplained_correlation('AKT1')
    print(res)
    res = db.find_next_unexplained_correlation('AKT1')
    print(res)
    res = db.find_next_unexplained_correlation('AKT1')
    print(res)


# TODO: Implement tests for the cases below
# db.find_next_correlation('AKT1',print_result)
# db.find_next_correlation('AKT1',print_result)
# db.find_next_correlation('AKT1',print_result)
# db.find_next_correlation('AKT1',print_result)
# db.find_next_correlation('AKT1',print_result)
# # db.find_next_correlation('AKT1',print_result)
# db.find_correlation_between('AKT1', 'BRAF')
# db.find_all_correlations('AKT1')
# print(db.find_mut_sig('TP53'))
# db.find_common_upstreams('RAC1', 'RAC2')
# db.find_common_upstreams(['AKT1', 'BRAF', 'MAPK1'], print_result)
