from causality_sbgnviz_interface import _resource_dir
import causality_agent

ca = causality_agent.CausalityAgent(_resource_dir)

def test_find_causality_targets_akt():
    def check(res):
        assert res == []
    source = {'id': 'AKT1',
              'rel': 'is-phosphorylated-by'}
    res = ca.find_causality_targets(source)
    check(res)


def test_find_causality_targets_braf():
    def check(res):
        res1 = res[0]
        print(res1)
        assert res1.get('rel') == 'is-phosphorylated-by'
        assert res1.get('id2') == 'MAPK1'
        assert res1.get('res2') == ['T']
        assert res1.get('pos2') == ['185']
    source = {'id': 'BRAF',
              'pSite': 'S365S',
              'rel': 'is-phosphorylated-by'}
    res = ca.find_causality_targets(source)
    check(res)


def test_find_causality_targets_mapk1():
    def check(res):
        print(res)
    source = {'id': 'MAPK1',
              'rel': 'phosphorylates'}
    res = ca.find_causality_targets(source)
    check(res)


def test_find_next_correlation_akt():
    # TODO: Implement checking the content of the results here
    res = ca.find_next_unexplained_correlation('AKT1')
    print(res)
    res = ca.find_next_unexplained_correlation('AKT1')
    print(res)
    res = ca.find_next_unexplained_correlation('AKT1')
    print(res)


# TODO: Implement tests for the cases below
# ca.find_next_correlation('AKT1',print_result)
# ca.find_next_correlation('AKT1',print_result)
# ca.find_next_correlation('AKT1',print_result)
# ca.find_next_correlation('AKT1',print_result)
# ca.find_next_correlation('AKT1',print_result)
# # ca.find_next_correlation('AKT1',print_result)
# ca.find_correlation_between('AKT1', 'BRAF')
# ca.find_all_correlations('AKT1')
# print(ca.find_mut_sig('TP53'))
# ca.find_common_upstreams('RAC1', 'RAC2')
# ca.find_common_upstreams(['AKT1', 'BRAF', 'MAPK1'], print_result)
