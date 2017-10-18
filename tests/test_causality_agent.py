import json
from kqml import KQMLList
from indra.statements import stmts_from_json
from causality_sbgnviz_interface import _resource_dir
import causality_agent
from causality_module import CausalityModule
from tests.integration import _IntegrationTest
from tests.util import ekb_kstring_from_text, get_request

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
        assert res1.get('mods2') == [{'residue': 'T', 'position': '185',
                                      'mod_type': 'phosphorylation',
                                      'is_modified': True}]
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


class TestCausalPath(_IntegrationTest):
    def __init__(self, *args):
        super(TestCausalPath, self).__init__(CausalityModule)

    def create_message(self):
        source = ekb_kstring_from_text('MAPK1')
        target = ekb_kstring_from_text('JUND')
        content = KQMLList('FIND-CAUSAL-PATH')
        content.set('source', source)
        content.set('target', target)
        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS'
        paths = output.gets('paths')
        jd = json.loads(paths)
        stmts = stmts_from_json(jd)
        assert len(stmts) == 1
        assert stmts[0].enz.name == 'MAPK1'
        assert stmts[0].sub.name == 'JUND'
        assert stmts[0].residue == 'S'
        assert stmts[0].position == '100'


class TestCausalityTarget(_IntegrationTest):
    def __init__(self, *args):
        super(TestCausalityTarget, self).__init__(CausalityModule)

    def create_message(self):
        target = ekb_kstring_from_text('MAPK1')
        content = KQMLList('FIND-CAUSALITY-TARGET')
        content.set('target', target)
        content.sets('type', 'phosphorylates')
        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        paths = output.gets('paths')
        jd = json.loads(paths)
        stmts = stmts_from_json(jd)
        assert len(stmts) == 19
        assert stmts[0].sub.name == 'EIF4EBP1'
        assert stmts[0].residue == 'S'
        assert stmts[0].position == '65'


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
