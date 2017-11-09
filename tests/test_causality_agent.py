import json
from kqml import KQMLList, KQMLString
from indra.statements import stmts_from_json
from causality_sbgnviz_interface import _resource_dir
import causality_agent
from causality_module import CausalityModule
from tests.integration import _IntegrationTest
from tests.util import ekb_kstring_from_text, ekb_from_text, get_request

ca = causality_agent.CausalityAgent(_resource_dir)


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
        content.sets('type', 'phosphorylation')
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


class TestCausalitySource(_IntegrationTest):
    def __init__(self, *args):
        super(TestCausalitySource, self).__init__(CausalityModule)

    def create_message(self):
        source = ekb_kstring_from_text('BRAF')
        content = KQMLList('FIND-CAUSALITY-SOURCE')
        content.set('source', source)
        content.sets('type', 'phosphorylation')
        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        paths = output.gets('paths')
        jd = json.loads(paths)
        stmts = stmts_from_json(jd)
        assert len(stmts) == 2
        assert stmts[0].sub.name == 'BRAF'
        assert stmts[0].enz.name == 'MAPK1'
        assert stmts[0].residue == 'S'
        assert stmts[0].position == '151'


class TestNextCorrelation(_IntegrationTest):
    def __init__(self, *args):
        super(TestNextCorrelation, self).__init__(CausalityModule)

    def create_message(self):
        source = ekb_kstring_from_text('AKT1')
        content = KQMLList('DATASET-CORRELATED-ENTITY')
        content.set('source', source)
        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        target = output.gets('target')
        correlation = output.gets('correlation')
        explainable = output.gets('explainable')

        assert target == 'BRAF'
        assert correlation == str(0.7610843243760473)
        assert explainable == '\"explainable\"'


class TestCommonUpstreams(_IntegrationTest):
    def __init__(self, *args):
        super(TestCommonUpstreams, self).__init__(CausalityModule)

    def create_message(self):
        content = KQMLList('FIND-COMMON-UPSTREAMS')
        genes = ekb_from_text('AKT1, BRAF and MAPK1')

        content.sets('genes', str(genes))

        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        upstreams = output.gets('upstreams')

        assert 'EGF' in upstreams


class TestMutex(_IntegrationTest):
    def __init__(self, *args):
        super(TestMutex, self).__init__(CausalityModule)

    def create_message(self):
        content = KQMLList('FIND-MUTEX')
        # gene = ekb_from_text('TP53')
        gene = ekb_from_text('TP53')
        content.set('gene', gene)

        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        mutex = output.gets('mutex')
        assert mutex == "[{'score': '0.0', 'group': ['CDH1', 'TP53']}, {'score': '0.0', 'group': ['TP53', 'CDH1']}, {'score': '0.0', 'group': ['GATA3', 'TP53']}, {'score': '0.05', 'group': ['FOXA1', 'TP53', 'GATA3']}]"


class TestMutSigOV(_IntegrationTest):
    def __init__(self, *args):
        super(TestMutSigOV, self).__init__(CausalityModule)

    def create_message(self):
        content = KQMLList('FIND-MUTATION-SIGNIFICANCE')
        gene = ekb_kstring_from_text('TP53')
        disease = ekb_from_text('Ovarian serous cystadenocarcinoma')
        content.set('gene', gene)
        content.set('disease', disease)

        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        mut_sig = output.gets('mutsig')
        assert mut_sig == "highly significant"

class TestMutSigPAAD(_IntegrationTest):
    def __init__(self, *args):
        super(TestMutSigPAAD, self).__init__(CausalityModule)

    def create_message(self):
        content = KQMLList('FIND-MUTATION-SIGNIFICANCE')
        gene = ekb_kstring_from_text('ACTN4')
        disease = ekb_from_text('pancreatic adenocarcinoma')
        content.set('gene', gene)
        content.set('disease', disease)

        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        mut_sig = output.gets('mutsig')
        assert mut_sig == "not significant"


class TestMutSigInvalid(_IntegrationTest):
    def __init__(self, *args):
        super(TestMutSigInvalid, self).__init__(CausalityModule)

    def create_message(self):
        content = KQMLList('FIND-MUTATION-SIGNIFICANCE')
        gene = ekb_kstring_from_text('ACTN4')
        disease = ekb_from_text('lung cancer')
        content.set('gene', gene)
        content.set('disease', disease)

        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        reason = output.gets('reason')
        assert reason == "INVALID_DISEASE"

