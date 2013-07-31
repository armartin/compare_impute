"""
analysis.py
author: Gerard Tse (gerardtse@gmail.com)

This file defines all sets of parameters that are to generate reference and test sets.
"""

DEPTH = 1
CHROMOSOMES = [22]

SAMPLED_DATA_BASE = '/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/workspace'

T_REF = "ref"
T_TEST_TRUTH = "test_truth"
T_TEST_AFFY = "test_affy"
T_TEST_ILLUMINA = "test_illumina"
T_TEST_OMNI = "test_omni"
T_TEST_EXOME = "test_exome"
T_TEST_AFFY_EXOME = "test_affy_exome" #am

class Analysis(object):
    
    ALL = []
    DATASET = [
        T_REF,
        T_TEST_TRUTH,
        T_TEST_AFFY,
        T_TEST_ILLUMINA,
        T_TEST_OMNI,
        T_TEST_EXOME,
        T_TEST_AFFY_EXOME #am
        ]

    def display_name(self):
        return self.path()

    def path(self):
        raise AssertionError("Subclass must implement path")

    def dirname(self, iteration):
        # No need to an extra layer if only one level
        if DEPTH > 1:
            return "%s/%s/%02d" % (SAMPLED_DATA_BASE, self.path(), iteration + 1)
        else:
            return "%s/%s" % (SAMPLED_DATA_BASE, self.path())

    def vcf_filename(self, vcf_type, chrom, iteration = 1):
        assert vcf_type in Analysis.DATASET
        return "%s/ch%02d-%s.vcf.gz" % (self.dirname(iteration), chrom, vcf_type.lower())

    def inds_filename(self, vcf_type, iteration = 1):
        assert vcf_type in Analysis.DATASET
        if vcf_type == T_REF:
            return "%s/ref.inds" % self.dirname(iteration)
        else:
            return "%s/test.inds" % self.dirname(iteration)

    def iterate(self):
        for chrom in CHROMOSOMES:
            for iteration in xrange(DEPTH):
                yield (chrom, iteration + 1)

class SampleSizeAnalysis(Analysis):
    def __init__(self, test_size, ref_size):
        self.testsize_ = test_size
        self.refsize_ = ref_size

    def path(self):
        return "samplesize/test=%d/ref=%d" % (self.testsize_, self.refsize_)

    def testsize(self):
        return self.testsize_

    def refsize(self):
        return self.refsize_

    ALL = []

SampleSizeAnalysis.ALL = [
        SampleSizeAnalysis(5, 5),   # Toy data for testing pipelines
        SampleSizeAnalysis(92, 63),
        SampleSizeAnalysis(92, 125),
        SampleSizeAnalysis(92, 250),
        SampleSizeAnalysis(92, 500),
        SampleSizeAnalysis(92, 1000),
        SampleSizeAnalysis(300, 62),
        SampleSizeAnalysis(300, 125),
        SampleSizeAnalysis(300, 250),
        SampleSizeAnalysis(300, 500),
        SampleSizeAnalysis(500, 62),
        SampleSizeAnalysis(500, 125),
        SampleSizeAnalysis(500, 250),
        SampleSizeAnalysis(500, 500),
    ]

_AF = 'Africa'
_AS = 'Asia'
_EU = 'Europe'
_AM = 'Americas'
POPULATIONS = {
    'ASW': _AF,
    'LWK': _AF,
    'YRI': _AF,
    'CHB': _AS,
    'CHS': _AS,
    'JPT': _AS,
    'CEU': _EU,
    'FIN': _EU,
    'GBR': _EU,
    'IBS': _EU,
    'TSI': _EU,
    'CLM': _AM,
    'MXL': _AM,
    'PUR': _AM
}

class BasePopulationAnalysis(Analysis):
    def __init__(self, pop, type_name):
        self.pop_ = pop
        self.type_name = type_name

    def path(self):
        return "population/%s/%s" % (self.type_name, self.pop_)

    def pop(self):
        return self.pop_

class HoldOnePopAnalysis(BasePopulationAnalysis):
    ALL = []
    def __init__(self, pop):
        BasePopulationAnalysis.__init__(self, pop, "ref_excludes_pop")
HoldOnePopAnalysis.ALL = [HoldOnePopAnalysis(pop) for pop in POPULATIONS.keys()]

class HoldNonePopAnalysis(BasePopulationAnalysis):
    ALL = []
    def __init__(self, pop):
        BasePopulationAnalysis.__init__(self, pop, "ref_includes_pop")
HoldNonePopAnalysis.ALL = [HoldNonePopAnalysis(pop) for pop in POPULATIONS.keys()]

class RegionBiasedControlAnalysis(BasePopulationAnalysis):
    ALL = []
    def __init__(self, pop):
        BasePopulationAnalysis.__init__(self, pop, "ref_unbiased")
RegionBiasedControlAnalysis.ALL = [RegionBiasedControlAnalysis(pop) for pop in POPULATIONS.keys()]

class RegionBiasedAnalysis(BasePopulationAnalysis):
    ALL = []
    def __init__(self, pop):
        BasePopulationAnalysis.__init__(self, pop, "ref_biased")

    def regional_pops(self):
        return [p for p in POPULATIONS.keys() if p != self.pop() and POPULATIONS[p] == POPULATIONS[self.pop()]]
RegionBiasedAnalysis.ALL = [RegionBiasedAnalysis(pop) for pop in POPULATIONS.keys()]

Analysis.ALL = SampleSizeAnalysis.ALL + HoldOnePopAnalysis.ALL + HoldNonePopAnalysis.ALL + RegionBiasedAnalysis.ALL + RegionBiasedControlAnalysis.ALL
