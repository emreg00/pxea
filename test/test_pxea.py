import unittest
import networkx # required for generating mock data
import sys # required to check python version as random generator results are inconsistent across 2 and 3
 
from pxea.utilities.network import calculate_proximity
from pxea.utilities.set_enrichment import get_enrichment_score_and_pval 

class TestUtilityMethods(unittest.TestCase):

    def setUp(self):
        # Mock data: 5 ranked pathways (w.r.t. proximity to a drug), 2 pathways (common to diseases), a path graph of 10 nodes
        self.ranked_list = ["d", "a", "e", "c", "b"]
        self.candidates = ["b", "c"]
        self.nodes_from = [1,2]
        self.nodes_to = [8,9]
        #self.pathway_to_genes = { "a": [1,2,8], "b": [0,2,3,4], "c": [6], "d":[2, 4], "e":[0,8,9]} # 7 in none of them
        self.network = networkx.path_graph(10)
        self.n_random = 1000


    def test_calculate_proximity(self):
        z, d = calculate_proximity(self.network, self.nodes_from, self.nodes_to, n_random=self.n_random)
        if sys.version_info[0] < 3:
            self.assertAlmostEqual(z, 72.5470234714173)
        else:
            self.assertAlmostEqual(z, 79.3296787911)
        self.assertEqual(d, 6.5)


    def test_get_enrichment_score_and_pval(self):
        score, pval = get_enrichment_score_and_pval(self.ranked_list, self.candidates, n_random=self.n_random, alternative="greater")
        self.assertAlmostEqual(score, -2.44948974278)
        self.assertEqual(pval, 1.0)


if __name__ == "__main__":
    unittest.main()

