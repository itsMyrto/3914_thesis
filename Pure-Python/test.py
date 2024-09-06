import unittest
from qsieve import quadratic_sieve

class TestQuadraticSieve(unittest.TestCase):
    def check_factorization(self, n, expected_factors):
        factors = quadratic_sieve(n)
        self.assertEqual(sorted(factors), sorted(expected_factors))

    def test_cases(self):
        numbers_and_factors = [
            (15347, (103, 149)),
            (87463, (587, 149)),
            (84923, (521, 163)),
            (832730084101, (830477, 1002713)),
            (851821581119671, (30613469, 27825059)),
            (7875168790028311, (61731809, 127570679)),
            (3207054426926827, (58216163, 55088729)),
            (12096819068999101, (110248973, 109722737)),
            (92092615464081619, (216708257, 424961267)),
            (64157244473449123, (248085329, 258609587)),
            (26408936706025597, (113935733, 231788009)),
            (959125210334783077, (978138209, 980562053)),
            (434686773884327407, (405561743, 1071814049)),
            (210491451967849183, (508788779, 413710877)),
            (2064846507704311861, (1049433389, 1967582249)),
            (6357994389398958601, (1853661119, 3429965879)),
            (3191071089482212003, (1979309237, 1612214519)),
            (16921456439215439701, (5915587277, 2860486313)),
            (93496418013679648963, (7192461209, 12999224507)),
            (54570430399383971173, (6650433887, 8205544379)),
            (27419891463310753159, (6910574183, 3967816673)),
            (14128513504013581789, (3661889003, 3858258263)),
            (861256316295598761961, (32781489479, 26272641359)),
            (533595842543374012417, (32025512693, 16661586269)),
            (200903802201060018373, (14684007509, 13681810097)),
            (7574625114799379190481, (132960332249, 56969059769)),
            (3541904643519702945253, (61052360297, 58014213149)),
            (1753044930908746416511, (26278498889, 66710238599)),
            (32621041168941237031687, (240021603299, 135908771213)),
            (14838142262537816848201, (126097405199, 117672066599)),
            (71346986589122957051491, (271572474749, 262718033759)),
            (854626116323831524991473, (980725424879, 871422413087)),
            (399066586857431823726709, (481294956713, 829151814893)),
            (235311326942746619548591, (467628321869, 503201615339)),
            (113576732865342496692451, (452466318749, 251016988799)),
            (3025577890235798683543591, (1656024390653, 1827012879347)),
            (1822370728996458306753277, (1819935140393, 1001338283189)),
            (46839566299936919234246726809, (100000000105583, 468395662504823)),
            (6172835808641975203638304919691358469663, (11111111111111111011, 555555222777777773333)),
            (982374584994591973035454918323883152991, (16842321112684836809, 58327743451863892199)),
            (3744843080529615909019181510330554205500926021947, (1123456667890987666543211, 3333322225555555777777777))
        ]

        for n, expected_factors in numbers_and_factors:
            with self.subTest(n=n):
                self.check_factorization(n, expected_factors)


if __name__ == '__main__':
    unittest.main()

# 41