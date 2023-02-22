import unittest

def fun(x):
	return x + 1

class MyTest(unittest.TestCase):
	def test(self):
		self.assertEqual(fun(3), 4)

		def test_whole_genome_pipeline:
			make_genome_bitvector('../../../test_files/vicWTct.txt', '../../../test_files/180122Rou_D18-476vicWTctL12align_small.sam')