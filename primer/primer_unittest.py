
import unittest

from primer import extract_coding_region, extract_seq_from_fasta, calculate_annealing_temp, design_primer

class NamesTestCase(unittest.TestCase):

  def test_coding_region(self):
      gene = extract_coding_region("AGTACATGAAAGGGTTTCCCGAACCGTAAGTTTAAGCT")
      self.assertEqual(gene, "ATGAAAGGGTTTCCCGAACCGTAA")
      gene = extract_coding_region("AGTACATGAAAGGGTTTCCCGAACCGTGAGTTTAAGCT") # Testing different stop codons
      self.assertEqual(gene, "ATGAAAGGGTTTCCCGAACCGTGA")
      gene = extract_coding_region("AGTACATGAAAGGGTTTCCCGAACCGTAGGTTTAAGCT") # Testing different stop codons
      self.assertEqual(gene, "ATGAAAGGGTTTCCCGAACCGTAG")

      gene = extract_coding_region("AGTACaTgaaaGgGTTTCCCGAACCGTaagTTTAAGCT") # Testing case sensitivity
      self.assertEqual(gene, "ATGAAAGGGTTTCCCGAACCGTAA")

      gene = extract_coding_region("AGTACATGTGAGTTTAAGCT")
      self.assertEqual(gene, "ATGTGA")
      gene = extract_coding_region("ATGAAAGGGTTTCCCGAACCGTGA")
      self.assertEqual(gene, "ATGAAAGGGTTTCCCGAACCGTGA")

      gene = extract_coding_region("AGTACAGGAAAGGGTTTCCCGAACCGTGAGTTTAAGCT") # Testing absence of start codon
      self.assertEqual(gene, "")
      gene = extract_coding_region("AGTACATGAAAGGGTTTCCCGAACCGTCACTTTCAGCT")  # Testing absence of stop codon
      self.assertEqual(gene, "ATGAAAGGGTTTCCCGAACCGTCACTTTCAGCT") # HIER FUNKTION ÜBERDENKEN




unittest.main()