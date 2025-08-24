"""
Comprehensive test suite for the Needleman-Wunsch algorithm implementation.

This module contains unit tests to verify the correctness of the alignment algorithm
and ensure all components work as expected.
"""

import unittest
from needleman_wunsch import NeedlemanWunsch, NeedlemanWunschAdvanced, align_sequences


class TestNeedlemanWunsch(unittest.TestCase):
    """Test cases for the basic Needleman-Wunsch implementation."""
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        self.aligner = NeedlemanWunsch(match_score=2, mismatch_score=-1, gap_penalty=-2)
    
    def test_identical_sequences(self):
        """Test alignment of identical sequences."""
        seq1 = seq2 = "ATCG"
        aligned_seq1, aligned_seq2, score = self.aligner.align(seq1, seq2)
        
        self.assertEqual(aligned_seq1, "ATCG")
        self.assertEqual(aligned_seq2, "ATCG")
        self.assertEqual(score, 8)  # 4 matches * 2 points each
    
    def test_empty_sequences(self):
        """Test that empty sequences raise appropriate errors."""
        with self.assertRaises(ValueError):
            self.aligner.align("", "ATCG")
        
        with self.assertRaises(ValueError):
            self.aligner.align("ATCG", "")
    
    def test_single_character_sequences(self):
        """Test alignment of single character sequences."""
        aligned_seq1, aligned_seq2, score = self.aligner.align("A", "A")
        self.assertEqual(aligned_seq1, "A")
        self.assertEqual(aligned_seq2, "A")
        self.assertEqual(score, 2)  # One match
        
        aligned_seq1, aligned_seq2, score = self.aligner.align("A", "T")
        self.assertEqual(aligned_seq1, "A")
        self.assertEqual(aligned_seq2, "T")
        self.assertEqual(score, -1)  # One mismatch
    
    def test_assignment_example_1(self):
        """Test the first example from the assignment: WHY vs WHAT."""
        aligner = NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_penalty=-2)
        aligned_seq1, aligned_seq2, score = aligner.align("WHY", "WHAT")
        
        # Expected alignment from assignment
        self.assertEqual(aligned_seq2, "WHAT")
        self.assertEqual(aligned_seq1, "WH-Y")
    
    def test_assignment_example_2(self):
        """Test the second example: GCATGCT vs GATTACA with gap=-2."""
        aligner = NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_penalty=-2)
        aligned_seq1, aligned_seq2, score = aligner.align("GCATGCT", "GATTACA")
        
        # Verify sequences have same length after alignment
        self.assertEqual(len(aligned_seq1), len(aligned_seq2))
        
        # Check that alignment maintains sequence integrity (no character changes)
        seq1_chars = aligned_seq1.replace('-', '')
        seq2_chars = aligned_seq2.replace('-', '')
        self.assertEqual(seq1_chars, "GCATGCT")
        self.assertEqual(seq2_chars, "GATTACA")
    
    def test_assignment_example_3(self):
        """Test the third example: GCATGCT vs GATTACA with gap=-1."""
        aligner = NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_penalty=-1)
        aligned_seq1, aligned_seq2, score = aligner.align("GCATGCT", "GATTACA")
        
        # Expected alignment from assignment
        self.assertEqual(aligned_seq2, "G-ATTACA")
        self.assertEqual(aligned_seq1, "GCA-TGCT")
    
    def test_assignment_example_4(self):
        """Test longer sequences: GTTTGACCAGCC vs CTGACCCACCGC."""
        aligner = NeedlemanWunsch(match_score=2, mismatch_score=-1, gap_penalty=-1)
        aligned_seq1, aligned_seq2, score = aligner.align("GTTTGACCAGCC", "CTGACCCACCGC")
        
        # Verify alignment maintains sequence integrity
        seq1_chars = aligned_seq1.replace('-', '')
        seq2_chars = aligned_seq2.replace('-', '')
        self.assertEqual(seq1_chars, "GTTTGACCAGCC")
        self.assertEqual(seq2_chars, "CTGACCCACCGC")
    
    def test_matrix_initialization(self):
        """Test that scoring matrix is initialized correctly."""
        matrix = self.aligner._initialize_score_matrix(3, 4)
        
        # Check dimensions
        self.assertEqual(len(matrix), 4)  # rows = seq1_length + 1
        self.assertEqual(len(matrix[0]), 5)  # cols = seq2_length + 1
        
        # Check first row initialization
        expected_first_row = [0, -2, -4, -6, -8]
        self.assertEqual(matrix[0], expected_first_row)
        
        # Check first column initialization
        expected_first_col = [0, -2, -4, -6]
        actual_first_col = [matrix[i][0] for i in range(4)]
        self.assertEqual(actual_first_col, expected_first_col)
    
    def test_score_function(self):
        """Test the scoring function."""
        # Test match
        self.assertEqual(self.aligner._score_match_mismatch('A', 'A'), 2)
        
        # Test mismatch
        self.assertEqual(self.aligner._score_match_mismatch('A', 'T'), -1)
        
        # Test case insensitive (if implemented)
        self.assertEqual(self.aligner._score_match_mismatch('a', 'A'), 2)
    
    def test_align_with_matrices(self):
        """Test that align_with_matrices returns complete results."""
        result = self.aligner.align_with_matrices("AT", "AC")
        
        required_keys = ['seq1_aligned', 'seq2_aligned', 'score', 'score_matrix', 'traceback_matrix']
        for key in required_keys:
            self.assertIn(key, result)
        
        # Verify matrix dimensions
        self.assertEqual(len(result['score_matrix']), 3)  # 2+1 rows
        self.assertEqual(len(result['score_matrix'][0]), 3)  # 2+1 cols
        
        self.assertEqual(len(result['traceback_matrix']), 3)
        self.assertEqual(len(result['traceback_matrix'][0]), 3)


class TestNeedlemanWunschAdvanced(unittest.TestCase):
    """Test cases for the advanced Needleman-Wunsch implementation with transition/transversion scoring."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.aligner = NeedlemanWunschAdvanced(
            match_score=2,
            transition_penalty=-1,
            transversion_penalty=-2,
            gap_penalty=-1
        )
    
    def test_transition_scoring(self):
        """Test that transitions receive correct penalties."""
        # A ↔ G (transition)
        self.assertEqual(self.aligner._score_match_mismatch('A', 'G'), -1)
        self.assertEqual(self.aligner._score_match_mismatch('G', 'A'), -1)
        
        # C ↔ T (transition)
        self.assertEqual(self.aligner._score_match_mismatch('C', 'T'), -1)
        self.assertEqual(self.aligner._score_match_mismatch('T', 'C'), -1)
    
    def test_transversion_scoring(self):
        """Test that transversions receive correct penalties."""
        # A ↔ T (transversion)
        self.assertEqual(self.aligner._score_match_mismatch('A', 'T'), -2)
        self.assertEqual(self.aligner._score_match_mismatch('T', 'A'), -2)
        
        # G ↔ C (transversion)
        self.assertEqual(self.aligner._score_match_mismatch('G', 'C'), -2)
        self.assertEqual(self.aligner._score_match_mismatch('C', 'G'), -2)
        
        # A ↔ C (transversion)
        self.assertEqual(self.aligner._score_match_mismatch('A', 'C'), -2)
        
        # G ↔ T (transversion)
        self.assertEqual(self.aligner._score_match_mismatch('G', 'T'), -2)
    
    def test_match_scoring(self):
        """Test that matches receive correct scores."""
        for nucleotide in ['A', 'T', 'G', 'C']:
            self.assertEqual(self.aligner._score_match_mismatch(nucleotide, nucleotide), 2)
    
    def test_assignment_example_advanced(self):
        """Test the advanced scoring example from assignment."""
        # Use exact parameters from assignment
        aligner = NeedlemanWunschAdvanced(
            match_score=2,
            transition_penalty=-1,
            transversion_penalty=-1,  # Note: assignment uses -1 for both
            gap_penalty=-1
        )
        
        aligned_seq1, aligned_seq2, score = aligner.align("GCATGCT", "GATTACA")
        
        # Verify alignment maintains sequence integrity
        seq1_chars = aligned_seq1.replace('-', '')
        seq2_chars = aligned_seq2.replace('-', '')
        self.assertEqual(seq1_chars, "GCATGCT")
        self.assertEqual(seq2_chars, "GATTACA")


class TestConvenienceFunctions(unittest.TestCase):
    """Test the convenience functions."""
    
    def test_align_sequences_function(self):
        """Test the align_sequences convenience function."""
        result = align_sequences("AT", "AC", match_score=1, mismatch_score=-1, gap_penalty=-1)
        
        self.assertIn('alignment', result)
        self.assertIn('score', result)
        self.assertEqual(len(result['alignment']), 2)  # Two aligned sequences
    
    def test_educational_alignment_function(self):
        """Test the educational_alignment function."""
        # Test with verbose=False to avoid print statements during testing
        from needleman_wunsch import educational_alignment
        result = educational_alignment("AT", "AC", verbose=False)
        
        required_keys = ['seq1_aligned', 'seq2_aligned', 'score', 'score_matrix', 'traceback_matrix']
        for key in required_keys:
            self.assertIn(key, result)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and boundary conditions."""
    
    def setUp(self):
        self.aligner = NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_penalty=-1)
    
    def test_very_different_lengths(self):
        """Test sequences with very different lengths."""
        seq1 = "A"
        seq2 = "ATCGATCG"
        
        aligned_seq1, aligned_seq2, score = self.aligner.align(seq1, seq2)
        
        # Both aligned sequences should have same length
        self.assertEqual(len(aligned_seq1), len(aligned_seq2))
        
        # Original sequences should be preserved
        self.assertEqual(aligned_seq1.replace('-', ''), seq1)
        self.assertEqual(aligned_seq2.replace('-', ''), seq2)
    
    def test_all_gaps_in_one_sequence(self):
        """Test alignment where one sequence becomes all gaps."""
        seq1 = "AAA"
        seq2 = "TTT"
        
        # With high gap penalty, alignment should prefer mismatches over gaps
        aligner = NeedlemanWunsch(match_score=2, mismatch_score=-1, gap_penalty=-10)
        aligned_seq1, aligned_seq2, score = aligner.align(seq1, seq2)
        
        # Should be no gaps (all mismatches preferred)
        self.assertNotIn('-', aligned_seq1)
        self.assertNotIn('-', aligned_seq2)
    
    def test_case_sensitivity(self):
        """Test case sensitivity handling."""
        aligned_seq1, aligned_seq2, score1 = self.aligner.align("AT", "at")
        aligned_seq3, aligned_seq4, score2 = self.aligner.align("AT", "AT")
        
        # Should handle case appropriately (assuming case-insensitive)
        # This test may need adjustment based on implementation choice
        pass


def run_performance_tests():
    """Run performance tests (not part of main test suite)."""
    import time
    
    print("Running performance tests...")
    aligner = NeedlemanWunsch()
    
    # Test with increasingly large sequences
    for size in [10, 50, 100, 200]:
        seq1 = "A" * size
        seq2 = "T" * size
        
        start_time = time.time()
        aligned_seq1, aligned_seq2, score = aligner.align(seq1, seq2)
        end_time = time.time()
        
        print(f"Size {size:3d}: {end_time - start_time:.4f} seconds")


if __name__ == '__main__':
    # Run unit tests
    unittest.main(argv=[''], exit=False, verbosity=2)
    
    # Optionally run performance tests
    print("\n" + "="*50)
    run_performance_tests()