"""
Needleman-Wunsch Algorithm for Global Sequence Alignment

This module provides a comprehensive implementation of the Needleman-Wunsch algorithm,
a fundamental dynamic programming approach for global sequence alignment in bioinformatics.

Author: [Your Name]
Date: August 2025
Version: 1.0.0
"""

from typing import Tuple, List, Union, Dict, Any
import sys


class NeedlemanWunsch:
    """
    Implementation of the Needleman-Wunsch algorithm for global sequence alignment.
    
    The algorithm uses dynamic programming to find the optimal alignment between
    two sequences by considering matches, mismatches, and gaps.
    
    Attributes:
        match_score (int): Score awarded for matching characters
        mismatch_score (int): Penalty for mismatching characters
        gap_penalty (int): Penalty for gaps/insertions/deletions
    """
    
    def __init__(self, match_score: int = 2, mismatch_score: int = -1, gap_penalty: int = -2):
        """
        Initialize the Needleman-Wunsch aligner.
        
        Args:
            match_score: Score for matching characters (default: 2)
            mismatch_score: Score for mismatching characters (default: -1)
            gap_penalty: Penalty for gaps (default: -2)
        """
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
    
    def _score_match_mismatch(self, char1: str, char2: str) -> int:
        """
        Calculate the score between two characters.
        
        Args:
            char1: Character from first sequence
            char2: Character from second sequence
            
        Returns:
            Score for the character pair
        """
        if char1 == char2:
            return self.match_score
        else:
            return self.mismatch_score
    
    def _initialize_score_matrix(self, seq1_length: int, seq2_length: int) -> List[List[int]]:
        """
        Initialize the scoring matrix for the Needleman-Wunsch algorithm.
        
        Args:
            seq1_length: Length of the first sequence
            seq2_length: Length of the second sequence
            
        Returns:
            Initialized scoring matrix
        """
        # Create matrix with dimensions (seq1_length + 1) x (seq2_length + 1)
        matrix = [[0 for _ in range(seq2_length + 1)] for _ in range(seq1_length + 1)]
        
        # Initialize first row (gaps in seq1)
        for j in range(1, seq2_length + 1):
            matrix[0][j] = j * self.gap_penalty
        
        # Initialize first column (gaps in seq2)
        for i in range(1, seq1_length + 1):
            matrix[i][0] = i * self.gap_penalty
            
        return matrix
    
    def _run_needleman_wunsch(self, seq1: str, seq2: str) -> Tuple[List[List[str]], List[List[int]]]:
        """
        Execute the main Needleman-Wunsch algorithm.
        
        Args:
            seq1: First sequence to align
            seq2: Second sequence to align
            
        Returns:
            Tuple containing (traceback_matrix, score_matrix)
        """
        seq1_length = len(seq1)
        seq2_length = len(seq2)
        
        # Initialize matrices
        score_matrix = self._initialize_score_matrix(seq1_length, seq2_length)
        traceback_matrix = [[None for _ in range(seq2_length + 1)] for _ in range(seq1_length + 1)]
        
        # Initialize traceback for first row and column
        for i in range(1, seq1_length + 1):
            traceback_matrix[i][0] = 'U'  # Up - gap in seq2
        for j in range(1, seq2_length + 1):
            traceback_matrix[0][j] = 'L'  # Left - gap in seq1
        
        # Fill the matrices using dynamic programming
        for i in range(1, seq1_length + 1):
            for j in range(1, seq2_length + 1):
                char1 = seq1[i-1]
                char2 = seq2[j-1]
                
                # Calculate scores for three possible moves
                diagonal_score = score_matrix[i-1][j-1] + self._score_match_mismatch(char1, char2)
                up_score = score_matrix[i-1][j] + self.gap_penalty  # Gap in seq2
                left_score = score_matrix[i][j-1] + self.gap_penalty  # Gap in seq1
                
                # Choose the maximum score
                max_score = max(diagonal_score, up_score, left_score)
                score_matrix[i][j] = max_score
                
                # Record the direction in traceback matrix
                if max_score == diagonal_score:
                    traceback_matrix[i][j] = 'D'  # Diagonal
                elif max_score == up_score:
                    traceback_matrix[i][j] = 'U'  # Up
                else:
                    traceback_matrix[i][j] = 'L'  # Left
        
        return traceback_matrix, score_matrix
    
    def _get_aligned_sequences(self, traceback_matrix: List[List[str]], seq1: str, seq2: str) -> Tuple[str, str]:
        """
        Reconstruct the aligned sequences using the traceback matrix.
        
        Args:
            traceback_matrix: Matrix containing traceback directions
            seq1: First original sequence
            seq2: Second original sequence
            
        Returns:
            Tuple of aligned sequences (seq1_aligned, seq2_aligned)
        """
        i = len(seq1)
        j = len(seq2)
        seq1_aligned = []
        seq2_aligned = []
        
        # Traceback from bottom-right to top-left
        while i > 0 or j > 0:
            if i == 0:
                # Must move left (gap in seq1)
                seq1_aligned.append('-')
                seq2_aligned.append(seq2[j-1])
                j -= 1
            elif j == 0:
                # Must move up (gap in seq2)
                seq1_aligned.append(seq1[i-1])
                seq2_aligned.append('-')
                i -= 1
            else:
                direction = traceback_matrix[i][j]
                if direction == 'D':  # Diagonal (match/mismatch)
                    seq1_aligned.append(seq1[i-1])
                    seq2_aligned.append(seq2[j-1])
                    i -= 1
                    j -= 1
                elif direction == 'U':  # Up (gap in seq2)
                    seq1_aligned.append(seq1[i-1])
                    seq2_aligned.append('-')
                    i -= 1
                elif direction == 'L':  # Left (gap in seq1)
                    seq1_aligned.append('-')
                    seq2_aligned.append(seq2[j-1])
                    j -= 1
        
        # Reverse the sequences (built backwards during traceback)
        return ''.join(seq1_aligned[::-1]), ''.join(seq2_aligned[::-1])
    
    def align(self, seq1: str, seq2: str) -> Tuple[str, str, int]:
        """
        Perform global sequence alignment using Needleman-Wunsch algorithm.
        
        Args:
            seq1: First sequence to align
            seq2: Second sequence to align
            
        Returns:
            Tuple containing (aligned_seq1, aligned_seq2, alignment_score)
        """
        if not seq1 or not seq2:
            raise ValueError("Both sequences must be non-empty")
        
        traceback_matrix, score_matrix = self._run_needleman_wunsch(seq1, seq2)
        aligned_seq1, aligned_seq2 = self._get_aligned_sequences(traceback_matrix, seq1, seq2)
        alignment_score = score_matrix[len(seq1)][len(seq2)]
        
        return aligned_seq1, aligned_seq2, alignment_score
    
    def align_with_matrices(self, seq1: str, seq2: str) -> Dict[str, Any]:
        """
        Perform alignment and return detailed results including matrices.
        
        Args:
            seq1: First sequence to align
            seq2: Second sequence to align
            
        Returns:
            Dictionary containing alignment results and matrices
        """
        traceback_matrix, score_matrix = self._run_needleman_wunsch(seq1, seq2)
        aligned_seq1, aligned_seq2 = self._get_aligned_sequences(traceback_matrix, seq1, seq2)
        alignment_score = score_matrix[len(seq1)][len(seq2)]
        
        return {
            'seq1_aligned': aligned_seq1,
            'seq2_aligned': aligned_seq2,
            'score': alignment_score,
            'score_matrix': score_matrix,
            'traceback_matrix': traceback_matrix
        }
    
    def print_score_matrix(self, score_matrix: List[List[int]], seq1: str = "", seq2: str = "") -> None:
        """
        Pretty print the scoring matrix.
        
        Args:
            score_matrix: The scoring matrix to display
            seq1: First sequence (for row labels)
            seq2: Second sequence (for column labels)
        """
        print("\nScoring Matrix:")
        print("=" * 50)
        
        # Print column headers
        print("     ", end="")
        if seq2:
            print("   ε", end="")  # Empty string symbol
            for char in seq2:
                print(f"   {char}", end="")
        print()
        
        # Print matrix rows
        for i, row in enumerate(score_matrix):
            if i == 0:
                print("  ε  ", end="")
            elif seq1 and i-1 < len(seq1):
                print(f"  {seq1[i-1]}  ", end="")
            else:
                print("     ", end="")
            
            for val in row:
                print(f"{val:4d}", end="")
            print()
    
    def print_traceback_matrix(self, traceback_matrix: List[List[str]], seq1: str = "", seq2: str = "") -> None:
        """
        Pretty print the traceback matrix with directional arrows.
        
        Args:
            traceback_matrix: The traceback matrix to display
            seq1: First sequence (for row labels)
            seq2: Second sequence (for column labels)
        """
        print("\nTraceback Matrix:")
        print("=" * 50)
        
        # Unicode arrows for better visualization
        arrows = {
            'U': '↑',    # Up - gap in seq2
            'L': '←',    # Left - gap in seq1
            'D': '↖',    # Diagonal - match/mismatch
            None: '•'    # No direction
        }
        
        # Print column headers
        print("     ", end="")
        if seq2:
            print("  ε", end="")
            for char in seq2:
                print(f"  {char}", end="")
        print()
        
        # Print matrix rows
        for i, row in enumerate(traceback_matrix):
            if i == 0:
                print("  ε  ", end="")
            elif seq1 and i-1 < len(seq1):
                print(f"  {seq1[i-1]}  ", end="")
            else:
                print("     ", end="")
            
            for val in row:
                arrow = arrows.get(val, '•')
                print(f"  {arrow}", end="")
            print()


class NeedlemanWunschAdvanced(NeedlemanWunsch):
    """
    Advanced Needleman-Wunsch implementation with transition/transversion scoring.
    
    This class extends the basic implementation to handle different penalties
    for transitions (A↔G, C↔T) versus transversions (A↔T, G↔C, A↔C, G↔T)
    in DNA sequence alignment.
    """
    
    def __init__(self, match_score: int = 2, transition_penalty: int = -1, 
                 transversion_penalty: int = -2, gap_penalty: int = -1):
        """
        Initialize the advanced aligner with transition/transversion scoring.
        
        Args:
            match_score: Score for matching nucleotides
            transition_penalty: Penalty for transitions (A↔G, C↔T)
            transversion_penalty: Penalty for transversions (A↔T, G↔C, A↔C, G↔T)
            gap_penalty: Penalty for gaps
        """
        super().__init__(match_score, transition_penalty, gap_penalty)  # Use transition as default mismatch
        self.transition_penalty = transition_penalty
        self.transversion_penalty = transversion_penalty
    
    def _score_match_mismatch(self, char1: str, char2: str) -> int:
        """
        Calculate score with transition/transversion penalties.
        
        Transitions (purine ↔ purine, pyrimidine ↔ pyrimidine):
        - A ↔ G (purine to purine)
        - C ↔ T (pyrimidine to pyrimidine)
        
        Transversions (purine ↔ pyrimidine):
        - A ↔ T, A ↔ C, G ↔ T, G ↔ C
        
        Args:
            char1: First nucleotide
            char2: Second nucleotide
            
        Returns:
            Appropriate score based on nucleotide relationship
        """
        if char1 == char2:
            return self.match_score
        
        # Define transitions
        transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
        nucleotide_pair = (char1.upper(), char2.upper())
        
        if nucleotide_pair in transitions:
            return self.transition_penalty
        else:
            return self.transversion_penalty


# Convenience functions for quick usage
def align_sequences(seq1: str, seq2: str, match_score: int = 2, mismatch_score: int = -1, 
                   gap_penalty: int = -2) -> Dict[str, Union[Tuple[str, str], int]]:
    """
    Quick function to align two sequences with default parameters.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        match_score: Score for matches
        mismatch_score: Score for mismatches
        gap_penalty: Penalty for gaps
        
    Returns:
        Dictionary with alignment results
    """
    aligner = NeedlemanWunsch(match_score, mismatch_score, gap_penalty)
    aligned_seq1, aligned_seq2, score = aligner.align(seq1, seq2)
    
    return {
        'alignment': (aligned_seq1, aligned_seq2),
        'score': score
    }


def educational_alignment(seq1: str, seq2: str, match_score: int = 1, mismatch_score: int = -1,
                         gap_penalty: int = -2, verbose: bool = True) -> Dict[str, Any]:
    """
    Educational function that shows the alignment process step by step.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        match_score: Score for matches
        mismatch_score: Score for mismatches
        gap_penalty: Penalty for gaps
        verbose: Whether to print detailed steps
        
    Returns:
        Complete alignment results with matrices
    """
    aligner = NeedlemanWunsch(match_score, mismatch_score, gap_penalty)
    
    if verbose:
        print(f"Aligning sequences: '{seq1}' and '{seq2}'")
        print(f"Scoring: Match={match_score}, Mismatch={mismatch_score}, Gap={gap_penalty}")
        print("-" * 60)
    
    result = aligner.align_with_matrices(seq1, seq2)
    
    if verbose:
        print(f"\nFinal Alignment:")
        print(f"Sequence 1: {result['seq1_aligned']}")
        print(f"Sequence 2: {result['seq2_aligned']}")
        print(f"Score: {result['score']}")
        
        aligner.print_score_matrix(result['score_matrix'], seq1, seq2)
        aligner.print_traceback_matrix(result['traceback_matrix'], seq1, seq2)
    
    return result


if __name__ == "__main__":
    # Example usage and testing
    print("Needleman-Wunsch Algorithm Demo")
    print("=" * 40)
    
    # Basic alignment
    aligner = NeedlemanWunsch(match_score=2, mismatch_score=-1, gap_penalty=-2)
    seq1, seq2 = "GCATGCT", "GATTACA"
    
    aligned_seq1, aligned_seq2, score = aligner.align(seq1, seq2)
    print(f"Sequences: {seq1} vs {seq2}")
    print(f"Alignment: {aligned_seq1}")
    print(f"           {aligned_seq2}")
    print(f"Score: {score}")
    
    # Advanced alignment with transition/transversion
    print("\n" + "=" * 40)
    print("Advanced Alignment (Transition/Transversion)")
    
    advanced_aligner = NeedlemanWunschAdvanced(
        match_score=2, transition_penalty=-1, transversion_penalty=-2, gap_penalty=-1
    )
    
    adv_aligned_seq1, adv_aligned_seq2, adv_score = advanced_aligner.align(seq1, seq2)
    print(f"Sequences: {seq1} vs {seq2}")
    print(f"Alignment: {adv_aligned_seq1}")
    print(f"           {adv_aligned_seq2}")
    print(f"Score: {adv_score}")