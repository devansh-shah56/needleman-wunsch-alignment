"""
Comprehensive examples demonstrating the Needleman-Wunsch algorithm implementation.

This module contains various examples from the original assignment and additional
test cases to showcase the algorithm's capabilities.
"""

from needleman_wunsch import NeedlemanWunsch, NeedlemanWunschAdvanced, align_sequences, educational_alignment


def example_1_basic_alignment():
    """Example 1: Basic alignment with WHAT vs WHY."""
    print("=" * 60)
    print("EXAMPLE 1: Basic Alignment - WHAT vs WHY")
    print("=" * 60)
    
    aligner = NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_penalty=-2)
    result = aligner.align_with_matrices("WHY", "WHAT")
    
    print(f"Sequence 1: WHY")
    print(f"Sequence 2: WHAT")
    print(f"Parameters: Match=1, Mismatch=-1, Gap=-2")
    print(f"\nAlignment Result:")
    print(f"Sequence 1 aligned: {result['seq1_aligned']}")
    print(f"Sequence 2 aligned: {result['seq2_aligned']}")
    print(f"Alignment Score: {result['score']}")
    
    # Display matrices
    aligner.print_score_matrix(result['score_matrix'], "WHY", "WHAT")
    aligner.print_traceback_matrix(result['traceback_matrix'], "WHY", "WHAT")
    print()


def example_2_dna_sequences():
    """Example 2: DNA sequence alignment with gap penalty -2."""
    print("=" * 60)
    print("EXAMPLE 2: DNA Sequence Alignment - GCATGCT vs GATTACA")
    print("=" * 60)
    
    seq1, seq2 = "GCATGCT", "GATTACA"
    aligner = NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_penalty=-2)
    
    aligned_seq1, aligned_seq2, score = aligner.align(seq1, seq2)
    
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"Parameters: Match=1, Mismatch=-1, Gap=-2")
    print(f"\nAlignment Result:")
    print(f"Sequence 1 aligned: {aligned_seq1}")
    print(f"Sequence 2 aligned: {aligned_seq2}")
    print(f"Alignment Score: {score}")
    print()


def example_3_different_gap_penalty():
    """Example 3: Same sequences with different gap penalty."""
    print("=" * 60)
    print("EXAMPLE 3: Same Sequences with Different Gap Penalty")
    print("=" * 60)
    
    seq1, seq2 = "GCATGCT", "GATTACA"
    aligner = NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_penalty=-1)
    
    aligned_seq1, aligned_seq2, score = aligner.align(seq1, seq2)
    
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"Parameters: Match=1, Mismatch=-1, Gap=-1 (less penalty)")
    print(f"\nAlignment Result:")
    print(f"Sequence 1 aligned: {aligned_seq1}")
    print(f"Sequence 2 aligned: {aligned_seq2}")
    print(f"Alignment Score: {score}")
    print()


def example_4_longer_sequences():
    """Example 4: Longer DNA sequences."""
    print("=" * 60)
    print("EXAMPLE 4: Longer DNA Sequences")
    print("=" * 60)
    
    seq1 = "GTTTGACCAGCC"
    seq2 = "CTGACCCACCGC"
    
    aligner = NeedlemanWunsch(match_score=2, mismatch_score=-1, gap_penalty=-1)
    aligned_seq1, aligned_seq2, score = aligner.align(seq1, seq2)
    
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"Parameters: Match=2, Mismatch=-1, Gap=-1")
    print(f"\nAlignment Result:")
    print(f"Sequence 1 aligned: {aligned_seq1}")
    print(f"Sequence 2 aligned: {aligned_seq2}")
    print(f"Alignment Score: {score}")
    print()


def example_5_transition_transversion():
    """Example 5: Advanced scoring with transitions and transversions."""
    print("=" * 60)
    print("EXAMPLE 5: Transition/Transversion Scoring")
    print("=" * 60)
    
    seq1, seq2 = "GCATGCT", "GATTACA"
    
    # Advanced aligner with different penalties for transitions vs transversions
    aligner = NeedlemanWunschAdvanced(
        match_score=2,
        transition_penalty=-1,      # A↔G, C↔T (less penalty)
        transversion_penalty=-2,    # A↔T, G↔C, A↔C, G↔T (more penalty)
        gap_penalty=-1
    )
    
    aligned_seq1, aligned_seq2, score = aligner.align(seq1, seq2)
    
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"Parameters:")
    print(f"  - Match: 2")
    print(f"  - Transition (A↔G, C↔T): -1")
    print(f"  - Transversion (A↔T, G↔C, etc.): -2")
    print(f"  - Gap: -1")
    print(f"\nAlignment Result:")
    print(f"Sequence 1 aligned: {aligned_seq1}")
    print(f"Sequence 2 aligned: {aligned_seq2}")
    print(f"Alignment Score: {score}")
    print()


def example_6_protein_sequences():
    """Example 6: Protein sequence alignment."""
    print("=" * 60)
    print("EXAMPLE 6: Protein Sequence Alignment")
    print("=" * 60)
    
    # Short protein sequences (amino acid single-letter codes)
    seq1 = "ACDEFGHIK"
    seq2 = "ACDGHIK"
    
    aligner = NeedlemanWunsch(match_score=3, mismatch_score=-2, gap_penalty=-2)
    aligned_seq1, aligned_seq2, score = aligner.align(seq1, seq2)
    
    print(f"Protein 1: {seq1}")
    print(f"Protein 2: {seq2}")
    print(f"Parameters: Match=3, Mismatch=-2, Gap=-2")
    print(f"\nAlignment Result:")
    print(f"Protein 1 aligned: {aligned_seq1}")
    print(f"Protein 2 aligned: {aligned_seq2}")
    print(f"Alignment Score: {score}")
    print()


def example_7_educational_walkthrough():
    """Example 7: Educational walkthrough showing all steps."""
    print("=" * 60)
    print("EXAMPLE 7: Educational Walkthrough")
    print("=" * 60)
    
    # Use the educational alignment function
    result = educational_alignment("CAT", "DOG", match_score=2, mismatch_score=-1, gap_penalty=-1)
    print("\nThis example shows the complete process including matrices.")
    print()


def benchmark_performance():
    """Benchmark the algorithm with sequences of different lengths."""
    print("=" * 60)
    print("PERFORMANCE BENCHMARK")
    print("=" * 60)
    
    import time
    import random
    
    def generate_random_sequence(length):
        """Generate a random DNA sequence."""
        nucleotides = ['A', 'T', 'G', 'C']
        return ''.join(random.choice(nucleotides) for _ in range(length))
    
    aligner = NeedlemanWunsch(match_score=2, mismatch_score=-1, gap_penalty=-1)
    
    lengths = [10, 25, 50, 100]
    
    print("Performance test results:")
    print("Length | Time (seconds) | Score")
    print("-" * 35)
    
    for length in lengths:
        seq1 = generate_random_sequence(length)
        seq2 = generate_random_sequence(length)
        
        start_time = time.time()
        _, _, score = aligner.align(seq1, seq2)
        end_time = time.time()
        
        print(f"{length:6d} | {end_time - start_time:12.6f} | {score:5d}")
    
    print()


def comparison_standard_vs_advanced():
    """Compare standard vs advanced (transition/transversion) scoring."""
    print("=" * 60)
    print("COMPARISON: Standard vs Advanced Scoring")
    print("=" * 60)
    
    seq1, seq2 = "ATCGATCG", "ATGCATGC"
    
    # Standard scoring
    standard_aligner = NeedlemanWunsch(match_score=2, mismatch_score=-1, gap_penalty=-1)
    std_seq1, std_seq2, std_score = standard_aligner.align(seq1, seq2)
    
    # Advanced scoring
    advanced_aligner = NeedlemanWunschAdvanced(
        match_score=2, transition_penalty=-1, transversion_penalty=-2, gap_penalty=-1
    )
    adv_seq1, adv_seq2, adv_score = advanced_aligner.align(seq1, seq2)
    
    print(f"Original sequences:")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print()
    
    print("Standard Scoring (all mismatches = -1):")
    print(f"Aligned 1: {std_seq1}")
    print(f"Aligned 2: {std_seq2}")
    print(f"Score: {std_score}")
    print()
    
    print("Advanced Scoring (transitions=-1, transversions=-2):")
    print(f"Aligned 1: {adv_seq1}")
    print(f"Aligned 2: {adv_seq2}")
    print(f"Score: {adv_score}")
    print()


def run_all_examples():
    """Run all examples sequentially."""
    print("🧬 NEEDLEMAN-WUNSCH ALGORITHM EXAMPLES")
    print("🔬 Comprehensive demonstration of sequence alignment")
    print()
    
    examples = [
        example_1_basic_alignment,
        example_2_dna_sequences,
        example_3_different_gap_penalty,
        example_4_longer_sequences,
        example_5_transition_transversion,
        example_6_protein_sequences,
        example_7_educational_walkthrough,
        comparison_standard_vs_advanced,
        benchmark_performance
    ]
    
    for i, example in enumerate(examples, 1):
        try:
            example()
            if i < len(examples):
                input("Press Enter to continue to the next example...")
                print("\n" * 2)
        except KeyboardInterrupt:
            print("\nExiting examples...")
            break
        except Exception as e:
            print(f"Error in example {i}: {e}")
            continue


if __name__ == "__main__":
    run_all_examples()