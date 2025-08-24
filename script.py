# Let me extract and analyze the key components from the notebook
# I'll create a clean Python implementation based on the notebook code

def analyze_notebook_structure():
    """
    Analyze the structure and identify the main components for the repository
    """
    
    components = {
        "core_functions": [
            "score_match_mismatch",
            "initialize_score_matrix_nw", 
            "run_needleman_wunsch",
            "get_aligned_seqs"
        ],
        "enhanced_functions": [
            "score_match_mismatch2",
            "run_needleman_wunsch2"
        ],
        "visualization_functions": [
            "pretty_print_score_matrix",
            "pretty_print_traceback_matrix"
        ],
        "test_cases": [
            "WHAT vs WHY",
            "GCATGCT vs GATTACA", 
            "GTTTGACCAGCC vs CTGACCCACCGC",
            "Transition/transversion scoring examples"
        ],
        "key_features": [
            "Global sequence alignment",
            "Dynamic programming implementation",
            "Customizable scoring matrices",
            "Traceback functionality", 
            "Matrix visualization",
            "Transition/transversion scoring"
        ]
    }
    
    return components

structure = analyze_notebook_structure()
print("Repository Structure Analysis:")
print("="*50)
for category, items in structure.items():
    print(f"\n{category.upper()}:")
    for item in items:
        print(f"  - {item}")