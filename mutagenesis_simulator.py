# Simple amino acid hydrophobicity scale (Kyte-Doolittle scale)
HYDROPHOBICITY = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}


def simulate_point_mutation(sequence, position, new_aa):
  """Simulates a point mutation in a protein sequence and compares the

  average hydrophobicity.
  """
  seq_list = list(sequence)
  index = position - 1  # 1-based indexing for convenience

  if index < 0 or index >= len(seq_list):
    return "Error: Position is out of sequence bounds."

  old_aa = seq_list[index]
  seq_list[index] = new_aa
  mutated_sequence = "".join(seq_list)

  # Calculate total hydrophobicity score
  old_score = sum(HYDROPHOBICITY.get(aa, 0) for aa in sequence) / len(sequence)
  new_score = sum(HYDROPHOBICITY.get(aa, 0) for aa in mutated_sequence) / len(
      mutated_sequence
  )

  print("\n--- IN SILICO MUTAGENESIS REPORT ---")
  print(f"Original sequence:  {sequence}")
  print(f"Mutated sequence:   {mutated_sequence}")
  print(f"Mutation: {old_aa} -> {new_aa} at position {position}")
  print(f"Average hydrophobicity (Before): {old_score:.2f}")
  print(f"Average hydrophobicity (After):  {new_score:.2f}")
  print("------------------------------------\n")

  return mutated_sequence


if __name__ == "__main__":
  # Sample protein sequence
  sample_protein = "MKTAYIAKQRQISF"

  # Let's try changing A (alanine) at position 4 to I (isoleucine)
  simulate_point_mutation(sample_protein, position=4, new_aa="I")
