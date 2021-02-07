from random import random

# Cumulative insertion probabilities based on previous nucleotide
# TODO: unknown data, equal probability currently
custom_synthesis = {
  0: { # A
    0.25: 0,
    0.5: 1,
    0.75: 2,
    1: 3
  },
  1: { # G
    0.25: 0,
    0.5: 1,
    0.75: 2,
    1: 3
  },
  2: { # T
    0.25: 0,
    0.5: 1,
    0.75: 2,
    1: 3
  },
  3: { # C
    0.25: 0,
    0.5: 1,
    0.75: 2,
    1: 3
  }
}

# Based on empirical sequencing data gathered from Louis's project.
custom_sequencing = {
  0: { # A
    0.1667: 0,
    0.4537: 1,
    0.7483: 2,
    1: 3
  },
  1: { # G
    0.2533: 0,
    0.387: 1,
    0.679: 2,
    1: 3
  },
  2: { # T
    0.2064: 0,
    0.4365: 1,
    0.8041: 2,
    1: 3
  },
  3: { # C
    0.3568: 0,
    0.5712: 1,
    0.8055: 2,
    1: 3
  }
}

method_lookup = {
  'custom_synthesis': custom_synthesis,
  'custom_sequencing': custom_sequencing
}

# Function to get the inserted nucleotide based on probabilities.
def getInsertedNucleotide(prev_nuc, method):
  rand = random()
  method_probs = method_lookup[method]
  probabilities = method_probs[prev_nuc]

  for p in probabilities:
    if rand < p:
      return probabilities[p]
