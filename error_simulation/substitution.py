from random import random

# Cumulative substitution probabilities
high_pr = {
  0: {
    0.35: 3,
    0.9: 1,
    1: 2
  },
  1: {
    0.18: 2,
    0.25: 3,
    1: 0
  },
  2: {
    0.12: 1,
    0.88: 3,
    1: 0
  },
  3: {
    0.14: 0,
    0.2: 1,
    1: 2
  }
}

low_pr = {
  0: {
    0.21: 3,
    0.79: 1,
    1: 2
  },
  1: {
    0.26: 2,
    0.3: 3,
    1: 0
  },
  2: {
    0.14: 1,
    0.86: 3,
    1: 0
  },
  3: {
    0.18: 0,
    0.21: 1,
    1: 2
  }
}

erlich = {
  0: {
    0.23: 3,
    0.59: 1,
    1: 2
  },
  1: {
    0.71: 2,
    0.84: 3,
    1: 0
  },
  2: {
    0.46: 1,
    0.75: 3,
    1: 0
  },
  3: {
    0.27: 0,
    0.5: 1,
    1: 2
  }
}

goldman = {
  0: {
    0.11: 3,
    0.74: 1,
    1: 2
  },
  1: {
    0.17: 2,
    0.45: 3,
    1: 0
  },
  2: {
    0.45: 1,
    0.78: 3,
    1: 0
  },
  3: {
    0.13: 0,
    0.4: 1,
    1: 2
  }
}

high_pr4t = {
  0: {
    0.3: 3,
    0.8: 1,
    1: 2
  },
  1: {
    0.08: 2,
    0.13: 3,
    1: 0
  },
  2: {
    0.1: 1,
    0.8: 3,
    1: 0
  },
  3: {
    0.07: 0,
    0.12: 1,
    1: 2
  }
}

average = {
  0 : {
    0.25: 3,
    0.75: 1,
    1: 2
  },
  1 : {
    0.07: 3,
    0.74: 0,
    1: 2
  },
  2 : {
    0.28: 1,
    0.85: 3,
    1: 0
  },
  3 : {
    0.175: 0,
    0.25: 1,
    1: 2
  }
}

# Based on empirical sequencing data gathered from Louis's project.
custom_sequencing = {
  0 : { # A
    0.4238: 1,
    0.7068: 2,
    1: 3
  },
  1 : { # G
    0.4426: 0,
    0.7176: 2,
    1: 3
  },
  2 : { # T
    0.3415: 0,
    0.6082: 1,
    1: 3
  },
  3 : { # C
    0.3327: 0,
    0.6009: 1,
    1: 2
  }
}

method_lookup = {
  'high_pr': high_pr,
  'low_pr': low_pr,
  'erlich': erlich,
  'goldman': goldman,
  'high_pr4t': high_pr4t,
  'average': average,
  'custom_sequencing': custom_sequencing
}

# Function to get the substituting nucleotide based on probabilities.
# TODO: Check if correct from the paper.
def getSubNucleotide(nuc, method):
  rand = random()
  method_probs = method_lookup[method]
  probabilities = method_probs[nuc]

  for p in probabilities:
    if rand < p:
      return probabilities[p]


def getSubNucleotideForPCR(nuc):
  rand = random()

  # A or T
  if nuc % 2 == 0:
    if rand < 0.85:
      return (nuc + 1) % 4

    elif rand < 0.97:
      return (nuc + 2) % 4

    else:
      return (nuc + 3) % 4

  # G or C
  else:
    if rand < 0.84:
      return (nuc + 3) % 4

    elif rand < 0.91:
      return (nuc + 2) % 4

    else:
      return (nuc + 1) % 4
