BASES = ["A", "T", "G", "C"]
BASES_IUPAC = ["A", "T", "G", "C", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
IUPAC_TO_BASES = {
    "A": ["A"], "T": ["T"], "C": ["C"], "G": ["G"],
    "R": ["A", "G"], "Y": ["C", "T"], 
    "S": ["G", "C"], "W": ["A", "T"], 
    "K": ["G", "T"], "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "T", "C", "G"]
}
COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G", 
    "N": "N", "R": "Y", "Y": "R", "S": "S", 
    "W": "W", "K": "M", "M": "K", "B": "V", 
    "D": "H", "H": "D", "V": "B", ".": ".",
    "[": "]", "]": "[", "N": "N"
}
BASE_TO_ONE_HOT = {
    "A": [1, 0, 0, 0],
    "T": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "C": [0, 0, 0, 1],
    "N": [1, 1, 1, 1],
    ".": [1, 1, 1, 1]
}
BASE_TO_INT = {"A":1, "T":2, "G":3, "C":4, "N":5, "R":6, "Y":7, "S":8, "W":9, "K":10, "M":11, "B":12, "D":13, "H":14, "V":15}
INT_TO_BASE = {i: n for n, i in BASE_TO_INT.items()}
MOD_TYPE_TO_CANONICAL = {
    "m":"C",
    "a":"A",
    "21839":"C"
}
MOD_CODE_TO_PRETTY = {
    "m":"5mC",
    "a":"6mA",
    "21839":"4mC"
}
MOD_PRETTY_TO_CODE = {v: k for k, v in MOD_CODE_TO_PRETTY.items()}