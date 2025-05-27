from quaos.paulis import PauliSum, PauliString, Pauli


dimension = 2

X = Pauli(1, 0, dimension)
Z = Pauli(0, 1, dimension)
Id = Pauli(0, 0, dimension)

print(X * Id)
print(Id * X)
print(X * X)
print(X * Z)

print(X + Z)


