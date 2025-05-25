# Very basic Clifford state simulator using AP form

from dataclasses import dataclass, field
from itertools import product
import numpy as np

@dataclass(frozen=True)
class Proj:
    phase2: bool
    qubits: frozenset[int]

    def __post_init__(self):
        if len(self.qubits) == 0:
            raise ValueError("Projector must act on at least one qubit")
        if not isinstance(self.qubits, frozenset):
            object.__setattr__(self, 'qubits', frozenset(self.qubits))

    def __repr__(self):
        return f"Proj({self.phase2}, {set(self.qubits)})"


@dataclass(frozen=True)
class CZ:
    qubits: frozenset[int]

    def __post_init__(self):
        if len(self.qubits) != 2:
            raise ValueError("CZ gate must be between two qubits")
        if not isinstance(self.qubits, frozenset):
            object.__setattr__(self, 'qubits', frozenset(self.qubits))

    def __repr__(self):
        return f"CZ({set(self.qubits)})"


@dataclass
class APState:
    n_qubits: int
    projectors: set[Proj] = field(default_factory=set)
    phases4: list[int] = field(init=False)
    czs: set[CZ] = field(default_factory=set)

    def __post_init__(self):
        self.phases4 = [0] * self.n_qubits

    def h(self, qubit: int):
        """Apply Hadamard gate to a qubit."""
        # find all projectors acting on that qubit
        connected_projs = [p for p in self.projectors if qubit in p.qubits]
        connected_czs = [cz for cz in self.czs if qubit in cz.qubits]
        self.projectors -= set(connected_projs)
        self.czs -= set(connected_czs)

        # if there are no projectors, add a new one
        if not connected_projs:
            if self.phases4[qubit] % 2 == 1:
                self.phases4[qubit] = (self.phases4[qubit] + 2) % 4
            else:
                self.projectors.add(Proj(
                    phase2=self.phases4[qubit] == 2,
                    qubits={qubit}
                ))
        else:
            proj, projs = connected_projs[0], connected_projs[1:]
            other_qubits = proj.qubits - {qubit}
            # toggle the phases and connections of the other projectors
            self.projectors |= {
                Proj(
                    phase2=p.phase2 ^ proj.phase2,
                    qubits=p.qubits ^ proj.qubits
                ) for p in projs
            }
            if self.phases4[qubit] % 2 == 0:
                self.czs ^= {CZ({qubit, q}) for q in other_qubits}
            else:
                self.czs ^= {
                    CZ({q1, q2})
                    for q1 in proj.qubits
                    for q2 in proj.qubits if q1 < q2
                }
            for q in other_qubits:
                self.phases4[q] += self.phases4[qubit] * (-1) ** proj.phase2
                self.phases4[q] %= 4
            self.phases4[qubit] = int(proj.phase2) * 2

        # pivoting: cz become cx gates
        cz_qubits = {q for cz in connected_czs for q in cz.qubits} - {qubit}
        for q in cz_qubits:
            self.cx(q, qubit)
        return self

    def cx(self, qubit1: int, qubit2: int):
        """Apply CX gate between two qubits."""
        if qubit1 == qubit2:
            raise ValueError("CX gate cannot be applied to the same qubit")
        connected_projs = [p for p in self.projectors if qubit2 in p.qubits]
        connected_czs = [cz for cz in self.czs if qubit2 in cz.qubits]
        cz_qubits = {q for cz in connected_czs for q in cz.qubits} - {qubit2}

        self.projectors -= set(connected_projs)
        self.projectors |= {
            Proj(
                phase2=p.phase2,
                qubits=p.qubits ^ {qubit1}
            ) for p in connected_projs
        }
        self.czs ^= {
            CZ({qubit1, q}) for q in cz_qubits if q != qubit1
        }
        if qubit1 in cz_qubits:
            self.phases4[qubit1] = (self.phases4[qubit1] + 2) % 4
        if self.phases4[qubit2] % 2 == 1:
            self.phases4[qubit1] += self.phases4[qubit2]
            self.phases4[qubit1] %= 4
            self.czs ^= {CZ({qubit1, qubit2})}
        elif self.phases4[qubit2] == 2:
            self.phases4[qubit1] = (self.phases4[qubit1] + 2) % 4
        return self

    def cz(self, qubit1: int, qubit2: int):
        """Apply CZ gate between two qubits."""
        self.czs ^= {CZ({qubit1, qubit2})}
        return self

    def s(self, qubit: int):
        """Apply S gate to a qubit."""
        self.phases4[qubit] = (self.phases4[qubit] + 1) % 4
        return self

    def print(self):
        """Print the state in basis."""
        vec = self.vec()
        for bits in list(product([0, 1], repeat=self.n_qubits)):
            index = sum(b << i for i, b in enumerate(bits))
            if vec[index] != 0:
                print(f"{''.join(map(str, bits))}: {vec[index]}")

    def vec(self):
        """Convert the state to a vector."""
        vec = np.zeros(2 ** self.n_qubits, dtype=complex)
        w = [1, 1j, -1, -1j]

        for bits in product([0, 1], repeat=self.n_qubits):
            if all(sum(bits[q] for q in proj.qubits) % 2 == proj.phase2 for proj in self.projectors):
                phase4 = sum(self.phases4[q] for q in range(self.n_qubits) if bits[q])
                phase4 += sum(2 * all(bits[q] for q in cz.qubits) for cz in self.czs)
                vec[sum(b << i for i, b in enumerate(bits))] = w[phase4 % 4]

        return vec / np.linalg.norm(vec)  # Normalize the vector
