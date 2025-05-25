from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
import numpy as np
from random import randint
import pytest

from state import APState


def check(qc, ap):
    """Check if the quantum circuit state matches the APState."""
    qc_state = Statevector.from_instruction(qc).data
    ap_state = ap.vec()
    assert np.isclose(np.abs(np.vdot(qc_state, qc_state)), 1)
    assert np.isclose(np.abs(np.vdot(ap_state, ap_state)), 1)
    assert np.isclose(np.abs(np.vdot(qc_state, ap_state)), 1)


def apply_gate(qc, ap, gate, *qubits):
    """Apply a gate to the APState and QuantumCircuit."""
    getattr(ap, gate)(*qubits)
    getattr(qc, gate)(*qubits)
    check(qc, ap)


def test_ghz_state():
    qc = QuantumCircuit(3)
    ap = APState(3).h(0).h(1).h(2)

    apply_gate(qc, ap, "h", 0)
    apply_gate(qc, ap, "cx", 0, 1)
    apply_gate(qc, ap, "cx", 0, 2)


@pytest.mark.parametrize("n_qubits", range(2, 12))
@pytest.mark.parametrize("n_gates", range(200, 205))
def test_random_state(n_qubits, n_gates):
    """Test random state generation and gate application."""
    qc = QuantumCircuit(n_qubits)
    ap = APState(n_qubits)
    for qubit in range(n_qubits):
        ap.h(qubit)

    check(qc, ap)
    gates = ["h", "cx", "s", "cz"]
    for _ in range(n_gates):
        gate = gates[randint(0, 3)]
        if gate in ["h", "s"]:
            qubit = randint(0, n_qubits - 1)
            apply_gate(qc, ap, gate, qubit)
        else:
            control = randint(0, n_qubits - 1)
            target = randint(0, n_qubits - 1)
            while target == control:
                target = randint(0, n_qubits - 1)
            apply_gate(qc, ap, gate, control, target)
